#pragma once
//
// GOPack C++ -- the repack engine.
//
// This is a direct, faithful port of code/@GOPacker/GOPacker.m and its
// methods (Gerald Orick's sparse-matrix repacking engine, as used in
// GOPack). Field names and method names intentionally mirror the MATLAB
// source closely (rather than being renamed to "idiomatic C++") so the two
// can be read side by side and diffed for correctness.
//
// INDEXING: vertices are 1-indexed everywhere, matching MATLAB. All
// per-vertex arrays are sized nodeCount+1 with index 0 unused. This is a
// deliberate fidelity choice: it keeps every loop bound and offset
// identical to the original, which is what makes this translation checkable
// against the source rather than a reinterpretation of it.
//
// CURRENT PORT SCOPE: the "max pack" pipeline (mode 1: maximal packing of
// the disc/plane or the sphere) is fully ported end to end -- readpack (the
// *.p FLOWERS format), complex_count, indxMatrices, setMode, layoutBdry /
// setHoroCenters, the continueRiffle/layoutCenters/setEffective iteration,
// and writepack/writeEucl with hyperbolic/spherical output conversion.
// Polygonal and rectangle packing mode (mode 2: setMode's polygonal branch,
// setPolyCenters, setRectCenters, getAspect) is also ported: layoutBdry
// dispatches to setPolyCenters/setRectCenters instead of setHoroCenters when
// mode==2, and the shared continueRiffle/layoutCenters/setEffective
// iteration handles both modes without change (setEffective already
// branches on the sign of vAims, which setMode sets appropriately per
// mode). The bare triangle-list ingestion (parseTriangles, parse_triangles.m)
// and orphan-vertex removal (pruneComplex, pruneComplex.m) are also ported,
// along with the pure boundary-sampling helper rand_bdry_pts.m (see
// gopack::geom::randBdryPts in Geometry.h). The random-triangulation-
// *generation* family (randTriangulation/randomDisc/randomSphere/
// randomRectangle/randomSquare/randomTri) is also ported -- see the
// randomDisc()/randomSphere()/randomRectangle()/randomSquare()/randomTri()
// static factory methods below, built on gopack::geom::delaunayPlane()/
// convexHull3()/randTriangulationSphere()/randTriangulationPlane()
// (RandomGen.h), themselves wrapping the vendored Triangle/Qhull libraries
// (see third_party/*/README.md). These factory methods (and RandomGen.h
// itself) only exist when GOPACK_HAVE_RANDOM_GEN is defined (see
// core/CMakeLists.txt's GOPACK_BUILD_RANDOM_GEN option, default ON).
//
// INPUT: readpack() (the *.p FLOWERS format) is one way to get a complex
// into a Packer; loadComplex() is the other -- it has no MATLAB counterpart,
// and exists purely for embedding (JNI/library callers that already hold
// the triangulation in memory, e.g. CirclePack's own per-vertex flower
// data, and would otherwise have to serialize it to text just to have this
// library immediately parse that text back out again). Both funnel into the
// same ingestFlowers()/finalizeComplex() logic, so a complex loaded either
// way behaves identically from that point on.
#include <array>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "gopack/Geometry.h"
#include "gopack/SparseLinearSolver.h"
#include "gopack/Types.h"

#ifdef GOPACK_HAVE_RANDOM_GEN
#include "gopack/RandomGen.h"
#endif

namespace gopack {

class NotImplementedError : public std::logic_error {
public:
    explicit NotImplementedError(const std::string& what) : std::logic_error(what) {}
};

struct RiffleResult {
    int cycles = 0;
    Scalar visMon = 0.0;
    long elapsedMs = 0;
};

// Return type for parseTriangles(); mirrors parse_triangles.m's
// [facecount,newIndx,oldIndx] outputs.
struct ParseTrianglesResult {
    Index faceCount = 0;
    // 1-indexed, size top+1 (index 0 unused; 'top' = the maximal vertex
    // number encountered in the input tList). newIndx[oldVertex] = the
    // vertex's renumbered (contiguous, 1..nodeCount) index, or 0 if that old
    // vertex was never encountered in a face connected to the first face
    // (e.g. the triangulation was disconnected, or the old numbering simply
    // had gaps).
    std::vector<Index> newIndx;
    // 1-indexed, size nodeCount+1 (index 0 unused). oldIndx[newVertex] = the
    // original vertex number -- the inverse of newIndx restricted to the
    // vertices that made it into the result.
    std::vector<Index> oldIndx;
};

class Packer {
public:
    Packer() { cleanse(); }

    // ---- combinatorial / packing-mode string tables (GOPacker.PACKMODES / HES) ----
    static constexpr const char* kPackModes[] = {"max_pack", "polygonal", "frozen_bdry",
                                                  "orth", "fixed_corners"};
    static constexpr const char* kHesNames[] = {"hyp", "eucl", "sph"};

    // ==== fields, parallel to GOPacker.m properties ====
    std::string fileName;
    Geometry hes = Geometry::Euclidean;
    Index alpha = 0;
    Index gamma = 0;
    Index nodeCount = 0;
    Index edgeCount = 0;
    Index faceCount = 0;

    std::vector<Index> vNum;                 // vNum[v] = petal count (1-indexed)
    std::vector<std::vector<Index>> flowers; // flowers[v] closed list iff interior
    std::vector<int> bdryFlags;              // 1 = open/bdry, 0 = closed/interior
    std::vector<Scalar> vAims;               // target angle sums
    std::vector<Index> vlist;                // utility vertex list

    std::vector<Index> intVerts;
    std::vector<Index> bdryList;   // closed cclw list, size bdryCount+1
    std::vector<Index> orphanVerts;
    Index intCount = 0;
    Index bdryCount = 0;
    Index orphanCount = 0;

    std::vector<Index> layoutVerts;
    std::vector<Index> rimVerts;
    std::vector<Index> v2indx; // true index -> matrix index (0 if unset)
    std::vector<Index> indx2v; // matrix index -> true index (1-indexed access: indx2v[k])

    std::vector<std::pair<Index, Index>> orphanEdges;

    std::vector<Scalar> origRadii;
    std::vector<Complex> origCenters;
    std::vector<Scalar> radii;
    std::vector<Complex> centers;
    std::vector<Scalar> localradii;
    std::vector<Complex> localcenters;

    // sparse system assembly data (see indxMatrices/layoutCenters)
    std::vector<Index> tranI, tranJ, tranJindx;
    Index tranIJcount = 0;
    std::vector<Index> rhsI, rhsJ, rhsJindx;
    Index rhsIJcount = 0;

    SparseMatrix transMatrix; // layCount x layCount
    SparseMatrix rhsMatrix;   // layCount x colCount
    std::vector<Complex> rhs; // rhs = rhsMatrix * zb (boundary centers), 1-indexed by k

    std::vector<std::vector<Scalar>> inRadii; // inRadii[k] indexed by layout-index k
    std::vector<Scalar> conduct;              // conduct[k] indexed by layout-index k

    std::vector<Index> corners; // polygonal mode: corner vertices, in cclw order
    std::vector<std::vector<Index>> sides; // sides[i] = closed bdry-vertex run from corners[i] to corners[i+1]

    int mode = 1; // 1 = max pack, 2 = polygonal

    std::vector<Scalar> angsumMonitor;
    std::vector<Scalar> l2Monitor;
    std::vector<Scalar> visErrMonitor;
    std::vector<long> ticMonitor;

    SolverStrategy solverStrategy = SolverStrategy::Auto;

    // ==== methods, parallel to GOPacker.m methods ====

    // cleanse.m
    void cleanse();

    // complex_count.m -- returns nodeCount, or -1 on error.
    Index complexCount();

    // pruneComplex.m -- removes 'orphan' vertices (vertices outside the
    // interior component and its immediate boundary neighbors) after a raw
    // triangulation has been ingested via parseTriangles(). Only meaningful
    // after complexCount() has already run (parseTriangles() calls it
    // internally) and populated intVerts/bdryList/orphanCount -- a no-op
    // (returns 0) if orphanCount is already 0. Renumbers everything
    // (flowers, vAims, intVerts, bdryList, alpha, gamma, vlist,
    // origRadii/origCenters/radii/centers/localradii/localcenters) to the
    // trimmed, contiguous 1..newNodeCount range and calls complexCount()
    // again internally; callers must still call indxMatrices() themselves
    // afterward (matching pruneComplex.m, whose callers always do so
    // separately). Returns the number of vertices cut.
    //
    // Two deliberate deviations from the literal MATLAB source (see the .cpp
    // definition for the full rationale): (1) the boundary-flower trimming
    // loop reads from a saved local copy of the old vlist rather than from
    // obj.vlist after it's been cleared, which is real undefined behavior
    // (out-of-bounds access) in a literal C++ port, not just a MATLAB quirk;
    // (2) 'origCenters(nv)=obj.origCenters(nv)' is fixed to read
    // obj.origCenters(v) (the old-numbered vertex), matching the parallel
    // origRadii line right above it and every other line in that loop --
    // this reads as an unambiguous typo, not an intentional behavior.
    Index pruneComplex();

    // indxMatrices.m
    void indxMatrices(const std::vector<Index>& varlist = {});

    // setHoroCenters.m
    void setHoroCenters();

    // continueRiffle.m
    Index continueRiffle(int passNum);

    // riffle.m
    RiffleResult riffle(int passNum = 20);

    // layoutBdry.m
    void layoutBdry();

    // setRectCenters.m / setPolyCenters.m -- polygonal/rectangle boundary
    // layout. setPolyCenters dispatches to setRectCenters itself when there
    // are 4 corners all within 1e-5 of a right angle (the default for
    // sideN==4 with no custom corner angles); layoutBdry() is what actually
    // calls setPolyCenters when mode==2.
    void setRectCenters();
    void setPolyCenters();

    // setMode.m. mdIn: 1 = max pack, 2 = polygonal. For mode 2, 'crns' is an
    // optional list of corner vertices (by original vertex index, in cclw
    // order); if empty, corners are inferred from 'vlist' or else chosen
    // pseudo-randomly, mirroring GOPacker.m's nargin<3 behavior (a by-value
    // C++ vector can't distinguish "omitted" from "explicitly empty" the
    // way MATLAB's nargin can, so an empty 'crns' here always means
    // "figure out the corners for me", the more useful default). 'angs' is
    // an optional matching list of corner target angles. NOTE: entering
    // mode 2 always resets 'hes' to Euclidean (a deliberate deviation from
    // GOPacker.m -- see the comment in setMode's .cpp definition), even if
    // the packing was originally read as hyperbolic or spherical.
    int setMode(int mdIn, const std::vector<Index>& crns = {}, const std::vector<Scalar>& angs = {});

    // layoutCenters.m
    void layoutCenters();

    // setEffective.m
    void setEffective();

    // reapResults.m. NOTE: for Spherical packings, this also recenters
    // centers/radii (via the same affine normalization writepack() uses for
    // its spherical output) so the packing's tangency-point centroid sits
    // near the origin in 3D -- not part of GOPacker.m's reapResults, which
    // has no such step; see the .cpp definition for the full rationale.
    // Callers reading centers/radii right after riffle() (including the JNI
    // bridge) get this for free, without needing to go through writepack().
    void reapResults();

    // getAspect.m -- ratio (top+bot)/(left+right) side lengths; mode 2 with
    // exactly 4 corners only. Returns -1 (and prints a diagnostic) otherwise.
    Scalar getAspect();

    // angsumErrors.m
    std::vector<Scalar> angsumErrors() const { return angsumErrorsImpl(nullptr); }
    std::vector<Scalar> angsumErrors(const std::vector<Scalar>& radiiOverride) const {
        return angsumErrorsImpl(&radiiOverride);
    }

    // packStatus.m
    void packStatus() const;

    // visualErrors.m
    std::vector<Scalar> visualErrors() const;

    // updateVdata.m
    void updateVdata();

    // readpack.m -- reads the *.p (FLOWERS) format described in
    // docs/GO_Formats.txt. A bare triangle-list input (no NODECOUNT:/
    // CHECKCOUNT: header) is reported and rejected here rather than
    // dispatched to parseTriangles() automatically -- readpack()'s file
    // format and parseTriangles()'s in-memory Nx3 matrix are different
    // enough (and parseTriangles() needs its own tList/cents arguments,
    // which a bare vertex-triple text stream doesn't cleanly map to) that
    // callers with a raw triangle list should call parseTriangles()
    // directly instead. The OFF format fallback is still NOT YET PORTED.
    Index readpack(const std::string& fname);

    // loadComplex -- ingest an already-known combinatorial complex directly,
    // without going through the *.p text format at all. This is the same
    // logic readpack() uses after parsing (alpha resolution, complex_count,
    // default radii/centers/vAims, indxMatrices), factored out so a caller
    // that already holds the triangulation in memory -- e.g. a JNI caller
    // passing CirclePack's own per-vertex flower data -- never has to
    // serialize it to text and pay for readpack()'s parsing just to hand it
    // straight back. Intended as the library entry point for JNI/embedding
    // use; readpack() remains the entry point for CLI/file use and is
    // implemented in terms of this same code path, so the two can never
    // silently diverge in behavior.
    //
    // flowersIn: 1-indexed, size nodeCountIn+1 (index 0 unused/ignored).
    // flowersIn[v] is v's petal list in the *.p FLOWERS format convention:
    // CLOSED (front()==back()) iff v is an interior vertex, OPEN
    // (front()!=back()) iff v is a boundary vertex. vNum/bdryFlags/bdryCount
    // are derived from this, not taken as separate inputs, so they can never
    // disagree with the flowers themselves.
    // geometryIn: Euclidean/Hyperbolic/Spherical, matching the *.p file's
    // GEOMETRY: line.
    // alphaIn/gammaIn: pass 0 for "let GOPack pick" (matching an absent
    // ALPHA/GAMMA: line); a negative alphaIn forces auto-search the way
    // readpack()'s "alpha < 0" branch does. Most callers should just pass 0
    // for both.
    // initRadii/initCenters: optional (nullptr = omitted, matching an absent
    // RADII:/CENTERS: section) starting circle data, sized nodeCountIn+1,
    // 1-indexed. When omitted, defaults match readpack() (radii 0.5, centers
    // 0). NOTE: mirrors readpack()'s own (slightly surprising) behavior of
    // only applying initRadii/initCenters for Euclidean or Hyperbolic input
    // -- Spherical input falls back to the defaults even if provided, same
    // as reading a spherical *.p file with a RADII:/CENTERS: section.
    // vlistIn: optional (nullptr = omitted) utility vertex list, matching
    // the *.p file's VERT_LIST: section (used e.g. as a corner-selection
    // hint by setMode's mode-2 branch).
    // vAimsIn: optional (nullptr = omitted) target angle sums, matching the
    // *.p file's ANGLE_AIMS: section; when omitted, defaults to the usual
    // 2*pi interior / -1 boundary convention.
    // label: optional display name for the "packing is loaded" console
    // message (mirroring readpack()'s use of the file's base name); pass ""
    // to suppress that message entirely, which most in-process callers will
    // want.
    //
    // Returns nodeCount on success, 0 on error (same convention as
    // readpack()).
    Index loadComplex(Index nodeCountIn, const std::vector<std::vector<Index>>& flowersIn,
                       Geometry geometryIn, Index alphaIn = 0, Index gammaIn = 0,
                       const std::vector<Scalar>* initRadii = nullptr,
                       const std::vector<Complex>* initCenters = nullptr,
                       const std::vector<Index>* vlistIn = nullptr,
                       const std::vector<Scalar>* vAimsIn = nullptr,
                       const std::string& label = "");

    // parse_triangles.m -- ingest a bare triangle list (an Nx3 matrix of
    // vertex numbers, one row per face) directly, without going through the
    // *.p FLOWERS text format at all. This is GOPack's entry point for
    // triangulations produced by a Delaunay/convex-hull step (the
    // randTriangulation.m family) rather than read from a file; it is
    // functionally a third way to get a complex into a Packer, alongside
    // readpack() and loadComplex(), but with its own quirks preserved from
    // the source (see below) since real callers depend on them.
    //
    // tList: one entry per face, {v1,v2,v3} in any consistent winding (the
    // algorithm determines overall orientation from the first face and
    // re-orients any face it discovers as reversed relative to that -- this
    // is why tList is taken BY VALUE rather than by reference: the source
    // mutates the triangle list in place as part of this walk, and passing
    // by value keeps that mutation from leaking out to the caller's own
    // copy). Vertex numbers need not be contiguous from 1 or start at 1;
    // parseTriangles() figures out the actual range encountered and
    // contiguously renumbers on output (see ParseTrianglesResult).
    // cents: optional (nullptr = omitted) center for every old/original
    // vertex number from 1 to the maximal vertex number ('top') encountered
    // in tList -- sized top+1, index 0 unused, matching this project's usual
    // 1-indexed convention (note this differs from parse_triangles.m's own
    // plain length-top MATLAB vector -- there is no "index 0" to skip in
    // MATLAB). These may represent (theta,phi) polar coordinates in the
    // spherical case, exactly as read from randTriangulation.m's sphere
    // branch. When omitted, all centers default to 0 (matching
    // finalizeComplex()'s default), to be filled in by a later riffle().
    //
    // Sets alpha from whatever the 'alpha' member was already set to before
    // the call (0 meaning "let GOPack pick a deep interior vertex, matching
    // the source's isempty/==0 check -- NOT the same as loadComplex()'s
    // convention where a *negative* alpha also forces auto-search); this
    // mirrors parse_triangles.m's own holdalpha/cleanse()/restore pattern,
    // under which a caller sets obj.alpha before calling parse_triangles,
    // exactly as randomRectangle.m does with the alpha returned by
    // randTriangulation.m.
    //
    // Unlike loadComplex()/readpack(), this does NOT call indxMatrices()
    // itself -- parse_triangles.m's real callers (randomSphere.m,
    // randomDisc.m, randomRectangle.m) all call indxMatrices() separately
    // afterward, sometimes after also calling pruneComplex() first (for
    // triangulations of a bounded region, which can produce orphan
    // vertices cut off from the interior). Radii default to 0.5 (matching
    // finalizeComplex()'s convention) and are always set here directly
    // (parse_triangles.m has no analog of loadComplex()'s
    // hyperbolic-conversion path -- triangulations are always fresh,
    // never already-hyperbolic input).
    ParseTrianglesResult parseTriangles(std::vector<std::array<Index, 3>> tList,
                                         const std::vector<Complex>* cents = nullptr);

#ifdef GOPACK_HAVE_RANDOM_GEN
    // ==== "geometrically random" triangulation generators -- see
    // core/src/PackerRandom.cpp and randomDisc.m/randomSphere.m/
    // randomRectangle.m/randomSquare.m/randomTri.m. Each builds a fresh,
    // ready-to-riffle Packer (parseTriangles()/pruneComplex()/
    // indxMatrices() already called, alpha/hes/mode already set) from a
    // Delaunay triangulation (gopack::geom::delaunayPlane()) or 3D convex
    // hull (gopack::geom::convexHull3()) of randomly (Poisson Point
    // Process) chosen points. Only declared/defined when
    // GOPACK_HAVE_RANDOM_GEN is set (see core/CMakeLists.txt's
    // GOPACK_BUILD_RANDOM_GEN option, default ON).

    // randomDisc.m -- a random triangulation of the unit disc with 'N'
    // total (interior+boundary) points, set up for hyperbolic (Poincare
    // disc) maximal packing (hes=Hyperbolic, mode=1).
    static Packer randomDisc(Index N);

    // randomSphere.m -- a random triangulation of the sphere with 'intN'
    // points (spherical maximal packing, hes=Spherical, mode=1). Prints a
    // diagnostic and returns a default-constructed (empty) Packer if
    // intN<4.
    static Packer randomSphere(Index intN);

    // randomRectangle.m -- a random triangulation of the rectangle
    // [-aspect,aspect]x[-1,1] with 'intN' interior points, already set up
    // for rectangle packing (setMode(2, ...) with the 4 corners closest to
    // the rectangle's actual corners, hes=Euclidean per setMode's own
    // deviation -- see setMode()'s doc comment). 'bdryN'<=0 (the default)
    // means "compute the default boundary point count from 'aspect'",
    // matching randomRectangle.m's nargin<3 branch; a positive 'bdryN'
    // overrides it, matching nargin==3. Prints a diagnostic and returns a
    // default-constructed (empty) Packer if intN<1. Calls pruneComplex()
    // after parseTriangles() (matching randomRectangle.m) to clean up any
    // interior point that ends up outside the *sampled* boundary polygon
    // randTriangulationPlane() actually triangulates against -- see that
    // function's doc comment in RandomGen.h for why a sparse 'bdryN' can let
    // this happen even for a point genuinely inside the rectangle.
    static Packer randomRectangle(Index intN, Scalar aspect = 1.0, Index bdryN = -1);

    // randomSquare.m -- thin wrapper: a random triangulation of the unit
    // square (aspect 1) with 'N' total (interior+boundary) points. Equal to
    // randomRectangle() with intN/bdryN split by the same interior/boundary
    // formula randomDisc() uses.
    static Packer randomSquare(Index N);

    // randomTri.m, nargin==1 overload -- a random triangulation of the
    // sphere; equivalent to randomSphere(intN) with hes explicitly
    // (redundantly, matching the source) set to Spherical.
    static Packer randomTri(Index intN);

    // randomTri.m, nargin>=4 overload -- a random triangulation of a plane
    // region bounded by the closed polygonal path 'graph' (matching
    // randTriangulationPlane()'s 'graph' parameter), with 'intN' interior
    // and 'bdryN' boundary points, set up for euclidean maximal packing
    // (hes=Euclidean, mode=1) -- NOT polygonal/rectangle mode (unlike
    // randomRectangle(), this is a generic region with no notion of
    // "corners"; a caller wanting polygonal mode can call setMode(2, ...)
    // on the result itself). 'cent', if given (non-null) and inside
    // 'graph', is used as the packing's alpha vertex.
    //
    // Matching randomTri.m, this does NOT call pruneComplex() (unlike
    // randomRectangle()) -- so nodeCount can come out slightly below
    // intN+bdryN for a sparse 'bdryN' relative to a highly concave 'graph';
    // see randTriangulationPlane()'s doc comment in RandomGen.h for why.
    // Callers who want that cleanup can call pruneComplex()+indxMatrices()
    // themselves on the result.
    //
    // Bug fix vs. randomTri.m: when 'cent' is omitted or falls outside
    // 'graph', the source leaves its local GOPacker's alpha explicitly set
    // to -1 (from an unconditional 'gop.alpha=-1;' at the top of the
    // function) rather than 0 -- but parse_triangles.m's own "pick a deep
    // alpha for me" auto-selection only triggers on alpha==0, not on a
    // negative alpha, so that MATLAB code path hands back a Packer with an
    // invalid (never auto-resolved) alpha=-1, whereas every other caller of
    // parse_triangles.m either supplies a valid alpha or relies on exactly
    // that 0-triggered auto-selection. This reads as an unambiguous
    // oversight (randomRectangle.m, which has the same "cent may not land
    // inside the region" situation, correctly leaves its Packer's alpha at
    // its GOPacker()-default 0 in that case, letting the auto-selection or
    // its own closest-to-origin fallback run), not intentional behavior, so
    // here alpha is simply left at its Packer()-default 0 when 'cent' isn't
    // placed, letting parseTriangles()'s own auto-selection produce a valid
    // vertex instead.
    static Packer randomTri(Index intN, Index bdryN, const std::vector<Complex>& graph,
                             const Complex* cent = nullptr);
#endif // GOPACK_HAVE_RANDOM_GEN

    // loadTangency.m
    std::vector<Complex> loadTangency(const std::vector<Complex>& centersIn,
                                       const std::vector<Scalar>& radiiIn) const;

    // writepack.m / writeEucl.m
    Index writepack(const std::string& fname, bool euclFlag = false);
    Index writeEucl(const std::string& fname);

    // FarVert.m
    Index farVert(const std::vector<Index>& seeds) const;

    bool isInterior(Index v) const { return flowers[v].front() == flowers[v].back(); }

private:
    std::vector<Scalar> angsumErrorsImpl(const std::vector<Scalar>* radiiOverride) const;

    // Shared building blocks behind readpack()/loadComplex() -- see the
    // .cpp for the full rationale. ingestFlowers populates
    // flowers/vNum/bdryFlags/bdryCount/intCount/faceCount/edgeCount from a
    // flowers list and resolves alpha (reading/writing the `alpha`/`gamma`
    // members, which the caller must have already set to the desired
    // starting value, 0 for "let GOPack pick"). finalizeComplex does
    // everything readpack() does after that: complexCount(), default/given
    // radii+centers, default/given vAims, indxMatrices(), and the optional
    // console message.
    Index ingestFlowers(Index nodeCountIn, const std::vector<std::vector<Index>>& flowersIn);
    Index finalizeComplex(const std::vector<Scalar>* initRadii,
                           const std::vector<Complex>* initCenters,
                           const std::vector<Scalar>* vAimsIn, const std::string& label);
};

} // namespace gopack

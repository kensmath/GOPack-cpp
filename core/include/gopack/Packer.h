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
// mode). The triangle-list and OFF input readers (parse_triangles), orphan
// removal (pruneComplex), and random-triangulation generation are
// understood (see the corresponding .m files) but not yet ported; calling
// those methods throws NotImplementedError rather than guessing.
//
// INPUT: readpack() (the *.p FLOWERS format) is one way to get a complex
// into a Packer; loadComplex() is the other -- it has no MATLAB counterpart,
// and exists purely for embedding (JNI/library callers that already hold
// the triangulation in memory, e.g. CirclePack's own per-vertex flower
// data, and would otherwise have to serialize it to text just to have this
// library immediately parse that text back out again). Both funnel into the
// same ingestFlowers()/finalizeComplex() logic, so a complex loaded either
// way behaves identically from that point on.
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "gopack/Geometry.h"
#include "gopack/SparseLinearSolver.h"
#include "gopack/Types.h"

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

    // pruneComplex.m -- NOT YET PORTED (only used for random rect/square
    // generation in the original; not required for reading real packing
    // files, so deferred).
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

    // reapResults.m
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
    // docs/GO_Formats.txt. The triangle-list / OFF fallback paths
    // (parse_triangles.m) are NOT YET PORTED.
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

    // parse_triangles.m -- NOT YET PORTED.
    Index parseTriangles();

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

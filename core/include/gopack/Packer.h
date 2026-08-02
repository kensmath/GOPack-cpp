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
// Polygonal and rectangle packing modes (setPolyCenters/setRectCenters),
// the triangle-list and OFF input readers (parse_triangles), orphan removal
// (pruneComplex), and random-triangulation generation are understood (see
// the corresponding .m files) but not yet ported; calling those methods
// throws NotImplementedError rather than guessing.
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

    std::vector<Index> corners; // polygonal mode: corner vertices (deferred feature)
    std::vector<std::vector<Index>> sides;

    int mode = 1; // 1 = max pack, 2 = polygonal (not yet ported)

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

    // setRectCenters.m / setPolyCenters.m -- NOT YET PORTED (polygonal mode).
    void setRectCenters();
    void setPolyCenters();

    // setMode.m
    int setMode(int mdIn, const std::vector<Index>& crns = {}, const std::vector<Scalar>& angs = {});

    // layoutCenters.m
    void layoutCenters();

    // setEffective.m
    void setEffective();

    // reapResults.m
    void reapResults();

    // getAspect.m -- NOT YET PORTED (polygonal/rectangle mode only).
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
};

} // namespace gopack

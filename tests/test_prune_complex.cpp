// Regression test for Packer::pruneComplex() -- orphan-vertex removal
// (pruneComplex.m), used after ingesting a raw triangulation (via
// parseTriangles()) that may have vertices cut off from the interior
// component (e.g. random rectangle/square generation, or any Delaunay
// triangulation of a non-convex region -- see randomRectangle.m).
//
// Strategy: build (via the already-verified loadComplex() entry point) the
// classical "hex flower" complex (see test_loadcomplex.cpp) with one extra
// "flap" vertex (8) hanging off boundary edge (2,3) -- a single face (2,3,8)
// whose apex has no interior neighbor, making it exactly one orphan vertex
// by construction -- and confirm pruneComplex() removes it, leaving a valid,
// still-convergent hex-flower packing.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "gopack/Packer.h"

namespace {

void checkTrue(bool cond, const char* what) {
    if (!cond) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        std::exit(1);
    }
}

} // namespace

int main() {
    using gopack::Complex;
    using gopack::Geometry;
    using gopack::Index;
    using gopack::Packer;

    std::vector<std::vector<Index>> flowers(9);
    flowers[1] = {2, 3, 4, 5, 6, 7, 2}; // closed hub
    flowers[2] = {7, 1, 3, 8};          // open, extended with flap vertex 8
    flowers[3] = {8, 2, 1, 4};          // open, extended with flap vertex 8
    flowers[4] = {3, 1, 5};
    flowers[5] = {4, 1, 6};
    flowers[6] = {5, 1, 7};
    flowers[7] = {6, 1, 2};
    flowers[8] = {2, 3}; // flap tip: only face (2,3,8), open

    Packer p;
    Index nc = p.loadComplex(8, flowers, Geometry::Euclidean, /*alphaIn=*/1, /*gammaIn=*/0);
    checkTrue(nc == 8, "loadComplex should succeed with 8 vertices");
    checkTrue(p.orphanCount == 1, "expected exactly 1 orphan vertex (the flap tip)");
    checkTrue(p.orphanVerts.size() == 1 && p.orphanVerts[0] == 8, "the orphan should be vertex 8");
    checkTrue(p.intCount == 1, "intCount should still be 1 before pruning");
    checkTrue(p.bdryCount == 6, "bdryCount should still be 6 before pruning");

    Index cutCount = p.pruneComplex();
    checkTrue(cutCount == 1, "pruneComplex should report cutting 1 vertex");
    checkTrue(p.nodeCount == 7, "nodeCount should be 7 after pruning");
    checkTrue(p.orphanCount == 0, "orphanCount should be 0 after pruning");
    checkTrue(p.intCount == 1, "intCount should still be 1 after pruning");
    checkTrue(p.bdryCount == 6, "bdryCount should still be 6 after pruning");
    checkTrue(p.flowers.size() == static_cast<size_t>(p.nodeCount) + 1,
              "flowers should be resized to nodeCount+1");
    checkTrue(std::abs(p.radii[1] - 0.5) < 1e-12,
              "radii should have been translated (not left default-zero) for surviving vertices");

    // pruneComplex() doesn't call indxMatrices() itself (matches
    // pruneComplex.m -- its real callers, randomRectangle.m/randomSquare.m,
    // always call it separately afterward); do that, then confirm the
    // pruned complex is still a usable, convergent packing.
    int md = p.setMode(1);
    checkTrue(md == 1, "setMode(1) should succeed after pruning");
    p.indxMatrices();
    gopack::RiffleResult rr = p.riffle(100);
    checkTrue(rr.cycles >= 0, "riffle should not report an error after pruning");
    auto diffs = p.angsumErrors();
    checkTrue(diffs.size() == 1, "one interior angle-sum entry after pruning");
    checkTrue(std::abs(diffs[0]) < 1e-4, "interior angle sum should converge near 2*pi after pruning");

    // pruneComplex() on an already-orphan-free complex should be a
    // documented no-op (returns 0 immediately, matching pruneComplex.m's own
    // early-out when orphanCount==0).
    Index cutCount2 = p.pruneComplex();
    checkTrue(cutCount2 == 0, "pruneComplex should no-op when orphanCount==0");

    std::fprintf(stdout,
                 "prune_complex: OK (cutCount=%d, %d riffle passes, angsum diff=%.3e)\n", cutCount,
                 rr.cycles, diffs[0]);
    return 0;
}

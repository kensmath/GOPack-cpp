// Regression test for Packer::parseTriangles() -- the bare triangle-list
// ingestion path (parse_triangles.m), GOPack's entry point for
// triangulations produced by a Delaunay/convex-hull step (the
// randTriangulation.m family) rather than read from a *.p file.
//
// Two shapes are exercised: a hex-flower fan (one interior vertex, 6
// boundary vertices -- a plane/disc triangulation) and a tetrahedron (a
// closed, boundary-less triangulation -- exercises the Spherical/no-boundary
// path, including complex_count.m's pseudo-boundary anchor triangle
// behavior already covered for readpack()/loadComplex() by
// test_sphere_normalize.cpp).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "gopack/Geometry.h"
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
    using gopack::Geometry;
    using gopack::Index;
    using gopack::Packer;
    using gopack::ParseTrianglesResult;

    // ==== Hex-flower fan: hub vertex 1 surrounded by ring 2..7. ====
    {
        std::vector<std::array<Index, 3>> tList = {
            {1, 2, 3}, {1, 3, 4}, {1, 4, 5}, {1, 5, 6}, {1, 6, 7}, {1, 7, 2},
        };
        Packer p;
        ParseTrianglesResult r = p.parseTriangles(tList);
        checkTrue(p.nodeCount == 7, "nodeCount should be 7");
        checkTrue(r.faceCount == 6, "faceCount should be 6");
        checkTrue(p.intCount == 1, "intCount should be 1");
        checkTrue(p.bdryCount == 6, "bdryCount should be 6");
        checkTrue(p.orphanCount == 0, "orphanCount should be 0");
        checkTrue(p.hes == Geometry::Euclidean, "hes should be Euclidean (has boundary)");
        checkTrue(p.alpha == 1, "alpha should resolve to the hub vertex 1");
        checkTrue(p.isInterior(1), "vertex 1 should be interior (closed flower)");
        for (Index v = 1; v <= 7; ++v) {
            checkTrue(r.newIndx[v] == v, "newIndx should be identity for already-contiguous input");
            checkTrue(r.oldIndx[v] == v, "oldIndx should be identity for already-contiguous input");
        }

        // Don't assume a particular orientation (parseTriangles() derives
        // cw-vs-ccw from the winding of the *input* tList, not from any
        // externally-imposed convention) -- verify the combinatorial
        // structure directly instead: hub 1 is surrounded by exactly
        // {2..7}, and each ring vertex v's open flower is exactly
        // {1, prev(v), next(v)} around the hex ring.
        checkTrue(p.flowers[1].front() == p.flowers[1].back(), "hub flower should be closed");
        std::vector<Index> hubCore(p.flowers[1].begin(), p.flowers[1].end() - 1);
        std::sort(hubCore.begin(), hubCore.end());
        checkTrue((hubCore == std::vector<Index>{2, 3, 4, 5, 6, 7}),
                  "hub flower should surround exactly {2,3,4,5,6,7}");
        for (Index v = 2; v <= 7; ++v) {
            checkTrue(p.flowers[v].front() != p.flowers[v].back(), "ring vertex flower should be open");
            checkTrue(p.flowers[v].size() == 3, "ring vertex flower should have 3 entries");
            Index prevV = (v == 2) ? 7 : v - 1;
            Index nextV = (v == 7) ? 2 : v + 1;
            std::vector<Index> got = p.flowers[v];
            std::sort(got.begin(), got.end());
            std::vector<Index> want = {1, prevV, nextV};
            std::sort(want.begin(), want.end());
            checkTrue(got == want, "ring vertex flower should be exactly {1, prev, next}");
        }

        checkTrue(p.radii.size() == 8, "radii should be sized nodeCount+1");
        checkTrue(std::abs(p.radii[2] - 0.5) < 1e-12, "default radius should be 0.5");

        // parseTriangles() does NOT call indxMatrices() itself (matching
        // parse_triangles.m, whose real callers always call it separately
        // afterward) -- the caller must.
        int md = p.setMode(1);
        checkTrue(md == 1, "setMode(1) should succeed");
        p.indxMatrices();
        gopack::RiffleResult rr = p.riffle(100);
        checkTrue(rr.cycles >= 0, "riffle should not report an error");
        auto diffs = p.angsumErrors();
        checkTrue(diffs.size() == 1, "one interior angle-sum entry");
        checkTrue(std::abs(diffs[0]) < 1e-4, "interior angle sum should converge near 2*pi");

        std::fprintf(stdout, "hex-fan: OK (%d riffle passes, angsum diff=%.3e)\n", rr.cycles,
                     diffs[0]);
    }

    // ==== Tetrahedron: closed, no boundary -- Spherical path. ====
    {
        std::vector<std::array<Index, 3>> tList = {
            {1, 2, 3}, {1, 3, 4}, {1, 4, 2}, {2, 4, 3},
        };
        Packer p;
        ParseTrianglesResult r = p.parseTriangles(tList);
        checkTrue(p.nodeCount == 4, "tetrahedron: nodeCount should be 4");
        checkTrue(r.faceCount == 4, "tetrahedron: faceCount should be 4");
        checkTrue(p.hes == Geometry::Spherical, "tetrahedron: hes should be Spherical (no boundary)");
        // complex_count.m always gives a boundary-less complex a 3-vertex
        // pseudo-boundary anchor triangle (see test_sphere_normalize.cpp for
        // the full rationale) -- not 0.
        checkTrue(p.bdryCount == 3, "tetrahedron: bdryCount should be 3 (pseudo-anchor triangle)");
        checkTrue(p.intCount == 1, "tetrahedron: intCount should be 1 (4 - 3 anchor)");
        checkTrue(p.orphanCount == 0, "tetrahedron: orphanCount should be 0");
        for (Index v = 1; v <= 4; ++v) {
            checkTrue(p.isInterior(v), "tetrahedron: every vertex of a closed complex has a closed flower");
        }
        std::fprintf(stdout, "tetrahedron: OK\n");
    }

    // ==== Empty input: should fail cleanly, not crash. ====
    {
        Packer p;
        std::vector<std::array<Index, 3>> empty;
        ParseTrianglesResult r = p.parseTriangles(empty);
        checkTrue(r.faceCount == 0, "empty tList should yield faceCount 0");
        checkTrue(p.nodeCount == 0, "empty tList should leave nodeCount 0");
        std::fprintf(stdout, "empty input: OK\n");
    }

    std::fprintf(stdout, "parse_triangles: ALL OK\n");
    return 0;
}

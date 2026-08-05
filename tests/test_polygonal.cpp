// Independent sanity/regression test for polygonal/rectangle packing mode
// (setMode(2, ...), setPolyCenters/setRectCenters, getAspect). Reuses the
// same hex-flower combinatorics as test_hex_flower.cpp (1 interior vertex,
// 6 boundary vertices) but this time in mode 2, treating 4 of the 6
// boundary vertices as the corners of a rectangle. With the default equal
// corner angles setMode(2, ...) assigns (pi/2 for 4 corners), this should
// dispatch to setRectCenters(); after riffling, the four corners should
// still sit at right angles from the origin (by symmetry of the rectangle
// layout) and getAspect() should report a finite, positive aspect ratio.
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

    Packer packer;
    packer.nodeCount = 7;
    packer.flowers.assign(8, {});
    packer.flowers[1] = {2, 3, 4, 5, 6, 7, 2};
    packer.flowers[2] = {7, 1, 3};
    packer.flowers[3] = {2, 1, 4};
    packer.flowers[4] = {3, 1, 5};
    packer.flowers[5] = {4, 1, 6};
    packer.flowers[6] = {5, 1, 7};
    packer.flowers[7] = {6, 1, 2};
    packer.alpha = 1;
    packer.gamma = 0;
    packer.hes = Geometry::Euclidean;

    gopack::Index nc = packer.complexCount();
    checkTrue(nc == 7, "complexCount should return nodeCount (7)");
    checkTrue(packer.bdryCount == 6, "expected exactly 6 boundary vertices");

    packer.centers.assign(8, Complex(0.0, 0.0));
    packer.radii.assign(8, 0.5);
    packer.localcenters = packer.centers;
    packer.localradii = packer.radii;

    // readpack() always sets this up before a real caller could reach
    // setMode(2, ...); since this test builds the Packer by hand (bypassing
    // readpack entirely, like test_hex_flower.cpp does), it has to
    // replicate that precondition itself.
    packer.vAims.assign(8, gopack::kTwoPi);
    for (Index k = 1; k <= packer.nodeCount; ++k) {
        if (packer.bdryFlags[k] != 0) packer.vAims[k] = -1.0;
    }

    // 4 of the 6 boundary vertices, in counterclockwise order (matches
    // bdryList's own traversal order: 2,7,6,5,4,3).
    std::vector<Index> corners = {2, 3, 5, 6};
    int md = packer.setMode(2, corners);
    checkTrue(md == 2, "setMode(2, corners) should succeed with mode==2");
    checkTrue(packer.corners.size() == 4, "expected 4 corners to be recorded");
    checkTrue(packer.sides.size() == 4, "expected 4 sides to be recorded");

    packer.indxMatrices();
    gopack::RiffleResult result = packer.riffle(100);
    checkTrue(result.cycles >= 0, "riffle should not report an error");

    // Every center/radius should be finite (no NaN/Inf leaking out of the
    // polygonal boundary layout or the area-based radius updates).
    for (Index v = 1; v <= packer.nodeCount; ++v) {
        checkTrue(std::isfinite(packer.radii[v]), "radius should be finite");
        checkTrue(std::isfinite(packer.centers[v].real()) && std::isfinite(packer.centers[v].imag()),
                  "center should be finite");
        checkTrue(packer.radii[v] > 0.0, "radius should be strictly positive");
    }

    double aspect = packer.getAspect();
    checkTrue(std::isfinite(aspect) && aspect > 0.0, "getAspect() should be finite and positive");

    // The 4 corner vertices should still sit at (approximately) right
    // angles to their neighbors around the rectangle -- check that the
    // quadrilateral they form has two pairs of (nearly) equal, (nearly)
    // parallel opposite sides, i.e. it's still rectangle-shaped after
    // riffling, not degenerate.
    Complex c0 = packer.centers[packer.corners[0]];
    Complex c1 = packer.centers[packer.corners[1]];
    Complex c2 = packer.centers[packer.corners[2]];
    Complex c3 = packer.centers[packer.corners[3]];
    double side01 = std::abs(c1 - c0);
    double side23 = std::abs(c3 - c2);
    double side12 = std::abs(c2 - c1);
    double side30 = std::abs(c0 - c3);
    checkTrue(side01 > 1e-6 && side12 > 1e-6, "rectangle sides should be non-degenerate");
    checkTrue(std::abs(side01 - side23) / side01 < 0.05,
              "opposite rectangle sides should be nearly equal in length");
    checkTrue(std::abs(side12 - side30) / side12 < 0.05,
              "opposite rectangle sides should be nearly equal in length");

    std::fprintf(stdout,
                 "test_polygonal: OK (%d riffle passes, aspect=%.4f, corner side lengths "
                 "%.4f/%.4f/%.4f/%.4f)\n",
                 result.cycles, aspect, side01, side12, side23, side30);
    return 0;
}

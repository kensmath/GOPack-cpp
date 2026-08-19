// Independent sanity/regression test: the classical "hex flower" -- one
// interior vertex of degree 6 surrounded by a ring of 6 boundary vertices.
// By symmetry, GOPack's max-pack mode should converge to a packing where
// (a) the interior vertex's angle sum error goes to ~0, and (b) all six
// boundary radii end up equal to each other (the combinatorics are
// rotationally symmetric, so any correct solver should find a symmetric
// fixed point). This is not taken from GOPack's own test suite (it doesn't
// appear to have one in the repository) -- it's a standard example from
// circle-packing theory used here as an independent correctness check on
// the ported riffle/layoutCenters/setEffective pipeline.
#include <cmath>
#include <cstdio>
#include <cstdlib>

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
    using gopack::Packer;

    Packer packer;
    packer.nodeCount = 7;
    packer.flowers.assign(8, {});
    // Boundary flowers listed in correct counterclockwise order (matching
    // vertex 1's own 2,3,4,5,6,7 CCW ordering -- see the fixed-up
    // tests/data/small_disc.p, which hand-authored this exact same hex
    // flower and originally had these six rows backward; this test happened
    // to still pass either way since the hexagon's symmetry makes the
    // angle-sum/radii checks below orientation-agnostic, which is exactly
    // why the mistake went unnoticed).
    packer.flowers[1] = {2, 3, 4, 5, 6, 7, 2};
    packer.flowers[2] = {3, 1, 7};
    packer.flowers[3] = {4, 1, 2};
    packer.flowers[4] = {5, 1, 3};
    packer.flowers[5] = {6, 1, 4};
    packer.flowers[6] = {7, 1, 5};
    packer.flowers[7] = {2, 1, 6};
    packer.alpha = 1;
    packer.gamma = 0;
    packer.hes = Geometry::Euclidean;

    gopack::Index nc = packer.complexCount();
    checkTrue(nc == 7, "complexCount should return nodeCount (7)");
    checkTrue(packer.intCount == 1, "expected exactly 1 interior vertex");
    checkTrue(packer.bdryCount == 6, "expected exactly 6 boundary vertices");
    checkTrue(packer.orphanCount == 0, "expected no orphan vertices");

    packer.centers.assign(8, Complex(0.0, 0.0));
    packer.radii.assign(8, 0.5);
    packer.localcenters = packer.centers;
    packer.localradii = packer.radii;

    int md = packer.setMode(1);
    checkTrue(md == 1, "setMode(1) should succeed with mode==1");

    packer.indxMatrices();
    gopack::RiffleResult result = packer.riffle(100);
    checkTrue(result.cycles >= 0, "riffle should not report an error");

    auto diffs = packer.angsumErrors();
    checkTrue(diffs.size() == 1, "angsumErrors should have one entry (the interior vertex)");
    checkTrue(std::abs(diffs[0]) < 1e-4,
              "interior angle sum should converge close to 2*pi (diff ~ 0)");

    // All six boundary radii should be (very nearly) equal by symmetry.
    double r0 = packer.radii[2];
    for (gopack::Index v = 3; v <= 7; ++v) {
        double diff = std::abs(packer.radii[v] - r0);
        if (diff > 1e-3) {
            std::fprintf(stderr,
                         "FAIL: boundary radius at vertex %d (%.6f) differs from vertex 2 "
                         "(%.6f) by more than tolerance\n",
                         v, packer.radii[v], r0);
            return 1;
        }
    }

    std::fprintf(stdout, "hex_flower: OK (%d riffle passes, angsum diff=%.3e, radii[1]=%.6f, "
                          "boundary radius=%.6f)\n",
                 result.cycles, diffs[0], packer.radii[1], r0);
    return 0;
}

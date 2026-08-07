// Regression test for Packer::loadComplex() -- the new in-memory entry
// point meant for library/JNI callers that already hold a triangulation
// in memory (e.g. CirclePack's own per-vertex flower data) and want to hand
// it to GOPack directly, without round-tripping it through the *.p text
// format the way readpack() requires.
//
// Strategy: build the classical "hex flower" complex (see
// test_hex_flower.cpp) two different ways -- once via loadComplex(), once
// by hand-populating Packer fields directly (mirroring what readpack()
// itself effectively does) -- and check that they converge to the same
// packing. This is the sharpest test of loadComplex()/ingestFlowers()/
// finalizeComplex(): if they don't reproduce exactly what readpack() would
// have produced for the same combinatorics, this is where it would show up.
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

std::vector<std::vector<gopack::Index>> hexFlowerFlowers() {
    std::vector<std::vector<gopack::Index>> flowers(8);
    flowers[1] = {2, 3, 4, 5, 6, 7, 2};
    flowers[2] = {7, 1, 3};
    flowers[3] = {2, 1, 4};
    flowers[4] = {3, 1, 5};
    flowers[5] = {4, 1, 6};
    flowers[6] = {5, 1, 7};
    flowers[7] = {6, 1, 2};
    return flowers;
}

} // namespace

int main() {
    using gopack::Complex;
    using gopack::Geometry;
    using gopack::Index;
    using gopack::Packer;

    // ---- Packer A: via loadComplex(), alpha auto-detected (pass 0) ----
    Packer a;
    Index nc = a.loadComplex(7, hexFlowerFlowers(), Geometry::Euclidean,
                              /*alphaIn=*/0, /*gammaIn=*/0);
    checkTrue(nc == 7, "loadComplex should return nodeCount (7)");
    checkTrue(a.intCount == 1, "expected exactly 1 interior vertex");
    checkTrue(a.bdryCount == 6, "expected exactly 6 boundary vertices");
    checkTrue(a.orphanCount == 0, "expected no orphan vertices");
    checkTrue(a.alpha == 1, "alpha should auto-resolve to the interior vertex (1)");
    checkTrue(a.radii.size() == 8, "radii should be sized nodeCount+1 with defaults applied");
    checkTrue(std::abs(a.radii[2] - 0.5) < 1e-12, "default radius should be 0.5, matching readpack()");
    checkTrue(!a.vAims.empty(), "vAims should get readpack()'s usual default (2*pi/-1)");
    checkTrue(std::abs(a.vAims[1] - gopack::kTwoPi) < 1e-12, "interior vAims default should be 2*pi");
    checkTrue(a.vAims[2] < 0, "boundary vAims default should be -1 (free)");

    int md = a.setMode(1);
    checkTrue(md == 1, "setMode(1) should succeed with mode==1");
    gopack::RiffleResult resultA = a.riffle(100);
    checkTrue(resultA.cycles >= 0, "riffle should not report an error (loadComplex path)");

    // ---- Packer B: hand-populated fields directly (mirrors what
    // readpack() itself would leave in place, and what test_hex_flower.cpp
    // already exercises independently) ----
    Packer b;
    b.nodeCount = 7;
    b.flowers = hexFlowerFlowers();
    b.alpha = 1;
    b.gamma = 0;
    b.hes = Geometry::Euclidean;
    b.complexCount();
    b.centers.assign(8, Complex(0.0, 0.0));
    b.radii.assign(8, 0.5);
    b.localcenters = b.centers;
    b.localradii = b.radii;
    b.setMode(1);
    b.indxMatrices();
    gopack::RiffleResult resultB = b.riffle(100);
    checkTrue(resultB.cycles >= 0, "riffle should not report an error (hand-built path)");

    // The two should converge to (very nearly) the same packing -- same
    // combinatorics, same starting radii/centers, same solver.
    checkTrue(resultA.cycles == resultB.cycles,
              "loadComplex() and hand-built packers should take the same number of "
              "riffle passes to converge");
    for (Index v = 1; v <= 7; ++v) {
        double dr = std::abs(a.radii[v] - b.radii[v]);
        if (dr > 1e-9) {
            std::fprintf(stderr,
                         "FAIL: radius at vertex %d differs between loadComplex (%.10f) and "
                         "hand-built (%.10f) packers by %.3e\n",
                         v, a.radii[v], b.radii[v], dr);
            return 1;
        }
    }

    // Sanity on the actual packing geometry, same checks as test_hex_flower.cpp.
    auto diffs = a.angsumErrors();
    checkTrue(diffs.size() == 1, "angsumErrors should have one entry (the interior vertex)");
    checkTrue(std::abs(diffs[0]) < 1e-4,
              "interior angle sum should converge close to 2*pi (diff ~ 0)");

    // ---- Packer C: loadComplex() with an explicitly-invalid alpha (a
    // boundary vertex) -- must be rejected and re-resolved, same as
    // readpack() does for a bad ALPHA/GAMMA: line. ----
    Packer c;
    Index ncC = c.loadComplex(7, hexFlowerFlowers(), Geometry::Euclidean,
                               /*alphaIn=*/3 /* a boundary vertex -- invalid */, 0);
    checkTrue(ncC == 7, "loadComplex should still succeed with an invalid explicit alpha");
    checkTrue(c.alpha == 1, "invalid explicit alpha should be re-resolved to the interior vertex");

    // ---- Packer D: loadComplex() with explicit initial radii/vAims,
    // exercising the optional-argument paths that mirror readpack()'s
    // optional RADII:/ANGLE_AIMS: sections. ----
    Packer d;
    std::vector<gopack::Scalar> initRadii(8, 0.25);
    std::vector<gopack::Scalar> customAims(8, gopack::kTwoPi);
    customAims[2] = -1.0; // still free -- just confirming an override is honored at all
    Index ncD = d.loadComplex(7, hexFlowerFlowers(), Geometry::Euclidean, 0, 0, &initRadii,
                               /*initCenters=*/nullptr, /*vlistIn=*/nullptr, &customAims);
    checkTrue(ncD == 7, "loadComplex with optional radii/vAims should succeed");
    checkTrue(std::abs(d.radii[1] - 0.25) < 1e-12,
              "explicit initRadii should override the 0.5 default");
    checkTrue(std::abs(d.vAims[2] - (-1.0)) < 1e-12, "explicit vAimsIn should be honored as-is");

    // ---- Packer E: malformed input (size mismatch) must fail cleanly
    // rather than reading out of bounds. ----
    Packer e;
    std::vector<std::vector<Index>> tooShort(5); // should be size 8 for nodeCount=7
    Index ncE = e.loadComplex(7, tooShort, Geometry::Euclidean);
    checkTrue(ncE == 0, "loadComplex should reject a flowersIn size mismatch instead of crashing");

    std::fprintf(stdout,
                 "loadComplex: OK (%d riffle passes, angsum diff=%.3e, radii[1]=%.6f, "
                 "loadComplex/hand-built agree to <1e-9)\n",
                 resultA.cycles, diffs[0], a.radii[1]);
    return 0;
}

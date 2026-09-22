// Independent sanity/regression test for orthopack mode (setMode mode 4,
// layoutBdry / setOrthoCenters) -- new in this port, no MATLAB counterpart;
// see Packer::setOrthoCenters()'s doc comment in Packer.h.
//
// Two things are checked:
//
//   1. setOrthoCenters() in isolation, on a hand-built boundary-only
//      "packer" with deliberately UNEQUAL boundary radii (bypassing
//      complexCount()/riffle() entirely -- the method only reads
//      nodeCount/bdryCount/bdryList/localradii and writes localcenters).
//      This directly checks the boundary-layout math itself: every boundary
//      circle should come out orthogonal to the unit circle
//      (|center|^2 == 1 + radius^2) and tangent to both cclw neighbors
//      (|center_j - center_{j+1}| == radius_j + radius_{j+1}, including the
//      wraparound pair), to near machine precision -- both identities are
//      exact by construction, independent of any riffle convergence.
//
//   2. The same hex-flower fixture test_polygonal.cpp/test_hex_flower.cpp
//      use (1 interior vertex, 6 boundary vertices), riffled in orthopack
//      mode via the normal setMode(4)/riffle() path, checking that mode/hes
//      end up as expected, every radius/center stays finite and positive,
//      and -- after one extra "settle" pass (continueRiffle(0), which reruns
//      layoutBdry()+layoutCenters() without an intervening setEffective(),
//      so the boundary layout is recomputed from whatever radii riffle()
//      already converged to) -- the same orthogonality/tangency identities
//      as above hold for the boundary vertices to tight tolerance. This
//      exercises the actual setMode(4)/layoutBdry() dispatch and its
//      interaction with the shared continueRiffle/layoutCenters/setEffective
//      iteration (Steps B/C), not just the boundary-layout formula alone.
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

bool approxEqual(double a, double b, double tol) { return std::abs(a - b) <= tol; }

} // namespace

int main() {
    using gopack::Complex;
    using gopack::Geometry;
    using gopack::Index;
    using gopack::Packer;

    // ---- Part 1: setOrthoCenters() in isolation, unequal radii ----
    {
        Packer packer;
        const Index n = 5;
        packer.nodeCount = n;
        packer.bdryCount = n;
        packer.bdryList = {1, 2, 3, 4, 5, 1}; // closed cclw list, size bdryCount+1
        packer.localradii.assign(static_cast<size_t>(n) + 1, 0.0);
        const double r[6] = {0.0, 0.3, 0.7, 0.2, 1.1, 0.5}; // index 0 unused
        for (Index v = 1; v <= n; ++v) packer.localradii[v] = r[v];
        packer.localcenters.assign(static_cast<size_t>(n) + 1, Complex(0.0, 0.0));

        packer.setOrthoCenters();

        constexpr double kTol = 1e-9;
        for (Index j = 1; j <= n; ++j) {
            Index v = packer.bdryList[static_cast<size_t>(j) - 1];
            checkTrue(std::isfinite(packer.localcenters[v].real()) &&
                          std::isfinite(packer.localcenters[v].imag()),
                      "isolated: boundary center should be finite");
            checkTrue(packer.localradii[v] > 0.0, "isolated: boundary radius should stay positive");

            double normSq = std::norm(packer.localcenters[v]);
            double want = 1.0 + packer.localradii[v] * packer.localradii[v];
            checkTrue(approxEqual(normSq, want, kTol),
                      "isolated: boundary circle should be orthogonal to the unit circle");
        }
        for (Index j = 1; j <= n; ++j) {
            Index vj = packer.bdryList[static_cast<size_t>(j) - 1];
            Index vk = packer.bdryList[static_cast<size_t>(j % n)]; // wraps j==n -> index 0 -> vertex 1
            double dist = std::abs(packer.localcenters[vj] - packer.localcenters[vk]);
            double want = packer.localradii[vj] + packer.localradii[vk];
            checkTrue(approxEqual(dist, want, kTol),
                      "isolated: consecutive boundary circles should be exactly tangent");
        }
        std::fprintf(stdout, "test_orthopack: part 1 (isolated setOrthoCenters, n=%d) OK\n", n);
    }

    // ---- Part 2: full setMode(4)/riffle() integration, hex-flower fixture ----
    {
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
        // Deliberately start as Hyperbolic: setMode(4, ...) should force this
        // back to Euclidean, exactly like setMode(2, ...) already does.
        packer.hes = Geometry::Hyperbolic;

        Index nc = packer.complexCount();
        checkTrue(nc == 7, "complexCount should return nodeCount (7)");
        checkTrue(packer.bdryCount == 6, "expected exactly 6 boundary vertices");

        packer.centers.assign(8, Complex(0.0, 0.0));
        packer.radii.assign(8, 0.5);
        packer.localcenters = packer.centers;
        packer.localradii = packer.radii;

        int md = packer.setMode(4);
        checkTrue(md == 4, "setMode(4) should succeed with mode==4");
        checkTrue(packer.mode == 4, "packer.mode should be 4 after setMode(4)");
        checkTrue(packer.hes == Geometry::Euclidean,
                  "setMode(4) should force hes back to Euclidean");

        packer.indxMatrices();
        gopack::RiffleResult result = packer.riffle(300);
        checkTrue(result.cycles >= 0, "riffle should not report an error");

        for (Index v = 1; v <= packer.nodeCount; ++v) {
            checkTrue(std::isfinite(packer.radii[v]), "radius should be finite");
            checkTrue(std::isfinite(packer.centers[v].real()) &&
                          std::isfinite(packer.centers[v].imag()),
                      "center should be finite");
            checkTrue(packer.radii[v] > 0.0, "radius should be strictly positive");
        }

        // One more "settle" pass (no setEffective()) so the boundary layout
        // is exactly consistent with whatever radii riffle() converged to,
        // then pull that snapshot into centers/radii via reapResults() --
        // see the file header comment for why this is needed to get an
        // exact (not just converged-to-tolerance) geometric snapshot.
        packer.continueRiffle(0);
        packer.reapResults();

        constexpr double kTol = 1e-6;
        for (Index k = 1; k <= packer.bdryCount; ++k) {
            Index v = packer.bdryList[static_cast<size_t>(k) - 1];
            double normSq = std::norm(packer.centers[v]);
            double want = 1.0 + packer.radii[v] * packer.radii[v];
            checkTrue(approxEqual(normSq, want, kTol),
                      "riffled: boundary circle should be orthogonal to the unit circle");
        }
        for (Index k = 1; k <= packer.bdryCount; ++k) {
            Index vj = packer.bdryList[static_cast<size_t>(k) - 1];
            Index vk = packer.bdryList[static_cast<size_t>(k % packer.bdryCount)];
            double dist = std::abs(packer.centers[vj] - packer.centers[vk]);
            double want = packer.radii[vj] + packer.radii[vk];
            checkTrue(approxEqual(dist, want, kTol),
                      "riffled: consecutive boundary circles should be tangent");
        }

        std::fprintf(stdout,
                     "test_orthopack: part 2 (setMode(4)/riffle integration, %d passes) OK\n",
                     result.cycles);
    }

    return 0;
}

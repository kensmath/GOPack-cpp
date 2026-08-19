// Regression test for the spherical centroid-normalization fix: previously,
// Packer::radii/Packer::centers were only recentered inside writepack()'s
// Spherical branch, so a caller reading them straight after riffle() (like
// the JNI bridge's computeMaximalPackingFromComplex, which never calls
// writepack()) got an un-normalized, possibly lopsided packing. The fix
// moved the normalization into reapResults(), which riffle() always calls.
//
// Reads a real 1000-vertex spherical triangulation (data/sphtest1000.p,
// GEOMETRY: spherical) via readpack() -- the same file-loading path a real
// user goes through -- rather than a small hand-built complex, so this
// isn't relying on a synthetic combinatorial structure whose correctness
// hasn't been independently verified.
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

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <path-to-sphtest1000.p>\n", argv[0]);
        return 2;
    }

    using gopack::Geometry;
    using gopack::Index;
    using gopack::Packer;

    Packer packer;
    Index nc = packer.readpack(argv[1]);
    checkTrue(nc > 0, "readpack should succeed");
    checkTrue(nc == 1000, "expected 1000 vertices");
    checkTrue(packer.hes == Geometry::Spherical,
              "input file's GEOMETRY: spherical should be read as Geometry::Spherical");

    int md = packer.setMode(1);
    checkTrue(md == 1, "setMode(1) should succeed");

    gopack::RiffleResult result = packer.riffle(50);
    checkTrue(result.cycles >= 0, "riffle should not report an error");

    for (Index v = 1; v <= packer.nodeCount; ++v) {
        checkTrue(std::isfinite(packer.radii[v]) && packer.radii[v] > 0.0,
                  "all radii should be finite and positive after riffle");
        checkTrue(std::isfinite(packer.centers[v].real()) && std::isfinite(packer.centers[v].imag()),
                  "all centers should be finite after riffle");
    }

    // The actual fix under test: Packer::centers/radii -- read directly,
    // the same way the JNI bridge does, with NO call to writepack() -- must
    // already be centroid-normalized on the sphere by the time riffle()
    // returns.
    std::vector<gopack::Complex> T = packer.loadTangency(packer.centers, packer.radii);
    gopack::geom::CentroidResult c = gopack::geom::centroid(T, {1.0, 0.0, 0.0});
    std::fprintf(stdout, "post-riffle centroid: (%.6f, %.6f, %.6f), normSq=%.6f\n", c.x, c.y, c.z,
                 c.normSq);
    checkTrue(c.normSq <= 0.001,
              "Packer::centers/radii should already be centroid-normalized after riffle(), "
              "without needing writepack()");

    std::fprintf(stdout,
                 "sphere_normalize: OK (%d riffle passes, %d vertices, centroid normSq=%.3e)\n",
                 result.cycles, packer.nodeCount, c.normSq);
    return 0;
}

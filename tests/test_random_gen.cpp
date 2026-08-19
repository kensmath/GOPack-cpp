// Regression test for gopack::geom::delaunayPlane()/convexHull3()
// (core/src/RandomGen.cpp) -- the bindings to the vendored Triangle/Qhull
// libraries (see third_party/*/README.md) that a future
// randTriangulation()-equivalent will build on. Only built when
// GOPACK_BUILD_RANDOM_GEN is on (see tests/CMakeLists.txt).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <vector>

#include "gopack/RandomGen.h"

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
    using gopack::Index;
    using gopack::Scalar;

    // ==== delaunayPlane: unconstrained -- 4 corners of a square + 1 center
    // point should triangulate into a 4-triangle fan, using every point. ====
    {
        std::vector<Complex> pts = {
            Complex(0, 0), Complex(2, 0), Complex(2, 2), Complex(0, 2), Complex(1, 1),
        };
        auto tris = gopack::geom::delaunayPlane(pts);
        checkTrue(tris.size() == 4, "unconstrained: 4 corners + center should give 4 triangles");
        std::set<Index> usedVerts;
        for (auto& t : tris) {
            for (Index v : t) {
                checkTrue(v >= 0 && v < static_cast<Index>(pts.size()), "vertex index in range");
                usedVerts.insert(v);
            }
            checkTrue(t[0] != t[1] && t[1] != t[2] && t[0] != t[2], "triangle has 3 distinct verts");
        }
        checkTrue(usedVerts.size() == 5, "all 5 points should appear in the triangulation");
        std::fprintf(stdout, "delaunayPlane (unconstrained square+center): OK, %zu triangles\n",
                     tris.size());
    }

    // ==== delaunayPlane: constrained -- an L-shaped (non-convex) boundary
    // plus one interior point. The constrained triangulation's total area
    // must match the L-shape (3.0), not its bounding square (4.0) --
    // confirming Triangle actually honored the boundary constraint instead
    // of filling in the missing corner. ====
    {
        std::vector<Complex> bdry = {
            Complex(0, 0), Complex(2, 0), Complex(2, 1), Complex(1, 1), Complex(1, 2), Complex(0, 2),
        };
        std::vector<Complex> pts = bdry;
        pts.push_back(Complex(0.5, 0.5)); // interior point, safely inside the L
        std::vector<std::array<Index, 2>> segs;
        for (size_t i = 0; i < bdry.size(); ++i) {
            segs.push_back({static_cast<Index>(i), static_cast<Index>((i + 1) % bdry.size())});
        }
        auto tris = gopack::geom::delaunayPlane(pts, segs);
        checkTrue(!tris.empty(), "constrained L-shape: should produce triangles");
        Scalar totalArea = 0;
        for (auto& t : tris) {
            Complex a = pts[static_cast<size_t>(t[0])];
            Complex b = pts[static_cast<size_t>(t[1])];
            Complex c = pts[static_cast<size_t>(t[2])];
            totalArea += 0.5 * std::abs((b.real() - a.real()) * (c.imag() - a.imag()) -
                                         (c.real() - a.real()) * (b.imag() - a.imag()));
        }
        checkTrue(std::abs(totalArea - 3.0) < 1e-9,
                  "constrained triangulation area should match the L-shape (3.0), confirming "
                  "the boundary constraint was honored, not its bounding square (4.0)");
        std::fprintf(stdout, "delaunayPlane (constrained L-shape): OK, %zu triangles, area=%.6f\n",
                     tris.size(), totalArea);
    }

    // ==== delaunayPlane: degenerate input ====
    checkTrue(gopack::geom::delaunayPlane({Complex(0, 0), Complex(1, 0)}).empty(),
              "fewer than 3 points should return empty");

    // ==== convexHull3: a regular octahedron -- exactly 8 triangular
    // faces, using all 6 vertices. ====
    {
        std::vector<std::array<Scalar, 3>> pts = {
            {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1},
        };
        auto tris = gopack::geom::convexHull3(pts);
        checkTrue(tris.size() == 8, "octahedron convex hull should have exactly 8 faces");
        std::set<Index> usedVerts;
        for (auto& t : tris) {
            for (Index v : t) {
                checkTrue(v >= 0 && v < static_cast<Index>(pts.size()), "vertex index in range");
                usedVerts.insert(v);
            }
            checkTrue(t[0] != t[1] && t[1] != t[2] && t[0] != t[2], "triangle has 3 distinct verts");
        }
        checkTrue(usedVerts.size() == 6, "all 6 octahedron vertices should appear in the hull");
        std::fprintf(stdout, "convexHull3 (octahedron): OK, %zu triangular faces\n", tris.size());
    }

    // ==== convexHull3: degenerate input ====
    checkTrue(gopack::geom::convexHull3({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}).empty(),
              "fewer than 4 points should return empty");

    // ==== convexHull3: a larger, deterministic "spiral" point set on the
    // unit sphere (no <random> needed for a reproducible test), checking
    // Euler's formula (F = 2V-4 for a closed all-triangle hull) holds. ====
    {
        const int n = 200;
        std::vector<std::array<Scalar, 3>> pts;
        pts.reserve(static_cast<size_t>(n));
        Scalar phi = 0.0;
        for (int i = 0; i < n; ++i) {
            Scalar h = -1.0 + 2.0 * i / (n - 1);
            Scalar theta = std::acos(std::max(-1.0, std::min(1.0, h)));
            if (i != 0 && i != n - 1) {
                phi += 3.6 / std::sqrt(static_cast<Scalar>(n) * (1.0 - h * h));
            }
            Scalar s = std::sin(theta);
            pts.push_back({s * std::cos(phi), s * std::sin(phi), h});
        }
        auto tris = gopack::geom::convexHull3(pts);
        checkTrue(!tris.empty(), "spiral sphere: should produce a hull");
        std::set<Index> usedVerts;
        for (auto& t : tris) {
            for (Index v : t) usedVerts.insert(v);
        }
        Index V = static_cast<Index>(usedVerts.size());
        Index F = static_cast<Index>(tris.size());
        checkTrue(F == 2 * V - 4, "closed all-triangle hull should satisfy Euler's formula F = 2V-4");
        std::fprintf(stdout, "convexHull3 (spiral sphere, n=%d): OK, V=%d F=%d\n", n, V, F);
    }

    std::fprintf(stdout, "random_gen: ALL OK\n");
    return 0;
}

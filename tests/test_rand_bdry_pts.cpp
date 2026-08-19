// Regression test for gopack::geom::randBdryPts() (rand_bdry_pts.m) -- pick
// M points uniformly at random, by arc length, along a closed polygonal
// path. Used by randTriangulation.m's plane-region case (still not
// ported -- see Packer.h) to seed a constrained Delaunay triangulation's
// boundary points; tested standalone here since it has no dependency on
// Packer at all.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "gopack/Geometry.h"

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

    // Closed unit square boundary (already closed -- last point == first).
    std::vector<Complex> square = {
        Complex(1, 1), Complex(-1, 1), Complex(-1, -1), Complex(1, -1), Complex(1, 1),
    };
    const Index M = 37;
    std::vector<Complex> pts = gopack::geom::randBdryPts(square, M);
    checkTrue(static_cast<Index>(pts.size()) == M, "should return exactly M points");

    for (const Complex& z : pts) {
        // distance from z to the nearest of the 4 edges should be ~0
        Scalar best = 1e18;
        for (size_t k = 0; k + 1 < square.size(); ++k) {
            Complex a = square[k], b = square[k + 1];
            Complex ab = b - a;
            Scalar t = ((z - a).real() * ab.real() + (z - a).imag() * ab.imag()) /
                       (ab.real() * ab.real() + ab.imag() * ab.imag());
            t = std::max(0.0, std::min(1.0, t));
            Complex proj = a + t * ab;
            best = std::min(best, std::abs(z - proj));
        }
        checkTrue(best < 1e-9, "every returned point should lie on the square boundary");
    }

    // Auto-close: an explicitly-unclosed path should get closed
    // automatically (matching rand_bdry_pts.m's own check) and still return
    // M points, all on one of its (now 5, including the closing edge) sides.
    // NOTE: the source's auto-close condition is "diffX>0.001 AND
    // diffY>0.001" (not OR) -- ported literally (see randBdryPts's header
    // doc comment) -- so this needs a shape whose first/last points differ
    // in BOTH coordinates, unlike an axis-aligned rectangle's corners.
    std::vector<Complex> pentagonOpen = {
        Complex(0, 0), Complex(2, 0), Complex(3, 2), Complex(1, 3), Complex(-1, 1),
    };
    std::vector<Complex> pts2 = gopack::geom::randBdryPts(pentagonOpen, M);
    checkTrue(static_cast<Index>(pts2.size()) == M, "auto-closed path should also return M points");
    std::vector<Complex> pentagonClosed = pentagonOpen;
    pentagonClosed.push_back(pentagonOpen.front());
    for (const Complex& z : pts2) {
        Scalar best = 1e18;
        for (size_t k = 0; k + 1 < pentagonClosed.size(); ++k) {
            Complex a = pentagonClosed[k], b = pentagonClosed[k + 1];
            Complex ab = b - a;
            Scalar t = ((z - a).real() * ab.real() + (z - a).imag() * ab.imag()) /
                       (ab.real() * ab.real() + ab.imag() * ab.imag());
            t = std::max(0.0, std::min(1.0, t));
            Complex proj = a + t * ab;
            best = std::min(best, std::abs(z - proj));
        }
        checkTrue(best < 1e-9,
                  "auto-closed path: every returned point should lie on one of its 5 sides "
                  "(including the auto-added closing edge)");
    }

    // degenerate-input handling: too few graph points, or M<3.
    std::vector<Complex> tooShort = {Complex(0, 0), Complex(1, 0)};
    checkTrue(gopack::geom::randBdryPts(tooShort, 10).empty(),
              "graph with fewer than 3 points should return empty");
    checkTrue(gopack::geom::randBdryPts(square, 2).empty(), "M<3 should return empty");

    std::fprintf(stdout, "rand_bdry_pts: OK (%d points, all on the square boundary)\n",
                 static_cast<int>(pts.size()));
    return 0;
}

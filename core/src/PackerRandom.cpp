// GOPack C++ -- "geometrically random" triangulation generators.
//
// Direct translations of randomDisc.m, randomSphere.m, randomRectangle.m,
// randomSquare.m, and randomTri.m, built on top of
// gopack::geom::randTriangulationSphere()/randTriangulationPlane()
// (RandomGen.h/.cpp -- themselves built on the vendored Triangle/Qhull
// libraries, see third_party/*/README.md) and the already-ported
// Packer::parseTriangles()/pruneComplex()/indxMatrices()/setMode(). Only
// compiled when GOPACK_HAVE_RANDOM_GEN is defined -- see
// core/CMakeLists.txt's GOPACK_BUILD_RANDOM_GEN option (default ON) and the
// matching #ifdef around these methods' declarations in Packer.h.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

#include "gopack/Packer.h"
#include "gopack/RandomGen.h"

namespace gopack {

namespace {

// Shift a 0-indexed triangle list (gopack::geom::delaunayPlane()/
// convexHull3()'s own convention -- see RandomGen.h's module comment) to
// the 1-indexed convention Packer::parseTriangles() expects.
std::vector<std::array<Index, 3>> shiftTriTo1Indexed(
    const std::vector<std::array<Index, 3>>& tri0) {
    std::vector<std::array<Index, 3>> t(tri0.size());
    for (size_t i = 0; i < tri0.size(); ++i) {
        t[i] = {tri0[i][0] + 1, tri0[i][1] + 1, tri0[i][2] + 1};
    }
    return t;
}

// Build a parseTriangles()-compatible 'cents' vector (1-indexed, size
// Z.size()+1, index 0 an unused placeholder) from a 0-indexed point array.
std::vector<Complex> buildCents1Indexed(const std::vector<Complex>& Z) {
    std::vector<Complex> cents(Z.size() + 1);
    cents[0] = Complex(0.0, 0.0);
    for (size_t i = 0; i < Z.size(); ++i) {
        cents[i + 1] = Z[i];
    }
    return cents;
}

} // namespace

Packer Packer::randomDisc(Index N) {
    Packer gop;

    // N total points: for intN interior, want roughly pi*sqrt(intN) bdry.
    Index intN =
        1 + static_cast<Index>(std::floor(
                std::pow(std::sqrt(kPi * kPi + 4.0 * N) / 2.0 - kPi / 2.0, 2)));
    Index bdryN = N - intN;

    static thread_local std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<Scalar> unit(0.0, 1.0);

    // bdryN points distributed uniformly on the unit circle.
    std::vector<Complex> bdry(static_cast<size_t>(bdryN));
    for (Index j = 0; j < bdryN; ++j) {
        Scalar ang = 2.0 * kPi * unit(rng);
        bdry[static_cast<size_t>(j)] = Complex(std::cos(ang), std::sin(ang));
    }

    // intN points chosen interior to the unit disc: point 0 near the
    // origin (to act as alpha), the rest via rejection sampling in
    // [-1,1]x[-1,1]. See randTriangulationPlane()'s doc comment ("Bug fix
    // vs. randTriangulation.m") for why the safety counter below is
    // incremented on every attempt, unlike randomDisc.m's own identical
    // (broken) loop.
    std::vector<Complex> interior(static_cast<size_t>(intN));
    Scalar invN = 1.0 / static_cast<Scalar>(N);
    interior[0] = Complex(invN * (-1.0 + 2.0 * unit(rng)), invN * (-1.0 + 2.0 * unit(rng)));
    Index hits = 1;
    Index count = 0;
    const Index maxCount = 100 * intN;
    while (hits < intN && count < maxCount) {
        Scalar x = -1.0 + 2.0 * unit(rng);
        Scalar y = -1.0 + 2.0 * unit(rng);
        if (x * x + y * y < 1.0) {
            interior[static_cast<size_t>(hits)] = Complex(x, y);
            ++hits;
        }
        ++count;
    }
    if (count >= maxCount) {
        std::fprintf(stderr, "overran safety check in randomDisc\n");
    }

    std::vector<Complex> Z = interior;
    Z.insert(Z.end(), bdry.begin(), bdry.end());

    auto tri0 = geom::delaunayPlane(Z); // unconstrained: convex hull of the unit disc's points
    auto tList = shiftTriTo1Indexed(tri0);
    auto cents = buildCents1Indexed(Z);

    gop.alpha = 1;
    gop.parseTriangles(tList, &cents);
    gop.indxMatrices();
    gop.hes = Geometry::Hyperbolic;
    gop.mode = 1;
    return gop;
}

Packer Packer::randomSphere(Index intN) {
    Packer gop;
    if (intN < 4) {
        std::fprintf(stderr, "randomSphere should have at least 4 points\n");
        return gop;
    }

    // Note: Z has (theta,phi) values, but not useful in GOpacker -- omitted
    // from the parseTriangles() call below, matching randomSphere.m's
    // 'gop.parse_triangles(tri);' (no 'cents' argument).
    auto res = geom::randTriangulationSphere(intN);
    auto tList = shiftTriTo1Indexed(res.tri);

    gop.parseTriangles(tList);
    gop.indxMatrices();

    std::fprintf(stdout, "GOpacker started with random spherical triangulation, %d vertices\n",
                 gop.nodeCount);
    return gop;
}

Packer Packer::randomRectangle(Index intN, Scalar aspect, Index bdryN) {
    Packer gop;
    if (intN < 1) {
        std::fprintf(stderr, "random rectangle should have at least 1 interior point\n");
        return gop;
    }

    Scalar asp = std::abs(aspect);
    Index bN = static_cast<Index>(std::floor(4.0 * std::sqrt(static_cast<Scalar>(intN)) * asp));
    if (bdryN > 0) bN = bdryN; // explicit override, matching randomRectangle.m's nargin==3 branch
    if (bN < 4) bN = 4;

    // Rectangle boundary path [-asp,asp]x[-1,1], counterclockwise, closed
    // (matches randomRectangle.m's gX=[asp;-asp;-asp;asp;asp],
    // gY=[1;1;-1;-1;1]).
    std::vector<Complex> graph = {Complex(asp, 1), Complex(-asp, 1), Complex(-asp, -1),
                                   Complex(asp, -1), Complex(asp, 1)};

    static thread_local std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<Scalar> unit(0.0, 1.0);
    Scalar iN = 1.0 / static_cast<Scalar>(intN);
    Complex center(iN * (-1.0 + 2.0 * unit(rng)), iN * (-1.0 + 2.0 * unit(rng)));

    auto res = geom::randTriangulationPlane(intN, bN, graph, &center);
    if (res.alpha >= 0) {
        gop.alpha = res.alpha + 1; // 0-indexed -> 1-indexed
    }

    auto tList = shiftTriTo1Indexed(res.tri);
    auto cents = buildCents1Indexed(res.Z);
    gop.parseTriangles(tList, &cents);

    // The triangulation may have orphan vertices cut off from the interior
    // (see pruneComplex()'s own doc comment).
    gop.pruneComplex();
    gop.indxMatrices();

    // Find the boundary vertices closest to the rectangle's 4 actual
    // corners -- 'gamma' is the upper-right corner.
    Scalar ur = asp, ul = asp, ll = asp, lr = asp;
    Index urV = gop.bdryList[0], ulV = gop.bdryList[0], llV = gop.bdryList[0],
          lrV = gop.bdryList[0];
    for (size_t i = 0; i + 1 < gop.bdryList.size(); ++i) {
        Index v = gop.bdryList[i];
        Complex z = gop.centers[static_cast<size_t>(v)];
        Scalar distur = std::abs(z - Complex(asp, 1));
        Scalar distul = std::abs(z - Complex(-asp, 1));
        Scalar distll = std::abs(z - Complex(-asp, -1));
        Scalar distlr = std::abs(z - Complex(asp, -1));
        if (distur < ur) { urV = v; ur = distur; }
        if (distul < ul) { ulV = v; ul = distul; }
        if (distll < ll) { llV = v; ll = distll; }
        if (distlr < lr) { lrV = v; lr = distlr; }
    }

    // list of corners in cclw order
    gop.vlist = {urV, ulV, llV, lrV};
    gop.gamma = urV;

    // set 'alpha' if not already set (from 'center' above)
    if (res.alpha < 0) {
        Index a = gop.intVerts[0];
        Scalar mindist = std::abs(gop.centers[static_cast<size_t>(a)]);
        for (Index j = 0; j < gop.intCount; ++j) {
            Index v = gop.intVerts[static_cast<size_t>(j)];
            Scalar dist = std::abs(gop.centers[static_cast<size_t>(v)]);
            if (dist < mindist) {
                mindist = dist;
                a = v;
            }
        }
        gop.alpha = a;
    }

    gop.setMode(2, gop.vlist);
    std::fprintf(stdout, "GOpacker started with random rectangle, aspect %f, %d vertices\n", asp,
                 gop.nodeCount);
    return gop;
}

Packer Packer::randomSquare(Index N) {
    // with intN interior, want about 4*sqrt(intN) on bdry
    Index intN = 1 + static_cast<Index>(std::floor(std::pow(-2.0 + std::sqrt(4.0 + N), 2)));
    Index bdryN = N - intN;
    return randomRectangle(intN, 1.0, bdryN);
}

Packer Packer::randomTri(Index intN) {
    Packer gop = randomSphere(intN);
    gop.hes = Geometry::Spherical; // redundant but faithful (see randomTri.m)
    return gop;
}

Packer Packer::randomTri(Index intN, Index bdryN, const std::vector<Complex>& graph,
                          const Complex* cent) {
    Packer gop;

    auto res = geom::randTriangulationPlane(intN, bdryN, graph, cent);
    if (res.tri.empty() || res.Z.empty()) {
        std::fprintf(stderr, "Error in generating the random triangulation.\n");
        return gop;
    }
    // Bug fix vs. randomTri.m: when 'cent' wasn't placed, leave gop.alpha at
    // its Packer()-default 0 instead of the source's -1 -- see this
    // overload's doc comment in Packer.h for the full rationale.
    if (res.alpha >= 0) {
        gop.alpha = res.alpha + 1; // 0-indexed -> 1-indexed
    }

    // Deliberate simplification vs. randomTri.m: the source re-derives its
    // own contiguous 1..nodeCount renumbering here (a ~30-line manual
    // 'indxhits' walk over tList) before ever calling parse_triangles() --
    // but parse_triangles.m (Packer::parseTriangles()) already performs the
    // identical renumbering internally (see its own doc comment in
    // Packer.h), so that block would just be duplicated, dead work here;
    // tList/cents are passed to parseTriangles() directly instead.
    auto tList = shiftTriTo1Indexed(res.tri);
    auto cents = buildCents1Indexed(res.Z);
    gop.parseTriangles(tList, &cents);
    gop.indxMatrices();
    gop.hes = Geometry::Euclidean;
    gop.mode = 1;
    return gop;
}

} // namespace gopack

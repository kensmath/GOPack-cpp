#include "gopack/RandomGen.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

// ---------------------------------------------------------------------
// Triangle (Shewchuk) -- see third_party/triangle/README.md.
//
// triangle.h expects the *includer* to #define REAL and VOID before
// including it (triangle.c does this for its own compilation; a separate
// translation unit like this one that also includes triangle.h needs to
// do the same, matching what triangle.c was actually compiled with -- see
// third_party/triangle/CMakeLists.txt's TRILIBRARY/ANSI_DECLARATORS/
// NO_TIMER defines). REAL must be 'double' here because triangle.c was
// built WITHOUT -DSINGLE (so its own REAL is double, and the
// 'triangulateio' struct's REAL* fields must agree on element size between
// the two translation units). ANSI_DECLARATORS must also be (re)defined
// here so triangle.h declares triangulate()/trifree() with real parameter
// types rather than empty-parens K&R style -- in C++, an empty-parens
// declaration means "takes no arguments", which would make calling
// triangulate() with its actual 4 arguments a hard compile error.
// triangle.h itself has no extern "C" guard (it predates that convention),
// so this translation unit supplies one.
#define REAL double
#define VOID void
#ifndef ANSI_DECLARATORS
#define ANSI_DECLARATORS
#endif
extern "C" {
#include "triangle.h"
}

// ---------------------------------------------------------------------
// Qhull -- see third_party/qhull/README.md. qhull_ra.h (the "all headers"
// convenience header the upstream examples use) already wraps its
// declarations in '#ifdef __cplusplus extern "C" { ... }', so no manual
// wrapping is needed here the way Triangle's header needs above.
#include "qhull_ra.h"

#include "gopack/Geometry.h" // projVecToS(), randBdryPts() -- used below by
                              // randTriangulationSphere()/randTriangulationPlane()

namespace gopack::geom {

std::vector<std::array<Index, 3>> delaunayPlane(
    const std::vector<Complex>& points, const std::vector<std::array<Index, 2>>& segments) {
    if (points.size() < 3) {
        return {};
    }

    std::vector<REAL> pointlist(points.size() * 2);
    for (size_t i = 0; i < points.size(); ++i) {
        pointlist[2 * i] = points[i].real();
        pointlist[2 * i + 1] = points[i].imag();
    }

    std::vector<int> segmentlist;
    if (!segments.empty()) {
        segmentlist.reserve(segments.size() * 2);
        for (const auto& s : segments) {
            segmentlist.push_back(static_cast<int>(s[0]));
            segmentlist.push_back(static_cast<int>(s[1]));
        }
    }

    struct triangulateio in {};
    struct triangulateio out {};
    in.pointlist = pointlist.data();
    in.numberofpoints = static_cast<int>(points.size());
    if (!segments.empty()) {
        in.segmentlist = segmentlist.data();
        in.numberofsegments = static_cast<int>(segments.size());
    }

    // z: number output from 0 (matches this function's 0-indexed
    //    convention). Q: quiet (suppress Triangle's routine informational
    //    output; error/warning messages still print). N: don't bother
    //    writing out->pointlist -- we already have the points and only
    //    want the triangle connectivity back. p: use the PSLG (segmentlist)
    //    for a constrained triangulation honoring the boundary, only when
    //    segments were actually given -- matching
    //    randTriangulation.m's delaunayTriangulation(X,Y,C) (segments
    //    given) vs. delaunayTriangulation(X,Y) (no constraints) split.
    std::string switches = "zQN";
    if (!segments.empty()) {
        switches += "p";
    }
    // triangulate()'s first parameter is 'char *', not 'const char *' (a
    // pre-C89-const-correctness API) -- it never writes through the
    // pointer, but a real (non-literal) mutable buffer is the standard,
    // safe way to satisfy the signature.
    std::vector<char> switchesBuf(switches.begin(), switches.end());
    switchesBuf.push_back('\0');

    triangulate(switchesBuf.data(), &in, &out, nullptr);

    std::vector<std::array<Index, 3>> result;
    result.reserve(static_cast<size_t>(out.numberoftriangles));
    for (int t = 0; t < out.numberoftriangles; ++t) {
        result.push_back({static_cast<Index>(out.trianglelist[3 * t]),
                           static_cast<Index>(out.trianglelist[3 * t + 1]),
                           static_cast<Index>(out.trianglelist[3 * t + 2])});
    }

    // Triangle allocates every 'out' array itself (all of 'out's pointers
    // were left null/zeroed on input above); free them with trifree(),
    // Triangle's own allocator-paired deallocator -- not a bare free(), in
    // case a caller ever swaps in a custom trimalloc()/trifree() (see
    // triangle.h). Safe to call on the null pointers among these (the
    // fields we never asked Triangle to populate, e.g. edges/neighbors --
    // trifree() just wraps free(), and free(NULL) is a defined no-op).
    trifree(out.pointlist);
    trifree(out.pointmarkerlist);
    trifree(out.trianglelist);
    trifree(out.segmentlist);
    trifree(out.segmentmarkerlist);

    return result;
}

std::vector<std::array<Index, 3>> convexHull3(const std::vector<std::array<Scalar, 3>>& points) {
    constexpr int kDim = 3;
    if (points.size() < 4) {
        return {};
    }

    std::vector<coordT> coords(points.size() * kDim);
    for (size_t i = 0; i < points.size(); ++i) {
        coords[kDim * i + 0] = points[i][0];
        coords[kDim * i + 1] = points[i][1];
        coords[kDim * i + 2] = points[i][2];
    }

    qhT qh_qh;
    qhT* qh = &qh_qh;
    QHULL_LIB_CHECK
    qh_zero(qh, stderr);

    // "qhull " prefix is required by qh_new_qhull() itself (a vestige of
    // argv[0]-style command parsing, not a real qhull option); 'Qt' forces
    // triangulated output (every facet a simplex/triangle) even for
    // degenerate input with more than 3 coplanar points on a facet, which
    // matters here since we always expect a pure triangle list back.
    // outfile=nullptr is explicitly documented as valid by qh_new_qhull()
    // ("outfile may be null") and skips Qhull's own output-formatting pass
    // entirely, since we only want the facet list, not printed text.
    char flags[] = "qhull Qt";
    int exitcode = qh_new_qhull(qh, kDim, static_cast<int>(points.size()), coords.data(),
                                 /*ismalloc=*/False, flags, /*outfile=*/nullptr, stderr);

    std::vector<std::array<Index, 3>> result;
    if (!exitcode) {
        facetT* facet;
        vertexT* vertex;
        vertexT** vertexp;
        FORALLfacets {
            std::array<Index, 3> tri{};
            int n = 0;
            FOREACHvertex_(facet->vertices) {
                if (n < 3) {
                    tri[static_cast<size_t>(n)] = static_cast<Index>(qh_pointid(qh, vertex->point));
                }
                ++n;
            }
            if (n == 3) {
                result.push_back(tri);
            } else {
                // Shouldn't happen with 'Qt' (triangulated output) -- skip
                // rather than silently emit a bogus/truncated triangle from
                // a non-triangular facet.
                std::fprintf(stderr,
                             "gopack::geom::convexHull3: skipping a non-triangular facet "
                             "(%d vertices) despite 'Qt'\n", n);
            }
        }
    } else {
        std::fprintf(stderr, "gopack::geom::convexHull3: Qhull reported an error (exitcode %d)\n",
                     exitcode);
    }

    // Standard reentrant-Qhull cleanup (matches upstream's own
    // user_eg_r.c): free long memory, then short memory + the allocator
    // itself. Not wrapped in '#ifdef qh_NOmem' since this project never
    // defines that symbol.
    qh_freeqhull(qh, !qh_ALL);
    int curlong = 0, totlong = 0;
    qh_memfreeshort(qh, &curlong, &totlong);
    if (curlong || totlong) {
        std::fprintf(stderr,
                     "gopack::geom::convexHull3: qhull internal warning: did not free %d bytes "
                     "of long memory (%d pieces)\n", totlong, curlong);
    }

    return result;
}

bool pointInPolygon(Complex pt, const std::vector<Complex>& poly) {
    const size_t n = poly.size();
    if (n < 3) return false;
    const Scalar x = pt.real();
    const Scalar y = pt.imag();
    bool inside = false;
    // Standard even-odd (crossing-number) ray-casting test: cast a ray in
    // the +x direction from (x,y) and count how many polygon edges it
    // crosses; odd == inside. 'j' trails 'i' by one, wrapping from the last
    // vertex back to the first, so 'poly' need not be explicitly closed.
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        const Scalar xi = poly[i].real();
        const Scalar yi = poly[i].imag();
        const Scalar xj = poly[j].real();
        const Scalar yj = poly[j].imag();
        const bool straddles = (yi > y) != (yj > y);
        if (straddles && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
            inside = !inside;
        }
    }
    return inside;
}

RandTriangulationSphereResult randTriangulationSphere(Index intN) {
    RandTriangulationSphereResult result;
    if (intN < 1) return result;

    static thread_local std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<Scalar> unit(0.0, 1.0);

    std::vector<std::array<Scalar, 3>> sphPts(static_cast<size_t>(intN));
    result.Z.resize(static_cast<size_t>(intN));
    for (Index i = 0; i < intN; ++i) {
        Scalar a = unit(rng) * 2.0 * kPi;
        Scalar z = unit(rng) * 2.0 - 1.0;
        // std::max guards against sqrt() of a tiny negative value from
        // floating-point rounding when z rounds to just outside [-1,1] --
        // a defensive addition beyond randTriangulation.m's literal
        // 'sqrt(1-z(i)^2)', which has no such guard.
        Scalar r = std::sqrt(std::max(Scalar(0), 1.0 - z * z));
        std::array<Scalar, 3> pt = {r * std::cos(a), r * std::sin(a), z};
        sphPts[static_cast<size_t>(i)] = pt;
        result.Z[static_cast<size_t>(i)] = projVecToS(pt);
    }

    result.tri = convexHull3(sphPts);
    return result;
}

RandTriangulationPlaneResult randTriangulationPlane(Index intN, Index bdryN,
                                                     const std::vector<Complex>& graph,
                                                     const Complex* cent) {
    RandTriangulationPlaneResult result;

    const Index graphNum = static_cast<Index>(graph.size());
    if (graphNum < 3 || bdryN < 3 || intN < 1) {
        std::fprintf(stderr, "usage: randTriangulationPlane(intN, bdryN, graph, [cent])\n");
        return result;
    }

    // bdryN random points on the boundary path, in arc-length order.
    std::vector<Complex> bdryPts = randBdryPts(graph, bdryN);
    if (bdryPts.empty()) {
        return result; // randBdryPts() already printed a diagnostic
    }

    // bounding rectangle of the (original, un-closed) boundary path.
    Scalar minx = graph.front().real(), maxx = minx;
    Scalar miny = graph.front().imag(), maxy = miny;
    for (const Complex& g : graph) {
        minx = std::min(minx, g.real());
        maxx = std::max(maxx, g.real());
        miny = std::min(miny, g.imag());
        maxy = std::max(maxy, g.imag());
    }
    const Scalar xrange = maxx - minx;
    const Scalar yrange = maxy - miny;

    static thread_local std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<Scalar> unit(0.0, 1.0);

    // is 'cent' usable as the alpha vertex?
    bool placeCent = false;
    if (cent && pointInPolygon(*cent, graph)) {
        placeCent = true;
    }

    std::vector<Complex> interior(static_cast<size_t>(intN));
    Index hits = 0; // next interior[] slot to fill (0-indexed)
    if (placeCent) {
        interior[0] = *cent;
        hits = 1;
        result.alpha = 0;
    }

    // Rejection-sample the remaining interior points inside 'graph'. See
    // this function's header doc comment ("Bug fix vs. randTriangulation.m")
    // for why 'count' is incremented on every attempt here, unlike the
    // MATLAB source.
    Index count = 0;
    const Index maxCount = 100 * intN;
    while (hits < intN && count < maxCount) {
        Scalar x = minx + unit(rng) * xrange;
        Scalar y = miny + unit(rng) * yrange;
        if (pointInPolygon(Complex(x, y), graph)) {
            interior[static_cast<size_t>(hits)] = Complex(x, y);
            ++hits;
        }
        ++count;
    }
    if (count >= maxCount) {
        std::fprintf(stderr, "overran safety check in randTriangulationPlane\n");
    }

    // Z = [interior points; boundary points], matching
    // randTriangulation.m's X=[intX;bdryX] ordering.
    result.Z = interior;
    result.Z.insert(result.Z.end(), bdryPts.begin(), bdryPts.end());

    // Closed loop of boundary-edge constraints, 0-indexed into result.Z
    // (the boundary points start at offset intN).
    std::vector<std::array<Index, 2>> segments(static_cast<size_t>(bdryN));
    for (Index k = 0; k < bdryN; ++k) {
        segments[static_cast<size_t>(k)] = {intN + k, intN + ((k + 1) % bdryN)};
    }

    // Constrained Delaunay triangulation of the region bounded by
    // 'segments'. Unlike randTriangulation.m, no post-hoc trimming of
    // out-of-region triangles is needed here -- see this function's header
    // doc comment ("Deliberate simplification") for why.
    result.tri = delaunayPlane(result.Z, segments);

    return result;
}

} // namespace gopack::geom

#pragma once
//
// GOPack C++ -- bindings to vendored third-party computational-geometry
// libraries (Triangle, Qhull), plus the randTriangulation.m-equivalent
// orchestration built directly on top of them (point generation, boundary
// sampling, constrained/unconstrained triangulation selection). See
// Packer.h's randomDisc()/randomSphere()/randomRectangle()/randomSquare()/
// randomTri() for the next layer up -- the randomDisc.m/randomSphere.m/
// randomRectangle.m/randomSquare.m/randomTri.m equivalents that turn this
// file's output into a ready-to-riffle Packer. See
// third_party/triangle/README.md and third_party/qhull/README.md for
// exactly what's vendored, what license terms apply, and why this project
// vendors (commits) the source rather than fetching it over the network the
// way Eigen's CMakeLists.txt does -- short version: Ken wants a
// self-contained, offline-buildable package, not a build-time dependency on
// two old, low-traffic upstream projects staying reachable.
//
// Unlike the rest of gopack (1-indexed, size N+1 with index 0 unused -- see
// Packer.h's module comment), delaunayPlane()/convexHull3()/
// randTriangulationSphere()/randTriangulationPlane() below all return plain
// 0-indexed triangle lists and point arrays, matching Triangle's/Qhull's own
// native conventions and the fact that these operate on a bare point array
// with no natural "vertex 0" placeholder to reserve. Packer.h's
// randomDisc()/randomSphere()/randomRectangle()/randomTri() are responsible
// for the +1 shift into parseTriangles()'s 1-indexed tList convention.
#include <array>
#include <vector>

#include "gopack/Types.h"

namespace gopack::geom {

// Constrained Delaunay triangulation of a 2D point set, via Triangle
// (Shewchuk) -- what randTriangulation.m's plane-region case uses MATLAB's
// delaunayTriangulation(X,Y,C) for.
//
// points: the point set to triangulate (randTriangulation.m's own
// convention is interior points followed by boundary points, matching its
// X=[intX;bdryX], but nothing here actually requires that particular
// ordering -- it's simply what a caller building this array will naturally
// have on hand).
// segments: pairs of 0-indexed indices into 'points', each pair a
// constraint edge the triangulation must include as-is (typically a closed
// boundary polygon's edges in order, as randTriangulation.m's own C matrix
// encodes); pass empty for a plain (unconstrained) Delaunay triangulation
// of the convex hull of 'points'.
//
// Returns one {v0,v1,v2} entry per triangle, 0-indexed into 'points', in
// Triangle's standard (counterclockwise) winding. Returns empty if
// points.size() < 3, or if Triangle reports an error (e.g.
// self-intersecting segments) -- Triangle's own diagnostic is printed to
// stderr in that case (Triangle's informational chatter is otherwise
// suppressed via its 'Q' quiet switch; error/warning output is not).
std::vector<std::array<Index, 3>> delaunayPlane(
    const std::vector<Complex>& points,
    const std::vector<std::array<Index, 2>>& segments = {});

// 3D convex hull of a point set, via Qhull -- what randTriangulation.m's
// sphere case uses MATLAB's convhulln() for. For a point set already lying
// on the unit sphere (the only way this project uses it -- see
// s_pt_to_vec()/randTriangulation.m's own sphere-projection loop), the
// convex hull coincides exactly with the sphere's own Delaunay
// triangulation; no paraboloid-lifting step (the kind Qhull's own
// qh_setdelaunay() performs for planar Delaunay-via-convex-hull) is needed
// or performed here -- this is a bare "3D convex hull" primitive.
//
// Returns one {v0,v1,v2} entry per triangle, 0-indexed into 'points'.
// Returns empty if points.size() < 4 (the minimum for a nondegenerate 3D
// hull), or if Qhull reports an error (e.g. all points coplanar) -- Qhull's
// own diagnostic is printed to stderr in that case.
std::vector<std::array<Index, 3>> convexHull3(
    const std::vector<std::array<Scalar, 3>>& points);

// inpolygon.m (MATLAB builtin) equivalent -- strict interior test via the
// standard even-odd (crossing-number) rule. This project only ever uses it
// to decide whether a randomly sampled point should be kept (rejection
// sampling) or whether a given 'cent' falls inside a region, never to
// distinguish an exact-boundary hit from outside, so the crossing-number
// rule (which does not special-case boundary points, unlike MATLAB's
// three-way in/on/out result) is sufficient here. 'poly' need not be
// explicitly closed (first==last); the edge from the last vertex back to the
// first is always implicitly included.
bool pointInPolygon(Complex pt, const std::vector<Complex>& poly);

// randTriangulation.m's sphere case (nargin==1): choose 'intN' random points
// via a Poisson Point Process on the unit sphere S^2 and return their
// Delaunay triangulation, computed as the 3D convex hull (which coincides
// with the sphere's Delaunay triangulation for points already lying on the
// sphere -- see convexHull3()'s own doc comment).
struct RandTriangulationSphereResult {
    // 0-indexed, one {v0,v1,v2} per triangle, straight from convexHull3().
    std::vector<std::array<Index, 3>> tri;
    // One entry per input point (0-indexed, parallel to the points convexHull3()
    // was given), in (theta,phi) polar form (projVecToS()'s convention) --
    // matches randTriangulation.m's sphere-branch 'Z' output. Note
    // randTriangulation.m's own 'alpha' output is never actually set on this
    // code path (the function returns before reaching the plane-region code
    // that sets it) -- its real callers only ever destructure the first two
    // outputs here (e.g. randomSphere.m's '[tri,~]=randTriangulation(intN);'),
    // so there is nothing to return here either.
    std::vector<Complex> Z;
};
RandTriangulationSphereResult randTriangulationSphere(Index intN);

// randTriangulation.m's plane-region case (nargin>1): a constrained Delaunay
// triangulation of a plane region bounded by the closed polygonal path
// 'graph', using 'bdryN' points sampled uniformly (by arc length) along the
// boundary (via randBdryPts(), see Geometry.h) and 'intN' points chosen
// interior to the region via a Poisson Point Process (rejection sampling in
// the bounding rectangle). If 'cent' is given and falls inside 'graph', it
// is placed as interior point 0 (so the caller can use it as a packing's
// designated 'alpha' vertex); otherwise all interior points are randomly
// sampled.
//
// Deliberate simplification vs. randTriangulation.m: the MATLAB source
// triangulates with plain (unconstrained-region) delaunayTriangulation(X,Y,C)
// and then runs ~90 lines of its own manual post-hoc trimming (checking
// corner convexity, then discarding boundary-only faces and faces with an
// out-of-region interior vertex) to remove triangles that spill outside the
// region bounded by the *sampled* boundary polygon (bdryX/bdryY -- see the
// paragraph below on what "the region" means here). MATLAB's
// delaunayTriangulation with constraint edges only enforces those edges as
// *present*; it does not by itself exclude the exterior of a non-convex
// region. Triangle (see delaunayPlane() above) does this natively: passing a
// closed segment loop without the '-c' ("enclose convex hull") switch makes
// Triangle itself discard any triangle outside the segment-bounded region,
// concavities included (see third_party/triangle/triangle.c's own switch
// documentation under '-c', and the constrained/L-shape case in
// tests/test_random_gen.cpp, which exercises exactly this directly against
// delaunayPlane(), using the true polygon's own vertices as the boundary,
// and checks the resulting triangulation's total area against the true
// L-shape area). So none of randTriangulation.m's manual corner-convexity/
// face-discarding trimming logic is ported here -- it would be dead weight
// duplicating what delaunayPlane()'s own constrained mode already
// guarantees, GIVEN a closed segment loop.
//
// IMPORTANT caveat this function inherits unchanged from randTriangulation.m
// (not something delaunayPlane()'s native carving fixes or could fix): "the
// region bounded by the boundary segments" means the region bounded by the
// bdryN *sampled* points chord-connected in arc-length order -- a polygon
// inscribed in 'graph', not 'graph' itself. Both here and in the MATLAB
// source, the boundary segments (and every inpolygon()/pointInPolygon()
// check) are built from randBdryPts()'s bdryN-point output, never from
// 'graph's own vertices directly, so the triangulated region only converges
// to 'graph's true shape as bdryN grows -- for a sparse bdryN relative to a
// concave shape's feature size (e.g. a whole edge of 'graph' going
// unsampled purely by chance), the inscribed chord polygon can visibly
// shortcut a concave notch, and any interior point that ends up outside that
// chord polygon (despite being inside the true 'graph') is correctly
// excluded by delaunayPlane() from the returned triangulation entirely (it
// simply never appears in any output triangle, so it is silently absent
// from Z's "used" set -- see Packer::randomRectangle()'s pruneComplex() call
// in PackerRandom.cpp, which exists precisely to clean up after this for a
// Packer built from this function's output; Packer::randomTri()'s generic
// plane overload does not call pruneComplex(), matching randomTri.m, so a
// caller there may see slightly fewer than intN+bdryN vertices survive into
// the final Packer for a sparse bdryN/highly concave 'graph').
//
// Bug fix vs. randTriangulation.m (and randomDisc.m, which has the identical
// bug in its own separate interior-sampling loop): the source's rejection-
// sampling safety counter ('count', capped at 100*intN) is only ever
// incremented on a *successful* hit (inside the 'if in==1' branch in
// randTriangulation.m; there is no 'count=count+1' anywhere in the loop
// body itself), not on every attempt -- so the safety cap can never actually
// trigger; a region with near-zero rejection-sampling acceptance probability
// (e.g. a very thin sliver polygon) would spin the real MATLAB loop
// indefinitely. This reads as an unambiguous oversight (the surrounding
// comment and the "overran safety check" message both describe a per-attempt
// cap), not intentional behavior, and a hang is a much worse failure mode in
// a batch C++ CLI than in an interactive MATLAB session -- so here 'count'
// is incremented on every attempt, successful or not, making the cap
// actually effective.
struct RandTriangulationPlaneResult {
    // 0-indexed, one {v0,v1,v2} per triangle, straight from delaunayPlane().
    std::vector<std::array<Index, 3>> tri;
    // One entry per generated point (0-indexed): the 'intN' interior points
    // first (index 0 is 'cent', if it was placed), then the 'bdryN' boundary
    // points, matching randTriangulation.m's own X=[intX;bdryX] ordering.
    std::vector<Complex> Z;
    // 0-indexed index into Z of the placed 'cent' point, or -1 if 'cent' was
    // not given (nullptr) or fell outside 'graph' -- the 0-indexed analog of
    // randTriangulation.m's 1-indexed 'alpha' output (1 if placed, else -1).
    Index alpha = -1;
};
RandTriangulationPlaneResult randTriangulationPlane(Index intN, Index bdryN,
                                                     const std::vector<Complex>& graph,
                                                     const Complex* cent = nullptr);

} // namespace gopack::geom

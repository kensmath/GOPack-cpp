# GOPack C++

A C++ port of [GOPack](https://github.com/kensmath/GOPack) (Collins, Orick,
Stephenson, 2017) -- MATLAB software for computing maximal circle packings
with very large numbers of circles -- structured to build as a native
library for Java (via JNI) and as standalone Windows/macOS executables.

**License:** GPL-3.0, same as the original GOPack (this is a derivative
work; see LICENSE.md).

## Status: max-pack and polygonal/rectangle modes ported

This is a **direct, line-by-line port**, not a reimplementation from the
paper's description. The translation was done by reading every `.m` file in
`code/@GOPacker/` and `code/` and translating each one into C++ against the
same field names, method names, and loop structure, specifically so the two
can be diffed and checked against each other rather than trusted on faith.
See the module-level comment in `core/include/gopack/Packer.h` for the
indexing convention (1-indexed, matching MATLAB, to keep every loop bound
identical to the source).

**Fully ported and covered by tests** (`gopack::Packer`):

- `readpack` -- reads the `*.p` FLOWERS format (docs/GO_Formats.txt)
- `loadComplex` -- **not part of the original MATLAB source**; a new entry
  point that ingests an already-known combinatorial complex directly from
  memory (nodeCount + per-vertex flower lists + geometry, matching the *.p
  FLOWERS format's own fields but as arrays instead of parsed text), for
  callers that already hold the triangulation in memory instead of a file
  (e.g. the JNI bridge's `computeMaximalPackingFromComplex`, see "Java
  usage" below). Implemented as the same post-parse logic readpack() itself
  uses (factored out so the two can't silently diverge), so a complex loaded
  either way behaves identically from that point on.
- `complex_count`, `indxMatrices`, `FarVert`
- `parse_triangles` -- ingest a bare `Nx3` triangle-list (as produced by a
  Delaunay/convex-hull step, e.g. the `randTriangulation` family) directly,
  building flowers/orientation/alpha from scratch, without going through the
  `*.p` text format. `pruneComplex` -- remove orphan vertices (cut off from
  the interior component) after such a triangulation, needed for bounded
  (non-convex) regions. `rand_bdry_pts` -- pick `M` points uniformly at
  random by arc length along a closed polygonal path (`gopack::geom::randBdryPts`
  in `Geometry.h`; the only piece of the `rand*`/`Triangulation` family with
  no external Delaunay/convex-hull dependency). See "Notes on
  `parse_triangles`/`pruneComplex`/`rand_bdry_pts`" below for the handful of
  deliberate deviations from the literal source.
- **Mode 1 (maximal packing)** of the disc/plane, or the sphere if the
  complex has no boundary: `setMode` mode 1, `layoutBdry` / `setHoroCenters`
- **Mode 2 (polygonal / rectangle packing)**: `setMode` mode 2 (corner
  selection -- explicit, from `vlist`, or pseudo-random), `layoutBdry`
  dispatching to `setPolyCenters` (general n-gon) / `setRectCenters`
  (automatic special case for 4 right-angle corners), and `getAspect`. The
  shared `continueRiffle`/`layoutCenters`/`setEffective` iteration below
  needed no changes to support this -- `setEffective` already branches on
  the sign of `vAims`, which `setMode` sets appropriately for either mode.
- the core iteration (used by both modes): `continueRiffle` / `layoutCenters`
  (the sparse linear solve) / `setEffective` / `updateVdata` / `visualErrors`
- `riffle`, `reapResults`, `angsumErrors`, `packStatus`
- `writepack` / `writeEucl`, including hyperbolic and spherical output
  conversion (`e_to_h_data`, `h_to_e_data`, `e_to_s_data`, `s_to_e_data`,
  `sph_tangent`, `affineNormalizer`, `Centroid`, `loadTangency`)
- `cosAngle`, `cosCorner`
- `gopack::geom::delaunayPlane` / `convexHull3` (`core/src/RandomGen.cpp`) --
  **not part of the original MATLAB source**; thin C++ bindings to the
  vendored Triangle (plane constrained Delaunay) and Qhull (3D convex hull)
  libraries (see "Vendored third-party libraries" below), doing for C++ what
  MATLAB's built-in `delaunayTriangulation`/`convhulln` do for
  `randTriangulation.m`. Themselves fully tested (`tests/test_random_gen.cpp`).
- Random triangulation *generation* -- `randTriangulation` (as
  `gopack::geom::randTriangulationSphere`/`randTriangulationPlane`,
  `core/src/RandomGen.cpp`, built directly on `delaunayPlane`/`convexHull3`
  above) and `randomDisc`/`randomSphere`/`randomRectangle`/`randomSquare`/
  `randomTri` (as `Packer::randomDisc`/`randomSphere`/`randomRectangle`/
  `randomSquare`/`randomTri`, `core/src/PackerRandom.cpp`, built on
  `parse_triangles`/`pruneComplex` above). Covered by
  `tests/test_random_packers.cpp`, and exposed from the standalone CLI (not
  just JNI, per Ken's request that the standalone executable have this
  capability too, independent of CirclePack) via `--random-disc`/
  `--random-sphere`/`--random-square`/`--random-rectangle`/`--random-tri` --
  see "CLI usage" below. `randomTri`'s arbitrary-polygon-region overload is
  also exposed via JNI as `computeRandomTri` -- see "Java usage" below. See
  "Notes on the `randTriangulation` family" below for the
  deliberate deviations from the literal source (two bug fixes, and one
  simplification enabled by Triangle's native boundary-carving).

**Understood but NOT yet ported**:

- The `OFF` file format fallback in `readpack` (calling it prints a
  diagnostic and returns 0, rather than guessing at behavior)
- The plotting method (`show`, which has no headless equivalent anyway)

If your workflow needs any of the above, say so -- the source for all of
them has already been read and understood (see the corresponding `.m`
files), so porting them is a bounded follow-up, not new research.

### Notes on the mode-2 (polygonal) port

- `setMode`'s C++ signature can't distinguish "corner list omitted" from
  "corner list explicitly passed as empty" the way MATLAB's `nargin` can; an
  empty `crns` here always falls through to the `vlist`-based or
  pseudo-random corner selection, which is the more useful default for a
  caller (see the doc comment on `Packer::setMode`).
- `setMode.m`'s "put `cornangs` in counterclockwise order" block recomputes
  a local `corners` variable using an unsorted `bdryIndx` (its
  `sort(bdryIndx);` call never captures a return value) and that
  recomputed value is never read again -- it's dead code in the source, and
  was omitted rather than ported literally, since porting genuinely dead
  code (with no observable effect either way) would just be noise.
- `setRectCenters.m` has a comment describing corners at
  "lowerleft (-aspect,-1), upper right (aspect,1)" that doesn't match the
  actual corner coordinates used in the same function (real part in
  `{+1,-1}`, imaginary part in `{+aspect,-aspect}`); ported literally as
  written, not as commented, with a note in the code.
- `setMode(2, ...)` always resets `hes` to Euclidean, even if the packing
  was originally read as hyperbolic or spherical -- `setMode.m` doesn't do
  this. GOPack always computes internally in euclidean coordinates
  regardless of `hes` (see `Geometry.h`); `hes` only controls whether
  `readpack()`/`writepack()` convert to/from hyperbolic or spherical circle
  data at the file boundary. A polygon/rectangle boundary is inherently a
  euclidean shape, so leaving `hes` at its original (hyperbolic/spherical)
  value would make `writepack()` apply a conversion to already-euclidean
  polygon output that was never intended for it.

### Notes on `parse_triangles`/`pruneComplex`/`rand_bdry_pts`

- `parse_triangles.m`'s own alpha-auto-selection fallback (used when no
  alpha was already set before the call) indexes its `utilFlag` boundary-seed
  array using the *new* (renumbered) vertex numbers as if they were the
  *old* ones -- only actually correct when the input's vertex numbering is
  already contiguous from 1 (`nodeCount == top`), which is true for this
  function's real callers (fresh Delaunay/hull output always numbers its
  points contiguously). Ported literally, since "fix" here would mean
  guessing at an unintended generalization rather than correcting a clear
  bug -- if a future caller passes genuinely sparse vertex numbers, alpha
  selection may pick a less-than-ideal (but still valid; `complex_count()`
  independently re-validates alpha regardless) vertex.
- `pruneComplex.m` reads `obj.vlist(j)` *after* `obj.vlist` was just cleared
  two lines earlier in the source -- real undefined behavior (out-of-bounds
  access) in a literal C++ port, not just a MATLAB quirk to preserve. Fixed
  to read from a saved local copy of the pre-clear `vlist` instead.
- `pruneComplex.m` has `origCenters(nv)=obj.origCenters(nv)`, almost
  certainly a typo for `obj.origCenters(v)` -- it should read from the OLD
  numbering, exactly like the `origRadii` line right above it and every
  other line in that loop. Ported as the evident intent (a literal port
  would silently read whichever old-numbered vertex happens to share `nv`'s
  *new* index -- wrong data, not just a crash).
- `pruneComplex.m`'s local `v2indx`/`indx2v` variables are renamed
  `oldToNew`/`newToOld` in the port: an unqualified `v2indx` inside a
  `Packer` member function would otherwise mean the real `Packer::v2indx`
  *member*, which serves an unrelated purpose (sparse-matrix layout
  indexing, set up by `indxMatrices()`) -- MATLAB's `obj.v2indx` vs. a
  same-named local are different namespaces, but C++ has no such
  distinction.
- `rand_bdry_pts.m`'s auto-close check is `abs(diffX)>0.001 AND
  abs(diffY)>0.001` (not OR) -- so a path whose first and last points
  coincide in exactly one coordinate (e.g. two corners of an axis-aligned
  rectangle) is *not* auto-closed. Ported literally as a preserved quirk,
  not a typo fix, since it's a plausible (if debatable) deliberate choice
  rather than an unambiguous slip; see `tests/test_rand_bdry_pts.cpp` for a
  worked example of what does and doesn't trigger it.
- `randBdryPts` adds one small defensive guard beyond the literal source:
  the interpolation ratio is taken as 0 (rather than dividing by zero) if
  two consecutive path points coincide exactly, and the forward arc-length
  scan is bounds-checked against running past the last segment -- both
  latent (if practically unreachable) issues in the MATLAB original too.

### Notes on the `randTriangulation` family

- **Simplification enabled by Triangle's native boundary carving**:
  `randTriangulation.m`'s plane-region case runs ~90 lines of its own manual
  post-hoc trimming (corner-convexity checks, then discarding boundary-only
  and out-of-region faces) after `delaunayTriangulation(X,Y,C)`, because
  MATLAB's triangulator only enforces constraint edges as *present* -- it
  doesn't exclude a non-convex region's exterior on its own. Triangle does
  this natively (passing a closed segment loop without the `-c` switch makes
  it discard everything outside the segment-bounded region, concavities
  included), so none of that manual trimming logic is ported --
  `gopack::geom::randTriangulationPlane` (`RandomGen.cpp`) relies on
  `delaunayPlane`'s constrained mode directly. Verified against a non-convex
  L-shaped region using its own true vertices as the boundary (checked by
  total area, confirming the missing corner isn't filled in) in both
  `tests/test_random_gen.cpp` (the raw `delaunayPlane` binding) and
  ASan/UBSan-instrumented sandbox validation of the full
  `randTriangulationPlane` orchestration.
- **Bug fix**: `randTriangulation.m`'s (and `randomDisc.m`'s own separate
  copy of the same loop's) rejection-sampling safety counter is only
  incremented on a *successful* hit, never on a rejected attempt -- so a
  region with low acceptance probability (e.g. a thin sliver) can spin the
  real MATLAB loop indefinitely, since the cap it's compared against can
  never actually be reached. Fixed to increment on every attempt; verified
  in sandbox validation that a deliberately thin sliver polygon with a large
  point request returns promptly rather than hanging.
- **Bug fix**: `randomTri.m` sets its local `GOPacker`'s `alpha` to `-1`
  before calling `randTriangulation`, and only overwrites it if
  `randTriangulation` returns a positive `alpha` (i.e. a valid `cent` was
  placed) -- but `parse_triangles.m`'s own alpha-auto-selection only
  triggers on `alpha==0`, not on a negative value, so a plain call without
  `cent` (or with a `cent` outside the region) hands back a `Packer` with an
  unresolved, invalid `alpha=-1`. `randomRectangle.m`, which has the same
  "`cent` may land outside the region" situation, avoids this because it
  never sets its `GOPacker`'s `alpha` away from its constructor default (0)
  in that case. `Packer::randomTri`'s plane-region overload leaves `alpha`
  at its `Packer()`-default 0 instead, letting the existing auto-selection
  produce a valid vertex.
- **Caveat inherited unchanged from `randTriangulation.m` (not a port bug,
  and not something Triangle's native carving above fixes or could fix)**:
  the boundary segments `randTriangulationPlane` triangulates against --
  and the source's own `inpolygon` rejection-sampling checks -- are built
  from `bdryN` *sampled* points chord-connected in arc-length order, not
  from `graph`'s own vertices directly, in both the C++ port and the
  original MATLAB (`rand_bdry_pts`/`randBdryPts`'s output feeds both). The
  triangulated region only converges to `graph`'s true shape as `bdryN`
  grows; for a sparse `bdryN` relative to a concave shape's feature size, the
  inscribed chord polygon can visibly shortcut a notch, and any point that
  ends up outside that chord polygon (despite being inside the true `graph`)
  is correctly excluded from the triangulation -- discovered and quantified
  during this port's sandbox validation (see `RandomGen.h`'s doc comment on
  `randTriangulationPlane` for the full explanation and a worked example).
  This only matters for `Packer::randomTri`'s generic plane-region overload
  with a caller-supplied concave `graph` and a small `bdryN`; `randomDisc`/
  `randomSphere`/`randomRectangle`/`randomSquare` (the CLI-exposed
  generators) never hit it, since a circle and a rectangle have no concave
  features to shortcut. `Packer::randomRectangle` also calls `pruneComplex()`
  (matching `randomRectangle.m`) specifically to clean up any point this
  does affect; `Packer::randomTri`'s plane overload does not (matching
  `randomTri.m`), so a caller who wants that cleanup there should call
  `pruneComplex()`/`indxMatrices()` on the result themselves.

### Notes on spherical packing normalization

- `reapResults()` now recenters `centers`/`radii` for Spherical packings
  (the same affine-normalization `affineNormalizer`/`centroid` in
  `Geometry.cpp` compute, moving the tangency-point centroid to the origin
  in 3D) -- **not part of `reapResults.m`**, which has no such step.
  `writepack()`'s Spherical branch used to be the only place this
  normalization happened, which meant any caller reading
  `Packer::centers`/`radii` directly after `riffle()` instead of going
  through `writepack()` -- notably the JNI bridge's
  `computeMaximalPackingFromComplex`, which never calls `writepack()` at
  all -- got an un-normalized, possibly lopsided packing. Every caller now
  gets a normalized packing for free as soon as `riffle()` returns.
- `writepack()` still does its own `affineNormalizer` call too, on purpose:
  it's the only safety net for a caller who calls `writepack()` on a
  Spherical packing that was just `readpack()`/`loadComplex()`-loaded but
  never `riffle()`d (so `reapResults()` never ran). `affineNormalizer` is
  idempotent on an already-centered input -- its very first centroid check
  is already within tolerance, so it returns the identity transform
  (`A=1, B=0`) immediately -- so this costs one cheap redundant pass over
  the tangency points on the normal riffle-then-write path, not a second
  real optimization.

## What has -- and hasn't -- been verified

- Ten regression tests (`tests/test_hex_flower.cpp`,
  `tests/test_readpack_roundtrip.cpp`, `tests/test_polygonal.cpp`,
  `tests/test_loadcomplex.cpp`, `tests/test_sphere_normalize.cpp`,
  `tests/test_parse_triangles.cpp`, `tests/test_prune_complex.cpp`,
  `tests/test_rand_bdry_pts.cpp`, `tests/test_random_gen.cpp`,
  `tests/test_random_packers.cpp`) exercise the pipeline on the classical
  "hex flower" complex (one interior vertex of degree 6 ringed by 6 boundary
  vertices) in both modes, checking angle-sum convergence and symmetric
  radii for mode 1, a read/riffle/write/re-read round trip, that a 4-corner
  rectangle layout (mode 2) stays finite, positive, and rectangle-shaped
  after riffling, that `loadComplex()` produces a packing that agrees with a
  hand-built (readpack()-equivalent) `Packer` to within 1e-9 -- including
  its alpha-resolution and optional-radii/vAims-override paths -- that
  `Packer::centers`/`radii` for a real 1000-vertex spherical triangulation
  (`tests/data/sphtest1000.p`) are already centroid-normalized right after
  `riffle()` with no call to `writepack()`, that `parse_triangles` correctly
  reconstructs a hex-flower fan and a closed tetrahedron (Spherical, with
  the 3-vertex pseudo-boundary anchor triangle `complex_count.m` always
  gives a boundary-less complex) from a bare triangle list, that
  `pruneComplex` removes an orphan vertex deliberately attached via a "flap"
  face and leaves a still-convergent packing behind, that `rand_bdry_pts`'s
  returned points all lie exactly on the source polygon's boundary
  (including its auto-close path), that `delaunayPlane`/`convexHull3` (the
  vendored Triangle/Qhull bindings) produce a correct constrained
  triangulation of a non-convex L-shaped region (checked by total area,
  confirming the missing corner wasn't filled in) and a correct 3D convex
  hull of points on a sphere (checked against Euler's formula, `F = 2V-4`,
  for a 200-point case), and that all five `Packer::random*` generators
  produce valid, riffle-convergent packings with the right `hes`/`mode`/
  `alpha`/corner metadata -- including a specific regression check that the
  `randomTri` `alpha=-1` bug fix above actually produces a valid vertex
  number, not the source's unresolved `-1`.
- **The vendored Triangle/Qhull integration, and the geometric-primitive
  layer built directly on it, was actually compiled and run in this
  sandbox** (unlike most of the rest of this port, which is syntax/logic-
  checked in isolation but only really build-tested on Ken's Windows
  machine): both libraries' real vendored source, plus
  `gopack::geom::delaunayPlane`/`convexHull3`/`randTriangulationSphere`/
  `randTriangulationPlane`/`pointInPolygon` (`RandomGen.cpp`) calling into
  them, were built with `-fsanitize=address,undefined` and exercised against
  real geometry (the L-shape and sphere cases above, plus an octahedron,
  unconstrained square+center-point case, and -- specifically for this
  session's new orchestration functions -- a dense-boundary-sampling
  convergence check, a `cent`-placement/rejection check, and a thin-sliver
  no-hang check for the rejection-sampling bug fix above) with no sanitizer
  reports. The `Packer`-level `randomDisc`/`randomSphere`/`randomRectangle`/
  `randomSquare`/`randomTri` generators (`PackerRandom.cpp`) were compile-
  checked for real against this sandbox's minimal Eigen stand-in (same as
  the rest of `Packer`, see below), and their Eigen-independent inner
  logic -- `parseTriangles()`/`pruneComplex()`, which is what these
  generators actually feed their triangulations through -- was additionally
  compiled, linked, and run for real against the new orchestration
  functions' actual output (not synthetic input), which is how the
  boundary-chord-approximation caveat documented above was discovered and
  quantified in the first place. This is a materially stronger verification
  bar than the "compiles cleanly against a minimal Eigen stand-in" level the
  rest of this port gets in-sandbox, since none of this layer depends on
  Eigen at all.
- **Mode 1 has been built and run successfully on real Windows hardware**,
  including on genuinely large inputs (up to ~500,000 vertices from
  `GOPack/data/lace500000_K.p`), producing packings that loaded correctly
  in CirclePack. Along the way this caught and fixed two real bugs that
  static review alone had missed: a Windows text-mode `ifstream`
  `tellg()`/`seekg()` desync in the `ALPHA/GAMMA:` parser (fixed by opening
  the file in binary mode -- see the comment on `Packer::readpack`), and an
  `Eigen::SparseLU`/`BiCGSTAB` solver-strategy default that assumed a
  symmetric system where GOPack's `transMatrix` isn't one.
- **Mode 2 (polygonal/rectangle) has been validated with a standalone
  ASan/UBSan-instrumented harness against the hex-flower fixture** (4-corner
  rectangle, 6-corner hexagon, 3-corner triangle, and auto-detected-corner
  cases all produced finite, geometrically sane output with no sanitizer
  reports), but **has not yet been exercised through a real build with a
  full riffle loop** the way mode 1 has -- that's the natural next step once
  this reaches your machine.
- **Numerical agreement with the original MATLAB has not been checked.**
  That requires running both implementations on the same input and diffing
  radii, which needs a MATLAB (or Octave) environment this session doesn't
  have. Recommended next step: run
  `GOPacker.readpack(...); GOPacker.setMode(1); GOPacker.riffle(50);` in
  MATLAB on one of the `data/*.p` files and compare the resulting
  `obj.radii` against this port's `gopack_cli` output on the same file.

## Vendored third-party libraries

GOPack-cpp vendors (commits directly into the repo) two small third-party C
libraries under `third_party/`, used by the random-triangulation
*generation* family (`gopack::geom::delaunayPlane`/`convexHull3`/
`randTriangulationSphere`/`randTriangulationPlane` and
`Packer::randomDisc`/`randomSphere`/`randomRectangle`/`randomSquare`/
`randomTri` -- see "Fully ported and covered by tests" above):

- **[Triangle](https://www.cs.cmu.edu/~quake/triangle.html)** (Jonathan
  Richard Shewchuk, v1.6, 2005) -- plane constrained Delaunay triangulation,
  what `randTriangulation.m`'s plane-region case uses MATLAB's
  `delaunayTriangulation(X,Y,C)` for.
- **[Qhull](https://github.com/qhull/qhull)** (v2020.2, `src/libqhull_r`
  only) -- 3D convex hull, what `randTriangulation.m`'s sphere case uses
  MATLAB's `convhulln` for.

See `third_party/triangle/README.md` and `third_party/qhull/README.md` for
exactly what's vendored (a deliberately minimal subset of each upstream
project -- not the CLI frontends, GUI tools, or docs), what license terms
apply, and how `core/src/RandomGen.cpp` calls into each.

**Why vendored (committed) rather than fetched at build time, unlike
Eigen above:** Eigen is a large, actively-maintained project with a fast git
host, so `FetchContent`-ing it on demand (falling back from a local
checkout) is a reasonable trade of repo size for staying current. Triangle
and Qhull are the opposite case -- Triangle hasn't been updated since 2005
and has no official git repository at all, and both are small enough
(under 2MB combined) that committing them costs nothing meaningful. Ken's
call: a working executable that builds identically offline, regardless of
whether `cs.cmu.edu` or `github.com` happen to be reachable (or whether
either project is still there at all) beats a network dependency on two
old, low-traffic upstream projects, and license terms aren't a concern for
this project's own (research/academic) use either way. Once vendored,
`GOPACK_BUILD_RANDOM_GEN=OFF` remains available for anyone who'd rather not
build them at all.

## Building

Requires CMake 3.18+ and a C compiler (for the vendored Triangle/Qhull
libraries) and C++17 compiler. Eigen is vendored via a local checkout at
`./eigen` if present (fully offline after that one-time clone), falling
back to CMake `FetchContent` (needs network access) if `./eigen` doesn't
exist -- Triangle and Qhull, by contrast, are committed directly into the
repo under `third_party/` and never touch the network at all (see "Vendored
third-party libraries" above). A JDK (`JAVA_HOME` set) is needed only for
the optional JNI target.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Build options (pass as `-D<OPTION>=OFF` to disable):

- `GOPACK_BUILD_CLI` (default ON) -- the `gopack` command-line executable
- `GOPACK_BUILD_JNI` (default ON, skipped automatically if no JDK is found)
  -- `gopack_jni.dll` / `libgopack_jni.dylib` for the Java bridge
  (`JNI.GOPackNative`, in `jni/java/`; package `JNI`, not
  `org.kensmath.gopack`, since 8/2026 -- see that class's doc comment)
- `GOPACK_BUILD_TESTS` (default ON)
- `GOPACK_BUILD_RANDOM_GEN` (default ON) -- build the vendored
  Triangle/Qhull libraries and the whole random-triangulation *generation*
  family built on them (`gopack::geom::delaunayPlane`/`convexHull3`/
  `randTriangulationSphere`/`randTriangulationPlane` and
  `Packer::randomDisc`/`randomSphere`/`randomRectangle`/`randomSquare`/
  `randomTri` -- see "Vendored third-party libraries" above). Disabling this
  also disables `tests/test_random_gen.cpp`/`tests/test_random_packers.cpp`
  and the CLI's `--random-*` flags (see "CLI usage" below);
  `GOPACK_HAVE_RANDOM_GEN` is defined for the rest of the codebase to
  `#ifdef` around when this is off.
- `GOPACK_USE_SUITESPARSE` (default OFF) -- swap Eigen's built-in sparse
  solvers for SuiteSparse/CHOLMOD; not wired up yet (Eigen was the chosen
  default for portability), but `core/include/gopack/SparseLinearSolver.h`
  is the seam where that would plug in

## CLI usage

```
gopack input.p -o output.p [--passes 200] [--eucl-out]
gopack input.p -o output.p --polygon [--corners v1,v2,v3,v4] [--angles a1,a2,a3,a4]
gopack --random-disc N -o output.p [--passes 200] [--eucl-out]
gopack --random-sphere N -o output.p [--passes 200]
gopack --random-square N -o output.p [--passes 200] [--eucl-out]
gopack --random-rectangle N[,aspect[,bdryN]] -o output.p [--passes 200] [--eucl-out]
gopack --random-tri intN,bdryN --graph x1,y1,x2,y2,... -o output.p [--cent cx,cy]
       [--passes 200] [--eucl-out]
```

`--polygon` switches to mode 2 (polygonal/rectangle packing). `--corners` is
a comma-separated list of 1-indexed boundary vertices, in counterclockwise
order, to use as polygon corners; if omitted, corners are inferred from the
input file's `VERT_LIST:` or chosen pseudo-randomly (mirroring
`setMode.m`'s own fallback behavior). `--angles` gives matching target
corner angles in radians; if omitted, all corners get equal angles (e.g.
exactly `pi/2` each for a 4-corner input, which is what triggers the
rectangle-specific layout in `setRectCenters`). With exactly 4 corners, the
CLI also prints the resulting aspect ratio (`getAspect`).

`--random-disc`/`--random-sphere`/`--random-square`/`--random-rectangle`
generate a fresh "geometrically random" triangulation (see
`Packer::randomDisc`/`randomSphere`/`randomSquare`/`randomRectangle` above)
instead of reading `<input.p>` -- these are mutually exclusive with the
positional `<input.p>` argument and with each other, and are only available
in a build with `GOPACK_BUILD_RANDOM_GEN` on (the default). `--random-square`
and `--random-rectangle` already leave the generated packing in polygonal
mode with 4 corners chosen automatically, so `--polygon`/`--corners`/
`--angles` are ignored for those two. `--random-rectangle` takes a single
comma-separated argument: `N` (interior point count) is required; `aspect`
(default 1) sets the rectangle to `[-aspect,aspect]x[-1,1]`; `bdryN` (default
computed from `N` and `aspect`, matching `randomRectangle.m`) overrides the
boundary point count.

`--random-tri intN,bdryN` generates a random triangulation of an **arbitrary**
closed polygonal region -- not just a disc/square/rectangle -- via
`Packer::randomTri(intN, bdryN, graph, cent)`. The boundary polygon is given
with the required `--graph x1,y1,x2,y2,...` flag: a flat x,y coordinate list,
at least 3 points, in order around the boundary (don't repeat the first
point at the end). `--cent cx,cy` optionally names a point inside `--graph`
to use as the packing's alpha (centering) vertex; if omitted, or the point
isn't actually inside the boundary, alpha is chosen automatically. Unlike
`--random-square`/`--random-rectangle`, this stays in max-pack mode (mode 1)
-- a generic region has no "corners" concept -- so `--polygon`/`--corners`/
`--angles` still apply afterward if you want polygonal mode on the result.

## Java usage

Two ways to compute a maximal packing (mode 1) from Java, depending on
whether your data starts out in a `*.p` file or already in memory:

```java
// From a file:
double[] radii = JNI.GOPackNative.computeMaximalPacking(
    "input.p", /* geometryHint (reserved) */ 0, /* tolerance (reserved) */ 0.0,
    /* maxPasses */ 200);

// From an in-memory complex (e.g. CirclePack's own per-vertex flower data):
int[][] flowers = new int[nodeCount + 1][]; // flowers[0] unused
// ... fill flowers[1..nodeCount], CLOSED (first==last) for interior
//     vertices, OPEN (first!=last) for boundary vertices ...
double[] radii2 = JNI.GOPackNative.computeMaximalPackingFromComplex(
    nodeCount, flowers, /* geometry: 0=eucl -1=hyp +1=sph */ 0,
    /* tolerance (reserved) */ 0.0, /* maxPasses */ 200);
```

A third method, `computeRandomTri`, generates a random triangulation of an
arbitrary closed polygonal region (the JNI counterpart to the CLI's
`--random-tri`) and computes its maximal packing in one call, and a fourth,
`computeRandomDisc`, does the same for the unit disc (the JNI counterpart to
`--random-disc`):

```java
double[] graphXY = { 0,0, 4,0, 4,3, 0,3 }; // a 4x3 rectangle, as an example
JNI.RandomComplexResult result =
    JNI.GOPackNative.computeRandomTri(
        /* intN */ 40, /* bdryN */ 20, graphXY,
        /* centX, centY, hasCent */ 0.0, 0.0, false,
        /* maxPasses */ 200);
// or: JNI.RandomComplexResult result = JNI.GOPackNative.computeRandomDisc(
//         /* n */ 200, /* maxPasses */ 200);
// result.nodeCount, result.flowers, result.radii, result.centersRe/centersIm,
// result.geometry, result.alpha, result.gamma are all populated (1-indexed,
// index 0 unused, same convention as radii/radii2 above) -- unlike
// computeMaximalPacking[FromComplex], the caller didn't supply the
// combinatorics, so the whole generated complex comes back, not just radii.
```

`radii`/`radii2` both have length `nodeCount+1`, and `radii[v]` is vertex
`v`'s euclidean radius for `v = 1..nodeCount` -- matching GOPack/CirclePack's
own 1-indexed vertex numbering (the same convention used throughout the C++
core), rather than shifting to a "natural" 0-indexed Java array. `radii[0]`
is unused. Load `gopack_jni.dll` / `libgopack_jni.dylib` via
`java.library.path`, or bundle both platform binaries and pick one at
runtime based on `os.name`/`os.arch` (see the class doc comment in
`GOPackNative.java`).

**Prefer `computeMaximalPackingFromComplex` over `computeMaximalPacking`
whenever the triangulation already exists in memory on the Java side** (as
it will for a CirclePack caller). `computeMaximalPacking` still has to open
and parse a `*.p` file on the C++ side (`Packer::readpack`), which for large
complexes can easily take longer than the packing computation itself --
that text parsing is real work regardless of whether it happens in a
subprocess or in-process. `computeMaximalPackingFromComplex` goes straight
to `Packer::loadComplex` (`core/src/PackerIO.cpp`), skipping both the
Java-side serialization to text and the C++-side parsing back out of it, so
its cost is close to the packing computation alone. `computeMaximalPacking`
remains the right choice when the data genuinely starts out as a file (a
`*.p` on disk with no in-memory representation yet).

Only mode 1 (maximal packing) is exposed through this JNI bridge so far;
mode 2 (polygonal/rectangle) is ported in the C++ core (both `readpack()`
and `loadComplex()` load a complex the same way regardless of which mode you
later select with `setMode`) but only reachable today via the CLI's
`--polygon` flag -- adding `computeMaximalPackingFromComplex`'s mode-2
counterpart is a small follow-up whenever you need it (same `loadComplex`
plumbing, just `setMode(2, corners, angles)` instead of `setMode(1)` before
`riffle`).

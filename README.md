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
- `complex_count`, `indxMatrices`, `FarVert`
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

**Understood but NOT yet ported** (calling these throws
`gopack::NotImplementedError` with a message pointing at the source file,
rather than guessing at behavior):

- The bare triangle-list and OFF file readers (`parse_triangles`)
- Orphan-vertex removal (`pruneComplex` -- only used by the random
  rectangle/square generators, not by loading real packing files)
- Random triangulation generation (`randTriangulation`, `randomDisc`,
  `randomSphere`, `randomRectangle`, etc.) and the plotting method (`show`,
  which has no headless equivalent anyway)

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

### One deliberate deviation from the literal source

`code/s_to_e_data.m` calls a function `proj_vec_to_sph` that does not exist
anywhere in the repository -- it's a typo for `proj_vec_to_s.m`. The C++
port (`geom::sToEData` in `core/src/Geometry.cpp`) calls the real function,
since that's unambiguously the intent; reproducing the bug would mean this
one rarely-hit code path (spherical circles enclosing the point at infinity)
throws in MATLAB and silently "works" here. This is the only place the port
knowingly differs from the literal source.

## What has -- and hasn't -- been verified

- Three regression tests (`tests/test_hex_flower.cpp`,
  `tests/test_readpack_roundtrip.cpp`, `tests/test_polygonal.cpp`) exercise
  the pipeline on the classical "hex flower" complex (one interior vertex of
  degree 6 ringed by 6 boundary vertices) in both modes, checking angle-sum
  convergence and symmetric radii for mode 1, a read/riffle/write/re-read
  round trip, and (for mode 2) that a 4-corner rectangle layout stays
  finite, positive, and rectangle-shaped after riffling.
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

## Building

Requires CMake 3.18+ and a C++17 compiler. Eigen is vendored via a local
checkout at `./eigen` if present (fully offline after that one-time clone),
falling back to CMake `FetchContent` (needs network access) if `./eigen`
doesn't exist. A JDK (`JAVA_HOME` set) is needed only for the optional JNI
target.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Build options (pass as `-D<OPTION>=OFF` to disable):

- `GOPACK_BUILD_CLI` (default ON) -- the `gopack` command-line executable
- `GOPACK_BUILD_JNI` (default ON, skipped automatically if no JDK is found)
  -- `gopack_jni.dll` / `libgopack_jni.dylib` for the Java bridge
  (`org.kensmath.gopack.GOPackNative`, in `jni/java/`)
- `GOPACK_BUILD_TESTS` (default ON)
- `GOPACK_USE_SUITESPARSE` (default OFF) -- swap Eigen's built-in sparse
  solvers for SuiteSparse/CHOLMOD; not wired up yet (Eigen was the chosen
  default for portability), but `core/include/gopack/SparseLinearSolver.h`
  is the seam where that would plug in

## CLI usage

```
gopack input.p -o output.p [--passes 20] [--eucl-out]
gopack input.p -o output.p --polygon [--corners v1,v2,v3,v4] [--angles a1,a2,a3,a4]
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

## Java usage

```java
double[] radii = org.kensmath.gopack.GOPackNative.computeMaximalPacking(
    "input.p", /* geometryHint (reserved) */ 0, /* tolerance (reserved) */ 0.0,
    /* maxPasses */ 20);
```

`radii` has length `nodeCount+1`, and `radii[v]` is vertex `v`'s euclidean
radius for `v = 1..nodeCount` -- matching GOPack/CirclePack's own 1-indexed
vertex numbering (the same convention used throughout the C++ core), rather
than shifting to a "natural" 0-indexed Java array. `radii[0]` is unused.
Load `gopack_jni.dll` / `libgopack_jni.dylib` via `java.library.path`, or
bundle both platform binaries and pick one at runtime based on
`os.name`/`os.arch` (see the class doc comment in `GOPackNative.java`).
Only mode 1 (maximal packing) is exposed through this JNI bridge so far;
mode 2 (polygonal/rectangle) is ported in the C++ core but only reachable
today via the CLI's `--polygon` flag.

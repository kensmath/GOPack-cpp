# GOPack C++

A C++ port of [GOPack](https://github.com/kensmath/GOPack) (Collins, Orick,
Stephenson, 2017) -- MATLAB software for computing maximal circle packings
with very large numbers of circles -- structured to build as a native
library for Java (via JNI) and as standalone Windows/macOS executables.

**License:** GPL-3.0, same as the original GOPack (this is a derivative
work; see LICENSE.md).

## Status: max-pack mode ported, everything else deferred

This is a **direct, line-by-line port**, not a reimplementation from the
paper's description. The translation was done by reading every `.m` file in
`code/@GOPacker/` and `code/` and translating each one into C++ against the
same field names, method names, and loop structure, specifically so the two
can be diffed and checked against each other rather than trusted on faith.
See the module-level comment in `core/include/gopack/Packer.h` for the
indexing convention (1-indexed, matching MATLAB, to keep every loop bound
identical to the source).

**Fully ported and covered by tests** (`gopack::Packer`, mode 1 = maximal
packing of the disc/plane, or the sphere if the complex has no boundary):

- `readpack` -- reads the `*.p` FLOWERS format (docs/GO_Formats.txt)
- `complex_count`, `indxMatrices`, `FarVert`
- `setMode` (mode 1), `layoutBdry` / `setHoroCenters`
- the core iteration: `continueRiffle` / `layoutCenters` (the sparse linear
  solve) / `setEffective` / `updateVdata` / `visualErrors`
- `riffle`, `reapResults`, `angsumErrors`, `packStatus`
- `writepack` / `writeEucl`, including hyperbolic and spherical output
  conversion (`e_to_h_data`, `h_to_e_data`, `e_to_s_data`, `s_to_e_data`,
  `sph_tangent`, `affineNormalizer`, `Centroid`, `loadTangency`)
- `cosAngle`, `cosCorner`

**Understood but NOT yet ported** (calling these throws
`gopack::NotImplementedError` with a message pointing at the source file,
rather than guessing at behavior):

- Polygonal / rectangle packing mode (`setMode` mode 2, `setPolyCenters`,
  `setRectCenters`, `getAspect`)
- The bare triangle-list and OFF file readers (`parse_triangles`)
- Orphan-vertex removal (`pruneComplex` -- only used by the random
  rectangle/square generators, not by loading real packing files)
- Random triangulation generation (`randTriangulation`, `randomDisc`,
  `randomSphere`, `randomRectangle`, etc.) and the plotting method (`show`,
  which has no headless equivalent anyway)

If your workflow needs any of the above, say so -- the source for all of
them has already been read and understood (see the corresponding `.m`
files), so porting them is a bounded follow-up, not new research.

### One deliberate deviation from the literal source

`code/s_to_e_data.m` calls a function `proj_vec_to_sph` that does not exist
anywhere in the repository -- it's a typo for `proj_vec_to_s.m`. The C++
port (`geom::sToEData` in `core/src/Geometry.cpp`) calls the real function,
since that's unambiguously the intent; reproducing the bug would mean this
one rarely-hit code path (spherical circles enclosing the point at infinity)
throws in MATLAB and silently "works" here. This is the only place the port
knowingly differs from the literal source.

## What has -- and hasn't -- been verified

- Two regression tests (`tests/test_hex_flower.cpp`,
  `tests/test_readpack_roundtrip.cpp`) exercise the full pipeline on the
  classical "hex flower" complex (one interior vertex of degree 6 ringed by
  6 boundary vertices) and check angle-sum convergence, symmetric radii, and
  a read/riffle/write/re-read round trip.
- **This code has not been compiled in the environment it was written in.**
  The sandbox this port was written in has no general internet access (no
  `apt`, `pip`, or `git clone` reaches outside a small allowlist), so Eigen
  could not be vendored locally to smoke-test the build. Every file was
  written and then re-read carefully for type errors and indexing mistakes,
  and the file-format assumptions in `readpack`/`writepack` were checked
  directly against real sample data (`data/Pinwheel_K.p`, 3081 vertices)
  line by line -- but a real compiler has not touched this code yet.
  `.github/workflows/build.yml` builds and runs the test suite on both
  Windows and macOS runners on every push; **treat the first CI run as part
  of this port, not as an afterthought**, and expect to fix a handful of
  compile errors.
- **Numerical agreement with the original MATLAB has not been checked.**
  That requires running both implementations on the same input and diffing
  radii, which needs a MATLAB (or Octave) environment this session doesn't
  have. Recommended next step: run
  `GOPacker.readpack(...); GOPacker.setMode(1); GOPacker.riffle(50);` in
  MATLAB on one of the `data/*.p` files and compare the resulting
  `obj.radii` against this port's `gopack_cli` output on the same file.

## Building

Requires CMake 3.18+, a C++17 compiler, and network access (to fetch Eigen
via `FetchContent` on first configure). A JDK (`JAVA_HOME` set) is needed
only for the optional JNI target.

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
```

## Java usage

```java
double[] radii = org.kensmath.gopack.GOPackNative.computeMaximalPacking(
    "input.p", /* geometryHint (reserved) */ 0, /* tolerance (reserved) */ 0.0,
    /* maxPasses */ 20);
```

`radii[i]` is the euclidean radius of vertex `i+1` (GOPack's 1-indexed
numbering). Load `gopack_jni.dll` / `libgopack_jni.dylib` via
`java.library.path`, or bundle both platform binaries and pick one at
runtime based on `os.name`/`os.arch` (see the class doc comment in
`GOPackNative.java`).

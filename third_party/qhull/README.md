# Vendored: Qhull (src/libqhull_r only)

**Version:** 2020.2 (8.1-alpha1), downloaded from the `master` branch of
<https://github.com/qhull/qhull> as a ZIP by hand (Ken) and committed here
directly -- not fetched by CMake, not a git submodule.

**License:** see `COPYING.txt` in this directory (kept verbatim). Summary: a
permissive attribution-style license -- free to copy, modify, and
redistribute provided copyright notices stay intact, `COPYING.txt` itself is
distributed alongside it, and any modifications are noted (which is why
GOPack-cpp's own logic lives entirely in `core/src/RandomGen.cpp` rather
than patching Qhull's vendored source -- there's nothing to note because
nothing here is modified). `REGISTER.txt` (optional voluntary registration
with the authors) is also kept, unmodified, per upstream's own
`README.txt`.

**What's vendored, and what's deliberately not:** only `src/libqhull_r/`
-- Qhull's reentrant (thread-safe) C library, the piece GOPack-cpp actually
links against (`qh_new_qhull()`, for `randTriangulation.m`'s sphere
convex-hull case -- the same underlying algorithm MATLAB's `convhulln`
wraps). The full upstream repository is roughly 4x this size once you
include:

- `src/libqhullcpp/`, `src/libqhullstatic*/` -- the C++ wrapper and the
  older non-reentrant C library. GOPack-cpp's binding
  (`core/src/RandomGen.cpp`) talks to `libqhull_r` directly in C, matching
  how it already talks to Triangle -- no need for either alternative.
- `src/qconvex/`, `src/qdelaunay/`, `src/qhalf/`, `src/qvoronoi/`,
  `src/rbox/`, `src/qhull/`, `src/testqset*/`, `src/user_eg*/` -- Qhull's
  own CLI frontends and test/example programs. Not needed; GOPack-cpp calls
  the library directly rather than shelling out to a `qconvex`/`qdelaunay`
  executable (which was the whole point of vendoring instead of the
  `ProcessBuilder`-based approach CirclePack already uses for `triangle`/
  `qhull` today -- see the root README).
- `html/`, `eg/`, `build/`, the Qt `.pro`/`.pri` project files, and the top
  level `CMakeLists.txt`/`Makefile` -- documentation and alternate build
  systems for the full project; GOPack-cpp builds the vendored subset via
  its own `CMakeLists.txt` in this directory instead.

Every file actually vendored (all 18 `.c` files and their headers under
`src/libqhull_r/`) is upstream's own, byte-for-byte, unmodified -- listed
explicitly in this directory's `CMakeLists.txt`, matching upstream's own
`libqhullr_SOURCES`/`libqhullr_HEADERS` lists in its top-level
`CMakeLists.txt` (as of the vendored version) so a future update just means
diffing that list against a fresh checkout, not re-deriving it from
scratch.

See `core/src/RandomGen.cpp` for how GOPack-cpp actually calls into this
(`qh_new_qhull()`, walking the resulting facet list, and the `qhT`
reentrant-context handle `libqhull_r`'s API threads through every call --
`libqhull_r.h`/`qhull_ra.h` document the reentrant API in more detail than
is worth repeating here).

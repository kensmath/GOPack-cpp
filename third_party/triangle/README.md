# Vendored: Triangle (Jonathan Richard Shewchuk)

**Version:** 1.6, July 28, 2005 (the last released version; this is genuinely
old, unmaintained code, which is exactly why it's vendored here rather than
fetched from the network at build time -- see the root README's "Notes on
vendoring vs. fetching" section).

**Source:** downloaded by hand from
<https://www.cs.cmu.edu/~quake/triangle.html> (Ken) and committed here
directly -- not fetched by CMake, not a git submodule. There is no official
git repository for Triangle; the CMU page is the canonical distribution
point.

**License:** see `README-upstream.txt` in this directory (the original
distribution's own README, kept verbatim) for the full text. Summary: free
for private, research, and institutional use; redistribution as part of a
commercial system requires direct arrangement with the author unless you're
just telling people where to get it for free. `triangle.c`'s own file header
carries the same notice and must not be removed (see the top of the file).

**What's vendored, and what's deliberately not:**

- `triangle.c`, `triangle.h` -- the actual library. Vendored as-is,
  byte-for-byte from the upstream download, **not modified** -- all of
  GOPack-cpp's own logic lives in `core/src/RandomGen.cpp`, which calls
  `triangulate()` rather than patching Triangle itself. This keeps future
  diffs against a fresh copy of Triangle trivial (there shouldn't be any).
- `README-upstream.txt` -- the original distribution's `README`, renamed to
  avoid colliding with this file; kept for the license text and file
  manifest.
- `tricall.c.reference` -- the original distribution's example program
  (`tricall.c`, renamed `.reference` so CMake never picks it up as a build
  source). Not compiled by GOPack-cpp; kept purely as a worked example of
  Triangle's calling convention, useful if `core/src/RandomGen.cpp`'s
  binding ever needs cross-checking against upstream's own usage.
- **Not vendored:** `showme.c` (a standalone X11 mesh-viewing GUI -- no
  headless equivalent needed, and pulling in an X11 dependency for a tool
  GOPack-cpp never calls made no sense), `A.poly` (a sample input file for
  Triangle's own CLI, not used by anything here), and the upstream
  `makefile` (GOPack-cpp builds Triangle via `CMakeLists.txt` in this
  directory instead -- see that file for the required `TRILIBRARY`/
  `ANSI_DECLARATORS`/`NO_TIMER` preprocessor switches and why each one is
  needed for a portable library build).

See `core/src/RandomGen.cpp` for how GOPack-cpp actually calls into this
(the `triangulateio` struct marshalling, and the `REAL`/`VOID` macro
convention `triangle.h` requires the includer to define -- undocumented
inside `triangle.h` itself beyond a comment, so worth reading that file's
binding code directly rather than guessing).

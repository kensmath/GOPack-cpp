# GOPack-cpp -- vendored copy of Jonathan Shewchuk's Triangle (v1.6, 2005),
# the plane Delaunay/constrained-Delaunay triangulator MATLAB's
# delaunayTriangulation ultimately wraps, and what randTriangulation.m's
# plane-region case needs. Only triangle.c/triangle.h are built -- showme.c
# (a standalone X11 mesh viewer, not something GOPack-cpp calls into) and
# the upstream makefile are intentionally not vendored at all; see
# README.md in this directory.
#
# Triangle must be compiled with three preprocessor switches to be usable
# (and portable) as a library rather than Triangle's own CLI:
#   TRILIBRARY       -- compiles triangulate()/trifree() instead of main();
#                        this is Triangle's own documented way to call it as
#                        a library (see the header comment in triangle.h).
#   ANSI_DECLARATORS -- modern (non-K&R) function declarations; also
#                        controls the REAL/VOID convention a caller of
#                        triangulate() must match, see gopack::geom's
#                        Triangle binding (core/src/RandomGen.cpp).
#   NO_TIMER         -- Triangle's un-guarded-by-TRILIBRARY internal timing
#                        code (inside the shared triangulate()/main() body)
#                        calls gettimeofday() from <sys/time.h> unless this
#                        is defined -- neither exists on Windows/MSVC, so
#                        this define is required for a portable build, not
#                        just an optimization.
add_library(gopack_triangle STATIC
    triangle.c
)
target_include_directories(gopack_triangle SYSTEM PUBLIC .)
target_compile_definitions(gopack_triangle PRIVATE
    TRILIBRARY
    ANSI_DECLARATORS
    NO_TIMER
)

set_target_properties(gopack_triangle PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)

if(MSVC)
    # Triangle's 1990s-era C uses plain malloc/free/printf and (in a few
    # spots even in TRILIBRARY mode) sprintf -- harmless, but MSVC's CRT
    # deprecates the unchecked forms by default; this is vendored code we
    # don't maintain, so suppress rather than "fix" it.
    target_compile_definitions(gopack_triangle PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _CRT_NONSTDC_NO_DEPRECATE)

    # Defensive hardening, not a confirmed fix for anything specific: this
    # file's own header comment (around the INEXACT/epsilon/splitter
    # definitions, ~line 260) documents that Triangle's robust orientation
    # predicates (exactinit()/the Shewchuk adaptive-precision routines
    # underlying every Delaunay decision it makes) are only correct under
    # strict, non-reassociated IEEE-754 double rounding -- the classic
    # failure mode is a compiler that's licensed to fuse/reorder/contract
    # floating-point expressions (MSVC's /fp:fast, or its default inherited
    # from some other CMAKE_CXX_FLAGS a toolchain/CI config might set)
    # silently producing a self-inconsistent triangulation (a face list that
    # looks structurally fine but doesn't actually satisfy the Delaunay
    # property/proper orientation), which is exactly the kind of malformed
    # input that could desync this project's own complexCount()/
    # indxMatrices() (see PackerCore.cpp's now-hardened boundary-walk
    # checks) or Triangle's own internal invariants further downstream.
    # cl.exe's default is already /fp:precise, so this is very likely a
    # no-op in practice -- but it removes the ambiguity of relying on
    # whatever the ambient toolchain default happens to be, for the one
    # translation unit in this whole project where getting it wrong would
    # matter, at essentially zero cost (this file is small and not
    # performance-critical).
    target_compile_options(gopack_triangle PRIVATE /fp:precise)
else()
    # Silence warnings from vendored code we don't maintain, matching the
    # SYSTEM include treatment above -- this target doesn't get GOPack-cpp's
    # own -Wall -Wextra (set globally in the root CMakeLists.txt).
    target_compile_options(gopack_triangle PRIVATE -w)

    # Same rationale as the MSVC /fp:precise above: GCC/Clang's equivalent
    # footguns for Shewchuk's exact predicates are -ffast-math (not used by
    # this project, but defensively excluded here too in case some ambient
    # CMAKE_C_FLAGS ever adds it) and FMA contraction reordering a chain of
    # double multiplies/adds into a single fused instruction with different
    # rounding -- explicit -ffp-contract=off removes that possibility
    # regardless of -O level or -march target.
    target_compile_options(gopack_triangle PRIVATE -ffp-contract=off)
endif()

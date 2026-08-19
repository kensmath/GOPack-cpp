# GOPack-cpp -- vendored subset of Qhull (src/libqhull_r only).
#
# This is deliberately NOT the full upstream qhull tree (docs, the C++
# wrapper libqhullcpp, the qconvex/qdelaunay/qhalf/qvoronoi/rbox CLI
# frontends, the build/ directory, and so on) -- just the reentrant C
# library (libqhull_r), which is all GOPack-cpp needs: it links directly
# against qh_new_qhull() for the sphere convex-hull case of
# randTriangulation.m, the same way MATLAB's own convhulln ultimately does.
# See README.md in this directory for what was omitted and why.
#
# Built as a plain STATIC library from the vendored source below -- no
# network access at configure or build time, matching GOPack-cpp's own
# Eigen-integration policy (see the root CMakeLists.txt comment).
add_library(gopack_qhull_r STATIC
    src/libqhull_r/accessors_r.c
    src/libqhull_r/geom2_r.c
    src/libqhull_r/geom_r.c
    src/libqhull_r/global_r.c
    src/libqhull_r/io_r.c
    src/libqhull_r/libqhull_r.c
    src/libqhull_r/mem_r.c
    src/libqhull_r/merge_r.c
    src/libqhull_r/poly2_r.c
    src/libqhull_r/poly_r.c
    src/libqhull_r/qset_r.c
    src/libqhull_r/random_r.c
    src/libqhull_r/rboxlib_r.c
    src/libqhull_r/stat_r.c
    src/libqhull_r/user_r.c
    src/libqhull_r/usermem_r.c
    src/libqhull_r/userprintf_r.c
    src/libqhull_r/userprintf_rbox_r.c
)

# SYSTEM: suppress -Wall/-Wextra/-W4 warnings originating from inside
# qhull's own headers, same rationale as the Eigen INTERFACE target in the
# root CMakeLists.txt -- we didn't write this code and don't want its
# warnings mixed in with ours.
target_include_directories(gopack_qhull_r SYSTEM PUBLIC src/libqhull_r)

set_target_properties(gopack_qhull_r PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)

if(MSVC)
    # Qhull's own I/O code uses plain sprintf/strcpy/fopen -- harmless (all
    # call sites use fixed-size internal buffers Qhull sizes itself), but
    # MSVC's CRT deprecates them by default and would otherwise flood the
    # build log with warnings from code we didn't write.
    target_compile_definitions(gopack_qhull_r PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _CRT_NONSTDC_NO_DEPRECATE)
else()
    # Silence warnings from vendored code we don't maintain, matching the
    # SYSTEM include treatment above -- these targets don't get GOPack-cpp's
    # own -Wall -Wextra (set globally in the root CMakeLists.txt).
    target_compile_options(gopack_qhull_r PRIVATE -w)
endif()

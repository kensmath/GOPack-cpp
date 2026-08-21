cmake_minimum_required(VERSION 3.18)
project(GOPackCpp
    VERSION 0.1.0
    DESCRIPTION "C++ port of GOPack (maximal circle packing) for large-scale packings"
    LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ---------------------------------------------------------------------------
# Build options
# ---------------------------------------------------------------------------
option(GOPACK_BUILD_CLI   "Build standalone command-line executable" ON)
option(GOPACK_BUILD_JNI   "Build JNI shared library for Java integration" ON)
option(GOPACK_BUILD_TESTS "Build unit / regression tests" ON)
option(GOPACK_USE_SUITESPARSE "Use SuiteSparse (CHOLMOD) instead of Eigen's built-in sparse solvers" OFF)
option(GOPACK_BUILD_RANDOM_GEN
       "Build the vendored Triangle/Qhull geometry libraries and gopack's random-triangulation generation (randTriangulation.m family)"
       ON)

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra)
endif()

# ---------------------------------------------------------------------------
# Static runtime linking
# ---------------------------------------------------------------------------
# gopack_jni is not just a build artifact for developers -- the prebuilt
# .dll/.so/.dylib is bundled straight into CirclePack's jar (cdeps/) and
# loaded via System.load() on end-user machines that may have nothing else
# installed. Relying on each platform's *default* (dynamic) C/C++ runtime
# linking would mean assuming the target machine already has a matching
# runtime present, which isn't safe on a genuinely clean install:
#   - Windows (MSVC): the CMake/MSVC default is /MD (dynamic CRT), which
#     makes gopack_jni.dll depend on the Visual C++ Redistributable being
#     installed on the target machine -- not guaranteed. Statically linking
#     the MSVC runtime (/MT) removes that dependency entirely.
#   - Linux (GCC/Clang, non-Apple): our CI builds inside a manylinux2014
#     (glibc 2.17) container using devtoolset-10 for a C++17-capable GCC, but
#     that newer GCC's own *default* libstdc++.so.6 can be newer than what's
#     actually present on the older target distros the manylinux2014
#     baseline is meant to support -- the classic manylinux/auditwheel
#     problem. -static-libgcc -static-libstdc++ statically links just those
#     two runtime libraries into gopack_jni.so, leaving only glibc/libc
#     itself as a runtime dependency (glibc keeps strong backward
#     compatibility, which is the whole point of the old baseline).
#   - macOS: no equivalent issue exists, so no fix is applied here --
#     Apple's system libc++/libSystem ship with the OS itself.
# This is set here (directory scope, before any add_subdirectory() calls
# below) so it's inherited by every target this project builds --
# third_party, core, cli, jni/cpp, and tests alike.
if(MSVC)
    # CMAKE_MSVC_RUNTIME_LIBRARY only takes effect under policy CMP0091=NEW;
    # that's already the default given cmake_minimum_required(VERSION 3.18)
    # above, but set it explicitly here so this block keeps working even if
    # the minimum-required version above is ever lowered.
    cmake_policy(SET CMP0091 NEW)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
elseif(NOT APPLE)
    add_link_options(-static-libgcc -static-libstdc++)
endif()

# ---------------------------------------------------------------------------
# Eigen (header-only sparse/dense linear algebra, MPL2 license)
# ---------------------------------------------------------------------------
# We only need the headers -- Eigen is header-only for our purposes, and we
# never touch its own tests/docs/benchmark targets -- so neither branch below
# runs Eigen's own CMakeLists.txt at all (no add_subdirectory(eigen), no
# FetchContent_MakeAvailable() which would do the same thing internally for
# the fetched copy). Both branches instead build a bare INTERFACE target,
# Eigen3::Eigen, that just points at Eigen's include directory. Two reasons:
#   1) It avoids a "Compatibility with CMake < 3.10 will be removed"
#      deprecation warning from Eigen's own cmake_minimum_required (as of the
#      3.4.0 tag pinned below) -- that warning is about Eigen's CMake
#      project, unrelated to this project's own CMAKE_MINIMUM_REQUIRED above,
#      and isn't fixable by editing anything in this repo.
#   2) It's a more direct match for "we only need the headers" than letting
#      Eigen's full CMake project configure (with EIGEN_BUILD_DOC/
#      BUILD_TESTING/etc. set OFF to suppress targets we'd never build).
if(EXISTS "${CMAKE_SOURCE_DIR}/eigen/Eigen/Core")
    # A local Eigen checkout (e.g. `git clone https://gitlab.com/libeigen/eigen.git`
    # into ./eigen) takes priority over fetching from the network -- this is
    # both faster on repeat configures and lets the project build fully
    # offline once Eigen has been cloned once.
    message(STATUS "GOPack: using local Eigen checkout at ${CMAKE_SOURCE_DIR}/eigen (headers only)")
    add_library(gopack_eigen_headers INTERFACE)
    target_include_directories(gopack_eigen_headers SYSTEM INTERFACE
        "${CMAKE_SOURCE_DIR}/eigen")
    add_library(Eigen3::Eigen ALIAS gopack_eigen_headers)
else()
    message(STATUS "GOPack: no local ./eigen checkout found; fetching Eigen 3.4.0 via CMake FetchContent (headers only)")
    include(FetchContent)
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.4.0
        GIT_SHALLOW    TRUE
    )
    # FetchContent_Populate (rather than FetchContent_MakeAvailable) only
    # downloads the source -- it deliberately does NOT add_subdirectory()
    # the result, for the same "headers only" reason as the local-checkout
    # branch above.
    FetchContent_GetProperties(Eigen3)
    if(NOT eigen3_POPULATED)
        FetchContent_Populate(Eigen3)
    endif()
    add_library(gopack_eigen_headers INTERFACE)
    target_include_directories(gopack_eigen_headers SYSTEM INTERFACE
        "${eigen3_SOURCE_DIR}")
    add_library(Eigen3::Eigen ALIAS gopack_eigen_headers)
endif()

# ---------------------------------------------------------------------------
# Vendored geometry libraries (Triangle, Qhull) -- see
# third_party/triangle/README.md and third_party/qhull/README.md.
# ---------------------------------------------------------------------------
# Unlike Eigen above, these are committed source, not fetched over the
# network: both are old, low-traffic upstream projects (Triangle in
# particular hasn't been updated since 2005 and has no official git repo at
# all), and this project's random-triangulation generation needs to build
# and work identically regardless of whether either upstream host is
# reachable at build time -- see the root README's "Notes on vendoring vs.
# fetching" section for the full rationale.
if(GOPACK_BUILD_RANDOM_GEN)
    add_subdirectory(third_party/triangle)
    add_subdirectory(third_party/qhull)
endif()

# ---------------------------------------------------------------------------
# Core engine (static library, no Java/CLI dependencies)
# ---------------------------------------------------------------------------
add_subdirectory(core)

if(GOPACK_BUILD_CLI)
    add_subdirectory(cli)
endif()

if(GOPACK_BUILD_JNI)
    add_subdirectory(jni/cpp)
endif()

if(GOPACK_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

install(FILES LICENSE.md DESTINATION share/gopack-cpp OPTIONAL)

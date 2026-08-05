cmake_minimum_required(VERSION 3.18)
project(GOPackCpp
    VERSION 0.1.0
    DESCRIPTION "C++ port of GOPack (maximal circle packing) for large-scale packings"
    LANGUAGES CXX)

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

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra)
endif()

# ---------------------------------------------------------------------------
# Eigen (header-only sparse/dense linear algebra, MPL2 license)
# ---------------------------------------------------------------------------
# We only need the headers; skip Eigen's own tests/docs/build targets.
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)

if(EXISTS "${CMAKE_SOURCE_DIR}/eigen/CMakeLists.txt")
    # A local Eigen checkout (e.g. `git clone https://gitlab.com/libeigen/eigen.git`
    # into ./eigen) takes priority over fetching from the network -- this is
    # both faster on repeat configures and lets the project build fully
    # offline once Eigen has been cloned once. add_subdirectory(eigen) defines
    # the same Eigen3::Eigen INTERFACE target that FetchContent_MakeAvailable
    # would below, so nothing downstream needs to know which path was taken.
    message(STATUS "GOPack: using local Eigen checkout at ${CMAKE_SOURCE_DIR}/eigen")
    add_subdirectory(eigen)
else()
    message(STATUS "GOPack: no local ./eigen checkout found; fetching Eigen 3.4.0 via CMake FetchContent")
    include(FetchContent)
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.4.0
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(Eigen3)
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

add_library(gopack_core STATIC
    src/Geometry.cpp
    src/PackerCore.cpp
    src/PackerIO.cpp
    src/PackerLayout.cpp
    src/SparseLinearSolver.cpp
)

target_include_directories(gopack_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(gopack_core PUBLIC Eigen3::Eigen)

target_compile_features(gopack_core PUBLIC cxx_std_17)

# RandomGen.cpp (gopack::geom::delaunayPlane/convexHull3/
# randTriangulationSphere/randTriangulationPlane) wraps the vendored
# Triangle/Qhull libraries (see third_party/*/README.md); PackerRandom.cpp
# (Packer::randomDisc/randomSphere/randomRectangle/randomSquare/randomTri)
# is built directly on top of it. Both are only built/linked when
# GOPACK_BUILD_RANDOM_GEN is on (the third_party subdirectories, and hence
# the gopack_triangle/gopack_qhull_r targets they define, only exist at all
# in that case -- see the root CMakeLists.txt). GOPACK_HAVE_RANDOM_GEN is a
# public compile definition (propagates to every consumer of gopack_core,
# including the CLI and JNI targets) so callers can #ifdef around
# RandomGen.h/the Packer.h randomTri()-family entry points when this
# feature is disabled, rather than getting a link error.
if(GOPACK_BUILD_RANDOM_GEN)
    target_sources(gopack_core PRIVATE src/RandomGen.cpp src/PackerRandom.cpp)
    target_link_libraries(gopack_core PRIVATE gopack_triangle gopack_qhull_r)
    target_compile_definitions(gopack_core PUBLIC GOPACK_HAVE_RANDOM_GEN=1)
endif()

if(GOPACK_USE_SUITESPARSE)
    find_package(CHOLMOD REQUIRED)
    target_link_libraries(gopack_core PUBLIC CHOLMOD::CHOLMOD)
    target_compile_definitions(gopack_core PUBLIC GOPACK_USE_SUITESPARSE=1)
endif()

install(TARGETS gopack_core ARCHIVE DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

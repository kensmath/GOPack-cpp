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

if(GOPACK_USE_SUITESPARSE)
    find_package(CHOLMOD REQUIRED)
    target_link_libraries(gopack_core PUBLIC CHOLMOD::CHOLMOD)
    target_compile_definitions(gopack_core PUBLIC GOPACK_USE_SUITESPARSE=1)
endif()

install(TARGETS gopack_core ARCHIVE DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

add_executable(test_hex_flower test_hex_flower.cpp)
target_link_libraries(test_hex_flower PRIVATE gopack_core)
add_test(NAME hex_flower COMMAND test_hex_flower)

add_executable(test_readpack_roundtrip test_readpack_roundtrip.cpp)
target_link_libraries(test_readpack_roundtrip PRIVATE gopack_core)
add_test(NAME readpack_roundtrip COMMAND test_readpack_roundtrip
         "${CMAKE_CURRENT_SOURCE_DIR}/data/small_disc.p")

add_executable(test_polygonal test_polygonal.cpp)
target_link_libraries(test_polygonal PRIVATE gopack_core)
add_test(NAME polygonal COMMAND test_polygonal)

add_executable(test_orthopack test_orthopack.cpp)
target_link_libraries(test_orthopack PRIVATE gopack_core)
add_test(NAME orthopack COMMAND test_orthopack)

add_executable(test_loadcomplex test_loadcomplex.cpp)
target_link_libraries(test_loadcomplex PRIVATE gopack_core)
add_test(NAME loadcomplex COMMAND test_loadcomplex)

add_executable(test_sphere_normalize test_sphere_normalize.cpp)
target_link_libraries(test_sphere_normalize PRIVATE gopack_core)
add_test(NAME sphere_normalize COMMAND test_sphere_normalize
         "${CMAKE_CURRENT_SOURCE_DIR}/data/sphtest1000.p")

add_executable(test_parse_triangles test_parse_triangles.cpp)
target_link_libraries(test_parse_triangles PRIVATE gopack_core)
add_test(NAME parse_triangles COMMAND test_parse_triangles)

add_executable(test_prune_complex test_prune_complex.cpp)
target_link_libraries(test_prune_complex PRIVATE gopack_core)
add_test(NAME prune_complex COMMAND test_prune_complex)

add_executable(test_rand_bdry_pts test_rand_bdry_pts.cpp)
target_link_libraries(test_rand_bdry_pts PRIVATE gopack_core)
add_test(NAME rand_bdry_pts COMMAND test_rand_bdry_pts)

if(GOPACK_BUILD_RANDOM_GEN)
    add_executable(test_random_gen test_random_gen.cpp)
    target_link_libraries(test_random_gen PRIVATE gopack_core)
    add_test(NAME random_gen COMMAND test_random_gen)

    add_executable(test_random_packers test_random_packers.cpp)
    target_link_libraries(test_random_packers PRIVATE gopack_core)
    add_test(NAME random_packers COMMAND test_random_packers)
endif()

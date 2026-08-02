add_executable(test_hex_flower test_hex_flower.cpp)
target_link_libraries(test_hex_flower PRIVATE gopack_core)
add_test(NAME hex_flower COMMAND test_hex_flower)

add_executable(test_readpack_roundtrip test_readpack_roundtrip.cpp)
target_link_libraries(test_readpack_roundtrip PRIVATE gopack_core)
add_test(NAME readpack_roundtrip COMMAND test_readpack_roundtrip
         "${CMAKE_CURRENT_SOURCE_DIR}/data/small_disc.p")

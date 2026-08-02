add_executable(gopack_cli main.cpp)
target_link_libraries(gopack_cli PRIVATE gopack_core)
set_target_properties(gopack_cli PROPERTIES OUTPUT_NAME "gopack")

install(TARGETS gopack_cli RUNTIME DESTINATION bin)

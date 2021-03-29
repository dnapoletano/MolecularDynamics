add_compile_options(${myCOMPILE_FLAGS})
AUX_SOURCE_DIRECTORY("${PROJECT_SOURCE_DIR}/src" source)
include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${YAML_PREFIX}/include")
link_directories ("${YAML_PREFIX}/lib")

link_libraries(yaml-cpp)

add_executable(${PROJECT_NAME} ${source})

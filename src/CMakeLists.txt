add_compile_options(${myCOMPILE_FLAGS})
AUX_SOURCE_DIRECTORY("${PROJECT_SOURCE_DIR}/src" source)
include_directories ("${PROJECT_SOURCE_DIR}/include")

add_executable(${PROJECT_NAME} ${source})

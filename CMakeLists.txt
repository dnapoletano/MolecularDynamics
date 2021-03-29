cmake_minimum_required(VERSION 3.0)
project(MolecularDynamics LANGUAGES CXX)
set(CMAKE_MAKE_PROGRAM "Makefile")

message("------------ Molecular Dynamics --------------")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(Yaml-cpp REQUIRED)

if(${Yaml-cpp_FOUND})
  message(STATUS "Yaml-correctly found")
  set(YAML_PREFIX "${YAML_CPP_INCLUDE_DIR}/..")
endif()

if(NOT YAML_PREFIX)
  message( FATAL_ERROR
    "yaml-cpp path not found automatically. If you have it installed, set it via -DYAML_PREFIX=<path-to-yaml-cpp>.
If you don't have yaml-cpp installed, please install it and try again" )
endif()

add_subdirectory(src)
message("------------   End Configure    --------------")
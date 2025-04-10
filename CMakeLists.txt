# CMAKE Version
cmake_minimum_required(VERSION 3.13)

# Project information
project("Merope")
set(VERSION "1.0.0")

option(MEROPE_WARNING "Decide whether warning are displayed at compilation" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts")

# Build types and coverage
set(CMAKE_BUILD_TYPES Release Debug)

if(MEROPE_COVERAGE)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
endif()


enable_testing()

add_subdirectory(modules)
add_subdirectory(tools)



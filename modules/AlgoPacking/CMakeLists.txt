cmake_minimum_required(VERSION 3.12)
project(rsa
	VERSION 1.1
	LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)

include(cmake/modules/compiler.cmake)
include(cmake/modules/common-compiler-flags.cmake)

add_subdirectory(src)
add_subdirectory(Interface_python)

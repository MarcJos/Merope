# CMAKE Version
cmake_minimum_required(VERSION 3.13)

# Project information
project("Merope")
set(VERSION "1.0.0")

option(MEROPE_WARNING "Decide whether warning are displayed at compilation" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CXX_STANDARD_REQUIRED ON)

add_subdirectory(modules)
add_subdirectory(tools)


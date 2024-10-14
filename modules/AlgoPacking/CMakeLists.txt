cmake_minimum_required(VERSION 3.12)
project(rsa
	VERSION 1.1
	LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)
add_subdirectory(Interface_python)

install(DIRECTORY src DESTINATION AlgoPacking)


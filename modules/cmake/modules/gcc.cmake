# Options de base
macro(add_cxx_flag flag)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag} -fPIC")
endmacro(add_cxx_flag)

add_cxx_flag("-std=c++2a -fconcepts -Wall -W -Wextra -pedantic -Wshadow")
add_cxx_flag("-Wpointer-arith -Wcast-qual -Wcast-align")
add_cxx_flag("-Wwrite-strings -Wctor-dtor-privacy -Wnon-virtual-dtor")
add_cxx_flag("-Woverloaded-virtual -Wreturn-type -Wfloat-equal")
add_cxx_flag("-Wno-endif-labels -Wsign-compare -Wmissing-format-attribute")
add_cxx_flag("-Wno-multichar -Wno-deprecated-declarations -Wpacked")
add_cxx_flag("-Wredundant-decls -Wlong-long")
add_cxx_flag("-Wunknown-pragmas -Wundef -Wreorder")


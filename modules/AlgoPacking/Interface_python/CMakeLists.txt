cmake_minimum_required(VERSION 3.4...3.18)
project(cmake_example)

add_subdirectory(pybind11)

# clean way
pybind11_add_module(sac_de_billes Interface_AlgoPacking.cxx)
target_link_libraries(sac_de_billes PRIVATE RSAalgolib)
target_compile_definitions(sac_de_billes PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_compile_options(sac_de_billes PRIVATE -fPIC)
install(TARGETS sac_de_billes DESTINATION lib)

# for retrocompatibility
pybind11_add_module(py_rsa Interface_AlgoPacking.cxx)
target_link_libraries(py_rsa PRIVATE RSAalgolib)
target_compile_definitions(py_rsa PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_compile_options(py_rsa PRIVATE -fPIC)
install(TARGETS py_rsa DESTINATION lib)


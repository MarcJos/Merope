add_executable(test_equivalenceClass equivalenceClass.cxx)
target_include_directories(test_equivalenceClass PUBLIC ../../include/)
target_link_libraries(test_equivalenceClass PUBLIC OpenMP::OpenMP_CXX)
install(TARGETS test_equivalenceClass DESTINATION bin)
### Voro++ : Cmake installer

add_library(vorolib cell.cc common.cc container.cc  unitcell.cc v_compute.cc c_loops.cc v_base.cc wall.cc pre_container.cc container_prd.cc)
target_include_directories(vorolib SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(voro++ cmd_line.cc)
target_link_libraries(voro++ vorolib)
target_include_directories(voro++ SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


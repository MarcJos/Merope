add_subdirectory(interface_amitex_fftp)
add_subdirectory(interface_tmfft)

FILE(GLOB PY_FILES_1 "*.py")
install(FILES ${PY_FILES_1} DESTINATION "tools/python")


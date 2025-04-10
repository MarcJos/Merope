# Define the function to add a test
function(add_my_test NAME)
    add_test(NAME ${NAME} COMMAND ${CMAKE_INSTALL_PREFIX}/bin/${NAME})
endfunction()

add_my_test(Laguerre_opti)
add_my_test(Main_pour_test)
add_my_test(Optimization_test)
add_my_test(Performance)
add_my_test(Performance_1)
add_my_test(RSAalgo)
add_my_test(Z_texture)
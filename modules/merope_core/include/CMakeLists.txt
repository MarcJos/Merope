FILE(GLOB HXX_FILES "*.hxx" "*/*.hxx")
FILE(GLOB IXX_FILES "*.ixx" "*/*.ixx")
install(FILES ${HXX_FILES} ${IXX_FILES} "localMFront.h" 
  DESTINATION "include")

FILE(GLOB HXX_FILES "FFTW3/*.hxx")
FILE(GLOB IXX_FILES "FFTW3/*.ixx")
install(FILES ${HXX_FILES} ${IXX_FILES} "FFTW3/NumericTypes.h" 
  DESTINATION "include/FFTW3")

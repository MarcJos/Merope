//! Copyright : see license.txt
//!
//! \briefNumeric type definitions

#ifndef _LIB_TMFFT_FFTW3NUMERICTYPES_H_
#define _LIB_TMFFT_FFTW3NUMERICTYPES_H_ 

#include <fftw/fftw3.h>

#ifdef	__cplusplus
extern "C" {
#endif /* __cplusplus */

  /*!
   * \brief numeric type used for all computations
   *
   * This type is defined in the global
   * namespace to be able to call external
   * C functions.
   */
#ifdef TMFFT_USE_FLOAT
  typedef float rfloat;
  typedef fftwf_complex cfloat;
  // Les différents plans
#define FFTW_PREF(name) fftwf_ ## name
#else /* TMFFT_USE_FLOAT */
#ifdef  TMFFT_USE_LONGDOUBLE
  typedef long double rfloat;
  typedef fftwl_complex cfloat;
  // Les différents plans
#define FFTW_PREF(name) fftwl_ ## name
#else
  typedef double rfloat;
  typedef fftw_complex cfloat;
  // Les différents plans
#define FFTW_PREF(name) fftw_ ## name
#endif /* TMFFT_USE_LONGDOUBLE */
#endif /* TMFFT_USE_FLOAT */

#ifdef	__cplusplus
}
#endif /* __cplusplus */

#endif /* _LIB_TMFFT_FFTW3NUMERICTYPES_H */


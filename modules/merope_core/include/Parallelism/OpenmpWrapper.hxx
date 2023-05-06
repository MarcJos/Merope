//! Copyright : see license.txt
//!
//! \brief OpenMP Wrapper
//!

#ifndef _OPENMP_WRAPPER_HXX
#define _OPENMP_WRAPPER_HXX 1

#include "../MeropeNamespace.hxx"
// #include <omp.h> -> bad interaction with python lambdas

namespace merope {
namespace localFFT {
//! \return the thread number
unsigned short get_num_threads();

//! Set the thread number
//! \param n Thread number
void set_num_threads(unsigned short n);
}

} // namespace merope


#ifdef TMFFT_USE_FFTW3
#include<fftw3.h>
#endif

#endif // __OPENMP_WRAPPER_HXX

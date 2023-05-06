//! Copyright : see license.txt
//!
//! \brief
//
#include <omp.h>

#include "Parallelism/OpenmpWrapper.hxx"
#include "Parallelism/SetNbOfThreads.hxx"


#include "MeropeNamespace.hxx"


namespace merope {

void setNbOfThreads(int n) {
  omp_set_num_threads(n);
  localFFT::set_num_threads(n);
#pragma omp parallel
  { // begin parallel section
#pragma omp master
    {
      cerr << "Number of openmp threads = " << omp_get_num_threads() << endl;
    }
  } // end parallel section
}

} // namespace merope

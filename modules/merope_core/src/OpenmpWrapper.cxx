//! Copyright : see license.txt
//!
//! \brief OpenMP Wrapper
//!


#include <stdexcept>
#include "FFTW3/NumericTypes.h"
#include <iostream>

#include "Parallelism/OpenmpWrapper.hxx"
#include "MeropeNamespace.hxx"


namespace merope {
namespace localFFT {

unsigned short get_num_threads() {
    return 1;
}

int estimateBestNbOfCores() {
#ifdef __APPLE__
    const char* sysinfo = "sysctl -n hw.physicalcpu";
#else
    const char* sysinfo = "LANG=C lscpu | grep ^Core | sed 's/.*: *//'";
#endif
    const unsigned BUFFER_SIZE = 512;
    char buffer[BUFFER_SIZE];
    // Retrieve number of cores from the system information
    FILE* pf = popen(sysinfo, "r");
    if (!pf) {
        fprintf(stderr, "Failed to retrieve number of cores.\n");
        return 1;
    }
    fgets(buffer, BUFFER_SIZE, pf);
    int nbcores = atoi(buffer);
    pclose(pf);

    // If OMP_NUM_THREADS is defined, impose number of cores
    // as OMP_NUM_THREADS value
    const char* envvar = "OMP_NUM_THREADS";

    // Get the value of the environment variable
    if (getenv(envvar)) {
        snprintf(buffer, BUFFER_SIZE, "%s", getenv(envvar));
        nbcores = atoi(buffer);
    }

    if (nbcores < 1) nbcores = 1;

    return nbcores;
}

void set_num_threads(unsigned short n) {
    if (!n) {
        n = estimateBestNbOfCores();
    }
#ifdef TMFFT_USE_FFTW3
    if (not FFTW_PREF(init_threads)()) {
        throw logic_error(" Init problem from multi-threaded implementation of FFTW");
    }
    FFTW_PREF(plan_with_nthreads)(static_cast<int>(n));
#endif
}

} // namespace localFFT
} // namespace merope


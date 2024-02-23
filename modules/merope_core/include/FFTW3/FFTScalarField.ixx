//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MEROPE_CORE_SRC_FFTW3_FFTSCALARFIELD_IXX_
#define MEROPE_CORE_SRC_FFTW3_FFTSCALARFIELD_IXX_

#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <type_traits>


#include "../MeropeNamespace.hxx"


namespace merope {

template<class C>
void FFTScalarField::setCov(const double Lx, const double Ly, const double Lz,
    const C& cs) {
    setSpatialField();
    size_t nx, ny;
    getDim(nx, ny);
    double dx, dy, dz;
    switch (d) {
    case 1:
        dx = dy = 0;
        dz = Lx / nz;
        break;
    case 2:
        dx = 0;
        dy = Lx / ny;
        dz = Ly / nz;
        break;
    default:
        dx = Lx / nx;
        dy = Ly / ny;
        dz = Lz / nz;
    }
    char Dnz = nzb - nz;
#pragma omp parallel
    { // begin parallel section
#pragma omp for schedule(static)
        for (size_t i = 0; i < nx; ++i) {
            rfloat* vi = f + i * ny * nzb;
            double hx;
            if (i < (nx + 1) / 2)
                hx = dx * i;
            else
                hx = ((int)i - (int)nx) * dx;
            for (size_t j = 0; j < ny; ++j) {
                double hy;
                if (j < (ny + 1) / 2)
                    hy = dy * j;
                else
                    hy = ((int)j - (int)ny) * dy;
                for (size_t k = 0; k < nz; ++k, ++vi) {
                    double hz;
                    if (k < (nz + 1) / 2)
                        hz = dz * k;
                    else
                        hz = ((int)k - (int)nz) * dz;
                    if constexpr (std::is_same<C, gaussianField::CovSum>::value) {
                        switch (d) {
                        case 1:
                            *vi = cs.cov(hz);
                            break;
                        case 2:
                            *vi = cs.cov(hy, hz);
                            break;
                        case 3:
                            *vi = cs.cov(hx, hy, hz);
                            break;
                        }
                    } else if constexpr (std::is_same<C, std::function<double(std::array<double, 3>)>>::value) {
                        *vi = cs(std::array<double, 3>{hx, hy, hz});
                    } else if constexpr (std::is_same<C, std::function<double(std::array<double, 2>)>>::value) {
                        *vi = cs(std::array<double, 2>{hy, hz});
                    } else {
                        cerr << __PRETTY_FUNCTION__ << endl;
                        throw runtime_error("Unexpected");
                    }
                }
                // In place shift
                vi += Dnz;
            }
        }
    } // end parallel section
}


template<class FUNCTION>
void FFTScalarField::loopOnFrequencies(const FUNCTION& function) {
    checkSpectral("loopOnFrequencies");
    size_t nx, ny;
    getDim(nx, ny);
#pragma omp parallel default(shared)
    {
        // Limits and indexes for loops
        int i2, j2, ncx, ncy, ncx2, ncy2;
        unsigned i, j, ncz, ncz2;
        // Dimension 1,2 and 3 constants
        ncx = nx / 2 + 1;
        if (nx % 2) {
            ncx2 = ncx;
        } else {
            ncx2 = ncx - 1;
        }

        ncy = ny / 2 + 1;
        if (ny % 2) {
            ncy2 = ncy;
        } else {
            ncy2 = ncy - 1;
        }

        ncz = nz / 2 + 1;
        if (nz % 2) {
            ncz2 = ncz;
        } else {
            ncz2 = ncz - 1;
        }
#pragma omp for collapse(2)
        // Loop on x frequencies
        for (i = 0; i < nx; ++i) {
            // Loop on y frequencies
            for (j = 0; j < ny; ++j) {

                cfloat* pF = &F[(i * ny + j) * ncz];
                i2 = i;
                if (i >= (unsigned)ncx) {
                    i2 = i - nx;
                }
                bool fx = (i2 != ncx2 && i2);
                j2 = j;
                if (j >= (unsigned)ncy) {
                    j2 = j - ny;
                }
                bool fy = (j2 != ncy2 && j2);

                // Loop on z frequencies
                for (unsigned k = 0; k < ncz; ++k) {
                    cfloat& rF = pF[k];
                    bool fz = (k != ncz2 && k);
                    function(fx, fy, fz, rF);
                }
            }
        }
    }
}

template<class COVARIANCE_TYPE>
void FFTScalarField::build(const Grid* grid, const COVARIANCE_TYPE& cs,
    bool IP, int seed, unsigned flags_, bool showCovariance) {
    if (flags_ != FFTW_ESTIMATE) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Unexpected");
    }
    // Variables number
    nv = 1;
    alloc(IP);
    // Fill the grid with a covariance function
    setCov(grid->getLx(), grid->getLy(), grid->getLz(), cs);
    forward();
    // Generate the random field
    if (showCovariance) {
        prepareCovarianceInFourier();
    } else {
        randFunc(seed);
    }
    backward();
}

} // namespace merope

#endif /* MEROPE_CORE_SRC_FFTW3_FFTSCALARFIELD_IXX_ */

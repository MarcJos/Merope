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
                    }
                    else if constexpr (std::is_same<C, std::function<double(std::array<double, 3>)>>::value) {
                        *vi = cs(std::array<double, 3>{hx, hy, hz});
                    }
                    else if constexpr (std::is_same<C, std::function<double(std::array<double, 2>)>>::value) {
                        *vi = cs(std::array<double, 2>{hy, hz});
                    }
                    else {
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

} // namespace merope

#endif /* MEROPE_CORE_SRC_FFTW3_FFTSCALARFIELD_IXX_ */

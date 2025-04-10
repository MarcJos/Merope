//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

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
    {  // begin parallel section
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
                    if constexpr (std::is_same_v<C, gaussianField::CovSum>) {
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
                    } else if constexpr (std::is_same_v<C, std::function<double(std::array<double, 3>)>>) {
                        *vi = cs(std::array<double, 3>{hx, hy, hz});
                    } else if constexpr (std::is_same_v<C, std::function<double(std::array<double, 2>)>>) {
                        *vi = cs(std::array<double, 2>{hy, hz});
                    } else {
                        Merope_error_impossible();
                    }
                }
                // In place shift
                vi += Dnz;
            }
        }
    }  // end parallel section
}


template<class FUNCTION>
void FFTScalarField::loopOnFrequencies(const FUNCTION& function, bool use_omp) {
    checkSpectral("loopOnFrequencies");
    size_t nx, ny;
    getDim(nx, ny);
#pragma omp parallel if(use_omp) default(shared)
    {
        // Limits and indexes for loops
        auto get_nc = [](size_t nx_) {
            size_t ncx = nx_ / 2 + 1;
            size_t ncx2 = (nx_ % 2) ? ncx : ncx - 1;
            return tuple<size_t, size_t>(ncx, ncx2);
            };
        auto get_f = [](auto i_, auto ncx2_) {
            return (i_ != ncx2_ && i_);
            };
        // Dimension 1,2 and 3 constants
        auto [ncx, ncx2] = get_nc(nx);
        auto [ncy, ncy2] = get_nc(ny);
        auto [ncz, ncz2] = get_nc(nz);

#pragma omp if(use_omp) for collapse(2)
        // Loop on x frequencies
        for (size_t i = 0; i < nx; ++i) {
            // Loop on y frequencies
            for (size_t j = 0; j < ny; ++j) {
                bool fx = get_f(i, ncx2);
                auto fy = get_f(j, ncy2);

                cfloat* pF = &F[(i * ny + j) * ncz];
                // Loop on z frequencies
                for (size_t k = 0; k < ncz; ++k) {
                    cfloat& rF = pF[k];
                    bool fz = get_f(k, ncz2);
                    function(fx, fy, fz, rF);
                }
            }
        }
    }
}

template<class COVARIANCE_TYPE>
void FFTScalarField::build(const Grid* grid, const COVARIANCE_TYPE& cs,
    bool IP, int seed, unsigned flags_, bool showCovariance) {
    Merope_assert(flags_ == FFTW_ESTIMATE, "Unexpected");
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

}  // namespace merope



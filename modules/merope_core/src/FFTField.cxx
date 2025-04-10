//! Copyright : see license.txt
//!
//! \briefFFT fields based on FFTW library
//!

#include "FFTW3/FFTField.hxx"

#include <stdexcept>  // use assert for all checks
#include <iomanip>
#include <cstring>


#include "Grid/Grid.hxx"


namespace merope {

// Number of active planes to handle memory freeing of FFTW
unsigned short activePlanesNb_g = 0;


//----------------------------------------------
// FFT field mother class
//----------------------------------------------

FFTField::FFTField(const Grid& grid, unsigned flags_i) :
  isSpatial(true), f(NULL), F(NULL), forwardPlan(NULL), backwardPlan(NULL), flags(flags_i) {
  d = grid.getDim();  // Grid dimensions
  n[0] = grid.getNx();
  n[1] = grid.getNy();
  n[2] = grid.getNz();

  switch (d) {
  case 1:
    nz = n[0];
    fSize = n[0];
    FSize = n[0] / 2 + 1;
    break;
  case 2:
    nz = n[1];
    fSize = (size_t)n[0] * n[1];
    FSize = (size_t)n[0] * (n[1] / 2 + 1);
    break;
  case 3:
    nz = n[2];
    fSize = (size_t)n[0] * n[1] * n[2];
    FSize = (size_t)n[0] * n[1] * (n[2] / 2 + 1);
    break;
  }
}


void FFTField::alloc(const bool IP) {
  F = (cfloat*)(FFTW_PREF(malloc)(sizeof(cfloat) * FSize * nv));
  if (IP) {
    f = (rfloat*)F;
    nzb = 2 * (nz / 2 + 1);
  } else {
    f = (rfloat*)(FFTW_PREF(malloc)(sizeof(rfloat) * fSize * nv));
    nzb = nz;
  }
}


void FFTField::getDim(size_t& nx, size_t& ny) const {
  switch (d) {
  case 1:
    // nz: FFT dimensions for 1D
    nx = 1;
    ny = 1;
    break;
  case 2:
    // nz,ny: FFT dimensions for 2D
    nx = 1;
    ny = n[0];
    break;
  case 3:
    // 3D: Immediate correspondance
    nx = n[0];
    ny = n[1];
    break;
  }
}


FFTField& FFTField::operator=(const FFTField& field) {
  if (this != &field) {
    if ((n[0] != field.n[0]) ||
      (n[1] != field.n[1]) ||
      (n[2] != field.n[2]) ||
      (nv != field.nv) ||
      (fSize != field.fSize)) {
      throw logic_error("Invalid affectation: FFTfields must have the same sizes");
    }
    isSpatial = field.isSpatial;
#ifdef _OPENMP
#pragma omp  parallel for simd schedule(static)
    for (size_t i = 0; i < FSize * nv; i++) {
      F[i][0] = field.F[i][0];
      F[i][1] = field.F[i][1];
    }
#else
    memcpy(F, field.F, sizeof(cfloat) * FSize * nv);
#endif
    }
  return *this;
  }


FFTField::~FFTField() {
  if (f != (rfloat*)F) { FFTW_PREF(free)(f); }
  if (F) { FFTW_PREF(free)(F); }
  if (forwardPlan) {
    FFTW_PREF(destroy_plan)(forwardPlan);
    --activePlanesNb_g;
  }
  if (backwardPlan) {
    FFTW_PREF(destroy_plan)(backwardPlan);
    --activePlanesNb_g;
  }

  // free all memory
  if (!activePlanesNb_g) { FFTW_PREF(cleanup)(); }
}


void FFTField::checkSpatial(const std::string& mName) const {
  if (!isSpatial) throw logic_error("FFTField::" + mName + ": field must be in spatial representation");
  if (!f) throw logic_error("FFTField::" + mName + ": spatial field must exist");
}


void FFTField::checkSpectral(const std::string& mName) const {
  if (isSpatial) throw logic_error("FFTField::" + mName + ": field must be in spectral representation");
  if (!F) throw logic_error("FFTField::" + mName + ": spectral field must exist");
}


void FFTField::buildForwardPlan() {
  forwardPlan = FFTW_PREF(plan_many_dft_r2c)(d, n, nv, f, NULL, nv, 1, F, NULL, nv, 1, flags);
  if (!forwardPlan) throw(logic_error("FFTField::buildForwardPlan: unable to build the plan"));
  ++activePlanesNb_g;
}


void FFTField::buildBackwardPlan() {
  backwardPlan = FFTW_PREF(plan_many_dft_c2r)(d, n, nv, F, NULL, nv, 1, f, NULL, nv, 1, flags);
  if (!backwardPlan) throw(logic_error("FFTField::buildBackwardPlan: unable to build the plan"));
  ++activePlanesNb_g;
}


void FFTField::forward() {
  checkSpatial("forward");
  if (!forwardPlan) { buildForwardPlan(); }
  FFTW_PREF(execute)(forwardPlan);
  isSpatial = false;

  double ifSize = 1. / fSize;
  // Normalization
#pragma omp parallel for simd schedule(static)
  for (size_t i = 0; i < FSize * nv; i++) {
    F[i][0] *= ifSize;
    F[i][1] *= ifSize;
  }
}


void FFTField::backward() {
  checkSpectral("backward");
  if (!backwardPlan) { buildBackwardPlan(); }
  FFTW_PREF(execute)(backwardPlan);
  isSpatial = true;
}


void FFTField::print(std::ostream& os) const {
  if (isSpatial) {
  } else {
    if (F) {
      int ncx = n[0] / 2 + 1;
      int ncz = n[2] / 2 + 1;
      for (unsigned char l = 0; l < nv; ++l) {
        // for (size_t i=0; i<FSize; ++i) {
        //   os << F[l+nv*i][0] <<","<< F[l+nv*i][1] << endl;
        // }

        // For comparison with 2DECOMP&FFT
        // for (int k=0; k<ncz; ++k) {
        //   for (int j=0; j<n[1]; ++j) {
        //     for (int i=0; i<n[0]; ++i) {
        //       os << F[l+nv*(k+ncz*(j+i*n[1]))][0] <<","<< F[l+nv*(k+ncz*(j+i*n[1]))][1] << endl;
        //     }
        //   }
        // }

        for (int k = 0; k < ncz; ++k) {
          for (int j = 0; j < n[1]; ++j) {
            for (int i = 0; i < ncx; ++i) {
              os << F[l + nv * (k + ncz * (j + i * n[1]))][0] << "," << F[l + nv * (k + ncz * (j + i * n[1]))][1] << endl;
            }
          }
        }
        for (int k = ncz; k < n[2]; ++k) {
          for (int j = 0; j < n[1]; ++j) {
            int j2;
            if (j)
              j2 = n[1] - j;
            else
              j2 = j;
            for (int i = 0; i < ncx; ++i) {
              int i2;
              if (i)
                i2 = n[0] - i;
              else
                i2 = i;
              os << F[l + nv * ((n[2] - k) + ncz * (j2 + i2 * n[1]))][0] << ","
                << -F[l + nv * ((n[2] - k) + ncz * (j2 + i2 * n[1]))][1] << endl;
            }
          }
        }

        os << endl;
      }
    }
  }
}


void FFTField::setSpatialField() {
  isSpatial = true;
}


void FFTField::setSpectralField() {
  isSpatial = false;
}


void FFTField::setSpectralField(const CastemReal* const values) {
  cfloat* F2 = F;

  // Null frequency
  for (unsigned char i = 0; i < nv; ++i, ++F2) {
    (*F2)[0] = values[i];
    (*F2)[1] = 0.;
  }

  // remaining
  for (size_t i = nv; i < nv * FSize; ++i, ++F2) {
    (*F2)[0] = 0.;
    (*F2)[1] = 0.;
  }
  isSpatial = false;
}


void FFTField::SpectralZeroMean() {
  checkSpectral("SpectralZeroMean");

  int nx = 1, ny = 1;
  switch (d) {
  case 2:
    ny = n[0];
    break;
  case 3:
    nx = n[0];
    ny = n[1];
    break;
  }

  unsigned ncz = nz / 2 + 1;
  bool xp = !(nx % 2);
  bool yp = !(ny % 2);
  bool zp = !(nz % 2);
  for (unsigned char i = 0; i < 2; ++i) {
    if (!i or xp) {
      unsigned i2;
      if (i)
        i2 = nx / 2;
      else
        i2 = 0;
      for (unsigned char j = 0; j < 2; ++j) {
        if (!j or yp) {
          unsigned j2;
          if (j)
            j2 = ny / 2;
          else
            j2 = 0;
          for (unsigned char k = 0; k < 2; ++k) {
            if (!k or zp) {
              unsigned k2;
              if (k) k2 = nz / 2;
              else k2 = 0;
              for (unsigned char l = 0; l < nv; ++l) {
                cfloat* F2 = F + l + nv * (k2 + ncz * (j2 + ny * (size_t)i2));
                (*F2)[0] = (*F2)[1] = 0.;
              }
            }
          }
        }
      }
    }
  }
}


void FFTField::resetSpectralField(const CastemReal* const values, const double k) {
  checkSpectral("resetSpectralField");

  cfloat* F2 = F;

  // Null frequency
  for (unsigned char i = 0; i < nv; ++i, ++F2) {
    (*F2)[0] = values[i];
    (*F2)[1] = 0.;
  }

  // remaining
  for (size_t i = nv; i < nv * FSize; ++i, ++F2) {
    (*F2)[0] *= k;
    (*F2)[1] *= k;
  }
}


unsigned char FFTField::getNv() const {
  return nv;
}


size_t FFTField::getfSize() const {
  return fSize;
}


rfloat* FFTField::getSpatialField(unsigned i) {
  assertSpatial();
  if (1 != d) {
    i = i % nz + nzb * (i / nz);
  }
  return f + nv * i;
}


const rfloat* FFTField::getSpatialField(unsigned i) const {
  assertSpatial();
  if (1 != d) {
    i = i % nz + nzb * (i / nz);
  }
  return f + nv * i;
}


cfloat* FFTField::getSpectralField() {
  assertSpectral();
  return F;
}


const cfloat* FFTField::getSpectralField() const {
  assertSpectral();
  return F;
}


void FFTField::minus(const FFTField& ff1, const FFTField& ff2) {
  if (ff1.isSpatial) {
    // Spatial representation
    if (!ff2.isSpatial) throw(logic_error("FFTField::minus: ff1 Spatial and ff2 Spectral"));
    throw(logic_error("FFTField::minus: Spatial difference not yet implemanted"));
    setSpatialField();
  } else {
    // Spectral representation
    if (ff2.isSpatial) throw(logic_error("FFTField::minus: ff1 Spectral and ff2 Spatial"));
    setSpectralField();
#pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < FSize * nv; i++) {
      F[i][0] = ff1.F[i][0] - ff2.F[i][0];
      F[i][1] = ff1.F[i][1] - ff2.F[i][1];
    }
  }
}


FFTField& FFTField::operator-=(const FFTField& ff2) {
  minus(*this, ff2);
  return*this;
}


void FFTField::plus(const FFTField& ff1, const FFTField& ff2) {
  if (ff1.isSpatial) {
    // Spatial representation
    if (!ff2.isSpatial) throw(logic_error("FFTField::minus: ff1 Spatial and ff2 Spectral"));
    throw(logic_error("FFTField::minus: Spatial difference not yet implemanted"));
    setSpatialField();
  } else {
    // Spectral representation
    if (ff2.isSpatial) throw(logic_error("FFTField::minus: ff1 Spectral and ff2 Spatial"));
    setSpectralField();
#pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < FSize * nv; i++) {
      F[i][0] = ff1.F[i][0] + ff2.F[i][0];
      F[i][1] = ff1.F[i][1] + ff2.F[i][1];
    }
  }
}


FFTField& FFTField::operator+=(const FFTField& ff2) {
  plus(*this, ff2);
  return*this;
}


long double FFTField::compare(const FFTField& ff2) const {
  size_t nx, ny;
  getDim(nx, ny);
  if (isSpatial) {
    // Spatial representation
    if (!ff2.isSpatial) throw(logic_error("FFTField::Compare: ff1 Spatial and ff2 Spectral"));
    long double sum = 0;
    size_t N12 = (size_t)fSize / nz, Nzb = (size_t)nv * nzb, Nzv = (size_t)nv * nz, i, k;
    for (i = 0; i < N12; ++i) {
      const rfloat* f1 = f + i * Nzb;
      const rfloat* f2 = ff2.f + i * Nzb;
      for (k = 0; k < Nzv; ++k, ++f1, ++f2) {
        rfloat df = *f1 - *f2;
        sum += df * df;
      }
    }
    return sum / fSize;
  } else {
    // Spectral representation
    if (ff2.isSpatial) throw(logic_error("FFTField::Compare: ff1 Spectral and ff2 Spatial"));

    unsigned short ncz = nz / 2 + 1;
    unsigned short ncz2;
    if (nz % 2) {
      ncz2 = ncz;
    } else {
      ncz2 = ncz - 1;
    }

    long double sum = 0, s;
    const cfloat* F1 = F, * F2 = ff2.F;
    for (size_t i = 0; i < (size_t)nx * ny; ++i) {
      for (unsigned short k = 0; k < ncz; k++) {
        s = 0;
        for (unsigned char l = 0; l < nv; ++l, ++F1, ++F2) {
          rfloat dFr = F1[0][0] - F2[0][0];
          rfloat dFi = F1[0][1] - F2[0][1];
          s += dFr * dFr + dFi * dFi;
        }
        if ((k != ncz2) && k) {
          sum += 2 * s;
        } else {
          sum += s;
        }
      }
    }
    return sum;
  }
}


inline long double FFTField::SpatialScalarProduct(const FFTField& ff2) const {
  long double sum = 0;
  size_t N12 = (size_t)fSize / nz, Nzb = (size_t)nv * nzb, Nzv = (size_t)nv * nz;
#pragma omp parallel
  {  // begin parallel section
#pragma omp for simd schedule(static) collapse(2) reduction(+:sum)
    for (size_t i = 0; i < N12; ++i) {
      for (size_t k = 0; k < Nzv; ++k) {
        sum += f[i * Nzb + k] * ff2.f[i * Nzb + k];
      }
    }
  }
  return sum / fSize;
}


inline long double FFTField::SpectralScalarProduct(const FFTField& ff2) const {
  size_t nx, ny;
  getDim(nx, ny);
  unsigned short ncz = nz / 2 + 1;
  unsigned short ncz2;
  if (nz % 2) {
    ncz2 = ncz;
  } else {
    ncz2 = ncz - 1;
  }
  long double sum = 0;
  const cfloat* F1 = F, * F2 = ff2.F;
  size_t NT = (size_t)nx * ny, NT2 = nv * ncz;
#pragma omp parallel reduction(+: sum)
  {  // begin parallel section
#pragma omp for schedule(static) collapse(2)
    for (size_t i = 0; i < NT; ++i) {
      for (unsigned short k = 0; k < ncz; k++) {
        F1 = F + i * NT2;
        F2 = ff2.F + i * NT2;
        long double s = 0;
        unsigned id = (unsigned)k * nv;
        const cfloat* lF1 = F1 + id, * lF2 = F2 + id;
        for (unsigned char l = 0; l < nv; ++l, ++lF1, ++lF2) {
          s += lF1[0][0] * lF2[0][0] + lF1[0][1] * lF2[0][1];
        }
        if ((k != ncz2) && k) {
          sum += 2 * s;
        } else {
          sum += s;
        }
      }
    }  // end parallel section
  }
  return sum;
}


long double FFTField::scalarProduct(const FFTField& ff2) const {
  if (isSpatial) {
    // Spatial representation
    if (!ff2.isSpatial) throw(logic_error("FFTField::scalarProduct: ff1 Spatial and ff2 Spectral"));
    return SpatialScalarProduct(ff2);
  } else {
    // Spectral representation
    if (ff2.isSpatial) throw(logic_error("FFTField::scalarProduct: ff1 Spectral and ff2 Spatial"));
    return SpectralScalarProduct(ff2);
  }
}


long double FFTField::norm2() const {
  if (isSpatial) {
    // Spatial representation
    return SpatialScalarProduct(*this);
  } else {
    // Spectral representation
    return SpectralScalarProduct(*this);
  }
}


void FFTField::linComb(const FFTField* const* V, const long double* l, const unsigned char N) {
  if (V[0]->isSpatial) {
    // Spatial representation
    throw(logic_error("linComb: Not yet defined in spatial representation"));
  }
  setSpectralField();

#pragma omp parallel default(shared)
  {  // begin parallel section
#pragma omp for schedule(static)
    for (size_t i = 0; i < FSize; ++i) {
      for (size_t k = 0; k < nv; ++k) {
        cfloat sum = { 0, 0 };
        const size_t ivk = i * nv + k;
        for (size_t jj = 0; jj < N; jj++) {
          sum[0] += (V[jj]->F)[ivk][0] * l[jj];
          sum[1] += (V[jj]->F)[ivk][1] * l[jj];
        }
        F[ivk][0] = sum[0];
        F[ivk][1] = sum[1];
      }
    }
  }  // end parallel section
}

}  // namespace merope

//! Copyright : see license.txt
//!
//! \brief Grid discretisation of a medium
//!


#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Grid/Grid.hxx"
#include "Parallelism/OpenmpWrapper.hxx"
#include "VTKinout/TIFFRead.hxx"
#include "VTKinout/VTKRead.hxx"
#include "VTKinout/VTKStream.hxx"


#include "MeropeNamespace.hxx"


namespace merope {

//---------------------
Grid::Grid() {}

void Grid::initCommon() {
    // Corrections for 1 or 2D
    switch (d) {
    case 1:
        // 1-Dimension
        ly = 0.1 * lx;
        ny = nz = 1;
        lz = 0;
        break;
    case 2:
        // 2-Dimensions
        nz = 1;
        lz = 0;
    }
    ng = nx * ny * nz;
}

void Grid::oneMat() {
    phases.resize(1);
    phases[0].resize(ng);
#pragma omp parallel
    { // begin parallel section
#pragma omp for schedule(static)
        for (size_t i = 0; i < ng; ++i) {
            phases[0][i] = i;
        }
    } // end parallel section
}

Grid::Grid(const double lx_i, const size_t nx_i):
    nx(nx_i), lx(lx_i) {
    d = 1;
    initCommon();
    oneMat();
}

Grid::Grid(const double lx_i, const double ly_i, const size_t nx_i,
    const size_t ny_i):
    nx(nx_i), ny(ny_i), lx(lx_i), ly(ly_i) {
    d = 2;
    initCommon();
    oneMat();
}

Grid::Grid(const double lx_i, const double ly_i, const double lz_i,
    const size_t nx_i, const size_t ny_i, const size_t nz_i):
    nx(nx_i), ny(ny_i), nz(nz_i), lx(lx_i), ly(ly_i), lz(lz_i) {
    d = 3;
    initCommon();
    oneMat();
}


// Grid constructor by reading a VTK File Grid
// vf: VTK File Grid
Grid::Grid(VTKRead& vf) {
    VTKRead::MType mt;
    vf.getParamG(d, nx, ny, nz, lx, ly, lz, mt);
    ng = nx * ny * nz;
    switch (mt) {
    case VTKRead::UC:
        readVTK<unsigned char>(vf);
        break;
    case VTKRead::US:
        readVTK<unsigned short>(vf);
        break;
    case VTKRead::ST:
        readVTK<short>(vf);
        break;
    case VTKRead::CH:
        readVTK<char>(vf);
        break;
    }
}

Grid::Grid(TIFFRead& tf, const double dx) {
    d = 2;
    nx = tf.TIFFGetField(TIFFTAG_IMAGEWIDTH);
    ny = tf.TIFFGetField(TIFFTAG_IMAGELENGTH);
    nz = 1;
    lx = nx * dx;
    ly = ny * dx;
    lz = 0;
    ng = nx * ny;

    // Number of phases
    phases.resize(2);

    tf.seekg(tf.TIFFGetField(TIFFTAG_STRIPOFFSETS), ios_base::beg);
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            char c;
            tf.read(&c, 1);
            if (c)
                phases[1].push_back(ny - 1 - j + i * ny);
            else phases[0].push_back(ny - 1 - j + i * ny);
        }
    }
}

//----------------------------------------------
// Getters and Setters
unsigned char Grid::getDim() const {
    return d;
}
size_t Grid::getNx() const {
    return nx;
}
size_t Grid::getNy() const {
    return ny;
}
size_t Grid::getNz() const {
    return nz;
}
size_t Grid::getNg() const {
    return ng;
}

double Grid::getLx() const {
    return lx;
}
double Grid::getLy() const {
    return ly;
}
double Grid::getLz() const {
    return lz;
}

double Grid::DivNormFactor() const {
    double sum = nx / (lx * lx);
    switch (d) {
    case 3:
        sum += nz / (lz * lz);
        /* fall through */
    case 2:
        sum += ny / (ly * ly);
    }
    return sqrt(sum);
}

void Grid::vtkReorderMaterialIdx(unsigned short* const ids) const {
    NodesList::const_iterator pIt;

    // For each phase index id, save its value at each point in this phase
    unsigned short id;
    vector<NodesList>::const_iterator l;
    for (l = phases.begin(), id = 0; l != phases.end(); ++l, ++id) {
        for (pIt = l->begin(); pIt != l->end(); ++pIt) {
            ids[*pIt] = id;
        }
    }
}

void Grid::VTKheaderCELL2D(VTKstream& fvtk) const {
    const double dx = lx / nx, dy = ly / ny;
    fvtk.STRUCTURED_POINTS(nx, ny, dx, dy);

    // Data type and associated color table
    fvtk.setCELL(ng);
}

void Grid::VTKheaderCELL3D(VTKstream& fvtk) const {
    const double dx = lx / nx, dy = ly / ny, dz = lz / nz;
    fvtk.STRUCTURED_POINTS(nx, ny, nz, dx, dy, dz);

    // Data type and associated color table
    fvtk.setCELL(ng);
}

void Grid::VTKheaderCELL(VTKstream& fvtk) const {
    switch (d) {
    case 1:
    case 2:
        VTKheaderCELL2D(fvtk);
        break;
    case 3:
        VTKheaderCELL3D(fvtk);
    }
}

void Grid::toVTKCELL2D(VTKstream& fvtk) const {
    // File header
    VTKheaderCELL2D(fvtk);
    fvtk << "SCALARS MaterialId unsigned_short" << endl;
    fvtk << "LOOKUP_TABLE default" << endl;

    // Geometry reconstruction
    vector<unsigned short> ids(ng);
    vtkReorderMaterialIdx(&ids[0]);

    // Write phase indices values in file
    size_t i, j;
    for (j = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            fvtk.write(ids[j + i * ny]);
        }
    }
}

void Grid::toVTKCELL3D(VTKstream& fvtk) const {
    // File header
    VTKheaderCELL3D(fvtk);
    fvtk << "SCALARS MaterialId unsigned_short" << endl;
    fvtk << "LOOKUP_TABLE default" << endl;

    // Geometry reconstruction
    vector<unsigned short> ids(ng);
    vtkReorderMaterialIdx(&ids[0]);

    // Write phase indices values in file
    size_t i, j, k;
    for (k = 0; k < nz; ++k) {
        for (j = 0; j < ny; ++j) {
            for (i = 0; i < nx; ++i) {
                fvtk.write(ids[k + nz * (j + i * ny)]);
            }
        }
    }
}

void Grid::toVTKCELL(VTKstream& fvtk) const {
    switch (d) {
    case 1:
    case 2:
        toVTKCELL2D(fvtk);
        break;
    case 3:
        toVTKCELL3D(fvtk);
    }
}

unsigned short Grid::getNbOfPhases() const {
    return (unsigned short)phases.size();
}

const std::vector<NodesList>& Grid::getPhases() const {
    return phases;
}

const NodesList* Grid::getPhase(const unsigned short alpha) const {
    if (alpha >= phases.size()) {
        stringstream msg;
        msg << "Grid::getPhase(): invalid phase number (" << alpha << ">"
            << phases.size() << ")";
        throw out_of_range(msg.str());
    }
    return &phases[alpha];
}

double Grid::phaseFraction(unsigned long i) const {
    unsigned long n = phases.size();
    if (i > n - 1)
        throw runtime_error(
            "Grid::phaseFraction : phase index must be less than "
            + to_string(n));
    return (double)phases[i].size() / ng;
}

void Grid::print(std::ostream& os) const {
    os << "====== Grid info: " << (int)d << "D [" << nx << "x" << ny << "x"
        << nz << "]\n";
    os << "Number of nodes  : " << ng << endl;
    os << "Number of phases : " << phases.size() << endl;
    vector<NodesList>::const_iterator i;
    unsigned short i2 = 0;
    for (i = phases.begin(); i != phases.end(); ++i, ++i2) {
        os << "  Phase[" << i2 << "]: " << i->size() << " nodes = "
            << (double)(i->size()) / (0.01 * ng) << " %" << endl;
    }
}

unsigned Grid::Go2VTKo(const size_t Ig) const {
    unsigned k = Ig % nz;
    unsigned I2 = Ig / nz;
    unsigned j = I2 % ny;
    unsigned i = I2 / ny;
    return i + nx * (j + k * ny);
}

void Grid::covX(const unsigned short a, const unsigned short b,
    ostream& fout) const {
    vector<double> var(nx / 2 + 1, 0); // Covariance function

    // Geometry reconstruction
    vector<unsigned short> ids(ng);
    vtkReorderMaterialIdx(&ids[0]);

    vector<double>::iterator vi;
    size_t nyz = nz * ny;
    size_t h;
    for (vi = var.begin(), h = 0; vi != var.end(); ++vi, ++h) {
        *vi = 0;
        for (size_t k = 0; k < nz; k++) {
            size_t cj = k;
            for (size_t j = 0; j < ny; ++j, cj += nz) {
                double sum = 0;
#pragma omp parallel
                { // begin parallel section
#pragma omp for schedule(static) reduction(+: sum)
                    for (size_t i = 0; i < nx; ++i) {
                        size_t i2 = i + h;
                        if (i2 >= nx) i2 -= nx;
                        sum += (a == ids[cj + nyz * i])
                            * (b == ids[cj + nyz * i2]);
                    }
                } // end parallel section
                *vi += sum;
            }
        }
    }

    // Phase proportions
    double fa = ((double)phases[a].size()) / ng;
    double fb = ((double)phases[b].size()) / ng;
    // Normalisation and output
    double dx = lx / nx, x;
    for (vi = var.begin(), x = 0; vi != var.end(); ++vi, x += dx) {
        fout << x << " " << *vi / ng - fa * fb << endl;
    }
}

void Grid::covY(const unsigned short a, const unsigned short b,
    ostream& fout) const {
    vector<double> var(ny / 2 + 1, 0); // Covariance function

    // Geometry reconstruction
    vector<unsigned short> ids(ng);
    vtkReorderMaterialIdx(&ids[0]);

    vector<double>::iterator vi;
    size_t nyz = nz * ny;
    size_t h;
    for (vi = var.begin(), h = 0; vi != var.end(); ++vi, ++h) {
        *vi = 0;
        for (size_t k = 0; k < nz; ++k) {
            size_t cj1 = k;
            size_t cj2 = k + h * nz;
            for (size_t j = 0, j2 = h; j < ny;
                ++j, ++j2, cj1 += nz, cj2 += nz) {
                if (j2 >= ny) {
                    j2 = 0;
                    cj2 = k;
                }
                double sum = 0;
#pragma omp parallel
                { // begin parallel section
#pragma omp for schedule(static) reduction(+: sum)
                    for (size_t i = 0; i < nx; ++i) {
                        sum += (a == ids[cj1 + i * nyz])
                            * (b == ids[cj2 + i * nyz]);
                    } // end parallel section
                }
                *vi += sum;
            }
        }
    }

    // Phase proportions
    double fa = ((double)phases[a].size()) / ng;
    double fb = ((double)phases[b].size()) / ng;
    // Normalisation and output
    double dy = ly / ny, y;
    for (vi = var.begin(), y = 0; vi != var.end(); ++vi, y += dy) {
        fout << y << " " << *vi / ng - fa * fb << endl;
    }
}

void Grid::covZ(const unsigned short a, const unsigned short b,
    ostream& fout) const {
    vector<double> var(nz / 2 + 1, 0); // Covariance function

    // Geometry reconstruction
    vector<unsigned short> ids(ng);
    vtkReorderMaterialIdx(&ids[0]);

    vector<double>::iterator vi;
    size_t nyz = nz * ny;
    size_t h;
    for (vi = var.begin(), h = 0; vi != var.end(); ++vi, ++h) {
        *vi = 0;
        for (size_t k = 0, k2 = h; k < nz; ++k, ++k2) {
            if (k2 >= nz) k2 = 0;
            size_t cj1 = k;
            size_t cj2 = k2;
            for (size_t j = 0; j < ny; ++j, cj1 += nz, cj2 += nz) {
                double sum = 0;
#pragma omp parallel
                { // begin parallel section
#pragma omp for schedule(static) reduction(+: sum)
                    for (size_t i = 0; i < nx; ++i) {
                        sum += (a == ids[cj1 + i * nyz])
                            * (b == ids[cj2 + i * nyz]);
                    }
                } // end parallel section
                *vi += sum;
            }
        }
    }

    // Phase proportions
    double fa = ((double)phases[a].size()) / ng;
    double fb = ((double)phases[b].size()) / ng;
    // Normalisation and output
    double dz = lz / nz, z;
    for (vi = var.begin(), z = 0; vi != var.end(); ++vi, z += dz) {
        fout << z << " " << *vi / ng - fa * fb << endl;
    }
}

void Grid::varioX(const CastemReal* const v, ostream& fout) const {
    size_t nyz = nz * ny;
    double dx = lx / nx;
    for (size_t h = 1; h <= nx / 2; ++h) {
        double sum = 0;
        size_t NT = 0;
#pragma omp parallel default(shared) reduction(+:sum,NT)
        {
#pragma omp for collapse(2) schedule(static) 
            for (size_t jk = 0; jk < nyz; jk++) {
                for (size_t i = 0; i < nx; ++i) {
                    size_t i2 = i + h;
                    if (i2 >= nx) i2 -= nx;
                    const CastemReal& v1 = v[jk + nyz * i];
                    const CastemReal& v2 = v[jk + nyz * i2];
                    if (isnormal(v1) and isnormal(v2)) {
                        double dv = v1 - v2;
                        sum += dv * dv;
                        ++NT;
                    }
                }
            }
        }
        fout << h * dx << " " << sum / (2 * NT) << endl;
    }
}

void Grid::varioY(const CastemReal* const v, ostream& fout) const {
    if (d < 2) throw(logic_error("Grid::varioY: No Y variogramm in 1D"));

    double dy = ly / ny;
    for (size_t h = 1; h <= ny / 2; ++h) {
        double sum = 0;
        size_t NT = 0;
#pragma omp parallel
        { // begin parallel section
#pragma omp for schedule(static) reduction(+: sum,NT)
            for (size_t k = 0; k < nz; ++k) {
                for (size_t j = 0, j2 = h; j < ny; ++j, ++j2) {
                    if (j2 >= ny) {
                        j2 = 0;
                    }
                    for (size_t i = 0; i < nx; ++i) {
                        const CastemReal& v1 = v[k + (j + i * ny) * nz];
                        const CastemReal& v2 = v[k + (j2 + i * ny) * nz];
                        if (isnormal(v1) and isnormal(v2)) {
                            double dv = v1 - v2;
                            sum += dv * dv;
                            ++NT;
                        }
                    }
                }
            }
        } // end parallel section
        fout << h * dy << " " << sum / (2 * NT) << endl;
    }
}

void Grid::varioZ(const CastemReal* const v, ostream& fout) const {
    if (d < 3) throw(logic_error("Grid::varioZ: No Z variogramm in 1 or 2D"));
    size_t nxy = nx * ny;
    double dz = lz / nz;
    for (size_t h = 1; h <= nz / 2; ++h) {
        double sum = 0;
        size_t NT = 0;
#pragma omp parallel
        { // begin parallel section
#pragma omp for schedule(static) reduction(+: sum,NT)
            for (size_t k = 0; k < nz; ++k) {
                size_t k2 = k + h;
                if (k2 >= nz) k2 -= nz;
                for (size_t ij = 0; ij < nxy; ++ij) {
                    const CastemReal& v1 = v[k + ij * nz];
                    const CastemReal& v2 = v[k2 + ij * nz];
                    if (isnormal(v1) and isnormal(v2)) {
                        double dv = v1 - v2;
                        sum += dv * dv;
                        ++NT;
                    }
                }
            }
        } // end parallel section
        fout << h * dz << " " << sum / (2 * NT) << endl;
    }
}

} // namespace merope

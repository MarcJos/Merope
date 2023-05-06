//! Copyright : see license.txt
//!
//! \brief 
//!


#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Grid/Grid_VER.hxx"
#include "Voxellation/Voxellation.hxx"


#include "MeropeNamespace.hxx"


namespace merope {

Grid_VER::Grid_VER():
    Grid() {}

Grid_VER::Grid_VER(VTKRead& geometry):
    Grid(geometry) {
    initCommon();
}

Grid_VER::Grid_VER(TIFFRead& tf, double dx_):
    Grid(tf, dx_) {
    initCommon();
}

Grid_VER::Grid_VER(const vox::Voxellation<3>& voxellisation) {
    this->phases = voxellisation.getGrid().phases;
    this->set_L<3>(voxellisation.getL());
    this->set_Nb<3>(voxellisation.getGrid().nbNodes);
    initCommon();
}

size_t Grid_VER::get_linear_index(size_t i, size_t j, size_t k) const {
    return k + (j + i * ny) * nz;
}

void Grid_VER::set_NbPhase(size_t s) {
    phases.resize(s);
}

void Grid_VER::set_Voxel(const array<size_t, 3>& ijk, size_t indexPhase) {
    set_Voxel(get_linear_index<3>(ijk), indexPhase);
}

void Grid_VER::set_Voxel(size_t i, size_t j, size_t k, size_t indexPhase) {
    set_Voxel(get_linear_index(i, j, k), indexPhase);
}

void Grid_VER::set_Voxel(size_t linear_index, size_t indexPhase) {
    assert(indexPhase < phases.size());
    phases[indexPhase].push_back(linear_index);
}

void Grid_VER::initCommon() {
    Grid::initCommon();
    L = vector<double>{ lx, ly, lz };
    nbNodes = array<size_t, 3> { nx, ny, nz };
    dx = array<double, 3> { 0., 0., 0. };
    inverse_dx = array<double, 3> { 0., 0., 0. };
    for (size_t i = 0; i < d; i++) {
        dx[i] = L[i] / nbNodes[i];
        inverse_dx[i] = 1. / dx[i];
    }
}

size_t vox::aux::symmetriseAuxi(size_t N1, size_t i, size_t N2) {
    assert(N1 < N2);
    if (i % 2 == 0) return N1 + i * N2;
    else return (i + 1) * N2 - N1 - 1;
}

void Grid_VER::removeUnusedPhase(vector<double>& coefficients) {
    if (coefficients.size() < this->phases.size()) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Problem");
    }
    ///
    size_t currentPhase = 0;
    vector<double> newCoeff{};
    for (size_t i = 0; i < this->phases.size(); i++) {
        if (this->phases[i].size() > 0) {
            newCoeff.push_back(coefficients[i]);
            if (i > currentPhase) {
                this->phases[currentPhase] = std::move(this->phases[i]);
            }
            currentPhase++;
        }
    }
    coefficients = newCoeff;
    this->phases.resize(currentPhase);
}

} // namespace merope

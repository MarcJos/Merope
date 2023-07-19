//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef FFTW3_VOXELLATIONGREEN_IXX_
#define FFTW3_VOXELLATIONGREEN_IXX_

#include "../Grid/AreCompatible.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<unsigned short DIM>
inline vox::VoxSimpleGauss<DIM>::VoxSimpleGauss(const CartesianField<DIM>& cartesianField_, const PreSubGrid<DIM>& gridParameters_) :
    VoxGrid<DIM, double>(gridParameters_),
    fieldGenerator{ &cartesianField_ } {
}

template<unsigned short DIM>
inline void vox::VoxSimpleGauss<DIM>::build() {
    for (size_t i = 0; i < DIM; i++) {
        if (this->getVoxGrid().getNMin()[i] != 0 or this->getVoxGrid().getNMax()[i] != this->getVoxGrid().getNbNodeBigGrid()[i]) {
            auxi_function::writeVectorToString(this->getVoxGrid().getNbNodeBigGrid(), cerr); cerr << endl;
            auxi_function::writeVectorToString(this->getVoxGrid().getNMin(), cerr); cerr << endl;
            auxi_function::writeVectorToString(this->getVoxGrid().getNMax(), cerr); cerr << endl;
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Unexpected");
        }
    }
    ////////////////////////////////
    if (fieldGenerator->getTypeField() == TypeField::Gaussian) { // case : gaussian
        Grid_VER grid{};
        grid.set_L<DIM>(this->getGridParameters().getL());
        grid.set_Nb<DIM>(this->getGridParameters().getNbNodeBigGrid());
        auto gaussianField_ = fieldGenerator->getGaussianField();
        auto randomField = gaussianField::createField(gaussianField_, grid, this->getVoxGrid().size());
        std::function<double(double)> nonLin = gaussianField_.nonlinearFunction;
        { // begin parallel section
        //#pragma omp parallel for firstprivate(nonLin) //inefficient when using python function
            for (size_t i = 0; i < this->getVoxGrid().size(); ++i) {
                this->getVoxGrid()[i] = nonLin(randomField[i]);
            }
        } // end parallel section
    } else if (fieldGenerator->getTypeField() == TypeField::NumericalCovariance) { // case : gaussian
        Grid_VER grid{};
        grid.set_L<DIM>(this->getGridParameters().getL());
        grid.set_Nb<DIM>(this->getGridParameters().getNbNodeBigGrid());
        auto covariance = fieldGenerator->getCovariance();
        auto covarianceField = gaussianField::createField(covariance, grid, this->getVoxGrid().size());
        { // begin parallel section
        //#pragma omp parallel for firstprivate(nonLin) //inefficient when using python function
            for (size_t i = 0; i < this->getVoxGrid().size(); ++i) {
                this->getVoxGrid()[i] = covarianceField[i];
            }
        } // end parallel section
    } else if (fieldGenerator->getTypeField() == TypeField::Scalar) { // case : scalar field
        auto nbVoxels = this->getVoxGrid().getNbNodeBigGrid();
        auto fonction = fieldGenerator->getScalarField().fieldFunction;
        auto dx = this->getVoxGrid().getDx();
        loop(nbVoxels, [this, &nbVoxels, &fonction, &dx](auto ijk) {
            this->getVoxGrid()[vox::auxi::get_linear_index<DIM>(ijk, nbVoxels)] = fonction(vox::auxi::origin<DIM>(ijk, dx));
            ;
            });
    } else if (fieldGenerator->getTypeField() == TypeField::Discretized) { // case discretized
        auto discretizedField = fieldGenerator->getDiscretizedField();
        if (areCompatible(discretizedField.getGridParameters(), this->getGridParameters())) {
            this->getVoxGrid() = discretizedField;
        }
    } else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw invalid_argument("TypeField");
    }
}

} //namespace vox
} // namespace merope


#endif /* FFTW3_VOXELLATIONGREEN_IXX_ */

//! Copyright : see license.txt
//!
//! \brief

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Field/GaussianField.hxx"
#include "FFTW3/FFTScalarField.hxx"


#include "MeropeNamespace.hxx"


namespace merope {

vector<double> gaussianField::createField(const NumericalCovariance<3>& covariance,
    const Grid& grid, size_t size_) {
    return createField_T<3, FFTScalarField>(covariance, grid, size_);
}

vector<double> gaussianField::createField(const NumericalCovariance<2>& covariance,
    const Grid& grid, size_t size_) {
    return createField_T<2, FFTScalarField>(covariance, grid, size_);
}

vector<double> gaussianField::createField(const SimpleGaussianField<3>& gaussianField,
    const Grid& grid, size_t size_) {
    return createField_T<3, FFTScalarField>(gaussianField, grid, size_);
}

vector<double> gaussianField::createField(const SimpleGaussianField<2>& gaussianField,
    const Grid& grid, size_t size_) {
    return createField_T<2, FFTScalarField>(gaussianField, grid, size_);
}

} // namespace merope

//! Copyright : see license.txt
//!
//! \brief
//
#pragma once



#include "../MeropeNamespace.hxx"


namespace merope {

namespace gaussianField {

template<unsigned short DIM>
template<class GaussianRV>
inline void SimpleGaussianField<DIM>::setNonlinearTransform(const GaussianRV& nonlinearTransform_) {
    auto local_copy = nonlinearTransform_;
    nonlinearFunction = [local_copy](double x) {return local_copy.inverseA(x);};
}

template<unsigned short DIM, class FFTScal>
vector<double> createField_T(const SimpleGaussianField<DIM>& gaussianField,
    const Grid& grid, size_t size_) {
    FFTScal Ff(grid, gaussianField.covariance, true, gaussianField.seed);
    vector<double> result(size_);
    for (size_t i = 0; i < size_; i++) {
        result[i] = *(Ff.getSpatialField(i));
    }
    return result;
}

template<unsigned short DIM, class FFTScal>
vector<double> createField_T(const NumericalCovariance<DIM>& covariance,
    const Grid& grid, size_t size_) {
    FFTScal Ff(grid, covariance.covariance, true, 0, true);
    vector<double> result(size_);
    for (size_t i = 0; i < size_; i++) {
        result[i] = *(Ff.getSpatialField(i));
    }
    return result;
}

}  // namespace gaussianField
}  // namespace merope



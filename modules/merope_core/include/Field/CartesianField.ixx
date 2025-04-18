//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace merope {

template<unsigned short DIM>
inline CartesianField<DIM>::CartesianField(
    const gaussianField::SimpleGaussianField<DIM>& gaussianField,
    const Point<DIM>& L) : CartesianField<DIM>(L) {
    set(gaussianField);
}


template<unsigned short DIM>
inline CartesianField<DIM>::CartesianField(const gaussianField::NumericalCovariance<DIM>& covariance,
    const Point<DIM>& L) : CartesianField<DIM>(L) {
    set(covariance);
}


template<unsigned short DIM>
inline CartesianField<DIM>::CartesianField(
    const realScalarField::Field<DIM>& scalarField, const Point<DIM>& L) : CartesianField<DIM>(L) {
    set(scalarField);
}

template<unsigned short DIM>
inline CartesianField<DIM>::CartesianField(const vox::GridField<DIM>& gridField,
    const Point<DIM>& L) : CartesianField<DIM>(L) {
    if (geomTools::distanceCarre<DIM>(gridField.getL(), L) > geomTools::EPSILON) {
        throw invalid_argument("Incompatible lengths");
    }
    set(gridField);
}

template<unsigned short DIM>
void CartesianField<DIM>::set(const gaussianField::NumericalCovariance<DIM>& covariance) {
    this->typeField = TypeField::NumericalCovariance;
    this->localField = covariance;
}

template<unsigned short DIM>
void CartesianField<DIM>::set(const gaussianField::SimpleGaussianField<DIM>& gaussianField) {
    this->typeField = TypeField::Gaussian;
    this->localField = gaussianField;
}

template<unsigned short DIM>
void CartesianField<DIM>::set(const realScalarField::Field<DIM>& scalarField) {
    this->typeField = TypeField::Scalar;
    this->localField = scalarField;
}

template<unsigned short DIM>
void CartesianField<DIM>::set(const vox::GridField<DIM>& discretizeField) {
    this->typeField = TypeField::Discretized;
    this->localField = discretizeField;
}

template<unsigned short DIM>
inline const gaussianField::SimpleGaussianField<DIM>& CartesianField<DIM>::getGaussianField() const {
    return std::get<gaussianField::SimpleGaussianField<DIM>>(localField);
}

template<unsigned short DIM>
inline const gaussianField::NumericalCovariance<DIM>& CartesianField<DIM>::getCovariance() const {
    return std::get<gaussianField::NumericalCovariance<DIM>>(localField);
}

template<unsigned short DIM>
inline const realScalarField::Field<DIM>& CartesianField<DIM>::getScalarField() const {
    return std::get<realScalarField::Field<DIM>>(localField);
}

template<unsigned short DIM>
inline const vox::GridField<DIM>& CartesianField<DIM>::getDiscretizedField() const {
    return std::get<vox::GridField<DIM>>(localField);
}

}  // namespace merope



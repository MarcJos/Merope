//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MEROPE_CORE_SRC_FIELD_CARTESIANFIELD_HXX_
#define MEROPE_CORE_SRC_FIELD_CARTESIANFIELD_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"


#include "../Field/GaussianField.hxx"
#include "../Field/RealScalarField.hxx"
#include "../Grid/PreGrid.hxx"
#include "../Grid/GridTypes.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {

enum class TypeField {
    Gaussian,
    NumericalCovariance,
    Scalar,
    Discretized,
    Void
};

template<unsigned short DIM>
class CartesianField : public InsideTorus<DIM> {
public:
    //! constructors
    CartesianField(const gaussianField::SimpleGaussianField<DIM>& gaussianField, const Point<DIM>& L);
    CartesianField(const gaussianField::NumericalCovariance<DIM>& covariance, const Point<DIM>& L);
    CartesianField(const realScalarField::Field<DIM>& scalarField, const Point<DIM>& L);
    CartesianField(const vox::GridField<DIM>& gridField, const Point<DIM>& L);
    //! getter
    const gaussianField::SimpleGaussianField<DIM>& getGaussianField() const;
    const gaussianField::NumericalCovariance<DIM>& getCovariance() const;
    const realScalarField::Field<DIM>& getScalarField() const;
    const vox::GridField<DIM>& getDiscretizedField() const;
    TypeField getTypeField() const { return typeField; }

protected:
    //! constructors
    CartesianField(const Point<DIM>& L) : InsideTorus<DIM>(L), localField{ nullptr }, typeField{ TypeField::Void } {}
    //! setter
    void set(const gaussianField::SimpleGaussianField<DIM>& gaussianField);
    void set(const gaussianField::NumericalCovariance<DIM>& covariance);
    void set(const realScalarField::Field<DIM>& scalarField);
    void set(const vox::GridField<DIM>& trueField);

private:
    //! holds either a gaussian field or a scalarField
    //! default is void*
    std::variant<gaussianField::SimpleGaussianField<DIM>,
        realScalarField::Field<DIM>,
        vox::GridField<DIM>,
        gaussianField::NumericalCovariance<DIM>,
        void*> localField;
    //! stores the type of field
    TypeField typeField;
};

} // namespace merope

#include "../Field/CartesianField.ixx"

#endif /* MEROPE_CORE_SRC_FIELD_CARTESIANFIELD_HXX_ */

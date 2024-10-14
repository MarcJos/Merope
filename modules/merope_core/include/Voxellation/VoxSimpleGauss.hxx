//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Field/CartesianField.hxx"
#include "../Grid/Grid_VER.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<unsigned short DIM>
class VoxSimpleGauss : public VoxGrid<DIM, double> {
public:
    //! constructor
    VoxSimpleGauss(const CartesianField<DIM>& cartesianField_, const GridParameters<DIM>& gridParameters);
protected:
    //! fills the vector CartesianGrid<DIM, double>
    void build() override;
    //! inner representation
    const CartesianField<DIM>* const fieldGenerator;
};

} /* namespace vox */
}  // namespace merope


#include "VoxSimpleGauss.ixx"



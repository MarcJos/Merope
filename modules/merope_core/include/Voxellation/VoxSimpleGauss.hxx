//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef VOXELLATIONGREEN_HXX_
#define VOXELLATIONGREEN_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Field/CartesianField.hxx"
#include "../Grid/Grid_VER.hxx"
// #include "../Voxellation/VoxRecurStructure.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<unsigned short DIM>
class VoxSimpleGauss: public VoxGrid<DIM, double> {
public:
    //! constructor
    VoxSimpleGauss(const CartesianField<DIM>& cartesianField_, const PreSubGrid<DIM>& gridParameters);
protected:
    //! fills the vector CartesianGrid<DIM, double>
    void build() override;
    //! inner representation
    const CartesianField<DIM>* const fieldGenerator;
};

} /* namespace vox */
} // namespace merope


#include "VoxSimpleGauss.ixx"

#endif /* VOXELLATIONGREEN_HXX_ */

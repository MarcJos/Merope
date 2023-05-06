//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef PREVOXELLATION_IXX_
#define PREVOXELLATION_IXX_


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

// PreVoxellation<DIM>
template<unsigned short DIM>
inline void PreVoxellation<DIM>::setGridNL(array<size_t, DIM> nbNodes_,
    array<double, DIM> L_) {
    this->set_nbNodes_L(nbNodes_, L_);
    grid.set_L<DIM>(L_);
    grid.set_Nb<DIM>(nbNodes_);
}

} // namespace vox
} // namespace merope


#endif /* PREVOXELLATION_IXX_ */

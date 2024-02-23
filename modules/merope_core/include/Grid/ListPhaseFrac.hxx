//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GRID_LISTPHASEFRAC_HXX_
#define GRID_LISTPHASEFRAC_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/SphereContainer.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {
namespace gridAuxi {

template<class TYPE_PHASE>
class ListPhaseFrac : public vector<auxi_SphereCollection::PhaseFrac<TYPE_PHASE>> {
    //! Class used to represent on a given point the different volume fractions of each phase
    //! Ideally, these volume fractions sum up to 100%
public:
    //! constructor
    ListPhaseFrac() : vector<auxi_SphereCollection::PhaseFrac<TYPE_PHASE>>({}) {}
    //! constructor
    ListPhaseFrac(initializer_list<auxi_SphereCollection::PhaseFrac<TYPE_PHASE>> lst) :vector<auxi_SphereCollection::PhaseFrac<TYPE_PHASE>>(lst) {}
    //! adds another phasefrac, but maintaining minimal memory requirements
    void add(const auxi_SphereCollection::PhaseFrac<TYPE_PHASE>& pf) { this->push_back(pf); }
    //! merge phases that are present twice (adding the volume fractions)
    void merge(double merge_criterion = 1e-6);
    //! renormalizes the volume fractions to 100%, eventually add the matrixPhase
    //! \param matrixPhase : phase to be eventually added to get 1 volume fraction
    //! \param is_there_matrix : if true, there is a matrix
    void renormalize(bool is_there_matrix = true, TYPE_PHASE matrixPhase = 0);
private:
    //! @brief renormalizes the volume fractions to 100%
    //! by dividing all volume fractions by the given inverseTotalVolumeFraction
    void renormalize_by_multiply(double inverseTotalVolumeFraction);
};

} // namespace gridAuxi
} // namespace vox
} // namespace merope

#include "../Grid/ListPhaseFrac.ixx"

#endif /* GRID_LISTPHASEFRAC_HXX_ */

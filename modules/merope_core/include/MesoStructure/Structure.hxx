//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/GeomTools.hxx"

#include "../Field/GaussianField.hxx"
#include "../Field/RealScalarField.hxx"
#include "../MultiInclusions/MultiInclusions.hxx"
#include "RecurMesoStructure.hxx"


namespace merope {

template<unsigned short DIM>
//! combines many MultiInclusions
//! according to rules auxiMicroStructure::TypeOfCombination
class Structure final : public RecurGeom<DIM> {
public:
    //! initializes a structure with a multiInclusion
    //! \param mi : multiInclusions
    Structure(const MultiInclusions<DIM>& mi) { this->setMultiInclusions(mi); }
    //! Constructors : Combination2
    //! \param s1, s2 : structures to be combined
    //! \param map_ : if phase of s2 == 1, then turn the phase of s1 into another phase through the map
    Structure(const Structure<DIM>& s1, const Structure<DIM>& s2, map<PhaseType, PhaseType> map_) { this->setCombination2(s1, s2, map_); }
    //! Constructors : Combination2
    //! \param s1, s2 : structures to be combined
    //! \param func : general function that combines phases of s1 and s2 to produce a new phase s3
    //! \param newPhases : collects all the novel phases produced by the function
    Structure(const Structure<DIM>& s1, const Structure<DIM>& s2, std::function<PhaseType(PhaseType, PhaseType)> func, vector<PhaseType> newPhases) { this->setCombination2(s1, s2, func, newPhases); }
    //! Constructor : Mask
    //! \param s1, s2 : structures to be combined
    //! \param mask : if phase of mask == 1, replace phase of s1 by phase of s2
    Structure(const Structure<DIM>& s1, const Structure<DIM>& s2, const Structure<DIM>& mask_) { this->setMask(s1, s2, mask_); }
    //! Structures = RecurGeoms in terms of internal data
    Structure(const RecurGeom<DIM>& preMicro) : RecurGeom<DIM>(preMicro) {}
};

template<unsigned short DIM>
//! combine many CartesianField
//! according to rules auxiMicroStructure::TypeOfCombination
class FieldStructure : public RecurField<DIM> {
public:
    //! initializes a structure with a multiInclusion
    //! \param cartesianField :cartesianField
    FieldStructure(const CartesianField<DIM>& cartesianField) { this->setBasicStruct(cartesianField); }
    //! Constructors : Combination2
    //! \param s1, s2 : structures to be combined
    //! \param func : general function that combines phases of s1 and s2 to produce a new phase s3
    FieldStructure(const FieldStructure<DIM>& s1, const FieldStructure<DIM>& s2, std::function<double(double, double)> func) { this->setCombination2(s1, s2, func); }
    //! Constructor : Mask
    //! \param s1, s2 : structures to be combined
    //! \param mask : if phase of mask > 0, replace phase of s1 by phase of s2
    //! \warning : should be the same with previous mask
    FieldStructure(const FieldStructure<DIM>& s1, const FieldStructure<DIM>& s2, const FieldStructure<DIM>& mask_) { this->setMask(s1, s2, mask_); }
    //! Structures = RecurGeoms in terms of internal data
    FieldStructure(const RecurField<DIM>& recurStructure) : RecurField<DIM>(recurStructure) {}
};

template<class STRUCTURE>
using Phase_Type_From_Structure_Type = std::conditional_t<std::is_same_v<STRUCTURE, Structure<2>> or std::is_same_v<STRUCTURE, Structure<3>>
    or is_same_v<STRUCTURE, MultiInclusions<2>> or is_same_v<STRUCTURE, MultiInclusions<3>>, PhaseType, double>;


}  // namespace merope



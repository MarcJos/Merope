//! Copyright : see license.txt
//!
//! \brief
//!

#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"
#include "../../../Geometry/include/AspRatio.hxx"

#include "../../../AlgoPacking/include/AlgoNames.hxx"

#include "../MesoStructure/InsideTorus.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"
#include "ObjectInclusions.hxx"


namespace merope {

template<unsigned short DIM>
class LaguerreTess final : public  PolyInclusions<DIM>, public WithAspratio<DIM> {
    // class implementing Laguerre tessellation, with possibility to impose aspect ratio
public:
    //! \param L_ : the torus dimensions
    //! \param seeds_  : a list of seeds parametrizing the tessellation (=centers of tessels + weight + phase Id)
    LaguerreTess(array<double, DIM> L_, vector<Sphere<DIM>> seeds_);
    //! compute the tessels appealing to voro++
    void computeTessels();
    //! @brief : getter
    const vector<Sphere<DIM>>& getSeeds() const { return seeds; }
    //! @return only the polyInclusions inside it, forgetting the paving nature of Laguerre tessellations
    PolyInclusions<DIM> toPolyInclusions() { this->computeTessels(); return static_cast<PolyInclusions<DIM>>(*this); };

protected:
    //! seeds giving rise to the tessels
    const vector<Sphere<DIM>> seeds;

private:
    //! impossible to add inclusions to LaguerreTess
    using PolyInclusions<DIM>::addInclusion;
    //! for aspect ratio = {1,1,1}
    void no_aspratio_computeTessels();
    //! for aspect ratio != {1,1,1}
    void with_aspratio_computeTessels();
};

}  // namespace merope


#include "../MultiInclusions/LaguerreTess.ixx"



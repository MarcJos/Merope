//! Copyright : see license.txt
//!
//! \brief
//!

#pragma once

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AlgoNames.hxx"
#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../Geometry/AspRatio.hxx"
#include "../MesoStructure/InsideTorus.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"
#include "PolyInclusions.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM>
class LaguerreTess final : public  PolyInclusions<DIM> {
    // class implementing Laguerre tessellation, with possibility to impose aspect ratio
public:
    //! \param L_ : the torus dimensions
    //! \param seeds_  : a list of seeds parametrizing the tessellation (=centers of tessels + weight + phase Id)
    LaguerreTess(array<double, DIM> L_, vector<Sphere<DIM>> seeds_);
    //! compute the tessels appealing to voro++
    void computeTessels();

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



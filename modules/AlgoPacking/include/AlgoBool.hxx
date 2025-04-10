//! Copyright : see license.txt
//!
//! \briefImplements random boolean sphere placement
//
#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

#include "../../Geometry/include/AmbiantSpace.hxx"
#include "../../Geometry/include/GlobalShape.hxx"

#include "RadiusGenerator.hxx"
#include "AlgoPacking.hxx"
#include "AlgoNames.hxx"

namespace sac_de_billes {

template<unsigned short DIM>
class AlgoBool_Template {
public:
    AlgoBool_Template(AmbiantSpace::BigShape<DIM>* T_,
        algoRSA_aux::RadiusGenerator<DIM>* radiusGen_, unsigned seed_);
    map<string, string> proceed(int method);
    //! vector containing the Spheres, refers to the motherGrid
    vector<Sphere<DIM>>* placedSpheres;
    //! Ambiant space
    AmbiantSpace::BigShape<DIM>* bigShape;
    //! destructor
    ~AlgoBool_Template() {
        delete placedSpheres;
    }
    //! \return the type of the algorithm
    algoSpheres::TypeAlgo getTypeAlgo() const { return algoSpheres::TypeAlgo::BOOL; }

private:
    //! radius generator
    algoRSA_aux::RadiusGenerator<DIM>* radiusGen;
    //! pick a point inside the torus
    Point<DIM> pickPoint();
    //! random generator
    mt19937 randGenerator;
};

template<unsigned short DIM>
using AlgoBool = AlgoInterface<DIM, AlgoBool_Template<DIM>>;

}  // namespace sac_de_billes

#include "AlgoBool.ixx"


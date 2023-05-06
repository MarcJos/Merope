//! Copyright : see license.txt
//!
//! \brief 
//!

#include "AlgoPacking.hxx"

#include "AlgoRSA.hxx"
#include "AmbiantSpace.hxx"

namespace sac_de_billes {

using namespace algoRSA_aux;
using namespace AmbiantSpace;

AlgoRSA2D::AlgoRSA2D(Point<2> L, vector<array<double, 2> > desiredRPhi_,
    double exclusionDistance_, unsigned seed, unsigned short method,
    string nameShape_) {
    setExclusionDistance(exclusionDistance_);
    setBigShape(vector<double> { L[0], L[1] }, nameShape_);
    setRadiusGenerator(desiredRPhi_);
    proceed(seed, method);
}

AlgoRSA3D::AlgoRSA3D(Point<3> L, vector<array<double, 2> > desiredRPhi_,
    double exclusionDistance_, unsigned seed, unsigned short method,
    string nameShape_): AlgoRSA3D(L, desiredRPhi_, exclusionDistance_, seed, method, AmbiantSpace::readShape(nameShape_)) {}


AlgoRSA3D::AlgoRSA3D(Point<3> L, vector<array<double, 2> > desiredRPhi_,
    double exclusionDistance_, unsigned seed, unsigned short method,
    AmbiantSpace::NameShape nameShape_) {
    setExclusionDistance(exclusionDistance_);
    setBigShape(vector<double> { L[0], L[1], L[2] }, nameShape_);
    setRadiusGenerator(desiredRPhi_);
    proceed(seed, method);
}


double fractionVolThMaxRSA(short unsigned d) {
    switch (d) {
    case 1:
        return 0.74750; // Reference : Remyi, Publ. Math. Inst. Hung. Acad. Sci. 3, 109 (1958)
    case 2:
        return 0.54689; // Reference : Torquato, Physical Review E 74, 061308 2006
    case 3:
        return 0.385; // Reference : Cooper, Physical Review A Vol38(1) 1988
    default:
        throw runtime_error("la dimension est 1, 2 ou 3.");
    }
}

} // namespace sac_de_billes

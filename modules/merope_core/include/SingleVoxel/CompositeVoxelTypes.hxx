//! Copyright : see license.txt
//!
//! \brief
//

#include "../../../Geometry/include/Point.hxx"
#include "../../../Geometry/include/BasicGeometricOperations.hxx"

#pragma once

namespace merope {
namespace vox {


template<class PHASE_TYPE_>
class PhaseFrac {
public:
    static constexpr unsigned short DIM = 0;
    using PHASE_TYPE = PHASE_TYPE_;
    PHASE_TYPE_ phase;
    double fracVol;  // volume fraction
    PhaseFrac(PHASE_TYPE_ p = 0, double f = 1) :
        phase{ p }, fracVol{ f } {}
    bool operator<(PhaseFrac<PHASE_TYPE_>& phfv2) {
        return this->phase < phfv2.phase;
    }
};

template<unsigned short DIM_, class PHASE_TYPE_>
class PhaseFracNormal {
public:
    static constexpr unsigned short DIM = DIM_;
    using PHASE_TYPE = PHASE_TYPE_;
    PHASE_TYPE phase;
    double fracVol;  // volume fraction
    Point<DIM> normal;  // implicitly of norm 1

    PhaseFracNormal(PHASE_TYPE p, double f, Point<DIM> normal_) :
        phase{ p }, fracVol{ f }, normal(normal_) {
        geomTools::renormalize<DIM>(normal);
    }

    PhaseFracNormal(PHASE_TYPE p = 0, double f = 1) :
        PhaseFracNormal<DIM, PHASE_TYPE>(p, f, create_array<DIM>(0.1)) {}


    bool operator<(PhaseFracNormal<DIM, PHASE_TYPE>& phfv2) {
        return this->phase < phfv2.phase;
    }
};


template<class C>
constexpr bool is_PhaseFrac = false;
template<class C>
constexpr bool is_PhaseFrac<PhaseFrac<C>> = true;

template<class C>
constexpr bool is_PhaseFracNormal = false;
template<class C, unsigned short DIM>
constexpr bool is_PhaseFracNormal<PhaseFracNormal<DIM, C>> = true;

}  // namespace vox
}  // namespace merope
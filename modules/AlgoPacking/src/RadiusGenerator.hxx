//! Copyright : see license.txt
//!
//! \brief Class for generating radii of spheres to be placed by algoRSA
//!	Builds at the very beginning the vector of desired radii and then gives it to algoRSA when appealing proceed.
//! As a class invariant, this vector is always decreasing.
#ifndef RADIUSGENERATOR_HXX_
#define RADIUSGENERATOR_HXX_

#include "StdHeaders.hxx"

#include "AmbiantSpace.hxx"
#include "GlobalShape.hxx"

namespace sac_de_billes {
using namespace std;

namespace algoRSA_aux {
class RadiusPhase {
//! represents the (radius,phase) of a sphere
public:
    RadiusPhase(double r, PhaseType p) :
            radius { r }, phase { p } {
    }
    double radius;
    PhaseType  phase;
};

class RPhiPhase: public RadiusPhase {
//! represents the (radius,phase,volumefraction) of a sphere
public:
    RPhiPhase(double r, double vf, PhaseType  p) :
            RadiusPhase(r, p), volFrac { vf } {
    }
    double volFrac;
};

template<unsigned short DIM>
class RadiusGenerator {
//! class used for generating successively radii of spheres to be placed
//! warning: the internal representation of radii takes into account a parameter called exclusionDistance, \see RadiusGenerator::setTabRadii
public:
    //! constructor from a list of desired R and phi
    //! \param desiredRPhi_ : vectors of desired {radius of spheres, volume fraction}
    //! \param exclusionDistance : minimal distance between spheres
    RadiusGenerator(AmbiantSpace::BigShape<DIM> *bigShape_,
            vector<array<double, 2>> desiredRPhi_, vector<PhaseType> phases_,
            double exclusionDistance);
    //! constructor from a vector of radii
    //! \param tabRadii_ : vector of radii
    //! \param exclusionDistance : minimal distance between spheres
    RadiusGenerator(AmbiantSpace::BigShape<DIM> *bigShape_,
            vector<double> tabRadii_, vector<PhaseType> phases_,
            double exclusionDistance);
    //! constructor from a file to read a vector of radii
    //! \param nameFile : file which radii are written
    //! 2 possible formats:
    //! 1) double radius on each line
    //! 2) (double radius, PhaseType phase) on each line
    //! \param exclusionDistance : minimal distance between spheres
    RadiusGenerator(AmbiantSpace::BigShape<DIM> *bigShape_, string nameFile,
            double exclusionDistance);

    //! generates a new radius.
    //! \return if a new sphere has to be placed
    bool nextRadius();
    //! \return the current radius
    double getRadius() const;
    //! \return the current phase
    PhaseType getPhase() const;
    //! \return the maximal radius
    double maxRadius() const;
    //! \return the minimal radius
    double minRadius() const;
    //! vector containing the tuples of achieved radii and volume fractions
    //! considers that the ball corresponding tabRadii[indexRadius] is the first not to be placed.
    vector<array<double, 2>> achievedRPhi() const;
    //! shrinks all radii by a factor
    void shrink(double factor);

private:
    //! constructor from a list of (radii,phase)
    //! \param tabRadPhase_ : vector of (radius,phase)
    //! \param exclusionDistance : minimal distance between spheres
    RadiusGenerator(AmbiantSpace::BigShape<DIM> *bigShape_,
            vector<RadiusPhase> tabRadPhase_, double exclusionDistance);
    //! tests whether there are radii to be placed and reduces the unused memory
    void emptyTabRadii();

    RadiusGenerator(AmbiantSpace::BigShape<DIM> *bigShape_) :
            tabRadii( { }), indexRadius(0) {
        bigShape = bigShape_;
    }
    //! Ambiant space
    AmbiantSpace::BigShape<DIM> *bigShape;
    //! stores all the radii and phases of the spheres to put
    vector<RadiusPhase> tabRadii;
    //! refers to the current radius
    size_t indexRadius;
    //! sets the desiredRPhi, and thus tabRadii
    void setDesiredRPhi(vector<RPhiPhase>);
    //! sets tabRadii
    //! \param tabRadii : computed table of radii
    //! \param exclusionDistance : the internal representation is radius -> radius + 0.5 exclusionDistance
    void setTabRadii(vector<RadiusPhase>, double exclusionDistance = 0);
};

//! defaut identifier for the phase
static constexpr PhaseType DEFAULT_PHASE = 0;

//! Sets default phases if possible
template<class T>
void checkPhase(const vector<T> &tabRadii_, vector<PhaseType> &phases);

//! create the vector of [R, phi] to insert into a boolean algorithm from a desired [R, phi]
template<class T>
vector<array<double, 2>> createRPhi_Bool(vector<array<double, 2>>);

} // namespace algoRSA_aux
} // namespace sac_de_billes

#include "RadiusGenerator.ixx"

#endif /* RADIUSGENERATOR_HXX_ */

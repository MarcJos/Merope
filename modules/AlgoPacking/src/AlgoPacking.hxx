//! Copyright : see license.txt
//!
//! \brief	RSA algorithm from
//! 				A Simple Algorithm for Maximal Poisson-Disk Sampling in High Dimensions
//! 				Mohamed S. Ebeida, Scott A. Mitchell, Anjul Patney, Andrew A. Davidson, and John D. Owens.
//!					Has been adapted in order to tackle as well different sphere sizes and
//!					objective volume fractions.
//!
//! This class defines interfaces for the library

#pragma once

#include "StdHeaders.hxx"

#include "AlgoRSA.hxx"
#include "AmbiantSpace.hxx"
#include "GlobalShape.hxx"
#include "Loops.hxx"
#include "RadiusGenerator.hxx"
#include "SphereContainer.hxx"
#include "AlgoNames.hxx"

namespace sac_de_billes {

class AlgoPacking_gen {
public:
        //! default constructor
        AlgoPacking_gen() {}
        //! default destructor
        virtual ~AlgoPacking_gen() {}

        //! OUTPUTS
        //! proceeds with the algorithm
        //! \param seed : seed for the stochastic process throwing the spheres
        //! \param method : can be 1 or 2, 1 corresponds to a slightly heavier method in memory, but faster, whereas 2 is lighter but slower. (1 is recommended; however the difference is small.)
        virtual map<string, string> proceed(unsigned seed, unsigned method) = 0;
        //! \return whether the spheres are actually intersecting (just for checking, very slow)
        virtual bool verifySphere() const = 0;

        //! INPUTS
        //! set the surrounding bigShape
        //! \param L_ : vector of size = dimension setting the dimensions of the bigCube
        //! \param nameShape: \see AmbiantSpace::NameShape
        virtual void setBigShape(vector<double> L, string nameShape_) = 0;
        //! set the exclusionDistance
        //! \param exclusionDistance : minimal distance between spheres
        virtual void setExclusionDistance(double exclusionDistance_) = 0;
        //! set the boundaryExclusionDistance
        //! \param boundaryExclusionDistance_ : minimal distance between the spheres and the boundary of the BigShape
        //! \warning : the boundaryExclusionDistance in bigShape representation takes the exclusionDistance into account. (Namely, it substracts it.)
        virtual void setBoundaryExclusionDistance(
                double boundaryExclusionDistance_) = 0;
        //! set the RadiusGenerator
        //! \see AlgoRSA_aux::RadiusGenerator<DIM>::RadiusGenerator
        virtual void setRadiusGenerator(vector<array<double, 2>> desiredRPhi_,
                vector<PhaseType> tabPhases_) = 0;
        virtual void setRadiusGenerator(vector<array<double, 2>> desiredRPhi_) = 0;
        //! set the RadiusGenerator
        //! \see AlgoRSA_aux::RadiusGenerator<DIM>::RadiusGenerator
        virtual void setRadiusGenerator(vector<double> tabRadii_,
                vector<PhaseType> tabPhases_ = { }) = 0;
        virtual void setRadiusGenerator(vector<double> tabRadii_) = 0;
        //! set the RadiusGenerator
        //! \see AlgoRSA_aux::RadiusGenerator<DIM>::RadiusGenerator
        virtual void setRadiusGenerator(string nameFile) = 0;
};

template<unsigned short DIM>
class AlgoPacking : public AlgoPacking_gen, public SphereContainer<DIM> {
public:
        AlgoPacking() :
                AlgoPacking_gen(), SphereContainer<DIM>(), exclusionDistance(0), boundaryExclusionDistance{
                        0 }, isSet_exclusionDistance(false), isSet_bigShape(false), isSet_radiusGen(
                                false), has_proceeded(false) {}

                        virtual ~AlgoPacking() {}

                        //! INPUTS
                        //! \see AlgoRSA_gen::setRadiusGenerator
                        void setRadiusGenerator(vector<array<double, 2>> desiredRPhi_,
                                vector<PhaseType> tabPhases_) override;
                        void setRadiusGenerator(vector<array<double, 2>> desiredRPhi_) override;
                        //! \see AlgoRSA_gen::setRadiusGenerator
                        void setRadiusGenerator(vector<double> tabRadii_, vector<PhaseType> tabPhases_)
                                override;
                        void setRadiusGenerator(vector<double> tabRadii_) override;
                        //! \see AlgoRSA_gen::setRadiusGenerator
                        void setRadiusGenerator(string nameFile) override;
                        void setRadiusGenerator(const algoRSA_aux::RadiusGenerator<DIM>&);
                        //! \see AlgoRSA_gen::setExclusionDistance
                        void setExclusionDistance(double exclusionDistance_) override;
                        //! \see AlgoRSA_gen::setExclusionDistance
                        void setBoundaryExclusionDistance(double boundaryExclusionDistance_)
                                override;
                        //! \see AlgoRSA_gen::setBigShape
                        void setBigShape(vector<double> L_, string nameShape_) override;
                        void setBigShape(vector<double> L_, AmbiantSpace::NameShape nameShape_);
                        void setBigShape(Point<DIM> L_, AmbiantSpace::NameShape nameShape_);

                        //! OUTPUTS
                        //! \see AlgoRSA_gen::verifySphere
                        bool verifySphere() const override;
                        //! returns the dimension of the box
                        Point<DIM> getLength() const override;
                        //! \return volume fraction
                        double volumeFraction() const;
                        //! \see AlgoPacking_gen::proceed
                        map<string, string> proceed(unsigned seed, unsigned method) override;
                        //! puts the output of the algorithm inside the sphere collection theSpheres
                        virtual void setSpheres() = 0;
                        //! get the type of the algorithm
                        virtual algoSpheres::TypeAlgo getTypeAlgo() const = 0;

protected:
        //! exclusion distance between spheres
        //! \warning because of this parameter, desiredRPhi differ from algoRSA to algoRSA_Template
        double exclusionDistance;
        //! exclusion distance between spheres and the boundary of the bigShape
        double boundaryExclusionDistance;
        //! pointer to the bigShape
        unique_ptr<AmbiantSpace::BigShape<DIM>> bigShape;
        //! pointer to the radius generator
        unique_ptr<algoRSA_aux::RadiusGenerator<DIM>> radiusGen;
        //! exclusionDistance shall be set once
        bool isSet_exclusionDistance;
        //! bigShape shall be set once
        bool isSet_bigShape;
        //! radiusGen shall be set once
        bool isSet_radiusGen;
        //! one can proceed only once
        bool has_proceeded;
        //! verify is exclusionDistance is not already set, and imposes to set the exclusionDistance before the boundary exclusionDistance
        void verifyExclusionDistance();
        //! verify if radiusGen is not set and sets it
        void verifyRadiusGenerator();
        //! verify if bigShape is not set and sets it
        void verifyBigShape();
        //!
        void printVolumeFraction();
        //! \see AlgoPacking<DIM>::proceed auxiliary function for proceeding, depend on the method employed
        virtual map<string, string> proceed_loc(unsigned seed, unsigned method) = 0;
};

template<unsigned short DIM, class T>  // T the type of algo used
class AlgoInterface : public AlgoPacking<DIM> {
public:
        AlgoInterface() :
                AlgoPacking<DIM>(), algo(nullptr) {}
        virtual ~AlgoInterface() {
                delete algo;
        }
        //! \see AlgoRSA_gen::proceed
        map<string, string> proceed_loc(unsigned seed, unsigned method) override;
        //! \return the type of the algorithm
        algoSpheres::TypeAlgo getTypeAlgo()  const override { return algo->getTypeAlgo(); };

private:
        //! \see AlgoPacking::setSpheres
        void setSpheres() override;
        //! pointer to the core algorithm
        T* algo;
};

template<unsigned short DIM>
using AlgoRSA = AlgoInterface<DIM, algoRSA_aux::AlgoRSA_Template<DIM>>;

class AlgoRSA2D final : public AlgoRSA<2> {
public:
        //! \param desiredRPhi_ : an array containing in 1st position the radii of the spheres, in 2nd position the desired volume fraction of them.
        //! \param exclusionDistance_ : exclusion distance between spheres
        //! \warning : because of this parameter, desiredRPhi differ from algoRSA to algoRSA_Template
        //! \param method : slightly different method chosen for checking intersection (!) does not change the output but only the time spent. Empirically, method 1 is 3x faster.
        AlgoRSA2D(Point<2> L, vector<array<double, 2> > desiredRPhi_,
                double exclusionDistance_, unsigned seed, unsigned short method = 1,
                string nameShape_ = "Tore");
};

class AlgoRSA3D final : public AlgoRSA<3> {
public:
        //! \param desiredRPhi_ : an array containing in 1st position the radii of the spheres, in 2nd position the desired volume fraction of them.
        //! \param exclusionDistance_ : exclusion distance between spheres
        //! \warning : because of this parameter, desiredRPhi differ from algoRSA to algoRSA_Template
        //! \param seed : seed for the random generator
        //! \param method : slightly different method chosen for checking intersection (!) does not change the output but only the time spent. Empirically, method 1 is 3x faster.
        AlgoRSA3D(Point<3> L, vector<array<double, 2> > desiredRPhi_,
                double exclusionDistance_, unsigned seed, unsigned short method = 1,
                string nameShape_ = "Tore");
        AlgoRSA3D(Point<3> L, vector<array<double, 2> > desiredRPhi_,
                double exclusionDistance_, unsigned seed, unsigned short method,
                AmbiantSpace::NameShape nameShape_);

        AlgoRSA3D() {}
};

double fractionVolThMaxRSA(short unsigned d);
}  // namespace sac_de_billes

#include "AlgoPacking.ixx"



//! Copyright : see license.txt
//!
//! \brief C++ easy interface for outer calls
//
#ifndef INTERFACE_HXX_
#define INTERFACE_HXX_

#include "StdHeaders.hxx"

#include "AmbiantSpace.hxx"
#include "AlgoNames.hxx"
#include "AlgoPacking.hxx"
#include "AlgoBool.hxx"
#include "AlgoRSA.hxx"
#include "AlgoWP.hxx"
#include "GlobalShape.hxx"

namespace sac_de_billes {
using namespace std;

namespace algoSpheres {

//! At the end, param algo points to the desired algorithm
template<unsigned short DIM>
unique_ptr<AlgoPacking<DIM>> newAlgo(const algoSpheres::TypeAlgo &typeAlgo);

//! \return a random distribution of spheres inside a prescribed shape
//!	\param typeAlgo : chooses the algorithm, RSA, WP or Bool
//! \param L : size of the box
//! \param nameShape : name of the global shape in which spheres are thrown
//! \param seed : seed for the random generator
//! \param desiredRPhi : desired (radius,volume_fraction) of spheres
//! \param tabPhases [optional] : phases relative to the desiredRPhi. If use {}, default phase = 0
//! \param mindist : repulsion distance of sphere, default=0. Inactive in case of bool
template<unsigned short DIM>
vector<Sphere<DIM>> throwSpheres(const algoSpheres::TypeAlgo &typeAlgo,
        AmbiantSpace::NameShape nameShape, Point<DIM> L, unsigned seed,
        vector<array<double, 2>> desiredRPhi, vector<PhaseType> tabPhases = { },
        double mindist = 0);

//! \return a random distribution of spheres inside a prescribed shape
//!	\param typeAlgo : chooses the algorithm, RSA, WP or Bool
//! \param L : size of the box
//! \param nameShape : name of the global shape in which spheres are thrown
//! \param seed : seed for the random generator
//! \param tabRadii : table of desired radii of spheres
//! \param tabPhases [optional] : phases relative to the tabRadii. If use {}, default phase = 0
//! \param mindist : repulsion distance of sphere, default=0. Inactive in case of bool
template<unsigned short DIM>
vector<Sphere<DIM>> throwSpheres(const algoSpheres::TypeAlgo &typeAlgo,
        AmbiantSpace::NameShape nameShape, Point<DIM> L, unsigned seed,
        vector<double> tabRadii, vector<PhaseType> tabPhases = { },
        double mindist = 0);

//! Evaluates the radius of monodisperse spheres in the torus, such that the volume fraction is prescribed (in expectation for boolean scheme).
//! \param L : dimensions of the torus
//! \param seedsNb : Seeds number
//! \param phi : prescribed volume fraction
//! \param typeAlgo : Spheres distribution method
template<unsigned short DIM>
double monoDispReferenceRadius(Point<DIM> L, size_t seedsNb, double phi, algoSpheres::TypeAlgo typeAlgo);
//! Evaluates the number of monodisperse spheres in the torus, such that the volume fraction is prescribed (in expectation for boolean scheme).
//! \param L : dimensions of the torus
//! \param radius : radius
//! \param phi : prescribed volume fraction
//! \param typeAlgo : Spheres distribution method
template<unsigned short DIM>
size_t monoDispNbSpheres(Point<DIM> L, double radius, double phi, algoSpheres::TypeAlgo typeAlgo);

//! return a list of spheres that are thrown according to the RSA in a near-maximal packing fraction
//! \param nameShape : name of the global shape in which spheres are thrown
//! \param L : size of the box
//! \param NbSpheres : number of spheres
//! \param seed : seed for the random generator
//! \param mindist : repulsion distance of sphere, default=0. Inactive in case of bool
template<unsigned short DIM>
vector<Sphere<DIM>> fillMaxRSA(AmbiantSpace::NameShape nameShape,
        Point<DIM> L, size_t NbSpheres, unsigned seed, double mindist = 0);


/*
//! return a vector of radii from a cumulative histogram
//! \param typeAlgo : Spheres distribution method
//! \param nameShape : name of the global shape in which spheres are thrown
//! \param L : size of the box
//! \param desiredRCumPhi : desired (radius,cumulated volume_fraction) of spheres
//! \param outFrac Output file name for saving initial and final fractions and numbers of spheres for each class
template<unsigned short DIM>
vector<double> fromCumHisto(algoSpheres::TypeAlgo typeAlgo,
        Point<DIM> L, AmbiantSpace::NameShape nameShape,
        vector<array<double, 2>> desiredRCumPhi, const char *const outFrac =
                nullptr);
*/
} // namespace algoSpheres
} // namespace sac_de_billes

#include "Interface.ixx"
#endif /* INTERFACE_HXX_ */

//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "StdHeaders.hxx"

#include "AlgoPacking.hxx"
#include "AlgoRSA.hxx"
#include "AmbiantSpace.hxx"
#include "GlobalShape.hxx"
#include "AlgoNames.hxx"

namespace sac_de_billes {
using namespace std;

template<unsigned short DIM>
class WPGrid;


namespace algoWP_aux {
class Pair : public array<size_t, 2> {
public:
    Pair(size_t i1, size_t i2) :
        array<size_t, 2>({ i1, i2 }) {}
    Pair() :
        Pair(0, 0) {}
    // lexicographical order
    bool operator<(const algoWP_aux::Pair& pair2) const;
};

template<unsigned short DIM>
class IntersectedSpheres {
public:
    IntersectedSpheres(WPGrid<DIM>* motherGrid_,
        AmbiantSpace::BigShape<DIM>* bigShape_) :
        motherGrid{ motherGrid_ }, bigShape{ bigShape_ } {}

    virtual ~IntersectedSpheres() {}

    //! index of spheres (wrt motherGrid->placedSpheres) that intersect other spheres
    vector<size_t> listOf;
    WPGrid<DIM>* motherGrid;
    //! computes the energy of distance = 1/4 \sum_{i != j} delta_ij^2.
    double energy_J() const;
    //! computes the new positions of the spheres after descending the gradient of J
    void updatePositions();
    vector<Sphere<DIM>> newPositions() const;
    //! add the pair of interacting spheres (index1,index2)
    void add(size_t index1, size_t index2);
    // setter
    void setListOfPairs(vector<algoWP_aux::Pair>&& a_listOfPairs);

private:
    //! for the building steps
    vector<algoWP_aux::Pair> listOfPairs;
    //! "time step" for removing overlapping spheres, cf [Schneider, Digital Microstructure ..., 2020; (2.4.3)]
    static constexpr double h_multiplier = 0.2505;
    //! Ambiant space
    AmbiantSpace::BigShape<DIM>* bigShape;
    //! derivative of the energy J
    //! BEWARE, there is an additional hidden parameter in this gradient, in the sense that very close spheres are slightly pushed
    vector<Point<DIM>> energy_NablaJ() const;
    //! build the lists listOf
    void buildLists();
};
}  // namespace algoWP_aux

template<unsigned short DIM>
class WPGrid : public algoRSA_aux::MotherGrid<DIM> {
public:
    WPGrid(DiscPoint<DIM> sizes_, AmbiantSpace::BigShape<DIM>* bigShape_,
        double voxelLength_, double minRadius_, double maxRadius_);
    //! pushes the spheres slightly further than necessary, see delta_ij
    static constexpr double MARGIN = 1 + 3e-3;

    //! shrinks the volume by multiplying the radii of spheres
    //! \param shrinkFactor is equal to 1/the factor by which the volume is shrinked = the factor by which the radius of sphere is increased
    void shrink(double shrinkFactor);
    //! removes overlapping spheres by moving them
    void removeOverlap();

private:
    //! designates the voxels in which a sphere is likely to overlap
    vector<algoRSA_aux::MotherVoxel<DIM>*> activatedVoxels;
    //! activate all voxels (at the beginning)
    void resetActivatedVoxels();
    //!	move the ovelapping spheres with two schemes
    //! sequentially NUM_METHOD_OVERLAP_= 1 -> slightly more efficient
    //! parallelly   NUM_METHOD_OVERLAP_= 2
    template<unsigned short NUM_METHOD_OVERLAP_>
    void moveOverlappingSpheres();
    void moveOverlappingSpheres_1();
    void moveOverlappingSpheres_2();
    //! if a sphere has been moved, computes the voxel where it has to be put
    void updateSphereIndex(size_t indexSphere,
        algoRSA_aux::MotherVoxel<DIM>* departureVoxel,
        algoRSA_aux::MotherVoxel<DIM>* arrivalVoxel);
    void updateSphereIndex(size_t indexSphere, const Sphere<DIM>& oldSphere,
        const Sphere<DIM>& newSphere);
    //! gathers all spheres that are near a given sphere
    vector<size_t> getNearSpheres(algoRSA_aux::MotherVoxel<DIM>* voxel);
    vector<size_t> getNearSpheres(algoRSA_aux::MotherVoxel<DIM>* voxel, std::vector<size_t>& theSpheres);
    //!
    algoWP_aux::IntersectedSpheres<DIM> findIntersectedSpheres();
    //! \see findIntersectedSpheres()
    //! auxiliary method for building a listOfPairs containing {index_1, index_2} of intersected spheres designated by index_1 and index_1.
    vector<algoWP_aux::Pair> findIntersectedSpheres_listOfPairs();
    //! maximum number of iterations for removeOverlap
    static constexpr size_t MAX_ITERATIONS = 10000;
    //! \see WPGrid::moveOverlappingSpheres
    static constexpr unsigned short NUM_METHOD_OVERLAP = 1;
};

namespace algoWP_aux {
template<unsigned short DIM>
class AlgoWP_Template {
public:
    //! constructor
    AlgoWP_Template(AmbiantSpace::BigShape<DIM>* T_,
        algoRSA_aux::RadiusGenerator<DIM>* radiusGen_, unsigned seed_);
    //! \see proceed_Template()
    map<string, string> proceed(unsigned short method);
    // Outputs
    //! vector containing the Spheres, refers to the grid
    vector<Sphere<DIM>>* placedSpheres;
    //! Ambiant space
    AmbiantSpace::BigShape<DIM>* bigShape;
    //! \return the type of the algorithm
    algoSpheres::TypeAlgo getTypeAlgo() const { return algoSpheres::TypeAlgo::WP; }

private:
    //! radius generator
    algoRSA_aux::RadiusGenerator<DIM>* radiusGen;
    //! global grid for localizing the spheres
    unique_ptr<WPGrid<DIM>> grid;
    //! gets a reasonable first guess for placing the spheres (by RSA)
    vector<Sphere<DIM>> initPlacedSpheres(unsigned short method);
    //! first steps before proceed() : first guess by RSA
    void initialize(unsigned short method);
    //! radius_multiplicator for RSA first pass
    static constexpr size_t NUMBER_OF_STEPS = 3;
    static constexpr double RADIUS_MULT = pow(0.4, 1. / DIM);
    static constexpr bool VERBOSE = false;
    //! says some informations
    void printFinalMessage();
    //! random seed
    unsigned seed;
};

template<unsigned short DIM>
inline double deltaij(const Sphere<DIM>& sph1, const Sphere<DIM>& sph2,
    double distance) {
    return max((sph1.radius + sph2.radius) * WPGrid<2>::MARGIN - distance, 0.);
}

template<unsigned short DIM>
inline double deltaij(const AmbiantSpace::BigShape<DIM>* bigShape,
    const Sphere<DIM>& sph1, const Sphere<DIM>& sph2) {
    return deltaij(sph1, sph2,
        sqrt(bigShape->distanceCarre(sph1.center, sph2.center)));
}
}  // namespace algoWP_aux

template<unsigned short DIM>
using AlgoWP = AlgoInterface<DIM, algoWP_aux::AlgoWP_Template<DIM>>;

}  // namespace sac_de_billes

#include "AlgoWP.ixx"



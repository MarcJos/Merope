//! Copyright : see license.txt
//!
//! \brief Main part of the code.
//! Description : 	RSA algorithm from
//! 				A Simple Algorithm for Maximal Poisson-Disk Sampling in High Dimensions
//! 				Mohamed S. Ebeida, Scott A. Mitchell, Anjul Patney, Andrew A. Davidson, and John D. Owens.
//!					Has been adapted in order to tackle as well different sphere sizes and
//!					objective volume fractions.

#ifndef ALGORSA_HXX_
#define ALGORSA_HXX_

#include "StdHeaders.hxx"

#include "AmbiantSpace.hxx"
#include "GlobalShape.hxx"
#include "Path.hxx"
#include "RadiusGenerator.hxx"
#include "Voxel.hxx"
#include "Loops.hxx"
#include "AlgoNames.hxx"
#include "MultiDArrayObject.hxx"
#include "ArrayDimensions.hxx"


namespace sac_de_billes {

using namespace std;

namespace algoRSA_aux {

template<unsigned short DIM>
class MotherVoxel;

template<unsigned short DIM>
class MotherGrid {
    //! Main grid, that contains the motherVoxels and informations necessary for find the spheres
public:
    MotherGrid(DiscPoint<DIM> sizes_, AmbiantSpace::BigShape<DIM>* bigShape_,
        double voxelLength_, double minRadius_, double maxRadius_);
    //! Numbers of nodes in x, y, z directions (spatial discretization)
    DiscPoint<DIM> sizes;
    //! Ambiant space
    AmbiantSpace::BigShape<DIM>* bigShape;
    //! 3D table of Voxels
    typedef typename merope::vox::MultiDArrayObject < DIM, MotherVoxel<DIM>, merope::vox::ArrayDimensions<DIM>> TYPE_TABVOXEL; //!< defines the type of tabVoxels
    TYPE_TABVOXEL tabVoxels;
    //! getter for tabVoxels
    template<typename T>
    MotherVoxel<DIM>* getVoxel(array<T, DIM> index);
    //! voxel length
    double voxelLength;
    //! inverse of the voxelLength
    double inverse_voxelLength;
    //! minimal radius of spheres
    double minRadius;
    //! maximal radius of spheres
    double maxRadius;
    //! vector containing the Spheres
    vector<Sphere<DIM>> placedSpheres;

    //! \return the number of voxels
    size_t getNbVoxels();
    //! \return true if we may place a given sphere
    template<unsigned short NUM_METHOD>
    bool isPlaceable(const Sphere<DIM>& sphere, const MotherVoxel<DIM>* mother);
    //! adds a sphere to the grid, and computes consequences
    template<unsigned short NUM_METHOD>
    void addSphere(const Sphere<DIM>& sphere, MotherVoxel<DIM>* motherVoxel);
    //! informs each surrounding voxel of a sphere that he may be intersected with it.
    void linkSphereToVoxel(MotherVoxel<DIM>* motherVoxel,
        const Sphere<DIM>& sphere, size_t iSphe);
    //! \return the spheres already placed that might intersect the point or voxel
    vector<Sphere<DIM>> neighborSpheres(const MotherVoxel<DIM>* voxel) const;
    //! \return the coordinates of the voxel corresponding to the given point (either in discrete or double coordinates)
    array<size_t, DIM> discCoordinate(const DiscPoint<DIM>& discCoord) const;
    array<size_t, DIM> discCoordinate(const Point<DIM>& doubleCoord) const;
    //! \warning :  this one does not that the modulo (the coordinates might be outside the grid)
    DiscPoint<DIM> discCoordinate_withoutModulo(const Point<DIM>& point) const;

protected:
    //! \return an unordered_map {relative_position, surrounding voxel}
    //! \param centralVoxel : discreteCoordinates of the central voxel
    //! \param iMin : minimal relativePosition
    //! \param iMax : maximal relativePosition
    vector<tuple<array<int, DIM>, algoRSA_aux::MotherVoxel<DIM>*>> surroundingVoxels(
        DiscPoint<DIM> centralVoxel, array<int, DIM> iMin,
        array<int, DIM> iMax);
    //! \return the list of voxels close to a given voxel, within a distance multiple of radius
    //! \param layerWidth : counts how many voxel from the central voxel to the surface of the cube
    vector<tuple<array<int, DIM>, algoRSA_aux::MotherVoxel<DIM>*>> surroundingVoxels(
        MotherVoxel<DIM>* centralVoxel, int layerWidth);

private:
    //! number of expected sphers intersecting the cube
    static constexpr long NB_EXPECTED_SPHERES_IN_CUBE = 25;

    //! \return true if we may place a given sphere
    //! \see isPlaceable
    bool isPlaceable_1(const Sphere<DIM>& sphere,
        const MotherVoxel<DIM>* mother) const;
    //! \return true if we may place a given sphere
    //! \see isPlaceable
    bool isPlaceable_2(const Sphere<DIM>& sphere,
        const MotherVoxel<DIM>* mother);
    //! auxiliary function for isPlaceable_2
    template<size_t taille>
    bool checksOtherSpheres_2(const DiscPoint<DIM>& discCoord,
        const Sphere<DIM>& sph1,
        const array<array<int, DIM>, nbSubcubes<DIM>(taille)>& tabPath);
    //! auxiliary function for MotherGrid::MotherGrid :  initializes the i-th voxel
    //! \param i : discrete coordinate of the voxel
    void initializeVoxel(const DiscPoint<DIM>& i);
    //! auxiliary function for MotherGrid::linkSphereToVoxel. Links a specific sphere to a specific voxel
    //! \param *currentVoxel : voxel the sphere should be linked to
    //! \param iSphe : index of the sphere
    //! \param i : index indicating the relative position of the voxel w.r.t. the sphere
    //! \param sphereOnBoundary : is the sphere close to boundary?
    void auxi_linkSphereToVoxel(MotherVoxel<DIM>* currentVoxel,
        const size_t& iSphe, const array<int, DIM>& i,
        const bool& sphereOnBoundary);
    //! enlarges the limits iMin & iMax if they cross the boundary of the cube
    //! \return whether the limits has been enlarged
    bool correctForBoundaries(algoRSA_aux::MotherVoxel<DIM>* motherVoxel,
        array<int, DIM>& iMin, array<int, DIM>& iMax);
};

template<unsigned short DIM>
class CurrentGrid {
    //! Grid containing uncovered voxels. Does not have direct access to the spheres.
    //! Decreases when algorithm goes.
public:
    //! constructor
    CurrentGrid(MotherGrid<DIM>* motherGrid_, mt19937* randGen_);

    //! list voxels potentially uncovered
    vector<Voxel<DIM>> tabVoxels;
    //! pointer to the mother grid
    MotherGrid<DIM>* motherGrid;
    //! uniform distribution to reach a voxel potentially uncovered
    uniform_int_distribution<> distribution;
    //! random generatpr
    mt19937* randGen;
    //! number of uncovered voxels ( =/= tabVoxels.size() in general)
    long nbUncoveredVoxels;

    //! set the tabVoxels (accordingly, nbUncoveredVoxels and distribution)
    void setTabVoxels(const vector<Voxel<DIM>>& tabVoxels_);
    //! \return a random activated voxel
    Voxel<DIM>* pickUncoveredVoxel();
    //! subdivides the grid
    template<unsigned short NUM_METHOD>
    void subdivide();
};

template<unsigned short DIM>
class AlgoRSA_Template {
    //! Main class for applying the RSA.
public:
    //! constructor
    AlgoRSA_Template(AmbiantSpace::BigShape<DIM>* T_,
        RadiusGenerator<DIM>* radiusGen_, unsigned seed);
    //! \see proceed_Template()
    map<string, string> proceed(unsigned short method);

    virtual ~AlgoRSA_Template();

    // Outputs
    //! vector containing the Spheres, refers to the motherGrid
    vector<Sphere<DIM>>* placedSpheres;
    //! Ambiant space
    AmbiantSpace::BigShape<DIM>* bigShape;
    //! \return the type of the algorithm
    algoSpheres::TypeAlgo getTypeAlgo() const { return algoSpheres::TypeAlgo::RSA; }

private:
    //! number of darts per uncovered voxel
    static constexpr double PARAM_A = 4;
    //! maximal number of subdivision
    static constexpr double NB_MAX_SUB_DIV = 30;
    //! is used to check whether something wrong happend when placing the balls.
    static constexpr long MULT_MAX = 100;
    //! maximal number of subvoxels before bugging
    static constexpr long NB_MAX_SUB_VOXELS = 9999999;
    //! Contains the global grid
    MotherGrid<DIM>* motherGrid;
    //! Contains the uncovered voxels
    CurrentGrid<DIM>* currentGrid;
    //! random generator
    mt19937 randGenerator;
    //! uniform law between [0,1)
    std::uniform_real_distribution<> randomReal;
    //! radius generator
    RadiusGenerator<DIM>* radiusGen;
    //! maximal radius of spheres
    double maxRadius;
    //! minimal radius of spheres
    double minRadius;
    //! maximal number of voxels (if more, the algorithm considers it cannot achieve a higher volume fraction)
    long nbMaxVoxels;

    //! throws a dart and gives the consequences
    //! \return whether this leads to a new ball or not
    template<unsigned short NUM_METHOD>
    bool throwDart();
    //! prints if not fully packed
    void printFinalMessage(bool noMoreSpheres = false) const;
    //! print information after each subdivision
    void printStep() const;
    //! proceeds the RSA algorithm
    //! \return if fully packed
    template<unsigned short NUM_METHOD>
    bool proceed_T();
};

} // namespace algoRSA_aux
} // namespace sac_de_billes

#include "AlgoRSA.ixx"

#endif /* ALGORSA_IXX_ */

//! Copyright : see license.txt
//!
//! \brief 2 classes for implementing voxel representation
//
#ifndef VOXEL_HXX_
#define VOXEL_HXX_

#include "StdHeaders.hxx"

#include "AmbiantSpace.hxx"
#include "GlobalShape.hxx"
#include "Path.hxx"


namespace sac_de_billes {
using namespace std;

namespace algoRSA_aux {

template<unsigned short DIM>
class MotherGrid;
template<unsigned short DIM>
class MotherVoxel;

template<unsigned short DIM>
class Voxel {
public:
    Voxel() : discreteCoordinates(create_array<DIM, long>(0)), length{ 0. }, covered{ false }, motherGrid{ nullptr }, mother{ nullptr } {}
    //! constructor
    Voxel(const DiscPoint<DIM>& discreteCoordinates_, double length_,
        MotherGrid<DIM>* motherGrid_);

    //! \return a random point picked uniformly inside the voxel
    //! checks that the point is indeed inside the bigShape
    Point<DIM> pickPoint(mt19937& randGenerator,
        uniform_real_distribution<>& randomReal);

    //! refreshes the status (covered or not) of the Voxel wrt the sphere that is added
    void updateCovered(const Sphere<DIM>& sphere, double minRadius);
    //!	\return whether it is totally covered by the sphere or not
    //bool isInSphere(const Sphere<DIM>& sphere, double minRadius) const;
    bool isInSphere(const Sphere<DIM>& sphere, const double& minRadius,
        const vector<array<unsigned short, DIM>>& pathForCorner =
        path::TabCorner<DIM>::get().getTab()) const;
    //! \return if covered by any of the speres in the sphereList
    bool isCovered(const vector<Sphere<DIM>>* sphereList,
        double minRadius) const;
    //! \return if it is outside the bigShape
    bool isOutside(AmbiantSpace::BigShape<DIM>*, double minRadius);

    //! \fill in the vector futureTabVoxel the new uncovered voxels when subdividing the voxel.
    //! BY METHOD 1: (copying the spheres into the motherVoxel it intersects)
    //! BY METHOD 2: (searching dynamically for near spheres)
    //! \param futureTabVoxel : voxel vector to be filled in
    //! \param minRadius : minimal radius
    template<unsigned short NUM_METHOD>
    void subdivide(vector<Voxel<DIM>>& futureTabVoxel) const;

    //! \return the corner coordinates
    //! chooses the corner w.r.t. path::TABCORNER
    //! \param index in [0,7],
    Point<DIM> corner(int index) const;
    //! \return the corner coordinates
    //! \param diCoord: indicates the position of the corner wrt the center of the cube
    Point<DIM> corner(const array<unsigned short, DIM>& discCoord) const;
    //! coordinates of the closest corner to 0
    Point<DIM> position() const;
    //! coordinate of the center of the voxel
    Point<DIM> center() const;

public:
    //! position of the corner of the voxel of minimal coordinates in terms of length
    DiscPoint<DIM> discreteCoordinates;
    //! length of the voxel
    double length;
    //! is it ENTIRELY covered by a SINGLE sphere?
    // (Should take into account also the (minimal) radius of the new spheres we wish to place.) ie Sphere of radius r + minRadius
    //! true = yes
    //! false = we do not know yet
    bool covered;
    //! refers to the motherGrid
    MotherGrid<DIM>* motherGrid;
    //! mother voxel
    MotherVoxel<DIM>* mother;
};

template<unsigned short DIM>
class MotherVoxel : public Voxel<DIM> {
    //! Special voxels of the motherGrid.
    //! Besides the normal properties, helps to find spheres.
public:
    MotherVoxel() : Voxel<DIM>{}, spheresInside{}, spheresRelativePosition{}, isClose2Boundary{ false } {}
    //! Constructor
    MotherVoxel(const DiscPoint<DIM>& discreteCoordinates_, double length_,
        MotherGrid<DIM>* motherGrid_);
    //! Gets the sphere through spheresInside, to placedSpheres
    inline Sphere<DIM> getSphere(int i) const;
    //! Can it feel influence of spheres through the boundary of the periodic cube? -> initializes the variable isClose2Boundary
    //! Checking if the cube cointaining any close voxel, in which a sphere can intersect a sphere of the original voxel is indeed in the BigShape
    bool isItClose2Boundary(AmbiantSpace::BigShape<DIM>* bigShape);

public:
    //! Spheres already placed inside the voxel. Refers to the list placedSpheres
    vector<size_t> spheresInside;
    //! Relative position of the spheres placed inside, wrt pathFor27Corners
    vector<int> spheresRelativePosition;
    //! Can it feel influence of spheres through the boundary of the periodic cube?
    bool isClose2Boundary;
};

template<unsigned short DIM>
class VoxelSubdivision {
    //! Auxiliary class to efficiently subdivide an uncovered voxel into 2^DIM uncovered subvoxels.
    //! Checks that subvoxels are indeed uncovered.
private:
    //! voxel containing the others
    const Voxel<DIM>* fatherVoxel;
    //! minimal radius
    double minRadius;
    //! Middle point of the subdivision
    Point<DIM> midPoint;
    //! Number of subvoxels
    static constexpr unsigned short SIZE = nbSubcubes<DIM>();
    //! Contains the candidate subvoxels
    array<Voxel<DIM>, SIZE> tabNewVoxels;
    //! tabUnCovered[i] says whether tabNewVoxels[j] is uncovered (at the end)
    array<bool, SIZE> tabUnCovered;
    //! Number of uncovered subvoxels in the table tabNewVoxels
    size_t NbUncoveredVoxels;
public:
    VoxelSubdivision(const Voxel<DIM>* voxel) :
        fatherVoxel(voxel), minRadius(voxel->motherGrid->minRadius), midPoint(
            voxel->center()), tabNewVoxels{ }, tabUnCovered{ }, NbUncoveredVoxels(
                SIZE), closeSpheres{ }, closeSpheresRelativePos{ } {
        initialize();
    }

    //! returns whether if all the subvoxels are covered
    bool isEmpty();
    //! checks whether subvoxels are covered or not by method 1
    void checkCoveredVoxel_1();
    //! checks whether subvoxels are covered or not by method 2
    void checkCoveredVoxel_2();
    //! appends futureTabVoxel with the remaining uncovered subvoxels
    void fillNewVoxels(vector<Voxel<DIM>>& futureTabVoxel);
    //! find and sets the spheres that may intersect the subvoxels
    void findSpheres();
    //! is the central point of the subdivision covered by the sphere?
    //! (if not, no chance that the sphere covers any of the subcubes)
    bool hitsMidPoint(Sphere<DIM>& sphere);
    //! tests whether the j-th subcube is covered by a sphere seen at relativePosition
    template<unsigned short NUM_METHOD, unsigned short NB_NGHB = 0>
    bool testSphere(const Sphere<DIM>& sphere, const size_t j,
        const int relativePosition = 0);
private:
    vector<Sphere<DIM>> closeSpheres;
    vector<int> closeSpheresRelativePos;
    //! Initialize the properties
    void initialize();
    //! Set that the jth subvoxel is covered.
    void setFalse(int j);
    //! checks whether the sphere covers the jth subvoxel.
    //! in case yes, it updates the subdivision
    bool checkSphere(const int& j, const Sphere<DIM>& sphere,
        const vector<array<short unsigned, DIM>>& pathForCorner);
    //! auxiliary function for checkCoveredVoxel_2
    template<unsigned short NB_NGHB>
    void testSubVoxels();
};

} // namespace algoRSA_aux
} // namespace sac_de_billes

#include "Voxel.ixx"

#endif /* VOXEL_HXX_ */

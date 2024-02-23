//! Copyright : see license.txt
//!
//! \brief 
//!
#ifndef VOROINTERFACE_HXX_
#define VOROINTERFACE_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include <concepts>

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../Geometry/GeomTools.hxx"
#include "../MesoStructure/InsideTorus.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"

#include "container.hh"
#include "pre_container.hh"

//! prepare the use of voro++
namespace voro {
class container_poly;
class voronoicell;
}

#include "../MeropeNamespace.hxx"


namespace merope {

namespace voroInterface {

struct SingleCell;

class PreparedVoroppContainer {
public:
    //! prepares the VoroppContainer
    //! \param centerTessels : list of the centers of the voronoi cells + weight for Laguerre
    //! \param L : lenghts of the torus
    //! \param voropp_container : pointer to be modified, shall contain at the end the tessellation
    //! \param periodicity : periodicity in all directions
    PreparedVoroppContainer(const vector<Sphere<3>>& centerTessels,
        array<double, 3> L, array<bool, 3> periodicity = { true, true, true });
    //! @brief : destructor
    virtual ~PreparedVoroppContainer() { delete voropp_container; }
    //! @brief : contains the interface to the voro++ code
    voro::container_poly* voropp_container;

private:
    // DELETE OPERATORS
    PreparedVoroppContainer(PreparedVoroppContainer&) = delete;
    PreparedVoroppContainer(const PreparedVoroppContainer&) = delete;
    PreparedVoroppContainer& operator=(PreparedVoroppContainer&) = delete;
    PreparedVoroppContainer& operator=(const PreparedVoroppContainer&) = delete;
    PreparedVoroppContainer operator=(PreparedVoroppContainer) = delete;
    //! parameter for voro++
    static constexpr int VOROPP_INIT_MEM = 8;
};

class ListOfVoroppCells : protected InsideTorus<3>, public PreparedVoroppContainer {
public:
    //! main constructor
    //! \param L : lengths of the torus
    //! \param centerTesssels : center and weights of the tessels
    ListOfVoroppCells(array<double, 3> L, const vector<Sphere<3>>& centerTessels);
    //! @brief  destructor
    ~ListOfVoroppCells();
    //! @brief appeals the voro++ code
    void build();
    //! @return whether the structure is coherent
    bool is_coherent() const { return vec_spheres.size() == vec_cells.size(); }

public:
    //! @brief : centers of the tessels with radius as weight
    vector<Sphere<3>> vec_spheres;
    //! @brief : contains the interface to the voro++ code
    vector<voro::voronoicell_neighbor*> vec_cells;
};


template<unsigned short DIM>
//! Interface for exposing some functionalities of Voro++
class VoroInterface : protected InsideTorus<DIM>, private PreparedVoroppContainer {
    static_assert(DIM == 2 or DIM == 3);
public:
    //! main constructor
    //! \param L : lengths of the torus
    //! \param centerTesssels : center and weights of the tessels
    VoroInterface(array<double, DIM> L, const vector<Sphere<DIM>>& centerTessels);
    //! \returns the tessel containing the point
    //! \param pt : a point of the torus
    int findTessel(const Point<DIM>& pt);
    //! draw the particles in gnuplot format
    void drawGnuPlot(string fileName);
    //! get all the cells
    vector<smallShape::ConvexPolyhedronInc<DIM>> getMicroInclusions();
    //! get all the cells
    //! \warning : will not work in dimension d = 2
    vector<SingleCell> getSingleCells();

protected:
    //! auxiliary function to add an inclusion to the list polyhedrons, used by getMicroInclusion()
    //! \param pp : pointer to the center of the cell
    //! \param index : index of the cell
    //! \param c : pointer to the cell
    //! \param polyhedrons : list of polyhedrons
    void addInclusion(double* pp, int index, voro::voronoicell_neighbor* c, vector<smallShape::ConvexPolyhedronInc<DIM>>& polyhedrons);
    //! auxiliary function to add a cell to the list polyhedrons, used by getSingleCells()
    //! \param pp : pointer to the center of the cell
    //! \param index : index of the cell
    //! \param c : pointer to the cell
    //! \param polyhedrons : list of cells
    void addCell(double* pp, int index, voro::voronoicell_neighbor* c, vector<SingleCell>& polyhedrons);
    //! if DIM=2, get the 3-DIM length
    array<double, 3> virtualLength;
};

class SingleCell {
public:
    //! constructor from a voro++ cell
    SingleCell(Identifier identifier_, const Point<3>& center_, voro::voronoicell_neighbor* c);
    //! trivial constructor
    SingleCell(Identifier identifier_,
        const Point<3>& center_,
        const vector<Point<3>>& vertices_,
        const vector<Point<3>>& faceNormal_,
        const vector<vector<Identifier>>& faceVertices_,
        const vector<Identifier>& neighbors_);
    //! for having an identifier
    Identifier identifier;
    //! center of the cell
    Point<3> center;
    //! vertices constituting the cell
    vector<Point<3>> vertices;
    //! normal of each face
    vector<Point<3>> faceNormal;
    //! vertices of each face
    vector<vector<Identifier>> faceVertices;
    //! to each face, the related neighbor
    vector<Identifier> neighbors;
    //! print into the stream f
    void print(std::ostream& f) const;
private:
    //! correct the normal so that none is (0, 0, 0)
    void correctNormals();
};

namespace voroInterface_aux {
//! transforms v into a list of Point<3>
//! \param v : a purely linear list of doubles
vector<Point<3>> breakList(const vector<double>& v);
//! if DIM = 2, the last virtual 3-DIM is VIRTUAL_LENGTH_MULTIPLIER smaller than the other lengths
constexpr double VIRTUAL_LENGTH_MULTIPLIER = 0.05;
//! if DIM=2, extends the length to a 3-DIM Length
template<unsigned short DIM>
Point<3> virtualLength(Point<DIM>);

//! \return a couple of list of normal, and of a list of vertices
//! \param cell : a voro++ cell
template<unsigned short DIM>
vector<HalfSpace<DIM>> getFaces(voro::voronoicell_neighbor* cell);
vector<HalfSpace<3>> getFaces3D(voro::voronoicell_neighbor* cell);
vector<HalfSpace<2>> getFaces2D(voro::voronoicell_neighbor* cell);
//! \return the cuboid in which the voronoi cell lie
//! \param cell : a voro++ cell
template<unsigned short DIM>
Cuboid<DIM> getCuboid(voro::voronoicell_neighbor* cell, const Point<DIM>& center);
//! \return the list of vertices of  the voronoi cell
//! the origin is the center given
//! \param cell : a voro++ cell
vector<Point<3>> getVertices(voro::voronoicell_neighbor* cell, const Point<3>& center);
//! \return the list of vertices of  the voronoi cell
//! the origin is the center of the polyhedron
//! \param cell : a voro++ cell
vector<Point<3>> getRenormalizedVertices(voro::voronoicell_neighbor* cell);
//! \return the list of indexes of vertices related to each face
//! \param cell : a voro++ cell
vector<vector<Identifier>> getFacesToVertices(voro::voronoicell_neighbor* cell);
//! get the list of normals of the voronoi cell
//! \param cell : a voro++ cell
vector<Point<3>> getNormals(voro::voronoicell_neighbor* cell);
//! \return the list of neighbors of a cell
//! \param cell : a voro++ cell
vector<long> getNeighbors(voro::voronoicell_neighbor* cell);

//! convert a point in 2 dimensions to a point in 3 dimensions by a adding a 3rd coordinate equal to 0.5 L[3]
//! \param L : length of the (virtual) 3-dimensional torus
Point<3> extendDimension(const Point<2>& oldPoint, const array<double, 3>& L);
//! for each sphere, convert its center in 2 dimensions to a center in 3 dimensions by a adding a 3rd coordinate equal to 0.5 L[3]
//! \param L : length of the (virtual) 3-dimensional torus
vector<Sphere<3>> extendDimension(const vector<Sphere<2>>& oldSpheres, const array<double, 3>& L);

//! \return the tessel containing the point
//! \param pt : a point of the torus
int findTessel(const Point<3>& pt, voro::container_poly* voropp_container);
} // namespace voroInterface_aux

//! for each sphere, get the closest neighbors
template<unsigned short DIM>
std::map<Identifier, vector<Identifier>> getClosestNeighbors(array<double, DIM> L, const vector<Sphere<DIM>>& centerTessels);

//! @return the list of the volume of each tessel
//! \param centerTessels : list of the centers of the voronoi cells + weight for Laguerre
//! \param L : length of the (virtual) 3-dimensional torus
vector<double> compute_volumes(const vector<Sphere<3>>& centerTessels, const Point<3>& L);
//! @return the list of centroid of each tessel
//! \param centerTessels : list of the centers of the voronoi cells + weight for Laguerre
//! \param L : length of the (virtual) 3-dimensional torus
vector<Point<3>> compute_relative_centroids(const vector<Sphere<3>>& centerTessels, const Point<3>& L);


//! @brief build the tessels and evaluate on each cell the function
//! @param my_function : lambda function to be applied to (double* center, size_t index_of_cell, voro::voronoicell_neighbor &c))
template<class Func>
void loop_on_voroppcontainer(voro::container_poly* voropp_container, Func my_function);

} // namespace voroInterface
} // namespace merope


#include "../Voronoi/VoroInterface.ixx"
#endif /* VOROINTERFACE_HXX_ */

//! Copyright : see license.txt
//!
//! \brief MicroInclusions
//
//
#ifndef MICROINCLUSION_HXX_
#define MICROINCLUSION_HXX_

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../Geometry/GeomTools.hxx"
#include "../Geometry/SpheroPolyhedron.hxx"
#include "../Geometry/AuxiConvexPolyhedron.hxx"

#include "../MeropeNamespace.hxx"

namespace merope {
namespace smallShape {
//! For describing a single small inclusion inside the material

//! create a cuboid covering a SOLID
//! fixme : case of ellipse not optimal
template <unsigned short DIM, class SOLID>
Cuboid<DIM> createCuboid(const SOLID& sph);
//! create a cuboid covering a SOLID, with an additional margin on each face
template <unsigned short DIM, class SOLID>
Cuboid<DIM> createCuboid(const SOLID& sph, double margin);

struct singleLayer {
    size_t index;
    array<double, 2> limits;
};

template <unsigned short DIM, class SOLID>
class MicroInclusion_ {
    //! MicroInclusion as a subset of the space. Is inside a cuboid (for voxellisation)
public:
    MicroInclusion_(Identifier i, Cuboid<DIM> cuboid_, Point<DIM> center_, SOLID solid): identifier{ i }, center(center_), cuboid{ cuboid_ }, layerPhases({ i }), layerIncrement({ 0 }), solids({ solid }) {}
    MicroInclusion_(const SOLID& solid);
    //! The space is dilated through (x_1,x_2,x_3) -> (lambda_1 x_1, lambda_2 x_2, lambda_3 x_3)
    void linearTransform(const Point<DIM>& linTransform);
    //! insert a new layer
    //! \param layerWidth : width of the new layer
    //! \param layerPhase : phase of the new layer
    void pushLayer(double layerWidth, PhaseType layerPhase_);
    //! get the phase of the ith component (0 = outer)
    PhaseType getPhaseForVoxellation(size_t i) const {
        assert(i < this->getNbOfLayers());
        return layerPhases[i];
    }
    //! get the number of layers
    size_t getNbOfLayers() const { return layerPhases.size(); }
    //! getter
    const vector<PhaseType>& getLayerPhases() const { return layerPhases; }
    //! getter for modification
    PhaseType& getPhaseGraphical(size_t i);
    //! getter
    const PhaseType& getPhaseGraphical(size_t i) const;
    //! getter
    const Cuboid<DIM>& getCuboid() const { return cuboid; }
    //! return all the nested inclusions as a list of faces
    const vector<SOLID>& getInnerInclusions() const { return solids; }

    //! BEWARE, in local coordinates!
    //! return (signed) the distance to the surface of the inclusion (negative if inside, positive is outside)
    //! \param vector_from_center_to_point : point in local coordinates
    //! \param layerIndex : index of the layer
    double distanceTo(const Point<DIM>& vector_from_center_to_point, size_t layerIndex = 0) const;
    //! BEWARE, in local coordinates!
    //! \return in which layer the point is, (-1 if not inside)
    //! \param vector_from_center_to_point : point in local coordinates
    long whichLayer(const Point<DIM>& point) const;
    //! BEWARE, in local coordinates!
    //! \param vector_from_center_to_point : point in local coordinates
    //! \param halfDiagonal : 0.5 * diagonal of the voxel
    //! \param layerIndex : index of the layer
    bool guaranteeInside(const Point<DIM>& vector_from_center_to_point, const double& halfDiagonal, const size_t& layerIndex = 0) const;
    //! BEWARE, in local coordinates!
    //! \param vector_from_center_to_point : point in local coordinates
    //! \param halfDiagonal : 0.5 * diagonal of the voxel
    //! \param layerIndex : index of the layer
    bool guaranteeOutside(const Point<DIM>& vector_from_center_to_point, const double& halfDiagonal, const size_t& layerIndex = 0) const;

public:
    //! designates the index of the Inclusion
    Identifier identifier;
    //! center of the polyhedron (seed of voronoi distribution)
    Point<DIM> center;

protected:
    //! a cube that contains the inclusion
    Cuboid<DIM> cuboid;
    //! contains the indices related to layerLimits
    //! the layers are order as follows :
    //! layerPhase[0] -> outer layer
    //! layerPhase[n] -> inner layer
    vector<PhaseType> layerPhases;
    //! contains [0, a_1, a_2, ..., a_n for layers [0,b_1], [b_{1},b_{2}], ... [b_n,\infty) from the surface to the center of the MicroInclusions
    //! for b_k - b_{k-1} = a_k
    vector<double> layerIncrement;
    //! contains the inner shapes
    vector<SOLID> solids;

private:
    //! get the layer width
    double getLayerIncrement(size_t i) const { return layerIncrement[i]; }
    //! BEWARE, in local coordinates!
    //! \param vector_from_center_to_point : point in local coordinates
    //! \param layerIndex : index of the layer
    bool isInside(const Point<DIM>& vector_from_center_to_point, const size_t& layerIndex = 0) const;
};

template <unsigned short DIM>
using SphereInc = MicroInclusion_<DIM, sac_de_billes::Sphere<DIM>>;

template <unsigned short DIM>
using EllipseInc = MicroInclusion_<DIM, sac_de_billes::Ellipse<DIM>>;

template <unsigned short DIM>
using ConvexPolyhedronInc = MicroInclusion_<DIM, ConvexPolyhedron<DIM>>;

template <unsigned short DIM>
using SpheroPolyhedronInc = MicroInclusion_<DIM, SpheroPolyhedron<DIM>>;

template <unsigned short DIM>
class Rectangle final: public ConvexPolyhedronInc<DIM> {
    // a cuboid MicroInclusion
public:
    Rectangle(Identifier i, Point<DIM> xmin, Point<DIM> xmax);
};

namespace auxi_MicroInclusions {
//! return the layerPhase index given a phaseIndex
size_t getIndexPhaseGraphical(size_t phaseIndex, size_t nbOfLayers);
} // namespace auxi_MicroInclusions

//! Standard type to add a layer
struct LayerInstructions {
    //! constructor
    LayerInstructions(Identifier identifier_, PhaseType phase_, double width_): identifier{ identifier_ }, phase{ phase_ }, width{ width_ } {};
    //! identifier of the inclusion particle
    Identifier identifier;
    //! phase to be added
    PhaseType phase;
    //! width of the layer
    double width;
    //! for sorting LayerInstructions
    bool operator<(const LayerInstructions& layer2) {
        return this->identifier < layer2.identifier;
    }
};


//! polyhedron factory
template <unsigned short DIM>
struct PolyhedronFactory {
    //! from the vertices
    //! \warning : the faces or segment should be oriented
    //!     - in 3D, the normal is such that the face is oriented in the trigonometric sense wrt it
    //!     - in 2D, the polygon is oriented in the trigonometric sense
    //! \param phase : phase of the polyhedron
    //! \param vertices : vertices, in absolute coordinates
    //! \param face_indices : array containing, for each face, the indices of the vertices contained in it
    //! \param minkowskiRadius : additional parameter for enlarging the polyhedron (due to coupling with rockable)
    ConvexPolyhedronInc<DIM> fromVertices(Identifier phase, const vector<Point<DIM>>& vertices, const vector<vector<long>>& face_indices);
};

//! spheroPolyhedron factory
template <unsigned short DIM>
struct SpheroPolyhedronFactory {
    //! from the vertices
    //! \warning : the faces or segment should be oriented
    //!     - in 3D, the normal is such that the face is oriented in the trigonometric sense wrt it
    //!     - in 2D, the polygon is oriented in the trigonometric sense
    //! \param phase : phase of the polyhedron
    //! \param vertices : vertices, in absolute coordinates
    //! \param face_indices : array containing, for each face, the indices of the vertices contained in it
    SpheroPolyhedronInc<DIM> fromVertices(Identifier phase, const vector<Point<DIM>>& vertices, const vector<vector<long>>& face_indices, double minkowskiRadius = 0);
};

namespace auxi_layerInstructions {
//! build a vector of layerInstruction
vector<LayerInstructions> buildInstructionVector(vector<Identifier> identifier, vector<PhaseType> phase, vector<double> width = {});

} // namespace auxi_layerInstructions

} // namespace smallShape
} // namespace merope

#include "../MicroInclusion/MicroInclusion.ixx"

#endif /* MICROINCLUSION_HXX_ */

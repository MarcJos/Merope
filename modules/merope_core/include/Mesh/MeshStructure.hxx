//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"

#include "../Mesh/GeoObjects.hxx"


namespace merope {

namespace mesh {
namespace meshStructure {
using namespace mesh::geoObjects;
using namespace mesh::sameThings;

template<unsigned short DIM>
class VoroMesh_Periodic;

//! raw data structures containing the basic information for building a mesh
template<unsigned short DIM>
struct VoroMesh_UnStructureData {
    // GeoObjects
    Point<DIM>              L;
    vector<GeoPoint<DIM>>   vecPoint;
    vector<Edge>            vecEdge;
    vector<CurveLoop>       vecCurveLoop;
    vector<Surface>         vecSurface;
    vector<SurfaceLoop>     vecSurfaceLoop;
    vector<Solid>           vecSolid;
    // periodic links
    vector<PerSurface<DIM>> vecPerSurface;
    vector<PerPoint> vecPerPoint;

    //! \return the maximal index of the objects of the structure
    Identifier getMaxIndex() const;
    //! shift the indices of all objects of the structure
    void shiftIndices(Identifier shift);
    //! apply a function on each vector
    template<class FUNCTION>
    void applyOnAllVectors(FUNCTION function);
    //! apply a function on each vector
    template<class FUNCTION>
    void applyOnAllVectors(FUNCTION function) const;
    //! copy from geoStructure
    void fromGeoPerStructure(const VoroMesh_Periodic<DIM>& geoStructure);
    //! reset all vectors of objects
    void reset();
    //! add the data of another structure
    void append(const VoroMesh_UnStructureData<DIM>& rawData);
    //! fixme
    void append_with_shift(VoroMesh_UnStructureData<DIM> rawData);
};


//! Contain the whole information of the mesh
//! std::map are used to expose the graph information
template<unsigned short DIM>
class VoroMesh_NotPeriodic {
public:
    // GeoObjects
    std::map<Identifier, GeoPoint<DIM>>     dictPoint;
    std::map<Identifier, Edge>              dictEdge;
    std::map<Identifier, CurveLoop>         dictCurveLoop;
    std::map<Identifier, Surface>           dictSurface;
    std::map<Identifier, SurfaceLoop>       dictSurfaceLoop;
    std::map<Identifier, Solid>             dictSolid;
    //! Ambiant space
    AmbiantSpace::Tore<DIM> torus;
public:
    //! constructor
    VoroMesh_NotPeriodic(const VoroMesh_UnStructureData<DIM>& rawData);
    //! destructor
    virtual ~VoroMesh_NotPeriodic() {}
    //! @brief add new raw mesh data
    void add_raw_mesh_data(const VoroMesh_UnStructureData<DIM>& rawData);
    //! verify if the structure is coherent
    bool isCoherent() const;
    //! prints into the stream f
    //! \param f: stream
    virtual void print(std::ostream& f) const;
    //! build the tree, namely connecting all the components
    void buildTree();
};


template<unsigned short DIM>
class VoroMesh_Periodic : public VoroMesh_NotPeriodic<DIM> {
public:
    //! to identify points copied by periodicity
    std::map<Identifier, PerPoint>          dictPerPoint;
    //! to identify surfaces copied by periodicity
    std::map<Identifier, PerSurface<DIM>>   dictPerSurface;
    //! constructor
    VoroMesh_Periodic(const VoroMesh_UnStructureData<DIM>& rawData, bool check_coherence);
    //! @brief : get the outer surface of the structure, assumed to be periodic
    geoObjects::PhysicalSurface getPeriodicOuterSurface(Identifier id) const;
    //! prints into the stream f
    //! \param f: stream
    void print(std::ostream& f) const override;
    //! restricts the polycrystal to its enveloppe only
    void restrictEnveloppe();
    //! \return the maximal index of the objects of the structure
    Identifier getMaxIndex() const;
    //! verify if the structure is strongly coherent
    //! in particular, each surface should be possessed by 2 volumes
    bool isStronglyCoherent() const;
protected:
    //! \return the identifiers of the surfaces the given point belong to
    //! \param pt_id : identifier of the given point
    std::set<Identifier> getSurfaces_from_Point(Identifier pt_id) const;
    //! \return the identifiers of the points delimitating a surface, in a order compatible with the orientation of the surface
    //! \param surf_id : identifier of a surface
    vector<Identifier> getPoints_from_Surface(Identifier surf_id) const;
private:
    //! criterion for identifier periodic surfaces
    //! magical constant!
    static constexpr double epsilon_per_surfaces = 1e-6;
    //! recover the periodic surfaces
    void buildPerSurface();
    //! put the correct order on periodic surfaces (compatible with gmsh, so that no "periodicity loop" exist)
    void orderPerSurface();
    //! remove all the singular components of the polyhedron mesh : (edges of 0 length, surface of 0 area, solids of 0 volume...)
    bool removeAllSingular();
    //! merge points, and compute the effects on the geometry
    bool mergeAll();
    //! Return true if surface_2 is the periodic copy of surface_1 with the given translation.
    //! the criterion is that they shall possess at least 3 points that are periodic copies
    //! \param surf_id1, surf_id2 : identifiers of the surface to compare
    //! \param translation : Expected translation in the Euclidean space to get surf_2 from sorf_1.
    bool comparePerSurface(Identifier surf_id1, Identifier surf_id2, const Point<DIM>& translation) const;
    //! verify the periodicity relationship
    //! fixme : very weak tests
    void verifyPeriodicity() const;
};

template<unsigned short DIM>
class VoroMesh_Periodic_Physical : public VoroMesh_Periodic<DIM> {
public:
    explicit VoroMesh_Periodic_Physical(VoroMesh_Periodic<DIM> voroMeshPer) :
        VoroMesh_Periodic<DIM>{ voroMeshPer }, dictPhysicalSurface({}), dictPhysicalVolume({}) {}
    //! to identify points copied by periodicity
    std::map<Identifier, PhysicalSurface>   dictPhysicalSurface;
    //! to identify surfaces copied by periodicity
    std::map<Identifier, PhysicalVolume>    dictPhysicalVolume;
};

//! translate a vector of GeoObjects into a map with identifier being its element
//! \param vec : a vector of GeoObjects
//! \param dict : the desired map of GeoObjects
template<class VEC, class DICT>
void translate(const VEC& vec, DICT& dict);
//! connect the leaves to their roots
//! \param dictThings : the map of leaves (to be modified)
//! \param dictRoot : the map of roots
template<class DICT_THINGS, class DICT_ROOT>
void connectRoot(DICT_THINGS& dictThings, const DICT_ROOT& dictRoots);
//! connect the leaves to their roots
//! \param dictThings : the map of leaves (to be modified)
//! \param dictPerRoot : the map of  periodic roots
template<class DICT_THINGS, class DICT_ROOT>
void connectPerRoot(DICT_THINGS& dictThings, const DICT_ROOT& dictPerRoots);
//! update the periodic components when merging elements
//! \param vecThings_id : well-posed vector containing the identifiers of things to be merged
//! \param dictThings : map of all things (to be modified)
//! \param dictPerThings : map of all corresponding periodic things (to be modified)
template<class DICT_THINGS, class DICT_PERIODIC_THINGS>
void updatePeriodicMerge(vector<SameThings<Identifier>> vecThings_id, const DICT_THINGS& dictThings, DICT_PERIODIC_THINGS& dictPerThings);
//! \return if a vector is a multiplier of the lengths of the torus
//! \param L : lengths of the torus
//! \param geomVector : vector
template<unsigned short DIM>
bool verifyTranslate(const AmbiantSpace::Tore<DIM>& torus, const Point<DIM>& translation, double epsilon);
//! erase all the leaves that are not connected to a given map of roots
//! \param leaf_map : the map of leaves (to be modified)
//! \param root_map : the map of roots
template<class DICT_LEAF, class DICT_ROOT>
void restrictTo_RootLeaves_withoutRootConnection(DICT_LEAF& dictThings, const DICT_ROOT& dictRoots);


}  // namespace meshStructure
}  // namespace mesh
}  // namespace merope

#include "../Mesh/MeshStructure.ixx"



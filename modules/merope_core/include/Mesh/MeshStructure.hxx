//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MESH_MESHSTRUCTURE_HXX_
#define MESH_MESHSTRUCTURE_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../Mesh/GeoObjects.hxx"

#include "../MeropeNamespace.hxx"


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
    vector<PerSurface<DIM>> vecPerSurface;

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
    void fromGeoPerStructure(VoroMesh_Periodic<DIM> geoStructure);
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
    VoroMesh_NotPeriodic(VoroMesh_UnStructureData<DIM> rawData);
    //! destructor
    virtual ~VoroMesh_NotPeriodic() {}
    //! verify if the structure is coherent
    bool isCoherent() const;
    //! prints into the stream f
    //! \param f: stream
    virtual void print(std::ostream& f) const;
protected:
    //! build the tree, namely connecting all the components
    void buildTree();
};

template<unsigned short DIM>
class VoroMesh_Periodic: public VoroMesh_NotPeriodic<DIM> {
public:
    //! to identify points copied by periodicity
    std::map<Identifier, PerPoint>          dictPerPoint;
    //! to identify surfaces copied by periodicity
    std::map<Identifier, PerSurface<DIM>>   dictPerSurface;
    //! constructor
    VoroMesh_Periodic(VoroMesh_UnStructureData<DIM> rawData, double adim_epsilon_0 = VoroMesh_Periodic::ADIM_EPSILON_0, double adim_epsilon_1 = VoroMesh_Periodic::ADIM_EPSILON_1);
    //! fixme not programmed yet
    geoObjects::PhysicalSurface getOuterSurface(Identifier id) const;
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
    //! adimensional distance to decide to merge two points, due to voro++ approximation errors
    static constexpr double ADIM_EPSILON_0 = 1.e-5;
    double epsilon_0;
    //! adimensional distance to decide to merge two points, for avoiding too small surfaces
    static constexpr double ADIM_EPSILON_1 = 1.e-5;
    double epsilon_1;
private:
    //! build the periodic structure
    //! fixme
    void buildPeriodicity(double epsilon_1);
    //! recovers the periodic points
    //! fixme
    void buildPerPoint(double epsilon_1_);
    //! recover the periodic surfaces
    void buildPerSurface();
    //! unify the Points that are too close to each others
    void unifyPoints();
    //! verify the periodicity relationship
    //! fixme : very weak tests
    //! fixme : strange way of swapping periodic surfaces
    void verifyPeriodicity();

    //! Merge points that are too close, but due to the Voronoi tessllation (no machine error)
    void removeClosePoints();
    //! remove all the singular components of the polyhedron mesh : (edges of 0 length, surface of 0 area, solids of 0 volume...)
    void removeAllSingular();
    //! merge two PerPoints
    //! \param samePoints : list of identifiers of periodic points that shall be merged. the merging process is supposed to be well-posed
    void mergePerPoints(const vector<SameThings<Identifier>>& samePoints);
    //! Merge two periodic points, namely perPt1 into perPt2
    //! For avoiding errors, a single reference point is chosen, and then, each of the physical point is aligned wrt to it
    //! \param perPt1, perPt2 : periodic points
    void mergeSinglePerPoint(const PerPoint& perPt1, const PerPoint& perPt2);
    //! align all the euclidean representants of a same periodic point
    //! \param perPoint_id : identifier of the periodic point
    void alignPerPoints(Identifier perPoint_id);
    //! merge points, and compute the effects on the geometry
    void mergeAll(vector<SameThings<Identifier>> listSamePoints);

    //! \param : closePointDoublePer = [[perPoint1, perPoint2, 1], ...] where perPoint1 and perPoint2 are close
    vector<SameThings<Identifier>> getClosePerPoints(const vector<SameThings<GeoPoint<DIM>>>& closePointDoublePer) const;
    //! \return points that should be merged
    vector<SameThings<Identifier>> tooClosePoints(double adimensionalDistance) const;
    //! \return : points that should be merged, taking the periodicity into account
    //! \warning : for some reasons, not exhaustive, unless the size of the result is 0
    tuple<vector<sameThings::SameThings<Identifier>>, vector<sameThings::SameThings<Identifier>>> periodicTooClosePoint(double adimensionalDistance) const;
    //! Return true if surface_2 is the periodic copy of surface_1 with the given translation.
    //! the criterion is that they shall possess at least 3 points that are periodic copies
    //! \param surf_id1, surf_id2 : identifiers of the surface to compare
    //! \param translation : Expected translation in the Euclidean space to get surf_2 from sorf_1.
    bool comparePerSurface(Identifier surf_id1, Identifier surf_id2, Point<DIM>& translation) const;
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


} // meshStructure
} // mesh
} // namespace merope

#include "../Mesh/MeshStructure.ixx"

#endif /* MESH_MESHSTRUCTURE_HXX_ */

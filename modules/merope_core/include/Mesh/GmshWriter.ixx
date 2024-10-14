//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../../AlgoPacking/src/AuxiFunctions.hxx"
#include "../Mesh/GeoObjects.hxx"



#include "../MeropeNamespace.hxx"


namespace merope {
namespace mesh {
namespace gmsh_writer {

template<class OBJ>
string name_of(const OBJ& object) {
    if constexpr (is_same_v<OBJ, geoObjects::Edge>) {
        if (object.typeEdge == geoObjects::TypeEdge::Segment)           return "Line";
        if (object.typeEdge == geoObjects::TypeEdge::Circle)            return "Circle";
        throw runtime_error("Unexpected!");
    } else if constexpr (is_same_v<OBJ, geoObjects::CurveLoop>)         return "Line Loop";
    else if constexpr (is_same_v<OBJ, geoObjects::Surface>) {
        if (object.typeSurface == geoObjects::TypeSurface::Plane)       return "Plane Surface";
        if (object.typeSurface == geoObjects::TypeSurface::Curved)      return "Surface";
        throw runtime_error("Unexpected!");
    } else if constexpr (is_same_v<OBJ, geoObjects::SurfaceLoop>)       return "Surface Loop";
    else if constexpr (is_same_v<OBJ, geoObjects::Solid>)               return "Volume";
    else if constexpr (is_same_v<OBJ, geoObjects::PhysicalSurface>)     return "Physical Surface";
    else if constexpr (is_same_v<OBJ, geoObjects::PhysicalVolume>)      return  "Physical Volume";
    throw runtime_error("Unexpected type!");
}

template<class OBJ>
inline void write(const OBJ& object, std::ostream& f) {
    if constexpr (is_same_v<OBJ, geoObjects::GeoPoint<3>>)                  auxi::write_point<3>(object, f);
    else if constexpr (is_same_v<OBJ, geoObjects::GeoPoint<2>>)             auxi::write_point<2>(object, f);
    else if constexpr (is_same_v<OBJ, geoObjects::PerPoint>)                auxi::write_simple(object, "PerPoint", f);
    else if constexpr (is_same_v<OBJ, geoObjects::Edge>
        or is_same_v<OBJ, geoObjects::CurveLoop>
        or is_same_v<OBJ, geoObjects::Surface>
        or is_same_v<OBJ, geoObjects::SurfaceLoop>
        or is_same_v<OBJ, geoObjects::Solid>
        or is_same_v<OBJ, geoObjects::PhysicalSurface>
        or is_same_v<OBJ, geoObjects::PhysicalVolume>)                      auxi::write_simple(object, gmsh_writer::name_of(object), f);
    else if constexpr (is_same_v<OBJ, geoObjects::PerSurface<3>>)           auxi::write_perSurface<3>(object, f);
    else if constexpr (is_same_v<OBJ, geoObjects::PerSurface<2>>)           auxi::write_perSurface<2>(object, f);
    else if constexpr (is_same_v<OBJ, meshStructure::VoroMesh_Periodic<3>>) auxi::write_geoPerStructure<3>(object, f);
    else if constexpr (is_same_v<OBJ, meshStructure::VoroMesh_Periodic<2>>) auxi::write_geoPerStructure<2>(object, f);
    else if constexpr (is_same_v<OBJ, Sphere<3>>)                           auxi::write_sphere<3>(object, f);
    else if constexpr (is_same_v<OBJ, Sphere<2>>)                           auxi::write_sphere<2>(object, f);
    else throw runtime_error("Unexpected type!");
}

/////////////////
// auxi
/////////////////

template<unsigned short DIM>
inline void auxi::write_geoPerStructure(const meshStructure::VoroMesh_Periodic<DIM>& geoPerStructure, std::ostream& f) {
    f.precision(17);
    auxi::writePreamble(f);
    auto write_loc = [&f](auto nameObject, const auto& dictThing) {
        auxi::writeDict(nameObject, dictThing, f);
        };
    write_loc("POINT", geoPerStructure.dictPoint);
    write_loc("LINE", geoPerStructure.dictEdge);
    write_loc("CURVE_LOOP", geoPerStructure.dictCurveLoop);
    write_loc("PLANE_SURFACE", geoPerStructure.dictSurface);
    write_loc("SURFACE_LOOP", geoPerStructure.dictSurfaceLoop);
    write_loc("VOLUME", geoPerStructure.dictSolid);
    write_loc("PERIODIC_SURFACE", geoPerStructure.dictPerSurface);
    write(geoPerStructure.getOuterSurface(1), f);
    auxi::writeEnd(f);
}

template<class DICT_THING>
inline void auxi::writeDict(string nameObject,
    const DICT_THING& dictThings, std::ostream& f) {
    std::vector<typename DICT_THING::mapped_type> vecThings;
    for (const auto& thing : dictThings) vecThings.push_back(thing.second);
    auxi::writeVect(nameObject, vecThings, f);
}

template<class VEC_THING>
inline void auxi::writeVect(string nameObject,
    VEC_THING vecThings, std::ostream& f) {
    std::sort(vecThings.begin(), vecThings.end(), [](const auto& e1, const auto& e2) {
        return e1.getId_forSort() < e2.getId_forSort();
        });
    f << "\n // BEGIN  " << nameObject << endl;
    for (const auto& x : vecThings) {
        write(x, f);
    }
    f << "// END  " << nameObject << endl << endl;
}

template<class OBJ>
inline void auxi::write_simple(const OBJ& object, std::string nameObject, std::ostream& f) {
    f << nameObject << "(" << object.identifier << ") = {";
    auxi_function::writeVectorToString(object.leaves, f);
    f << "};" << endl;
}

template<unsigned short DIM>
inline void auxi::write_point(const geoObjects::GeoPoint<DIM>& object, std::ostream& f) {
    f << "Point(" << object.identifier << ")={";
    auxi_function::writeVectorToString(object.coordinates, f);
    f << "," << mesh_size << "};" << endl;
}

template<unsigned short DIM>
inline void auxi::write_perSurface(const geoObjects::PerSurface<DIM>& object, std::ostream& f) {
    f << "Periodic Surface {" << object.leaves[0];
    f << " }={ " << object.leaves[1] << " } Translate {";
    auxi_function::writeVectorToString(object.translation, f);
    f << "};" << endl;
}

template<unsigned short DIM>
inline void auxi::write_sphere(Sphere<DIM> sphere, std::ostream& f, vector<Identifier> removedVolumes) {
    if constexpr (DIM == 2) throw runtime_error("mesh::gmsh_writer::auxi::write_sphere not programmed for DIM = 2");
    // inspired from https://sites.google.com/site/auxcapucins/maillage-3d-en-gmsh---maillage-d-une-sphere
    f << "//BEGIN OF SPHERE \n";
    Identifier index = sphere.phase;
    //
    auto translateIndex = [](auto decalage, auto i) {
        if (i < 0) return i - decalage;
        else return i + decalage;
        };
    //
    auto easyCoord = [&sphere](const auto& i, const auto& j, const auto& k) {
        return Point<DIM>({ sphere.center[0] + i * sphere.radius, sphere.center[1] + j * sphere.radius, sphere.center[2] + k * sphere.radius });
        };
    //
    auto writePt = [&index, &translateIndex, &easyCoord, &f](Identifier indexLoc, long i, long j, long k) {
        mesh::geoObjects::GeoPoint<DIM> pt(translateIndex(index, indexLoc), easyCoord(i, j, k));
        write(pt, f);
        };
    writePt(1, 0, 0, 0);
    writePt(2, 1, 0, 0);
    writePt(3, 0, 1, 0);
    writePt(4, 0, 0, 1);
    writePt(5, -1, 0, 0);
    writePt(6, 0, -1, 0);
    writePt(7, 0, 0, -1);
    //
    auto writeCircle = [&index, &translateIndex, &f](Identifier indexLoc, long i, long j, long k) {
        f << "Circle(" << translateIndex(index, indexLoc) << ")={"
            << translateIndex(index, i) << "," << translateIndex(index, j)
            << "," << translateIndex(index, k) << "};\n";
        };
    //
    writeCircle(1, 2, 1, 3);
    writeCircle(2, 3, 1, 5);
    writeCircle(3, 5, 1, 6);
    writeCircle(4, 6, 1, 2);
    writeCircle(5, 2, 1, 7);
    writeCircle(6, 7, 1, 5);
    writeCircle(7, 5, 1, 4);
    writeCircle(8, 4, 1, 2);
    writeCircle(9, 6, 1, 7);
    writeCircle(10, 7, 1, 3);
    writeCircle(11, 3, 1, 4);
    writeCircle(12, 4, 1, 6);
    //
    auto writeLineLoop = [&index, &translateIndex, &f](Identifier indexLoc, long i, long j, long k) {
        f << "Line Loop(" << translateIndex(index, indexLoc) << ") = {"
            << translateIndex(index, i) << "," << translateIndex(index, j)
            << "," << translateIndex(index, k) << "};\n";
        };
    //
    writeLineLoop(1, 1, 11, 8);
    writeLineLoop(2, 2, 7, -11);
    writeLineLoop(3, 3, -12, -7);
    writeLineLoop(4, 4, -8, 12);
    writeLineLoop(5, 5, 10, -1);
    writeLineLoop(6, -2, -10, 6);
    writeLineLoop(7, -3, -6, -9);
    writeLineLoop(8, -4, 9, -5);
    //
    auto writeRuleSurface = [&index, &translateIndex, &f](Identifier indexLoc) {
        f << "Surface(" << translateIndex(index, indexLoc) << ")={" << translateIndex(index, indexLoc) << "};\n";
        };
    //
    vector<long> identifiers{};
    for (long i = 1; i < 9; i++) {
        writeRuleSurface(i);
        identifiers.push_back(translateIndex(index, i));
    }
    //
    f << "Surface Loop(" << index << ")={";
    auxi_function::writeVectorToString(identifiers, f);
    f << "};\n";
    //
    if (removedVolumes.size() == 0) {
        f << "Volume(" << index << ")={" << index << "};\n";
    } else {
        f << "Volume(" << index << ")={" << index << ",";
        auxi_function::writeVectorToString(removedVolumes, f);
        f << "};\n";
    }
    //
    f << "//END OF SPHERE \n";
}

}  // namespace gmsh_writer
}  // namespace mesh
}  // namespace merope




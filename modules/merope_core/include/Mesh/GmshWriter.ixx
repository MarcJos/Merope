//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../GenericTools/AuxiFunctions.hxx"

#include "../Mesh/GeoObjects.hxx"


namespace merope {
namespace mesh {
namespace gmsh_writer {

template<class OBJ>
string name_of(const OBJ& object) {
    if constexpr (is_same_v<OBJ, geoObjects::Edge>) {
        if (object.typeEdge == geoObjects::TypeEdge::Segment)           return "Line";
        if (object.typeEdge == geoObjects::TypeEdge::Circle)            return "Circle";
        Merope_error_impossible();
    } else if constexpr (is_same_v<OBJ, geoObjects::CurveLoop>)         return "Line Loop";
    else if constexpr (is_same_v<OBJ, geoObjects::Surface>) {
        if (object.typeSurface == geoObjects::TypeSurface::Plane)       return "Plane Surface";
        if (object.typeSurface == geoObjects::TypeSurface::Curved)      return "Surface";
        Merope_error_impossible();
    } else if constexpr (is_same_v<OBJ, geoObjects::SurfaceLoop>)       return "Surface Loop";
    else if constexpr (is_same_v<OBJ, geoObjects::Solid>)               return "Volume";
    else if constexpr (is_same_v<OBJ, geoObjects::PhysicalSurface>)     return "Physical Surface";
    else if constexpr (is_same_v<OBJ, geoObjects::PhysicalVolume>)      return  "Physical Volume";
    else Merope_static_error(OBJ, "Unexpected type!");
}

template<class OBJ>
string name_new_(const OBJ& object) {
    if constexpr (is_same_v<OBJ, geoObjects::GeoPoint<3>>)              return "newp";
    else if constexpr (is_same_v<OBJ, geoObjects::Edge>)                return "newc";
    else if constexpr (is_same_v<OBJ, geoObjects::CurveLoop>)           return "newcl";
    else if constexpr (is_same_v<OBJ, geoObjects::Surface>)             return "news";
    else if constexpr (is_same_v<OBJ, geoObjects::SurfaceLoop>)         return "newsl";
    else if constexpr (is_same_v<OBJ, geoObjects::Solid>)               return "newv";
    else Merope_static_error(OBJ, "Unexpected type!");
}

template<class OBJ>
string letter_for(const OBJ& object) {
    if constexpr (is_same_v<OBJ, geoObjects::GeoPoint<3>>)               return "p";
    else if constexpr (is_same_v<OBJ, geoObjects::Edge>)                return "c";
    else if constexpr (is_same_v<OBJ, geoObjects::CurveLoop>)           return "cl";
    else if constexpr (is_same_v<OBJ, geoObjects::Surface>)             return "s";
    else if constexpr (is_same_v<OBJ, geoObjects::SurfaceLoop>)         return "sl";
    else if constexpr (is_same_v<OBJ, geoObjects::Solid>)               return "v";
    else Merope_static_error(OBJ, "Unexpected type!");
}

template<class OBJ>
string letter_for_leaf(const OBJ& object) {
    if constexpr (is_same_v<OBJ, geoObjects::Edge>)                     return "p";
    else if constexpr (is_same_v<OBJ, geoObjects::CurveLoop>)           return "c";
    else if constexpr (is_same_v<OBJ, geoObjects::Surface>)             return "cl";
    else if constexpr (is_same_v<OBJ, geoObjects::SurfaceLoop>)         return "s";
    else if constexpr (is_same_v<OBJ, geoObjects::Solid>)               return "sl";
    else if constexpr (is_same_v<OBJ, geoObjects::PhysicalSurface>)     return "s";
    else if constexpr (is_same_v<OBJ, geoObjects::PerSurface<3>>)       return "s";
    else if constexpr (is_same_v<OBJ, geoObjects::PhysicalVolume>)      return "v";
    else Merope_static_error(OBJ, "Unexpected type!");
}

template<class OBJ>
bool does_preserve_sign_for_leaf(const OBJ& object) {
    if constexpr (is_same_v<OBJ, geoObjects::SurfaceLoop>)          return false;
    else                                                            return true;
}

template<MeshMethod meshMethod, class OBJ>
void write(const OBJ& object, std::ostream& f) {
    if constexpr (is_same_v<OBJ, geoObjects::GeoPoint<3>>)                  auxi::write_point<3, meshMethod>(object, f);
    else if constexpr (is_same_v<OBJ, geoObjects::GeoPoint<2>>)             auxi::write_point<2, meshMethod>(object, f);
    else if constexpr (is_same_v<OBJ, geoObjects::PerPoint>)                auxi::write_simple<meshMethod>(object, "//PerPoint", f);
    else if constexpr (is_same_v<OBJ, geoObjects::Edge>
        or is_same_v<OBJ, geoObjects::CurveLoop>
        or is_same_v<OBJ, geoObjects::Surface>
        or is_same_v<OBJ, geoObjects::SurfaceLoop>
        or is_same_v<OBJ, geoObjects::Solid>
        or is_same_v<OBJ, geoObjects::PhysicalSurface>
        or is_same_v<OBJ, geoObjects::PhysicalVolume>)                      auxi::write_simple<meshMethod>(object, gmsh_writer::name_of(object), f);
    else if constexpr (is_same_v<OBJ, geoObjects::PerSurface<3>>)           auxi::write_perSurface<3, meshMethod>(object, f);
    else if constexpr (is_same_v<OBJ, geoObjects::PerSurface<2>>)           auxi::write_perSurface<2, meshMethod>(object, f);
    else if constexpr (is_same_v<OBJ, meshStructure::VoroMesh_Periodic<3>>) auxi::write_geoPerStructure<3, meshMethod>(object, f);
    else if constexpr (is_same_v<OBJ, meshStructure::VoroMesh_Periodic<2>>) auxi::write_geoPerStructure<2, meshMethod>(object, f);
    else if constexpr (is_same_v<OBJ, Sphere<3>>)                           auxi::write_sphere<3, meshMethod>(object, f);
    else if constexpr (is_same_v<OBJ, Sphere<2>>)                           auxi::write_sphere<2, meshMethod>(object, f);
    else Merope_static_error(OBJ, "Unexpected type!");
}

/////////////////
// auxi
/////////////////

template<unsigned short DIM, MeshMethod meshMethod>
void auxi::write_geoPerStructure(const meshStructure::VoroMesh_Periodic<DIM>& geoPerStructure, std::ostream& f) {
    f.precision(17);
    auxi::writePreamble<meshMethod>(f);
    auto write_loc = [&f](auto nameObject, const auto& dictThing) {
        auxi::writeDict<meshMethod>(nameObject, dictThing, f);
        };
    write_loc("POINT", geoPerStructure.dictPoint);
    write_loc("LINE", geoPerStructure.dictEdge);
    write_loc("CURVE_LOOP", geoPerStructure.dictCurveLoop);
    write_loc("PLANE_SURFACE", geoPerStructure.dictSurface);
    write_loc("SURFACE_LOOP", geoPerStructure.dictSurfaceLoop);
    write_loc("VOLUME", geoPerStructure.dictSolid);
    write_loc("PERIODIC_SURFACE", geoPerStructure.dictPerSurface);
    write<meshMethod>(geoPerStructure.getPeriodicOuterSurface(1), f);
    auxi::writeEnd(f);
}

template<MeshMethod meshMethod, class DICT_THING>
void auxi::writeDict(string nameObject,
    const DICT_THING& dictThings, std::ostream& f) {
    std::vector<typename DICT_THING::mapped_type> vecThings;
    for (const auto& thing : dictThings) vecThings.push_back(thing.second);
    auxi::writeVect<meshMethod>(nameObject, vecThings, f);
}

template<MeshMethod meshMethod, class VEC_THING>
void auxi::writeVect(string nameObject,
    VEC_THING vecThings, std::ostream& f) {
    std::sort(vecThings.begin(), vecThings.end(), [](const auto& e1, const auto& e2) {
        return e1.getId_forSort() < e2.getId_forSort();
        });
    f << "\n // BEGIN  " << nameObject << endl;
    for (const auto& x : vecThings) {
        write<meshMethod>(x, f);
    }
    f << "// END  " << nameObject << endl << endl;
}

template<MeshMethod meshMethod, class OBJ>
void auxi::write_simple(const OBJ& object, std::string nameObject, std::ostream& f) {
    if constexpr (meshMethod == MeshMethod::native_gmsh) {
        f << nameObject << "(" << object.identifier << ") = {";
        auxi_function::writeVectorToString(object.leaves, f);
        f << "};" << endl;
    } else if constexpr (meshMethod == MeshMethod::OpenCascade) {
        string name_p = to_string(object.identifier);
        if constexpr (is_same_v<OBJ, geoObjects::Edge>
            or is_same_v<OBJ, geoObjects::CurveLoop>
            or is_same_v<OBJ, geoObjects::Surface>
            or is_same_v<OBJ, geoObjects::SurfaceLoop>
            or is_same_v<OBJ, geoObjects::Solid>) {
            name_p = letter_for(object) + to_string(object.identifier);
            f << name_p << " = " << name_new_(object) << "; ";
        }
        //
        f << gmsh_writer::name_of(object) << "(" << name_p << ") = {";
        auxi_function::writeVectorToString_with_prefix(object.leaves, f, letter_for_leaf(object), does_preserve_sign_for_leaf(object));
        f << "};" << endl;
    } else {
        Merope_assert(false, "meshMethod incorrect");
    }
}

template<unsigned short DIM, MeshMethod meshMethod>
void auxi::write_point(const geoObjects::GeoPoint<DIM>& object, std::ostream& f) {
    if constexpr (meshMethod == MeshMethod::native_gmsh) {
        f << "Point(" << object.identifier << ")={";
        auxi_function::writeVectorToString(object.coordinates, f);
        f << "," << mesh_size() << "};" << endl;
    } else if constexpr (meshMethod == MeshMethod::OpenCascade) {
        auto name_p = letter_for(object) + to_string(object.identifier);
        f << name_p << " = " << name_new_(object) << "; ";
        f << "Point(" << name_p << ")={";
        auxi_function::writeVectorToString(object.coordinates, f);
        f << "," << mesh_size() << "};" << endl;
    } else {
        Merope_assert(false, "meshMethod incorrect");
    }
}

template<unsigned short DIM, MeshMethod meshMethod>
void auxi::write_perSurface(const geoObjects::PerSurface<DIM>& object, std::ostream& f) {
    f << "Periodic Surface";
    if constexpr (meshMethod == MeshMethod::native_gmsh) {
        f << " {" << object.leaves[0] << " }={ " << object.leaves[1] << " }";
    } else if constexpr (meshMethod == MeshMethod::OpenCascade) {
        f << " { " << letter_for_leaf(object) << object.leaves[0] << " }={ " << letter_for_leaf(object) << object.leaves[1] << " }";
    } else {
        Merope_assert(false, "meshMethod incorrect");
    }
    f << " Translate {";
    auxi_function::writeVectorToString(object.translation, f);
    f << "};" << endl;
}

template<unsigned short DIM, MeshMethod meshMethod>
void auxi::write_sphere(Sphere<DIM> sphere, std::ostream& f, vector<Identifier> removedVolumes) {
    if constexpr (DIM == 2) Merope_error_not_done();
    // inspired from https://sites.google.com/site/auxcapucins/maillage-3d-en-gmsh---maillage-d-une-sphere
    if constexpr (meshMethod == MeshMethod::native_gmsh) {
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
    } else if constexpr (meshMethod == MeshMethod::OpenCascade) {
        string name_s = "v" + to_string(sphere.phase);
        f << name_s << "= newv; ";
        f << "Sphere(" << name_s << ") = ";
        f << "{" << sphere.center[0] << ", " << sphere.center[1] << ", " << sphere.center[2] << ", " << sphere.radius << "};\n";

    }
}

template<MeshMethod meshMethod>
void auxi::writePreamble(std::ostream& f) {
    f.precision(17);
    f << "// File automatically generated by Merope \n";
    if constexpr (meshMethod == MeshMethod::OpenCascade) {
        f << R"(SetFactory("OpenCASCADE");)" << "\n";
        f << "// Ask OpenCASCADE to compute more accurate bounding boxes of entities using the \n";
        f << "// STL mesh: \n";
        f << "Geometry.OCCBoundsUseStl = 1;\n";
    } else if constexpr (meshMethod == MeshMethod::native_gmsh) {
        f << R"(//SetFactory("OpenCASCADE");)" << "\n";
        f << "// Ask OpenCASCADE to compute more accurate bounding boxes of entities using the \n";
        f << "// STL mesh: \n";
        f << "//Geometry.OCCBoundsUseStl = 1;\n";
    } else {
        Merope_assert(false, "Incorrect meshMethod");
    }
    f << mesh_size() << " = " << 0 << ";  // WARNING : Default value. Correct value is at the end of the script.\n";
}

}  // namespace gmsh_writer
}  // namespace mesh
}  // namespace merope




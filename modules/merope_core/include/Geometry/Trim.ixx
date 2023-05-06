//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MEROPE_CORE_SRC_GEOMETRY_TRIM_IXX_
#define MEROPE_CORE_SRC_GEOMETRY_TRIM_IXX_


template<unsigned short DIM, class SOLID>
inline SOLID sac_de_billes::geomTools::trim(const SOLID& solid,
    double layerWidth) {
    // case 1 : sphere
    if constexpr (std::is_same<SOLID, Sphere<DIM>>::value) {
        return Sphere<DIM>(solid.center, solid.radius - layerWidth, 0);
    }
    //
    else if constexpr (std::is_same<SOLID, ConvexPolyhedron<DIM>>::value) {
        vector<HalfSpace<DIM>> faces = solid.faces;
        for (auto& hs : faces) {
            hs.c() -= layerWidth;
        }
        return ConvexPolyhedron<DIM>(solid.center, faces);
    }
    //
    else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Not programmed yet");
    }
}


#endif /* MEROPE_CORE_SRC_GEOMETRY_TRIM_IXX_ */

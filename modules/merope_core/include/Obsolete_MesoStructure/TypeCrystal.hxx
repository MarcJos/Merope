//! Copyright : see license.txt
//!
//! \brief 
//!
//
#ifndef TYPECRYSTAL_HXX_
#define TYPECRYSTAL_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

//! This class enumerates the possible types of crystal
enum class TypeCrystal {
    Voronoi, Laguerre, JohnsonMehl, Spheres
};

inline TypeCrystal str_2_TypeCrystal(const string s) {
    TypeCrystal typeCrystal;
    if (s == "Voronoi" or s == "voronoi") {
        typeCrystal = TypeCrystal::Voronoi;
    }
    else if (s == "Laguerre" or s == "laguerre") {
        typeCrystal = TypeCrystal::Laguerre;
    }
    else if (s == "JohnsonMehl" or s == "johnsonmehl") {
        typeCrystal = TypeCrystal::JohnsonMehl;
    }
    else if (s == "Inclusions" or s == "inclusions" or s == "Inclusion"
        or s == "inclusion") {
        typeCrystal = TypeCrystal::Spheres;
    }
    else {
        throw runtime_error("PolyCrystal::setTypeCrystal, Unknown crystal");
    }
    return typeCrystal;
}

//! This class parametrizes the way of building the microstructure
enum class ColorMaterialID {
    //! each moncrystal/inclusion has it own color
    Poly,
    //! Poly + there is an interface between crystals, which has its own color
    Erode,
    //! ==Erode, but each monocrystal has the same color
    Erode2Mat,
    //! == Erode2Mat, but an additional collection of spheres intersects the interface between crystals, giving rise to a new MaterialId
    Erode3Mat
};

inline ColorMaterialID str_2_ColorMaterialID(const string s) {
    if (s == "Poly") {
        return ColorMaterialID::Poly;
    }
    else if (s == "Erode") {
        return ColorMaterialID::Erode;
    }
    else if (s == "Erode2Mat") {
        return ColorMaterialID::Erode2Mat;
    }
    else if (s == "Erode3Mat") {
        return ColorMaterialID::Erode3Mat;
    }
    else {
        throw invalid_argument("FromString_2_ColorMaterialID : invalid string");
    }
}

} // namespace merope


#endif /* TYPECRYSTAL_HXX_ */

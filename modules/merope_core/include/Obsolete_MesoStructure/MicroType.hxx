//! Copyright : see license.txt
//!
//! \brief
//!

#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "TypeCrystal.hxx"


namespace merope {

class MicroType {
public:
    MicroType() : typeCrystal{ TypeCrystal::Voronoi } {};
    //! type of crystal
    TypeCrystal typeCrystal;

    //! sets the crystal type
    inline void setTypeCrystal(TypeCrystal tc) {
        typeCrystal = tc;
    }
};

class Colorize {
public:
    Colorize() : colorization{ ColorMaterialID::Poly } {};
    //! Says how to set MaterialId
    ColorMaterialID colorization;
    //!
    inline void setColorization(ColorMaterialID color) {
        colorization = color;
    }
};

class MicroType_Ext : public Colorize, public MicroType {
public:
    using MicroType::setTypeCrystal;
    using Colorize::setColorization;

    inline void setTypeCrystal(const string& s) {
        setTypeCrystal(str_2_TypeCrystal(s));
    }
    inline void setColorization(const string& color) {
        setColorization(str_2_ColorMaterialID(color));
    }
};

}  // namespace merope




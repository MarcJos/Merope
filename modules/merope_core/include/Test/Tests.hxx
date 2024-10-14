//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

class Tests {
public:
    static void polyCrystal11();
    static void polyCrystal10();
    static void polyCrystal0();
    static void polyCrystal1();
    static void polyCrystal_fit_volumes_3D(string fileVTK);
    static void polyCrystal_fit_volumes_2D(string fileVTK);
    static void polyCrystal8();
    static void polyCrystal9();
    static void spheres1();
    static void extraction();
    static void testPerf0();
    static void testFields();
    static void testFields2();
    static void polyCrystalCentroidal(bool use_acceleration);
    static void outputLaminate();
};

}  // namespace merope



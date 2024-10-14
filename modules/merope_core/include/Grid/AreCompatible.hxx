//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../MesoStructure/InsideTorus.hxx"
#include "../Grid/PreGrid.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM>
bool areCompatible(const InsideTorus<DIM>& i1, const InsideTorus<DIM>& i2) {
    return (sqrt(geomTools::distanceCarre<DIM>(i1.getL(), i2.getL()) / geomTools::norme<DIM>(i1.getL())) < geomTools::EPSILON);
}

template<unsigned short DIM>
bool areCompatible(const vox::GridParameters<DIM>& pg1, const vox::GridParameters<DIM>& pg2) {
    bool result = (sqrt(geomTools::distanceCarre<DIM>(pg1.getL(), pg2.getL()) / geomTools::norme<DIM>(pg1.getL())) < geomTools::EPSILON);
    for (size_t i = 0; i < DIM; i++) {
        result = result and (pg1.getNbNodes()[i] == pg2.getNbNodes()[i]);
        result = result and (pg1.getNMax()[i] == pg2.getNMax()[i]);
        result = result and (pg1.getNMin()[i] == pg2.getNMin()[i]);
    }
    return result;
}

}  // namespace merope



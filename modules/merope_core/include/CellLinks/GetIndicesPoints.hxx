//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"


template<unsigned short DIM>
//! @brief class for computing indices from points inside a shape.
//! these indices are used for comparing nearby points
struct Helper_for_getIndicesPoints {
    //! @brief 
    //! @param torus_ 
    //! @param nbPoints 
    //! @param minDistance 
    Helper_for_getIndicesPoints(const AmbiantSpace::Tore<DIM> torus_, size_t nbPoints, double minDistance);
    vector<Identifier> getIndices(Point<DIM> x);
private:
    double l;
    DiscPoint<DIM> nbNodes;
    const AmbiantSpace::Tore<DIM> torus;
};

#include "GetIndicesPoints.ixx"
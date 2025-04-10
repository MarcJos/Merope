//! Copyright : see license.txt
//!
//! \briefFor Test: a 2 phases geometry, 
//! with a rectangle included in a matrix
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../MesoStructure/InsideTorus.hxx"


namespace merope {

//! Periodical 2 phases geometry with a rectangle included in a Matrix.
template<unsigned short DIM>
class Rectangle : public InsideTorus<DIM> {
public:
    //! Constructors
    //! \param L : Periodical geometry dimensions
    //! \param rec : 'rectangle' inclusion dimensions
    //! constructor.
    Rectangle(array<double, DIM> L, array<double, DIM> recL_) : InsideTorus<DIM>(L), recL{ recL_ }{};
    //! Number of phases in the geometry.
    size_t getNumberOfPhases() const;
    //! Rectangular inclusion dimension in x direction.
    array<double, DIM> recL;
};

template<unsigned short DIM>
size_t Rectangle<DIM>::getNumberOfPhases() const {
    return 2;
}

}  // namespace merope





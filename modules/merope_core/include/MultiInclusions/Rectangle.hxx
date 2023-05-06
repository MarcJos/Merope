//! Copyright : see license.txt
//!
//! \brief For Test: a 2 phases geometry, 
//! with a rectangle included in a matrix
//
#ifndef _RECTANGLE_HXX
#define _RECTANGLE_HXX 1


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../MesoStructure/InsideTorus.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {

//! Periodical 2 phases geometry with a rectangle included in a Matrix.
template<unsigned short DIM>
class Rectangle: public InsideTorus<DIM> {
public:
    //! Constructors
    //! \param L : Periodical geometry dimensions
    //! \param rec : 'rectangle' inclusion dimensions
    //! constructor.
    Rectangle(array<double, DIM> L, array<double, DIM> recL_): InsideTorus<DIM>(L), recL{ recL_ }{};
    //! Number of phases in the geometry.
    size_t getNumberOfPhases() const;
    //! Rectangular inclusion dimension in x direction.
    array<double, DIM> recL;
};

template<unsigned short DIM>
size_t Rectangle<DIM>::getNumberOfPhases() const {
    return 2;
}

} // namespace merope


#endif // _RECTANGLE_HXX


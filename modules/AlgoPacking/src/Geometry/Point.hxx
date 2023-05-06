//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef POINT_HXX_
#define POINT_HXX_

#include "../Geometry/GeomTypes.hxx"

namespace sac_de_billes {
using namespace std;


//! creates an array initialized with the same x
template<unsigned short DIM, class C = double>
constexpr array<C, DIM> create_array(C x);


/* TOO MUCH CHANGES FOR INTERFACING WITH PYTHON
template<unsigned short DIM>
class Point: public array<double, DIM>{
public:
    //! @brief default constructor
    Point(): array<double, DIM>(create_array<DIM>(0.)) {}
    //! @brief constructor from an array
    //! @param pt
    Point(const array<double, DIM>& pt):array<double, DIM>(pt) {}
    Point(array<double, DIM>&& pt):array<double, DIM>(std::move(pt)) {}
    //! @brief constructor from a normalized point
    Point(const RenormPoint<DIM>&);
    //! @brief constructor from list of double
    template<unsigned short DIM2 = DIM, typename = std::enable_if_t<DIM2 == 3>>
    Point(double x0, double x1, double x2): array<double, DIM>{x0, x1, x2} {}
    template<unsigned short DIM2 = DIM, typename = std::enable_if_t<DIM2 == 2>>
    Point(double x0, double x1) : array<double, DIM>{x0, x1} {}
    template<unsigned short DIM2 = DIM, typename = std::enable_if_t<DIM2 == 1>>
    Point(double x0) : array<double, DIM>{x0} {}
};
*/


// for defining all operators +=, -=, *=, /=, on Points
// parameter binop : binary operator +, -, *, /
#define MEROPE_IMPLEMENT_POINT_BINOPEG(binop) \
template<class T, typename std::enable_if<is_Point<T, 3>, bool>::type = true > \
inline Point<3>& operator binop ## = (Point<3>& x, const T& y) { x[0] binop ## = y[0]; x[1] binop ## = y[1]; x[2] binop ## = y[2]; return x; } \
template<class T, typename std::enable_if<is_Point<T, 2>, bool>::type = true > \
inline Point<2>& operator binop ## = (Point<2>& x, const T& y) { x[0] binop ## = y[0]; x[1] binop ## = y[1]; return x; }\
\
inline Point<3>& operator binop ## =(Point<3>& x, double a) { x[0] binop ## = a; x[1] binop ## = a; x[2] binop ## = a; return x; }\
inline Point<2>& operator binop ## =(Point<2>& x, double a) { x[0] binop ## = a; x[1] binop ## = a; return x; }\
\

// for defining all operators +, -, *, / for points
// parameter binop : binary operator +, -, *, /
// also takes advantage of move without using new memory
#define MEROPE_IMPLEMENT_POINT_BINOP(DIM, binop) \
template<class T1, class T2, typename std::enable_if<is_Point<T1, DIM> and is_Point<T2, DIM>, bool>::type = true > \
inline Point<DIM> operator binop(const T1& x, const T2& y) { Point<DIM> z = x; z binop ## = y; return z; } \
\
template<class T1, typename std::enable_if<is_Point<T1, 3>, bool>::type = true > \
inline Point<DIM> operator binop(Point<DIM>&& x, const T1& y) {Point<DIM> z = std::move(x); z binop ## = y; return z; } \
\

// for defining all interesting operators
#define MEROPE_IMPLEMENT_POINT_ALL_OP(binop) \
MEROPE_IMPLEMENT_POINT_BINOPEG(binop)\
MEROPE_IMPLEMENT_POINT_BINOP(3, binop)\
MEROPE_IMPLEMENT_POINT_BINOP(2, binop)\
\

// operator unary - and multiplication by real
#define MEROPE_IMPLEMENT_POINT_SPECOP(DIM)\
inline Point<DIM> operator-(Point<DIM>&& x) {Point<DIM> z = std::move(x); z *= (-1); return z; }\
\
template<class T, typename std::enable_if<is_Point<T, DIM>, bool>::type = true >\
inline Point<DIM> operator-(const T& x) { Point<DIM> z = x; z *= (-1); return z; }\
\
inline Point<DIM> operator*(double a, Point<DIM>&& x) {Point<DIM> z = std::move(x); z *= a; return z; }\
\
inline Point<DIM> operator*(double a, const Point<DIM>& x) {Point<DIM> z = x; z *= a; return z; }\
\
template<class T, typename std::enable_if<is_Point<T, DIM>, bool>::type = true >\
inline Point<DIM> operator*(double a, const T& x) { Point<DIM> z = x.getPoint(); z *= a; return z; }\
\

//
MEROPE_IMPLEMENT_POINT_ALL_OP(+)
MEROPE_IMPLEMENT_POINT_ALL_OP(-)
MEROPE_IMPLEMENT_POINT_ALL_OP(*)
MEROPE_IMPLEMENT_POINT_ALL_OP(/ )
MEROPE_IMPLEMENT_POINT_SPECOP(3)
MEROPE_IMPLEMENT_POINT_SPECOP(2)


//
template<unsigned short DIM>
Point<DIM> average(const vector<Point<DIM>>& vPoints);

// for defining all operators +=, -=, *=, /=, on Points
// parameter binop : binary operator +, -, *, /
#define MEROPE_IMPLEMENT_DISCPOINT_BINOPEG(binop) \
inline DiscPoint<3>& operator binop ##=(DiscPoint<3>& x, const DiscPoint<3>& y) { x[0] binop ##= y[0]; x[1] binop ##= y[1]; x[2] binop ##= y[2]; return x; } \
inline DiscPoint<2>& operator binop ##=(DiscPoint<2>& x, const DiscPoint<2>& y) { x[0] binop ##= y[0]; x[1] binop ##= y[1]; return x; }\
\
inline DiscPoint<3>& operator binop ## =(DiscPoint<3>& x, long a) { x[0] binop ## = a; x[1] binop ## = a; x[2] binop ## = a; return x; }\
inline DiscPoint<2>& operator binop ## =(DiscPoint<2>& x, long a) { x[0] binop ## = a; x[1] binop ## = a; return x; }\
\

// for defining all operators +, -, *, / for points
// parameter binop : binary operator +, -, *, /
// also takes advantage of move without using new memory
#define MEROPE_IMPLEMENT_DISCPOINT_BINOP(DIM, binop) \
inline DiscPoint<DIM> operator binop(const DiscPoint<DIM>& x, const DiscPoint<DIM>& y) { DiscPoint<DIM> z = x; z binop ## = y; return z; } \
\
inline DiscPoint<DIM> operator binop(DiscPoint<DIM>&& x, const DiscPoint<DIM>& y) { DiscPoint<DIM> z = std::move(x); z binop ## = y; return z; }

// for defining all interesting operators

#define MEROPE_IMPLEMENT_DISCPOINT_ALL_OP(binop) \
MEROPE_IMPLEMENT_DISCPOINT_BINOPEG(binop)\
MEROPE_IMPLEMENT_DISCPOINT_BINOP(3, binop)\
MEROPE_IMPLEMENT_DISCPOINT_BINOP(2, binop)\

// all operators +, -, *, +=, -=, *=. Beware, no operator /=
MEROPE_IMPLEMENT_DISCPOINT_ALL_OP(+)
MEROPE_IMPLEMENT_DISCPOINT_ALL_OP(-)
MEROPE_IMPLEMENT_DISCPOINT_ALL_OP(*)
// unary operators
inline DiscPoint<3> operator-(const DiscPoint<3>& x) { DiscPoint<3> z = x; z *= (-1); return z; }
inline DiscPoint<2> operator-(const DiscPoint<2>& x) { DiscPoint<2> z = x; z *= (-1); return z; }

//
template<unsigned short DIM>
static ostream& print(ostream& out, const Point<DIM>& p);

template<unsigned short DIM>
ostream& operator<<(ostream& out, const Point<DIM>& p);

} // namespace sac_de_billes

#include "../Geometry/Point.ixx"

#endif /* POINT_HXX_ */

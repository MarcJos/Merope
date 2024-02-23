//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef POINT_HXX_
#define POINT_HXX_

#include "../Geometry/GeomTypes.hxx"
#include "../AuxiFunctions.hxx"

namespace sac_de_billes {
using namespace std;


//! creates an array initialized with the same x
template<unsigned short DIM, class C = double>
constexpr array<C, DIM> create_array(C x);

// for defining all operators +=, -=, *=, /=, on Points
// parameter binop : binary operator +, -, *, /
#define MEROPE_IMPLEMENT_POINT_BINOPEG(DIM, binop) \
template<class T, typename std::enable_if<is_Point<T, DIM>, bool>::type = true > \
inline Point<DIM>& operator binop ## = (Point<DIM>& x, const T& y) \
{ sac_de_billes::auxi_function::for_constexpr<0,DIM>([&](int i){x[i] binop ## = y[i];}); return x; } \
\
inline Point<DIM>& operator binop ## =(Point<DIM>& x, double a)\
{ sac_de_billes::auxi_function::for_constexpr<0,DIM>([&](int i){x[i] binop ## = a;}); return x; }\
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

// for defining all interesting operators
#define MEROPE_IMPLEMENT_POINT_ALL_BINOP(DIM, binop) \
MEROPE_IMPLEMENT_POINT_BINOPEG(DIM, binop)\
MEROPE_IMPLEMENT_POINT_BINOP(DIM, binop)\
\

//
#define MEROPE_IMPLEMENT_POINT_ALL_OP(DIM) \
MEROPE_IMPLEMENT_POINT_ALL_BINOP(DIM, +) \
MEROPE_IMPLEMENT_POINT_ALL_BINOP(DIM, -) \
MEROPE_IMPLEMENT_POINT_ALL_BINOP(DIM, *) \
MEROPE_IMPLEMENT_POINT_ALL_BINOP(DIM, / ) \
MEROPE_IMPLEMENT_POINT_SPECOP(DIM) \
\

//
template<unsigned short DIM>
Point<DIM> average(const vector<Point<DIM>>& vPoints);

// for defining all operators +=, -=, *=, /=, on DiscPoints
// parameter binop : binary operator +, -, *, /
#define MEROPE_IMPLEMENT_DISCPOINT_BINOPEG(DIM, binop) \
inline DiscPoint<DIM>& operator binop ##=(DiscPoint<DIM>& x, const DiscPoint<DIM>& y)\
{ sac_de_billes::auxi_function::for_constexpr<0,DIM>([&](int i){x[i] binop ## = y[i];}); return x; } \
\
inline DiscPoint<DIM>& operator binop ## =(DiscPoint<DIM>& x, long a)\
{ sac_de_billes::auxi_function::for_constexpr<0,DIM>([&](int i){x[i] binop ## = a;}); return x; }\
\

// for defining all operators +, -, *, / for DiscPoints
// parameter binop : binary operator +, -, *, /
// also takes advantage of move without using new memory
#define MEROPE_IMPLEMENT_DISCPOINT_BINOP(DIM, binop) \
inline DiscPoint<DIM> operator binop(const DiscPoint<DIM>& x, const DiscPoint<DIM>& y) { DiscPoint<DIM> z = x; z binop ## = y; return z; } \
\
inline DiscPoint<DIM> operator binop(DiscPoint<DIM>&& x, const DiscPoint<DIM>& y) { DiscPoint<DIM> z = std::move(x); z binop ## = y; return z; }

// for defining all interesting operators

#define MEROPE_IMPLEMENT_DISCPOINT_ALL_BINOP(DIM, binop) \
MEROPE_IMPLEMENT_DISCPOINT_BINOPEG(DIM, binop)\
MEROPE_IMPLEMENT_DISCPOINT_BINOP(DIM, binop)\
\

// all operators +, -, *, +=, -=, *=. Beware, no operator /=
#define MEROPE_IMPLEMENT_DISCPOINT_ALL_OP(DIM) \
MEROPE_IMPLEMENT_DISCPOINT_ALL_BINOP(DIM, +)\
MEROPE_IMPLEMENT_DISCPOINT_ALL_BINOP(DIM, -)\
MEROPE_IMPLEMENT_DISCPOINT_ALL_BINOP(DIM, *)\
inline DiscPoint<DIM> operator-(const DiscPoint<DIM>& x) { DiscPoint<DIM> z = x; z *= (-1); return z; } /* unary operators */ \
\

#define MEROPE_IMPLEMENT_ALL_OP(DIM) \
MEROPE_IMPLEMENT_DISCPOINT_ALL_OP(DIM) \
MEROPE_IMPLEMENT_POINT_ALL_OP(DIM) \

MEROPE_IMPLEMENT_ALL_OP(1)
MEROPE_IMPLEMENT_ALL_OP(2)
MEROPE_IMPLEMENT_ALL_OP(3)
MEROPE_IMPLEMENT_ALL_OP(4)
MEROPE_IMPLEMENT_ALL_OP(5)
MEROPE_IMPLEMENT_ALL_OP(6)
MEROPE_IMPLEMENT_ALL_OP(7)
MEROPE_IMPLEMENT_ALL_OP(8)
MEROPE_IMPLEMENT_ALL_OP(9)



//
template<unsigned short DIM>
static ostream& print(ostream& out, const Point<DIM>& p);

template<unsigned short DIM>
ostream& operator<<(ostream& out, const Point<DIM>& p);

} // namespace sac_de_billes

#include "../Geometry/Point.ixx"

#endif /* POINT_HXX_ */

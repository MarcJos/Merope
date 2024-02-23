//! Copyright : see license.txt
//!
//! \brief
//
#ifndef AUXIFUNCTIONS_HXX_
#define AUXIFUNCTIONS_HXX_

#include "StdHeaders.hxx"

namespace sac_de_billes {
using namespace std;

namespace auxi_function {

// \return the product of all the values present in the table
//! \param C : an object one can iterate on
//! \param NB : size of the desired product
template<class NB, class C>
NB productOf(const C& table);

//! compte the result discoord modulo s
//! presumably |discoord| < 2 s
template<class T1, class T2>
T1 fast_modulo(T1 discCoord, T2 s);

//! \return x^EXPONENT
template<unsigned short EXPONENT, class C>
constexpr C puissance(const C& x);

//! \return x^exponent
template<class C>
C puissance(C x, size_t exponent);


//! in the std style
//! extracts from a first list={x1, x2, ..., xn} a list according to the criterion
//! xi is extracted if there exists yj in the second list such that __fun1(xi)==__fun2(yj)
//! \param __first1 : list1.begin()
//! \param __last1 : list1.end()
//! \param __first2 : list2.begin()
//! \param __last2 : list2.end()
//! \param __result : extractedList.begin()
//! \param __fun1 : function to be applied to elements of list1
//! \param __fun2 : function to be applied to elements of list2
template<typename _InputIterator1, typename _InputIterator2,
        typename _OutputIterator, typename _Function1, typename _Function2>
_OutputIterator extract_list(_InputIterator1 __first1, _InputIterator1 __last1,
        _InputIterator2 __first2, _InputIterator2 __last2,
        _OutputIterator __result, _Function1 __fun1, _Function2 __fun2);

//! generic function for converting a vector to a string
//! \param f : output flux
template<class VECTOR>
inline void writeVectorToString(const VECTOR& vect, std::ostream& f, std::string separator = ",");

//! erase from the container if the predicate is true
template<class MAP, typename PREDICATE>
size_t erase_if(MAP& c, PREDICATE pred);

template<size_t N, class A, class B>
array<B, N> convertArray(array<A, N> array_);

template<int Begin, int... I, typename Func>
void for_constexpr(Func f, std::index_sequence<I...>);

template<int Begin, int End, typename Func>
void for_constexpr(Func f);

//! @brief make the data in the vector circulate so that
//! new_v[i+step] <- old_v[i]
//! inplace modification
//! @tparam VECTOR : array-type data
//! @warning : inefficient (additional copy)
template<class VECTOR>
void circulate(VECTOR& vect, long step);

} // namespace auxi_function
} // namespace sac_de_billes


#include "AuxiFunctions.ixx"

#endif /* AUXIFUNCTIONS_HXX_ */

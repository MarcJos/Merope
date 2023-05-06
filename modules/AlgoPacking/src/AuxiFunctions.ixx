#include "AuxiFunctions.hxx"
//! Copyright : see license.txt
//!
//! \brief Implementation of AuxiFunctions.hxx
//
#ifndef AUXIFUNCTIONS_IXX_
#define AUXIFUNCTIONS_IXX_

namespace sac_de_billes {
namespace auxi_function {

template<class NB, class C>
NB productOf(const C& table) {
    NB res = 1;
    for (const auto& nb : table) {
        res *= nb;
    }
    return res;
}

template<class T1, class T2>
T1 fast_modulo(T1 discCoord, T2 s) {
    while (discCoord < 0) {
        discCoord += s;
    }
    while (discCoord >= s) {
        discCoord -= s;
    }
    return discCoord;
}

template<unsigned short EXPONENT, class C>
inline constexpr C puissance(const C& x) {
    if constexpr (EXPONENT == 2) {
        return x * x;
    }
    else if constexpr (EXPONENT == 3) {
        return x * x * x;
    }
    else {
        constexpr unsigned short new_expo = EXPONENT / 2;
        return puissance<new_expo>(x) * puissance<EXPONENT - new_expo>(x);
    }
}

template<class C>
C puissance(C x, size_t exponent) {
    if (exponent == 0) {
        return 1;
    }
    else if (exponent == 1) {
        return x;
    }
    else {
        size_t expo1 = exponent / 2;
        return puissance(x, expo1) * puissance(x, exponent - expo1);
    }
}

template<typename _InputIterator1, typename _InputIterator2,
    typename _OutputIterator, typename _Function1, typename _Function2>
_OutputIterator extract_list(_InputIterator1 __first1, _InputIterator1 __last1,
    _InputIterator2 __first2, _InputIterator2 __last2,
    _OutputIterator __result, _Function1 __fun1, _Function2 __fun2) {
    sort(__first1, __last1, [__fun1](const auto& i1, const auto& i2) {
        return __fun1(i1) < __fun1(i2);
        });
    sort(__first2, __last2, [__fun2](const auto& i1, const auto& i2) {
        return __fun2(i1) < __fun2(i2);
        });
    while (__first1 != __last1 and __first2 != __last2) {
        auto f1 = __fun1(*__first1);
        auto f2 = __fun2(*__first2);
        if (f1 == f2) {
            *__result = *__first1;
            __result++;
            __first1++;
        }
        else if (f1 < f2) {
            __first1++;
        }
        else {
            __first2++;
        }
    }
    return __result;
}

template<class VECTOR>
inline void writeVectorToString(const VECTOR& vect, std::ostream& f) {
    for (auto it = vect.begin(); it != vect.end(); it++) {
        if (it != vect.begin()) f << ",";
        f << *it;
    }
}

template<class MAP, typename PREDICATE>
size_t erase_if(MAP& c, PREDICATE pred) {
    auto old_size = c.size();
    for (auto i = c.begin(), last = c.end(); i != last; ) {
        if (pred(*i)) {
            i = c.erase(i);
        }
        else {
            ++i;
        }
    }
    return old_size - c.size();
}

template<size_t N, class A, class B>
array<B, N> convertArray(array<A, N> array_) {
    array<B, N> result{};
    for (size_t i = 0; i < N; i++) {
        result[i] = static_cast<B>(array_[i]);
    }
    return result;
}

} // namespace auxi_function
} // namespace sac_de_billes


#endif /* AUXIFUNCTIONS_IXX_ */

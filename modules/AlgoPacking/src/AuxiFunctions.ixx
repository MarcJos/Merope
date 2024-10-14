#include "AuxiFunctions.hxx"
//! Copyright : see license.txt
//!
//! \briefImplementation of AuxiFunctions.hxx
//
#pragma once

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

// Primary template
template<typename T, typename T2 = bool>
struct make_signed_gen {
    using type = T;
};

// Specialization for unsigned types
template<typename T>
struct make_signed_gen<T, std::enable_if_t<std::is_unsigned<T>::value>> {
    using type = std::make_signed_t<T>;
};


template<class T1, class T2>
T1 fast_modulo(T1 discCoord, T2 s) {
    using signedT1 = typename make_signed_gen<T1>::type;
    signedT1 discCoord_ = static_cast<signedT1>(discCoord);
    while (discCoord_ < 0) {
        discCoord_ += s;
    }
    while (discCoord_ >= static_cast<typename make_signed_gen<T2>::type>(s)) {
        discCoord_ -= s;
    }
    return discCoord_;
}

template<unsigned short EXPONENT, class C>
inline constexpr C puissance(const C& x) {
    if constexpr (EXPONENT == 0) {
        return 1;
    } else if constexpr (EXPONENT == 1) {
        return x;
    } else {
        constexpr unsigned short expo1 = EXPONENT / 2;
        return puissance<expo1>(x) * puissance<EXPONENT - expo1>(x);
    }
}

template<class C>
C puissance(C x, size_t exponent) {
    if (exponent == 0) {
        return 1;
    } else if (exponent == 1) {
        return x;
    } else {
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
        } else if (f1 < f2) {
            __first1++;
        } else {
            __first2++;
        }
    }
    return __result;
}

template<class VECTOR>
inline void writeVectorToString(const VECTOR& vect, std::ostream& f, std::string separator) {
    for (auto it = vect.begin(); it != vect.end(); it++) {
        if (it != vect.begin()) f << separator;
        f << *it;
    }
}

template<class MAP, typename PREDICATE>
size_t erase_if(MAP& c, PREDICATE pred) {
    auto old_size = c.size();
    for (auto i = c.begin(); i != c.end(); ) {
        if (pred(*i)) {
            i = c.erase(i);
        } else {
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

template<int Begin, size_t... I, typename Func>
void for_constexpr(Func f, std::index_sequence<I...>) {
    (f(I + Begin), ...);
}

template<int Begin, int End, typename Func>
void for_constexpr(Func f) {
    static_assert(End > Begin);
    for_constexpr<Begin>(f, std::make_index_sequence<static_cast<size_t>(End - Begin)>{});
}

template<class VECTOR>
void circulate(VECTOR& vect, long step) {
    if (vect.size() == 0) {
        return;
    }
    //
    auto vect_copy = vect;
    size_t vectSize = vect.size();
    for (size_t i = 0; i < vectSize; i++) {
        vect[auxi_function::fast_modulo(i + step, vectSize)] = vect_copy[i];
    }
}


}  // namespace auxi_function
}  // namespace sac_de_billes




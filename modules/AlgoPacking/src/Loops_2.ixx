//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef LOOPS_2_IXX_
#define LOOPS_2_IXX_

namespace sac_de_billes {
using namespace std;

template<size_t nI, typename Callable>
void loop_internals(array<long, nI> &indices, size_t, Callable &f) {
    f(indices);
} // end of loop_internals

template<long I_BEGIN, long I_END, long ... In, size_t nI, typename Callable>
void loop_internals(array<long, nI> &indices, size_t c, Callable &f) {
    for (long i = I_BEGIN; i != I_END; ++i) {
        indices[c] = i;
        loop_internals<In...>(indices, c + 1, f);
    }
} // end of loop_internals

template<long ... In, typename Callable> void loop(Callable f) {
    auto indices = array<long, sizeof...(In) / 2> { };
    loop_internals<In...>(indices, 0, f);
}

// version dynamique

template<size_t I, size_t N, typename Callable>
void loop_internals(array<size_t, N> &current_indices,
        const array<size_t, N> &bounds, Callable &f) {
    if constexpr (I == N) {
        f(current_indices);
    } else {
        for (size_t i = 0; i != bounds[I]; ++i) {
            current_indices[I] = i;
            loop_internals<I + 1, N>(current_indices, bounds, f);
        }
    }
}

template<size_t N, typename Callable>
void loop(const array<size_t, N> &indices, Callable f) {
    auto i = array<size_t, N> { };
    loop_internals<0, N>(i, indices, f);
} // end of loop

template<typename Callable>
void loop(const size_t i1, const size_t i2, Callable f) {
    auto bounds = array<size_t, 2u> { i1, i2 };
    loop(bounds, f);
} // end of loop

template<typename Callable>
void loop(const size_t i1, const size_t i2, const size_t i3, Callable f) {
    auto bounds = array<size_t, 3u> { i1, i2, i3 };
    loop(bounds, f);
} // end of loop

// version dynamique

template<size_t I, size_t N, typename Callable>
void loop_with_break_internals(bool &c, array<size_t, N> &current_indices,
        const array<size_t, N> &bounds, Callable &f) {
    if constexpr (I == N) {
        f(c, current_indices);
    } else {
        for (size_t i = 0; i != bounds[I] && c; ++i) {
            current_indices[I] = i;
            loop_with_break_internals<I + 1, N>(c, current_indices, bounds, f);
        }
    }
}

template<size_t N, typename Callable>
void loop_with_break(bool &c, const array<size_t, N> &indices, Callable f) {
    auto i = array<size_t, N> { };
    loop_with_break_internals<0, N>(c, i, indices, f);
} // end of loop_with_break

template<typename Callable>
void loop_with_break(const size_t i1, const size_t i2, Callable f) {
    auto bounds = array<size_t, 2u> { i1, i2 };
    auto c = true;
    loop_with_break(c, bounds, f);
} // end of loop_with_break

template<typename Callable>
void loop_with_break(const size_t i1, const size_t i2, const size_t i3,
        Callable f) {
    auto bounds = array<size_t, 3u> { i1, i2, i3 };
    auto c = true;
    loop_with_break(c, bounds, f);
} // end of loop_with_break


/*int main() {
 cout << "** statique:\n";
 loop<2, 3>([](const array<size_t, 2u> &i) {
 cout << "{" << i[0] << ", " << i[1] << "}\n";
 });
 cout << "** dynamique:\n";
 loop(array<size_t, 2u>{2, 3}, [](const array<size_t, 2u> &i) {
 cout << "{" << i[0] << ", " << i[1] << "}\n";
 });
 cout << "** dynamique2:\n";
 loop(2, 3, [](const array<size_t, 2u> &i) {
 cout << "{" << i[0] << ", " << i[1] << "}\n";
 });
 cout << "** loop with break:\n";
 loop_with_break(2, 3, [](bool &c, const array<size_t, 2u> &i) {
 c = !((i[0] == 1) && (i[1] == 1));
 cout << "{" << i[0] << ", " << i[1] << "}\n";
 });
 return EXIT_SUCCESS;
 }*/

} // namespace sac_de_billes

#endif /* LOOPS_2_IXX_ */

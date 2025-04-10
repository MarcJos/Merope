//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace merope {

template<bool use_omp, long I_BEGIN, long I_END, size_t nI, typename Callable>
void loop_internals(array<long, nI>& indices, size_t c, Callable& f) {
    if constexpr (use_omp) {
#pragma omp parallel for private (indices)
        for (long i = I_BEGIN; i != I_END; ++i) {
            indices[c] = i;
            f(indices);
        }
    } else {
        for (long i = I_BEGIN; i != I_END; ++i) {
            indices[c] = i;
            f(indices);
        }
    }
}  // end of loop_internals

template<bool use_omp, long I_BEGIN, long I_END, long I_BEGIN2, long I_END2, long ... In, size_t nI, typename Callable>
void loop_internals(array<long, nI>& indices, size_t c, Callable& f) {
    for (long i = I_BEGIN; i != I_END; ++i) {
        indices[c] = i;
        loop_internals<use_omp, I_BEGIN2, I_END2, In...>(indices, c + 1, f);
    }
}  // end of loop_internals

template<bool use_omp, long ... In, typename Callable> void loop(Callable f) {
    auto indices = array<long, sizeof...(In) / 2> { };
    loop_internals<use_omp, In...>(indices, 0, f);
}

// version dynamique

template<bool use_omp, size_t I, size_t N, typename Callable>
void loop_internals(array<size_t, N>& current_indices,
    const array<size_t, N>& bounds, Callable& f) {
    if constexpr (I == N) {
        f(current_indices);
    } else {
        for (size_t i = 0; i != bounds[I]; ++i) {
            current_indices[I] = i;
            loop_internals<use_omp, I + 1, N>(current_indices, bounds, f);
        }
    }
}

template<bool use_omp, size_t N, typename Callable>
void loop(const array<size_t, N>& indices, Callable f) {
    auto i = array<size_t, N> { };
    loop_internals<use_omp, 0, N>(i, indices, f);
}  // end of loop

template<bool use_omp, typename Callable>
void loop(const size_t i1, const size_t i2, Callable f) {
    auto bounds = array<size_t, 2u> { i1, i2 };
    loop<use_omp>(bounds, f);
}  // end of loop

template<bool use_omp, typename Callable>
void loop(const size_t i1, const size_t i2, const size_t i3, Callable f) {
    auto bounds = array<size_t, 3u> { i1, i2, i3 };
    loop<use_omp>(bounds, f);
}  // end of loop

// version dynamique

template<size_t I, size_t N, typename Callable>
void loop_with_break_internals(bool& c, array<size_t, N>& current_indices,
    const array<size_t, N>& bounds, Callable& f) {
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
void loop_with_break(bool& c, const array<size_t, N>& indices, Callable f) {
    auto i = array<size_t, N> { };
    loop_with_break_internals<0, N>(c, i, indices, f);
}  // end of loop_with_break

template<typename Callable>
void loop_with_break(const size_t i1, const size_t i2, Callable f) {
    auto bounds = array<size_t, 2u> { i1, i2 };
    auto c = true;
    loop_with_break(c, bounds, f);
}  // end of loop_with_break

template<typename Callable>
void loop_with_break(const size_t i1, const size_t i2, const size_t i3,
    Callable f) {
    auto bounds = array<size_t, 3u> { i1, i2, i3 };
    auto c = true;
    loop_with_break(c, bounds, f);
}  // end of loop_with_break


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

}  // namespace merope



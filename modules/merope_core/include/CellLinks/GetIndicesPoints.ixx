//! Copyright : see license.txt
//!
//! \brief

#pragma once


template<unsigned short DIM>
Helper_for_getIndicesPoints<DIM>::Helper_for_getIndicesPoints(
    const AmbiantSpace::Tore<DIM> torus_, size_t nbPoints, double minDistance) : l{ 0. }, nbNodes{}, torus{ torus_ } {
    l = max(pow(auxi_function::productOf<double>(torus.L) / nbPoints, 1. / DIM), 1.1 * minDistance);
    for (size_t d = 0; d < DIM; d++) {
        nbNodes[d] = static_cast<long>(torus.L[d] / l + 1);
    }
}

template<unsigned short DIM>
vector<Identifier> Helper_for_getIndicesPoints<DIM>::getIndices(Point<DIM> x) {
    vector<Identifier> result = {};
    torus.projection(x);
    DiscPoint<DIM> ijk_0;
    for (size_t d = 0; d < DIM; d++) {
        ijk_0[d] = static_cast<long>(x[d] / l);
    }
    array<size_t, DIM> tab_indices = create_array<DIM, size_t>(1);
    loop<false>(tab_indices, [&result, this, ijk_0](const auto& id_corner) {
        DiscPoint<DIM> ijk = ijk_0;
        for (size_t d = 0; d < DIM; d++) {
            ijk[d] += id_corner[d];
        }
        result.push_back(vox::auxi::get_linear_index_periodic<DIM>(ijk, this->nbNodes));
        });
    return result;
}

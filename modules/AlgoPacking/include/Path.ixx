#include "Path.hxx"
//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace sac_de_billes {

template<unsigned short DIM, size_t NB_NGHB>
inline array<array<int, DIM>, nbSubcubes<DIM>(NB_NGHB)> Path::createParcours() {
    array<array<int, DIM>, nbSubcubes<DIM>(NB_NGHB)> path_ = { };
    size_t compteur = 0;
    constexpr auto ix = IX;
    if constexpr (DIM == 3) {
        loop<false, 0, NB_NGHB, 0, NB_NGHB, 0, NB_NGHB>([&](const array<long, 3u>& i) {
            path_[compteur] = array<int, DIM> { ix[i[0]], ix[i[1]], ix[i[2]] };
            compteur++;
            });
    } else if constexpr (DIM == 2) {
        loop<false, 0, NB_NGHB, 0, NB_NGHB>([&](const array<long, 2u>& i) {
            path_[compteur] = array<int, DIM> { ix[i[0]], ix[i[1]] };
            compteur++;
            });
    } else {
        throw invalid_argument("path::createParcours(). Incorrect DIM.");
    }

    auto norme2 = [](auto x) {
        double res = 0;
        for (size_t j = 0; j < DIM; j++) {
            res += x[j] * x[j];
        }
        return res;
        };

    auto criterionSort = [norme2](auto x, auto y) {
        return norme2(x) < norme2(y);
        };
    sort(path_.begin(), path_.end(), criterionSort);
    return path_;
}

template<unsigned DIM, size_t NB_NGHB>
inline const array<array<int, DIM>, nbSubcubes<DIM>(NB_NGHB)>& Path::pathForNeighbors() const {
    if constexpr (DIM == 3) {
        if constexpr (NB_NGHB == 5) {
            return pathNghb_3D_5;
        } else if constexpr (NB_NGHB == 7) {
            return pathNghb_3D_7;
        }
    } else if constexpr (DIM == 2) {
        if constexpr (NB_NGHB == 5) {
            return pathNghb_2D_5;
        } else if constexpr (NB_NGHB == 7) {
            return pathNghb_2D_7;
        }
    }
    throw invalid_argument("Path::pathForNeighbors Incorrect DIM");
}

template<unsigned short DIM>
inline size_t Path::fromCorner2Index(array<int, DIM> posRelative) const {
    size_t result = 0;
    array<size_t, 3> exponents({ 9, 3, 1 });
    for (size_t j = 0; j < DIM; j++) {
        if (posRelative[j] > 0) {
            result += 2 * exponents[j];
        } else if (posRelative[j] == 0) {
            result += 1 * exponents[j];
        }
        // else if(posRelative[j]<0){
        //  ;
        // }
    }
    return result;
}

namespace path {

template<unsigned short DIM>
TabCorner<DIM>::TabCorner() : storage_tab_corner{} {
    if constexpr (DIM == 2) {
        storage_tab_corner = {
            array<unsigned short, 2> { 0, 0 },
                array<unsigned short, 2> { 1, 1 },
                array<unsigned short, 2> { 0, 1 },
                array<unsigned short, 2> { 1, 0 }
        };
    } else if constexpr (DIM == 3) {
        storage_tab_corner = {
            array<unsigned short, 3> { 0, 0, 0 },
                array<unsigned short, 3> { 1, 1, 1 },
                array<unsigned short, 3> { 0, 0, 1 },
                array<unsigned short, 3> { 1, 1, 0 },
                array<unsigned short, 3> { 0, 1, 0 },
                array<unsigned short, 3> { 1, 0, 1 },
                array<unsigned short, 3> { 0, 1, 1 },
                array<unsigned short, 3> { 1, 0, 0 }
        };
    } else {
        size_t nb_corners = auxi_function::puissance<DIM>(2);
        storage_tab_corner.resize(nb_corners);
        for (size_t i = 0; i < nb_corners; i++) {
            int k = i;
            for (size_t j = DIM; j > 0; --j) {
                storage_tab_corner[i][j] = k / auxi_function::puissance(2, j);
                k -= auxi_function::puissance(2, j) * storage_tab_corner[i][j];
            }
        }
    }
}

template<unsigned short DIM>
inline const TabCorner<DIM>& TabCorner<DIM>::get() {
    static TabCorner<DIM> tabC{};
    return tabC;
}

}  // namespace  path

}  // namespace sac_de_billes



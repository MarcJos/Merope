//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../../AlgoPacking/src/StdHeaders.hxx"
#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../Voronoi/VoroInterface.hxx"
#include "../../../Optimization/include/Anderson.hxx"
#include "../../../Optimization/include/BarzilaiBorwein.hxx"

namespace merope {
using namespace sac_de_billes;

namespace optimizeLaguerreTess {

template<unsigned short DIM>
void extract_center(const vector<Sphere<DIM>>& centerTessels_loc,
    vector<double>& centers_loc) {
    for (size_t i = 0; i < centerTessels_loc.size(); i++) {
        for (size_t k = 0; k < DIM; k++) {
            centers_loc[3 * i + k] = centerTessels_loc[i].center[k];
        }
    }
}

template<unsigned short DIM>
void put_center(Point<DIM> L, const vector<double>& centers_loc,
    vector<Sphere<DIM>>& centerTessels_loc) {
    for (size_t i = 0; i < centerTessels_loc.size(); i++) {
        for (size_t k = 0; k < DIM; k++) {
            centerTessels_loc[i].center[k] = centers_loc[3 * i + k];
            geomTools::projection_periodic_1D(centerTessels_loc[i].center[k], L[k]);
        }
    }
}


//! @brief : accelerated implementation of the Lloyd algorithm for getting centroidal Laguerre tessellations
//! @see : Fast methods for computing centroidal Laguerre tessellations for prescribed volume fractions with applications to microstructur generation of polycrystalline materials (Kuhn & al)
//! @param L : lengths of the cube
//! @param centerTessels : centers, weight = r_i^2
template<unsigned short DIM, unsigned short NB_ITER = 2>
void makeCentroidal(const Point<DIM>& L, vector<Sphere<DIM>>& centerTessels,
    size_t max_iter, double stopping_criterion,
    bool use_acceleration,
    bool verbose = false) {
    static_assert(DIM == 3);

    // local object
    auto centerTessels_loc = centerTessels;
    auto centers_loc = vector<double>(DIM * centerTessels.size());
    //
    auto fixed_point_application = [&](const vector<double>& original_centers) {
        put_center<DIM>(L, original_centers, centerTessels_loc);
        auto relative_centers = voroInterface::compute_relative_centroids(centerTessels_loc, L);
        for (size_t i = 0; i < centerTessels_loc.size(); i++) {
            for (size_t k = 0; k < DIM; k++) {
                centers_loc[3 * i + k] = relative_centers[i][k] + original_centers[3 * i + k];
            }
        }
        return centers_loc;
        };
    //
    extract_center<DIM>(centerTessels, centers_loc);
    optimization::Anderson::algorithm<NB_ITER, decltype(centers_loc), decltype(fixed_point_application)>
        algo(centers_loc, fixed_point_application);
    algo.proceed(max_iter, stopping_criterion, use_acceleration);
    put_center<DIM>(L, algo.get_x(), centerTessels);
    //
    if (verbose) {
        algo.show_message(std::cerr);
    }
}

//! @brief : accelerated implementation of the Lloyd algorithm for getting centroidal Laguerre tessellations
//! @see : Fast methods for computing centroidal Laguerre tessellations for prescribed volume fractions with applications to microstructur generation of polycrystalline materials (Kuhn & al)
//! @param L : lengths of the cube
//! @param centerTessels : centers, weight = r_i^2
template<unsigned short DIM>
void makeCentroidal__2(const Point<DIM>& L, vector<Sphere<DIM>>& centerTessels,
    size_t max_iter, double stopping_criterion, bool use_acceleration,
    bool verbose = false) {
    static_assert(DIM == 3);
    // local object
    auto centerTessels_loc = centerTessels;
    auto vec_concatenated = vector<double>(DIM * centerTessels.size());
    //
    auto nabla_f = [&](const vector<double>& original_centers) {
        put_center<DIM>(L, original_centers, centerTessels_loc);
        auto centers = voroInterface::compute_relative_centroids(centerTessels_loc, L);
        auto volumes = voroInterface::compute_volumes(centerTessels_loc, L);
        for (size_t i = 0; i < centerTessels_loc.size(); i++) {
            for (size_t k = 0; k < DIM; k++) {
                vec_concatenated[3 * i + k] = -centers[i][k];
                vec_concatenated[3 * i + k] *= 2 * volumes[i];
            }
        }
        //! for verbosity
        /*
        std::cerr << "delta = " << sqrt(std::inner_product(vec_concatenated.begin(), vec_concatenated.end(), vec_concatenated.begin(), 0.))
        / (2 * auxi_function::productOf<double>(L) / centerTessels.size()) << endl;
        */
        //! for verbosity
        return vec_concatenated;
        };
    //
    extract_center<DIM>(centerTessels_loc, vec_concatenated);
    optimization::Barzilai_Borwein::algorithm<vector<double>, decltype(nabla_f)> algo(vec_concatenated, nabla_f);
    algo.proceed(max_iter, stopping_criterion);
    put_center<DIM>(L, algo.get_x(), centerTessels);
    //
    if (verbose) {
        algo.show_message(std::cerr);
    }
}

} // namespace optimizeLaguerreTess

} // namespace merope
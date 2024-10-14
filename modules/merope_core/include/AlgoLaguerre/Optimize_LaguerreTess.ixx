//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../Voronoi/VoroInterface.hxx"

namespace merope {

namespace optimizeLaguerreTess {

template<unsigned short DIM>
void algo_fit_volumes<DIM>::proceed(double max_delta_Vol, size_t max_iter, bool verbose) {
    if constexpr (DIM == 2) {
        max_delta_Vol /= internal_algo.L[2];
    }
    internal_algo.proceed(max_delta_Vol, max_iter, verbose);
}

template<unsigned short DIM>
std::conditional_t<DIM == 3, const vector<Sphere<DIM>>&, vector<Sphere<DIM>>> algo_fit_volumes<DIM>::getCenterTessels() const {
    if constexpr (DIM == 3) {
        return internal_algo.getCenterTessels();
    }
    if constexpr (DIM == 2) {
        return transform_2D(internal_algo.getCenterTessels());
    }
}

template<unsigned short DIM>
double algo_fit_volumes<DIM>::maxDeltaVolumes() const {
    if constexpr (DIM == 3) {
        return internal_algo.maxDeltaVolumes();
    }
    if constexpr (DIM == 2) {
        return internal_algo.maxDeltaVolumes() / internal_algo.L[2];
    }
}

template<unsigned short DIM>
std::conditional_t<DIM == 3, const vector<double>&, vector<double>> algo_fit_volumes<DIM>::getCurrentVolumes() const {
    if constexpr (DIM == 3) {
        return internal_algo.getCurrentVolumes();
    }
    if constexpr (DIM == 2) {
        return transform_2D(internal_algo.getCurrentVolumes(), internal_algo.L);
    }
}

template<unsigned short DIM>
algo_fit_volumes<DIM>::algo_fit_volumes(const Point<DIM>& L, const vector<Sphere<DIM>>& centerTessels, const vector<double>& desiredVolumes)
    : internal_algo(transform_3D(L), transform_3D(centerTessels), transform_3D(desiredVolumes, transform_3D(L))) {}

template<unsigned short DIM>
Point<3> algo_fit_volumes<DIM>::transform_3D(const Point<2>& L) {
    Point<3> L_extended = { L[0], L[1],
        L_multiplier * min(L[0], L[1]) };
    return L_extended;
}

template<unsigned short DIM>
vector<Sphere<3>> algo_fit_volumes<DIM>::transform_3D(const vector<Sphere<2>>& centerTessels) {
    vector<Sphere<3>> theSpheres(centerTessels.size());
    for (size_t i = 0; i < centerTessels.size(); i++) {
        theSpheres[i].center = { centerTessels[i].center[0], centerTessels[i].center[1], 0. };
        theSpheres[i].radius = centerTessels[i].radius;
        theSpheres[i].phase = centerTessels[i].phase;
    }
    return theSpheres;
}

template<unsigned short DIM>
vector<double> algo_fit_volumes<DIM>::transform_3D(const vector<double>& desiredVolumes_2D, const Point<3>& L) {
    vector<double> rescaled_volumes = desiredVolumes_2D;
    if constexpr (DIM == 2) {
        for (auto& vol : rescaled_volumes) {
            vol *= L[2];
        }
    }
    return rescaled_volumes;
}

template<unsigned short DIM>
vector<Sphere<2>> algo_fit_volumes<DIM>::transform_2D(const vector<Sphere<3>>& centerTessels) {
    vector<Sphere<2>> theSpheres(centerTessels.size());
    for (size_t i = 0; i < centerTessels.size(); i++) {
        theSpheres[i].center = { centerTessels[i].center[0], centerTessels[i].center[1] };
        theSpheres[i].radius = centerTessels[i].radius;
        theSpheres[i].phase = centerTessels[i].phase;
    }
    return theSpheres;
}

template<unsigned short DIM>
vector<double> algo_fit_volumes<DIM>::transform_2D(const vector<double>& desiredVolumes_3D, const Point<3>& L) {
    vector<double> rescaled_volumes = desiredVolumes_3D;
    for (auto& vol : rescaled_volumes) {
        vol /= L[2];
    }
    return rescaled_volumes;
}


}  // namespace  optimizeLaguerreTess

}  // namespace  merope

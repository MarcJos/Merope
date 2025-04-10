//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"

namespace merope {
namespace optimizeLaguerreTess {

//! optimize the laguerre tessellation by changing the radii so that the volume fraction of each
//! tessel is equal to the prescribed volumeFraction
//! implements the Barzilai-Borwein algorithm described in
//! Fast methods for computing centroidal Laguerre tessellations for
//! prescribed volume fractions with applications to microstructure
//! generation of polycrystalline materials (Kuhn & al)
//! @return : new centerTessels, of unchanged center, but radius changed
class algo_fit_volumes_3D {
public:
    //! @brief : constructor
    //! @param L : lengths of the cube
    //! @param centerTessels : centers, weight = r_i^2
    //! @param desiredVolumes : prescribed final volumes
    algo_fit_volumes_3D(const Point<3>& L, const vector<Sphere<3>>& centerTessels, const vector<double>& desiredVolumes);
    //! @brief proceed with the algorithm
    //! @param max_delta_Vol : maximal error on the volume accepted
    //! @param max_iter : maximal number of iterations for Barzilai-Borwein
    void proceed(double max_delta_Vol = 1e-10, size_t max_iter = 3000, bool verbose = false);
    //! @return the current volumes
    //! @warning : should proceed first with computeCurrentVolumes
    const vector<double>& getCurrentVolumes() const { return volumes_current; }
    //! @return the maximal error between current volumes and objective volumes
    //! @warning : should proceed first with computeCurrentVolumes
    double maxDeltaVolumes() const;
    //! @brief : const getter
    const vector<Sphere<3>>& getCenterTessels() const { return centerTessels; }

    //! lengths of the cube
    const Point<3> L;
    //! aimed for objective
    const vector<double> volumes_objective;

private:
    //! centerTessels
    vector<Sphere<3>> centerTessels;
    //! storage for storing current volumes
    vector<double> volumes_current;
    //! gradient of the function
    //! warning :stores the w_i
    vector<double> nabla_g(const vector<double>& w_i);
    //! @brief set the radii of centerTessels according to the weights
    //! @param w_i : given weigths in a Laguerre tessel
    void set_radii_out_of_weights(const vector<double>& w_i);
    //! @brief write a verbose output
    void verbose_output(std::ostream& f, const auto& algo);
    //! @return the current volumes occupied by each tessel
    void computeCurrentVolumes();
};

//! @brief general interface template class for 2D and 3D
//! @see algo_fit_volumes_3D
template<unsigned short DIM>
class algo_fit_volumes {
    static_assert(DIM == 2 or DIM == 3);
public:
    //! @brief : constructor
    //! @param L : lengths of the cube
    //! @param centerTessels : centers, weight = r_i^2
    //! @param desiredVolumes : prescribed final volumes
    algo_fit_volumes(const Point<DIM>& L, const vector<Sphere<DIM>>& centerTessels, const vector<double>& desiredVolumes);
    //! @brief proceed with the algorithm
    //! @param max_delta_Vol : maximal error on the volume accepted
    //! @param max_iter : maximal number of iterations for Barzilai-Borwein
    void proceed(double max_delta_Vol = 1e-8, size_t max_iter = 3000, bool verbose = false);
    //! @brief : const getter
    std::conditional_t<DIM == 3, const vector<Sphere<DIM>>&, vector<Sphere<DIM>>> getCenterTessels() const;
    //! @return the maximal error between current volumes and objective volumes
    double maxDeltaVolumes() const;
    //! @return the current volumes
    std::conditional_t<DIM == 3, const vector<double>&, vector<double>> getCurrentVolumes() const;
private:
    //! @brief engine
    algo_fit_volumes_3D internal_algo;
    //! trivial transformations
    static const Point<3>& transform_3D(const Point<3>& L) { return L; }
    static const vector<Sphere<3>>& transform_3D(const vector<Sphere<3>>& centerTessels) { return centerTessels; }
    //!
    static constexpr double L_multiplier = 0.01;
    static Point<3> transform_3D(const Point<2>& L);
    static vector<Sphere<3>> transform_3D(const vector<Sphere<2>>& centerTessels);
    static vector<double> transform_3D(const vector<double>& desiredVolumes_2D, const Point<3>& L);
    static vector<Sphere<2>> transform_2D(const vector<Sphere<3>>& centerTessels);
    static vector<double> transform_2D(const vector<double>& desiredVolumes_3D, const Point<3>& L);
};

namespace auxi {
//! @return radii r_i corresponding to weigths. If all weights are equal r_i^2 = w_i.
//! If not, weigths should be first shifted.
vector<double> get_radii_out_of_weights(const vector<double>& w_i);

bool preconditions(const vector<Sphere<3>>& centerTessels, const Point<3>& L,
    const vector<double>& volumeFractions, double max_delta_Vol);

}  // namespace  auxi

}  // namespace  optimizeLaguerreTess

}  // namespace  merope

#include "Optimize_LaguerreTess.ixx"



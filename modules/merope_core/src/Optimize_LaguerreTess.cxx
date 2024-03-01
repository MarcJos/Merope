//! Copyright : see license.txt
//!
//! \brief

#include "AlgoLaguerre/Optimize_LaguerreTess.hxx"
#include "../../Optimization/include/BarzilaiBorwein.hxx"


namespace merope {
using namespace sac_de_billes;

namespace optimizeLaguerreTess {

algo_fit_volumes_3D::algo_fit_volumes_3D(const Point<3>& L_, const vector<Sphere<3>>& centerTessels_, const vector<double>& desiredVolumes)
    :L(L_), volumes_objective(desiredVolumes), centerTessels(centerTessels_), volumes_current{} {
    computeCurrentVolumes();
}

void algo_fit_volumes_3D::set_radii_out_of_weights(const vector<double>& w_i) {
    auto radii = auxi::get_radii_out_of_weights(w_i);
    for (size_t i = 0; i < centerTessels.size(); i++) {
        centerTessels[i].radius = radii[i];
    }
}

void algo_fit_volumes_3D::computeCurrentVolumes() {
    volumes_current = voroInterface::compute_volumes(centerTessels, L);
}

vector<double>  algo_fit_volumes_3D::nabla_g(const vector<double>& w_i) {
    auto res = vector<double>(w_i.size());
    // compute volumes related to weigths
    set_radii_out_of_weights(w_i);
    computeCurrentVolumes();
    //
    for (size_t i = 0; i < res.size(); i++) {
        res[i] = volumes_current[i] - volumes_objective[i];
    }
    return res;
}

void algo_fit_volumes_3D::proceed(double max_delta_Vol, size_t max_iter, bool verbose) {
    // warning! preconditions
    if (not auxi::preconditions(centerTessels, L, volumes_objective, max_delta_Vol)) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Preconditions are not met!");
    }
    //////////////////////////////////////////////////
    // initial stage
    vector<double> w_init(centerTessels.size());
    for (size_t i = 0; i < centerTessels.size(); i++) {
        w_init[i] = auxi_function::puissance<2>(centerTessels[i].radius);
    }
    // optimize
    int iter_number = 0;
    optimization::Barzilai_Borwein::algorithm algo(w_init, [&](const auto& w_i) {
        iter_number += 1;
        if (iter_number % 10 == 0) {
            std::cout << "Iteration : " << iter_number << std::endl;
        }
        return nabla_g(w_i);}
    );
    bool has_succeeded_optimization = algo.proceed(max_iter, [&](auto) {return  maxDeltaVolumes() < max_delta_Vol;});
    if (not has_succeeded_optimization) {
        throw runtime_error("Optimization is not a success");
    }
    // post-process
    set_radii_out_of_weights(algo.get_x());
    // outputs
    if (verbose) {
        verbose_output(cout, algo);
    }
}

double algo_fit_volumes_3D::maxDeltaVolumes() const {
    double max_deltaVolume = 0;
    for (size_t i = 0; i < volumes_current.size(); i++) {
        double loc_deltaV = abs(volumes_current[i] - volumes_objective[i]);
        max_deltaVolume = max(loc_deltaV, max_deltaVolume);
    }
    return max_deltaVolume;
}

void algo_fit_volumes_3D::verbose_output(std::ostream& f, const auto& algo) {
    algo.show_message(f);
}

vector<double> auxi::get_radii_out_of_weights(const vector<double>& w_i) {
    double minimum = *(std::min_element(w_i.begin(), w_i.end()));
    double shift_w = max(-minimum, 0.); // shift negative values for w_i
    vector<double> radii(w_i.size());
    for (size_t i = 0; i < radii.size(); i++) {
        radii[i] = sqrt(max(shift_w + w_i[i], 0.));
    }
    return radii;
}

bool auxi::preconditions(const vector<Sphere<3>>& centerTessels, const Point<3>& L,
    const vector<double>& volumes_objective, double max_delta_Vol) {
    if (centerTessels.size() != volumes_objective.size()) {
        cerr << "incompatible sizes" << endl;
        return false;
    }
    if (centerTessels.size() == 0) {
        cerr << "centerTessels.size() == 0" << endl;
        return false;
    } else if (abs(std::accumulate(volumes_objective.begin(), volumes_objective.end(), 0.)
        - auxi_function::productOf<double>(L)) > 1e-6) {
        cerr << "sum volumes not equal to total volume : "
            << std::accumulate(volumes_objective.begin(), volumes_objective.end(), 0.)
            << " vs " << auxi_function::productOf<double>(L) << endl;
        return false;
    } else if (max_delta_Vol <= 0.) {
        cerr << "negative max_delta_Vol" << endl;
        return false;
    }
    for (size_t i = 0; i < volumes_objective.size(); i++) {
        if (volumes_objective[i] <= 0) {
            cerr << "zero volume objective" << endl;
            return false;
        }
    }
    return true;
}

} // namespace  optimizeLaguerreTess

} // namespace  merope

//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../GenericMerope/StdHeaders.hxx"
#include "Anderson.hxx"

namespace optimization {
namespace tests {

template<bool USE_ANDERSON, int NB_ITER>
void test_Anderson() {//
    using namespace std;
    cerr << "################" << endl;
    cerr << "test_Anders()" << endl;
    cerr << std::boolalpha;
    cerr << "USE_ANDERSON : " << USE_ANDERSON << " ; " << "NB_ITER : " << NB_ITER << endl;
    cerr << "################" << endl;
    auto h = [](auto w) {
        return vector<double>{12 + 0.2 * w[1] + 0.01 * sqrt(abs(w[0] * w[1])),
            4 + 0.3 * w[0],
            19 + 0.2 * sin(w[3]),
            25 + 0.1 * abs(w[0] + w[1] + w[2] + w[3]),
            0.1 * abs(w[4]) + 0.05 * cos(w[2]),
            1 + 0.7 * w[5]};
        };
    vector<double> w_0 = { 0., 0., 0., 0., 0., 1. };
    optimization::Anderson::algorithm<NB_ITER, decltype(w_0), decltype(h)> algo{ w_0, h };
    algo.proceed(100, 1e-10, USE_ANDERSON);
    algo.show_message(std::cerr);
    auto w = algo.get_x();
    //
    auto print_loc = [](auto w_) {
        for (size_t i = 0; i < w_.size(); i++) {
            std::cerr << w_[i] << " , ";
        }
        cerr << endl;
        };
    //
    auto h_w = h(w);
    print_loc(w);
    std::cerr << endl;
    print_loc(h_w);
}

}  // namespace  tests
}  // namespace  optimization
//! Copyright : see license.txt
//!
//! \brief

#include "../include/Tests.hxx"
#include "../include/Anderson.hxx"
#include "../include/BarzilaiBorwein.hxx"

namespace optimization {
namespace tests {

void print(auto my_vector) {
    for (size_t i = 0; i < my_vector.size(); i++) {
        cout << my_vector[i] << ",";
    }
    cout << endl;
}

double f(const std::vector<double>& x) {
    return sin(std::inner_product(x.begin(), x.end(), x.begin(), 0.));
}

std::vector<double> nabla_f(const std::vector<double>& x) {
    double factor = 2 * cos(std::inner_product(x.begin(), x.end(), x.begin(), 0.));
    auto res = x;
    for (size_t i = 0; i < x.size(); i++) {
        res[i] *= factor;
    }
    return res;
}

void test1_BarzilaiBorwein() {
    cerr << "################" << endl;
    cerr << "test1_BarzilaiBorwein()" << endl;
    cerr << "################" << endl;
    std::vector<double> v_init = { 1., 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    auto nabf = [](const auto& x) {return nabla_f(x);};
    auto algo = Barzilai_Borwein::algorithm<std::vector<double>, decltype(nabf)>(v_init, nabf);
    for (size_t i = 0; i < 10; i++) {
        algo.iterate();
        print(algo.get_x());
        std::cerr << "## : " << f(algo.get_x()) << endl;
    }
}

void test2_BarzilaiBorwein() {
    cerr << "################" << endl;
    cerr << "test2_BarzilaiBorwein()" << endl;
    cerr << "################" << endl;
    std::vector<double> v_init = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    //! test with final criterion
    auto nabf = [](const auto& x) {return nabla_f(x);};
    auto algo2 = Barzilai_Borwein::algorithm<std::vector<double>, decltype(nabf)>(v_init, nabf);
    algo2.proceed(30, 0.0001);
    print(algo2.get_x());
    std::cerr << "## : " << f(algo2.get_x()) << endl;
}

}  // namespace  tests
}  // namespace  optimization


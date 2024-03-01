//! Copyright : see license.txt
//!
//! \brief

#pragma once

namespace optimization {

namespace Barzilai_Borwein {
template<class VECTOR, class NABLA_F>
algorithm<VECTOR, NABLA_F>::algorithm(const VECTOR& w_0, NABLA_F compute_nabla_f_, double first_alpha) :
    compteur{ 0 }, compute_nabla_f{ compute_nabla_f_ }, x_n(w_0), x_n_1(0), nabla_f(w_0), nabla_f_1(w_0) {
    nabla_f = compute_nabla_f_(x_n);
    iterate_gradient_descent(first_alpha);
}

template<class VECTOR, class NABLA_F>
template<class STOPPING_TEST>
bool algorithm<VECTOR, NABLA_F>::proceed(size_t max_iter, STOPPING_TEST stopping_criterion) {
    size_t compte = 0;
    auto w_temp = x_n;
    do {
        iterate();
        for (size_t i = 0; i < x_n.size(); i++) {
            w_temp[i] = x_n[i] - x_n_1[i];
        }
        compte++;
        compteur++;
    } while (not(stopping_criterion(w_temp)) and compte < max_iter);
    return compte < max_iter;
}

template<class VECTOR, class NABLA_F>
bool algorithm<VECTOR, NABLA_F>::proceed(size_t max_iter, double stopping_criterion) {
    auto stopping_criterion_func = [stopping_criterion](const auto& w_temp) {
        return std::inner_product(w_temp.begin(), w_temp.end(), w_temp.begin(), 0.)
            < stopping_criterion * stopping_criterion;
        };
    return proceed(max_iter, stopping_criterion_func);
}

template<class VECTOR, class NABLA_F>
void algorithm<VECTOR, NABLA_F>::show_message(std::ostream& os) const {
    os << "############\n";
    os << "Execution of Barzilai-Borwein\n";
    os << "Nb iterations : " << compteur << endl;
    os << "############\n";
}

template<class VECTOR, class NABLA_F>
void algorithm<VECTOR, NABLA_F>::iterate() {
    double alpha_n = auxi::compute_alpha_n(x_n, x_n_1, nabla_f, nabla_f_1);
    iterate_gradient_descent(alpha_n);
}

template<class VECTOR, class NABLA_F>
void algorithm<VECTOR, NABLA_F>::iterate_gradient_descent(double alpha_n) {
    // save
    x_n_1 = x_n;
    nabla_f_1 = nabla_f;
    // update
    for (size_t i = 0; i < x_n.size(); i++) {
        x_n[i] -= alpha_n * nabla_f[i];
    }
    nabla_f = compute_nabla_f(x_n);
}

template<class VECTOR>
double auxi::compute_alpha_n(const VECTOR& x_n, const VECTOR& x_n_1, const VECTOR& nabla_f, const VECTOR& nabla_f_1) {
    VECTOR delta_w = x_n;
    VECTOR delta_nabla_f = nabla_f;
    for (size_t i = 0; i < x_n.size(); i++) {
        delta_w[i] -= x_n_1[i];
        delta_nabla_f[i] -= nabla_f_1[i];
    }
    double alpha_n = std::inner_product(delta_w.begin(), delta_w.end(), delta_w.begin(), 0.)
        / std::inner_product(delta_nabla_f.begin(), delta_nabla_f.end(), delta_w.begin(), 0.);
    return alpha_n;
}

} // namespace Barzilai_Borwein

} // namespace  optimization
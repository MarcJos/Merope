#include "Anderson.hxx"
//! Copyright : see license.txt
//!
//! \brief

#pragma once

namespace optimization {

using namespace std;

namespace Anderson {

template<unsigned short NB_ITER, class VECTOR, class FIXED_POINT>
algorithm<NB_ITER, VECTOR, FIXED_POINT>::algorithm(VECTOR w_0, FIXED_POINT h_) : h{ h_ },
x_n{}, h_n{}, size{ w_0.size() } {
    //! compute the NB_ITER first with fixed point iterations
    x_n.reserve(NB_ITER + 1);
    h_n.reserve(NB_ITER + 1);
    x_n.push_back(w_0);
    h_n.push_back(h(w_0));
    for (size_t i = 1; i < NB_ITER + 1; i++) {
        x_n.push_back(h_n[i - 1]);
        h_n.push_back(h(x_n[i]));
    }
}


template<unsigned short NB_ITER, class VECTOR, class FIXED_POINT>
template<class STOPPING_TEST>
bool algorithm<NB_ITER, VECTOR, FIXED_POINT>::proceed(size_t max_iter, STOPPING_TEST stopping_criterion, bool use_acceleration) {
    size_t compte = 0;
    do {
        if (use_acceleration) {
            iterate();
        } else {
            iterate_fixed_point();
        };
        compte++;
    } while (not(stopping_criterion(get_x())) and compte < max_iter);
    return compte == max_iter;
}

template<unsigned short NB_ITER, class VECTOR, class FIXED_POINT>
bool algorithm<NB_ITER, VECTOR, FIXED_POINT>::proceed(size_t max_iter, double stopping_tolerance, bool use_acceleration) {
    auto stopping_criterion = [&](const auto&) {
        VECTOR Delta_X = compute_Delta_X();
        double delta_square = std::inner_product(Delta_X.begin(), Delta_X.end(), Delta_X.begin(), 0.);
        return delta_square < stopping_tolerance * stopping_tolerance;
        };
    return proceed(max_iter, stopping_criterion, use_acceleration);
}

template<unsigned short NB_ITER, class VECTOR, class FIXED_POINT>
VECTOR algorithm<NB_ITER, VECTOR, FIXED_POINT>::compute_Delta_X() const {
    VECTOR Delta_X(size);
    for (size_t k = 0; k < size; k++) {
        Delta_X[k] = h_n[NB_ITER][k] - x_n[NB_ITER][k];
    }
    return Delta_X;
}

template<unsigned short NB_ITER, class VECTOR, class FIXED_POINT>
void algorithm<NB_ITER, VECTOR, FIXED_POINT>::update_with(const VECTOR& new_x_n, bool use_acceleration) {
    for (size_t i = 0; i < NB_ITER; i++) {
        x_n[i] = x_n[i + 1];
        h_n[i] = h_n[i + 1];
    }
    x_n[NB_ITER] = new_x_n;
    h_n[NB_ITER] = h(new_x_n);
    //!
    VECTOR Delta_X = compute_Delta_X();
    double delta = sqrt(std::inner_product(Delta_X.begin(), Delta_X.end(), Delta_X.begin(), 0.));
    // keeping track of the way the algorithm behaves
    this->compteur++;
    this->sequence_of_residues.push_back(delta);
    this->sequence_of_successes.push_back(use_acceleration);
}

template<unsigned short NB_ITER, class VECTOR, class FIXED_POINT>
void algorithm<NB_ITER, VECTOR, FIXED_POINT>::iterate_worse() {
    VECTOR Delta_X = compute_Delta_X();
    vector<VECTOR> delta_R_k(NB_ITER);
    vector<VECTOR> delta_h_k(NB_ITER);
    //!
    //!
    for (size_t i = 0; i < NB_ITER; i++) {
        delta_h_k[i].resize(size);
        delta_R_k[i].resize(size);
        for (size_t k = 0; k < size; k++) {
            delta_h_k[i][k] = h_n[i + 1][k] - h_n[i][k];
            delta_R_k[i][k] = (h_n[i + 1][k] - x_n[i + 1][k]) - (h_n[i][k] - x_n[i][k]);
        }
    }

    //! evaluation of (40) in Iterative residual-based vector methods to accelerate fixed
    // point iterations, Ramiere & Helfer
    VECTOR new_x_n = h_n[NB_ITER];
    //!
    Eigen::MatrixXd matrix(NB_ITER, NB_ITER);
    for (size_t i = 0; i < NB_ITER; i++) {
        for (size_t j = 0; j < NB_ITER; j++) {
            matrix(i, j) = std::inner_product(delta_R_k[i].begin(), delta_R_k[i].end(), delta_R_k[j].begin(), 0.);
        }
    }
    Eigen::MatrixXd matrix_inverse = matrix.inverse();
    //! goes back to normal fixed-point iterations. Badly conditionned matrix
    if ((isnan(matrix_inverse.norm())) or matrix.norm() * matrix_inverse.norm() > 1e4) {
        /*
        std::cerr << std::boolalpha;
        std::cerr << "Is nan inverse : " << isnan(matrix_inverse.norm()) << endl;
        std::cerr << "precond number : " << matrix.norm() * matrix_inverse.norm() << endl;
        std::cerr << "Anderson disabled. Iteration nb : " << compteur << endl;
        */
        iterate_fixed_point();
        return;
    }
    //! END goes back to normal iterations
    //!
    Eigen::VectorXd rightmost_vector(NB_ITER);
    for (size_t i = 0; i < NB_ITER; i++) {
        rightmost_vector(i) = std::inner_product(delta_R_k[i].begin(), delta_R_k[i].end(), Delta_X.begin(), 0.);
    }
    //!
    auto rightmost_terms = matrix_inverse * rightmost_vector;
    for (size_t i = 0; i < NB_ITER; i++) {
        for (size_t k = 0; k < size; k++) {
            new_x_n[k] -= delta_h_k[i][k] * rightmost_terms(i);
        }
    }
    update_with(new_x_n, true);
}


template<unsigned short NB_ITER, class VECTOR, class FIXED_POINT>
void algorithm<NB_ITER, VECTOR, FIXED_POINT>::iterate() {
    //! see Schneider & Kuhn
    vector<VECTOR> R_k(NB_ITER + 1);
    for (size_t i = 0; i < NB_ITER + 1; i++) {
        R_k[i].resize(size);
        for (size_t k = 0; k < size; k++) {
            R_k[i][k] = h_n[i][k] - x_n[i][k];
        }
    }
    vector<VECTOR> R_k_minus_R_m(NB_ITER);
    for (size_t i = 0; i < NB_ITER; i++) {
        R_k_minus_R_m[i].resize(size);
        for (size_t k = 0; k < size; k++) {
            R_k_minus_R_m[i][k] = R_k[i][k] - R_k[NB_ITER][k];
        }
    }

    Eigen::MatrixXd my_matrix(NB_ITER, NB_ITER);
    for (size_t i = 0; i < NB_ITER; i++) {
        for (size_t j = 0; j < NB_ITER; j++) {
            my_matrix(i, j) = std::inner_product(R_k_minus_R_m[i].begin(), R_k_minus_R_m[i].end(),
                R_k_minus_R_m[j].begin(), 0.);
        }
    }
    Eigen::VectorXd my_vector(NB_ITER);
    for (size_t i = 0; i < NB_ITER; i++) {
        my_vector(i) = -std::inner_product(R_k_minus_R_m[i].begin(), R_k_minus_R_m[i].end(),
            R_k[NB_ITER].begin(), 0.);
    }
    Eigen::VectorXd alpha = my_matrix.colPivHouseholderQr().solve(my_vector);
    //!
    VECTOR new_x_n(size, 0.);
    for (size_t i = 0; i < NB_ITER; i++) {
        for (size_t k = 0; k < size; k++) {
            new_x_n[k] += alpha[i] * h_n[i][k];
        }
    }
    double alpha_NB_ITER = 1 - std::accumulate(alpha.begin(), alpha.end(), 0.);
    for (size_t k = 0; k < size; k++) {
        new_x_n[k] += alpha_NB_ITER * h_n[NB_ITER][k];
    }
    //!
    update_with(new_x_n, true);
}

template<unsigned short NB_ITER, class VECTOR, class FIXED_POINT>
void algorithm<NB_ITER, VECTOR, FIXED_POINT>::iterate_fixed_point() {
    update_with(h(get_x()), false);
}

template<unsigned short NB_ITER, class VECTOR, class FIXED_POINT>
void algorithm<NB_ITER, VECTOR, FIXED_POINT>::show_message(std::ostream& os) const {
    os << "############\n";
    os << "############\n";
    os << "Execution of Anderson\n";
    os << "Nb iterations : " << compteur << endl;
    //
    auto delta_X = compute_Delta_X();
    os << "Error : " << sqrt(std::inner_product(delta_X.begin(), delta_X.end(), delta_X.begin(), 0.)) << endl;
    //
    os << "############\n";
    os << "Sequence of fails | residues" << endl;
    os << std::boolalpha;
    for (size_t i = 0; i < sequence_of_residues.size(); i++) {
        os << "    " << sequence_of_successes[i] << "  " << sequence_of_residues[i] << endl;
    }
    os << "############\n";
    os << "############\n";
}

} // namespace Anderson
} // namespace optimization
//! Copyright : see license.txt
//!
//! \brief

#pragma once
#include <math.h>
#include <iostream>
#include <numeric>
#include <vector>

namespace optimization {

using namespace std;

namespace Barzilai_Borwein {

template<class VECTOR, class NABLA_F>
class algorithm {
    //! find the minimum of a given function f : \R^d -> \R, given its gradient
    //!
    //! implements the Barzilai-Borwein algorithm described in
    //! Fast methods for computing centroidal Laguerre tessellations for
    //! prescribed volume fractions with applications to microstructure
    //! generation of polycrystalline materials (Kuhn & al)

public:
    //! @brief : constructor
    //! @param x_0 : initial guess
    //! @param compute_nabla_f_ : lambda for computing nabla_f
    //! @param first_alpha : first parameter for gradient descent
    algorithm(const VECTOR& x_0, NABLA_F compute_nabla_f_, double first_alpha = 1);

    //! @brief : proceed with the algorith
    //! @param max_iter : maximum nb of iterations
    //! @param stopping_criterion : test algorithm success
    //! @return : if algorithm succeeded
    template<class STOPPING_TEST>
    bool proceed(size_t max_iter, STOPPING_TEST stopping_criterion);

    //! @brief : proceed with the algorith
    //! @param max_iter : maximum nb of iterations
    //! @param stopping_criterion : max distance between 2 subsequent values
    //! @return : if algorithm succeeded
    bool proceed(size_t max_iter, double stopping_criterion);
    //! @brief displays more on the algorithm
    void show_message(std::ostream& os) const;
    //! @return : the obtained minimizer
    const VECTOR& get_x() const { return x_n; }
    //! @brief : BB iteration of gradient descent
    void iterate();
    //! @brief : standard iteration of gradient descent with fixed parameter
    //! @param alpha_n : step
    void iterate_gradient_descent(double alpha_n);

private:
    size_t compteur;
    NABLA_F compute_nabla_f;
    VECTOR x_n;
    VECTOR x_n_1;
    VECTOR nabla_f;
    VECTOR nabla_f_1;
};

namespace auxi {
//! @return the coefficient in gradient descent computed by Barzilai-Borwein
//! @param x_n
//! @param x_n_1
//! @param nabla_f
//! @param nabla_f_1
template<class VECTOR>
double compute_alpha_n(const VECTOR& x_n, const VECTOR& x_n_1, const VECTOR& nabla_f, const VECTOR& nabla_f_1);
}  // namespace  auxi

}  // namespace Barzilai_Borwein

}  // namespace  optimization

#include "BarzilaiBorwein.ixx"

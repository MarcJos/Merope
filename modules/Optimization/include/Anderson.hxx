//! Copyright : see license.txt
//!
//! \brief

#pragma once


#include <iostream>
#include <numeric>
#include <vector>
#include <math.h>

#include "../../Eigen/Dense"


namespace optimization {

using namespace std;

namespace Anderson {

template<unsigned short NB_ITER, class VECTOR, class FIXED_POINT>
class algorithm {
    //! algorithm for finding x such that h(x) = x by means of the Anderson algorithm
    //! see Iterative residual-based vector methods to accelerate fixed
    // point iterations, Ramiere & Helfer, 2015
public:
    algorithm(VECTOR w_0, FIXED_POINT h_);
    //! @brief : const getter
    const VECTOR& get_x() const { return x_n[NB_ITER]; }
    //! @brief proceed with the Anderson algorithm
    //! @return wether the algorithm was successful
    template<class STOPPING_TEST>
    bool proceed(size_t max_iter, STOPPING_TEST stopping_criterion, bool use_acceleration = true);
    //! @brief proceed with the Anderson algorithm
    //! @return wether the algorithm was successful
    //! @param stopping_tolerance : the difference ||h(X)-X|| < stopping_tolerance
    bool proceed(size_t max_iter, double stopping_tolerance, bool use_acceleration = true);
    //! @brief displays a message
    void show_message(std::ostream& os) const;


private:
    //! @brief does 1 iteration of the Anderson method. Probably more precise and efficient than iterate_worse
    void iterate();
    //! @brief does 1 iteration of the Anderson method
    void iterate_worse();
    //! @brief does 1 iteration of fixed point
    void iterate_fixed_point();
    //! compute the vector Delta_X = current_h - current_w
    VECTOR compute_Delta_X() const;
    //! @brief update the memories x_n and h_n
    //! @param new_x_n : new x_n
    //! @param use_acceleration : has used acceleration or not for updating?
    void update_with(const VECTOR& new_x_n, bool use_acceleration);
    //! @brief counts the nb of applications of iterate
    size_t compteur;
    //! @brief application of which we attempt to find a fixed point
    FIXED_POINT h;
    //! @brief stores memory of x
    vector<VECTOR> x_n;
    //! @brief stores memory of h(x)
    vector<VECTOR> h_n;
    //! @brief size of the vector x
    size_t size;
    //! @brief : (optional) stores the sequence of residues
    vector<double> sequence_of_residues;
    //! @brief : (optional) stores the sequence of Anderson success (fails = false)
    vector<bool> sequence_of_successes;
};


} // namespace  Anderson

} // namespace optimization


#include "Anderson.ixx"
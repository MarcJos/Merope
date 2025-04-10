//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

namespace optimization {
namespace tests {

//! f = sin(|x|^2)
double f(const std::vector<double>& x);

// its gradient
std::vector<double> nabla_f(const std::vector<double>& x);

void test1_BarzilaiBorwein();

void test2_BarzilaiBorwein();

void print(auto my_vector);

template<bool USE_ANDERSON, int NB_ITER>
void test_Anderson();

}  // namespace  tests
}  // namespace  optimization

#include "Tests.ixx"

//! Copyright : see license.txt
//!
//! \brief

#include "../include/Tests.hxx"

using namespace std;

int main() {
    optimization::tests::test1_BarzilaiBorwein();
    optimization::tests::test2_BarzilaiBorwein();
    optimization::tests::test_Anderson<true, 2>();
    optimization::tests::test_Anderson<true, 3>();
    optimization::tests::test_Anderson<true, 4>();
    optimization::tests::test_Anderson<false, 2>();
    return 0;
}

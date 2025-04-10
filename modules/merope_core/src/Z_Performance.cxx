//! Copyright : see license.txt
//!
//! \brief


#include "Test/Tests.hxx"

using namespace std;
using namespace sac_de_billes;
using namespace merope;

int main() {
    clock_t t_0 = clock();
    Tests::polyCrystal9();
    cerr << "Time : " << (static_cast<double>(clock() - t_0)) / CLOCKS_PER_SEC << endl;
    t_0 = clock();
    Tests::polyCrystal8();

    cerr << "Time : " << (static_cast<double>(clock() - t_0)) / CLOCKS_PER_SEC << endl;
    return 0;
}

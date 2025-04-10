//! Copyright : see license.txt
//!
//! \brief

#include "../../GenericTools/CPP_Functions.hxx"

using namespace std;
using namespace merope::cppFunctions;


int main() {
    vector<pair<int, int>> pairs = { {1,2}, {1, 3}, {1, 4}, {3, 7}, {4, 8}, {5, 6} };
    auto result = merope::cppFunctions::computeEquivalenceClass<int>(pairs);
    for (size_t i = 0; i < result.size(); i++) {
        const auto& r = result[i];
        std::cout << "Class " << i << " : ";
        for (const auto& id : r) {
            std::cout << id << " ";
        }
        auto minimum = *(std::min_element(r.begin(), r.end()));
        std::cout << " represented by " << minimum;
        std::cout << std::endl;
    }
    // verify
    vector<vector<int>> result_expected = { {8,7,4,3,2,1}, {6,5} };
    for (size_t i = 0; i < result_expected.size(); i++) {
        const auto& vect_sets = result_expected[i];
        for (const auto& elem : vect_sets) {
            if (result[i].find(elem) == result[i].end()) {
                return EXIT_FAILURE;
            }
        }
    }
    return EXIT_SUCCESS;
}
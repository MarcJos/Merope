//! Copyright : see license.txt
//!
//! \briefPurely informatic function for doing some C++ manipulations
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {
namespace cppFunctions {

template<class C>
vector<const C*> fromVectorToConstVectorPointer(const vector<C>& input) {
    vector<const C*> result{};
    result.reserve(input.size());
    for (const C& element : input) {
        result.push_back(&element);
    }
    return result;
}

template<class FUNCTION, class NAMEFUNCTION>
double print_time_execution(FUNCTION f, NAMEFUNCTION nameFunction, vector<int> liste_omp_nb_threads = vector<int>({ 1, 2, 4, 8, 16, 32 })) {
    cerr << static_cast<string>(nameFunction);
    double result = 0;
    for (auto nbProc : liste_omp_nb_threads) {
        omp_set_num_threads(nbProc);
#pragma omp parallel
        {  // begin parallel section
#pragma omp master
            {
                cerr << " Number of threads " << omp_get_num_threads() << " ";
            }
        }  // end parallel section
        auto  t_0 = std::chrono::steady_clock::now();
        f();
        std::chrono::duration<double> elapsed_seconds = std::chrono::steady_clock::now() - t_0;
        //
        double result = static_cast<double>(elapsed_seconds.count());
        cerr << "Time " << result << endl;
    }
    return result;
}

template<class STRING>
void verify(bool test, STRING __pretty_function__, std::string message) {
    if (not test) {
        cerr << __pretty_function__ << endl;
        throw runtime_error(message);
    }
}

#define VERIFIE_SI(test, message) cppFunctions::verify(test, __PRETTY_FUNCTION__, message)
// French name to avoid problems

}  // namespace cppFunctions
}  // namespace merope




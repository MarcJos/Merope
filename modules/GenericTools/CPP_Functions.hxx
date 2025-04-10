//! Copyright : see license.txt
//!
//! \briefPurely informatic function for doing some C++ manipulations
//
#pragma once


#include "../../GenericMerope/StdHeaders.hxx"

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
        double elapsed_time = static_cast<double>(elapsed_seconds.count());
        cerr << "Time " << elapsed_time << endl;
    }
    return result;
}

class my_timer {
public:
    my_timer() : t_0(clock()) {}

    double get_time() {
        t_0 = clock();
        return static_cast<double>(clock() - t_0) / CLOCKS_PER_SEC;
    }

    void write_time(string s) { std::cerr << s << " " << get_time() << endl; }

private:
    clock_t t_0;
};


//! \brief Computes minimal equivalence classes from given pairs of elements.
//! This function template takes a vector of pairs representing equivalence relationships and computes the minimal set of equivalence classes using a Union-Find data structure.
//! \tparam T Type of elements in pairs (e.g., int, char).
//! \param equivalencePairs Vector of pairs representing equivalence relationships between elements.
//! \return A vector of unordered_set<T>, where each set contains the members of a single equivalence class.
template<typename T>
vector<unordered_set<T>> computeEquivalenceClass(const vector<pair<T, T>>& equivalencePairs);


}  // namespace cppFunctions
}  // namespace merope

#include "CPP_Functions.ixx"


//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../GenericMerope/StdHeaders.hxx"
#include "../../Geometry/include/Point.hxx"

using namespace sac_de_billes;

/////////////////////
//! One cannot trust pybind for binding lambda functions with C++ functions.
//! It does not properly interact with OpenMP.
//! It is allowed to do this only if OpenMP is disabled for Fields.
//! In the other case, use Mamba to get a pointer to a C function
//! which is then used at the proper place through a converter "Interf_FuncPointer".
/////////////////////

// prevent the use of raw python lambda if openmp is used
#ifdef USE_OPENMP_FOR_FFT
static constexpr bool USE_OPENMP_FOR_FFT_LOC = true;
#else
static constexpr bool USE_OPENMP_FOR_FFT_LOC = false;
#endif


inline void error_use_omp() {
    if constexpr (USE_OPENMP_FOR_FFT_LOC) {
        throw runtime_error("You shall : \n 1) compile with another option for using this function.See USE_OPENMP_FOR_FFT \n OR \n 2) use pointer to functions");
    }
}


//! @brief : class for storing function : \R^dim_from x \Z ->  \R^dim_to
struct Interf_TexturePointer {
    //! @brief  constructor
    //! @param pointer_address_ : pointer adress
    Interf_TexturePointer(std::size_t pointer_address_, std::array<std::size_t, 2>dim_from_dim_to) :
        pointer_address{ pointer_address_ }, dim_from{ dim_from_dim_to[0] }, dim_to{ dim_from_dim_to[1] } {}

    template<std::size_t DIM>
    using POINT_DOUBLE = std::conditional_t<DIM == 1, double, Point<DIM>>;

    template<std::size_t DIM>
    using INTERMEDIATE_POINT_FORMAT = std::conditional_t<DIM == 1, double, double*>;

    //! @brief turn the pointer into a function R -> R
    //! @return the nonlinear function
    template<std::size_t DIM_FROM, std::size_t DIM_TO>
    inline std::function<POINT_DOUBLE<DIM_TO>(POINT_DOUBLE<DIM_FROM>, long)> get() const {
        if (dim_from != DIM_FROM or dim_to != DIM_TO) {
            std::cerr << __PRETTY_FUNCTION__ << "\n";
            std::cerr << "Was requested a function from R**" << DIM_FROM << " to R**" << DIM_TO << "\n";
            std::cerr << "But can only  a function from R**" << dim_from << " to R**" << dim_to << "\n";
            throw runtime_error("Uncompatible dimensions!");
        }
        if (DIM_TO > 1) {
            cerr << "Requested a function with output in dimension>1.";
            Merope_error_not_done();
        }
        auto function_pointer = reinterpret_cast<INTERMEDIATE_POINT_FORMAT<DIM_TO>(*)(INTERMEDIATE_POINT_FORMAT<DIM_FROM>, long)>(
            reinterpret_cast<void*>(pointer_address));
        auto final_function = [function_pointer](POINT_DOUBLE<DIM_FROM> pt, long i) {
            if constexpr (DIM_FROM == 1) {
                return function_pointer(pt, i);
            } else {
                return function_pointer(pt.data(), i);
            }
            };
        return final_function;
    }

private:
    //! @brief adress of the pointer
    const std::size_t pointer_address;
    const std::size_t dim_from;
    const std::size_t dim_to;
};

//! @brief : class for storing function : \R^dim_from ->  \R^dim_to
struct Interf_FuncPointer {
    //! @brief  constructor
    //! @param pointer_address_ : pointer adress
    Interf_FuncPointer(std::size_t pointer_address_, std::array<std::size_t, 2>dim_from_dim_to) :
        pointer_address{ pointer_address_ }, dim_from{ dim_from_dim_to[0] }, dim_to{ dim_from_dim_to[1] } {}

    template<std::size_t DIM>
    using POINT_DOUBLE = std::conditional_t<DIM == 1, double, Point<DIM>>;

    template<std::size_t DIM>
    using INTERMEDIATE_POINT_FORMAT = std::conditional_t<DIM == 1, double, double*>;

    //! @brief turn the pointer into a function R -> R
    //! @return the nonlinear function
    template<std::size_t DIM_FROM, std::size_t DIM_TO>
    inline std::function<POINT_DOUBLE<DIM_TO>(POINT_DOUBLE<DIM_FROM>)> get() const {
        if (dim_from != DIM_FROM or dim_to != DIM_TO) {
            std::cerr << __PRETTY_FUNCTION__ << "\n";
            std::cerr << "Was requested a function from R**" << DIM_FROM << " to R**" << DIM_TO << "\n";
            std::cerr << "But can only  a function from R**" << dim_from << " to R**" << dim_to << "\n";
            throw runtime_error("Uncompatible dimensions!");
        }
        if (DIM_TO > 1) {
            cerr << "Requested a function with output in dimension>1.\n";
            Merope_error_not_done();
        }
        auto function_pointer = reinterpret_cast<INTERMEDIATE_POINT_FORMAT<DIM_TO>(*)(INTERMEDIATE_POINT_FORMAT<DIM_FROM>)>(
            reinterpret_cast<void*>(pointer_address));
        auto final_function = [function_pointer](POINT_DOUBLE<DIM_FROM> pt) {
            if constexpr (DIM_FROM == 1) {
                return function_pointer(pt);
            } else {
                return function_pointer(pt.data());
            }
            };
        return final_function;
    }

private:
    //! @brief adress of the pointer
    const std::size_t pointer_address;
    const std::size_t dim_from;
    const std::size_t dim_to;
};

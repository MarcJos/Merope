//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GAUSSIANFIELD_HXX_
#define GAUSSIANFIELD_HXX_

#include "../../../AlgoPacking/src/StdHeaders.hxx"


#include "../../../AlgoPacking/src/Geometry/GeomTypes.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {

class Grid;

namespace gaussianField {

struct Interf_FuncPointer{
    Interf_FuncPointer(size_t address_, bool):address{address_}{}

    const size_t address;
};

template<unsigned short DIM>
class CovarianceField;

//! simulate a nonlinear transform of a gaussian field A(g), with g a gaussian field of zero mean and given covariance
template<unsigned short DIM>
class SimpleGaussianField {
public:
    //! constructor
    SimpleGaussianField(std::function<double(Point<DIM>)> covariance_,
        std::function<double(double)> nonlinearFunction_) :
        covariance(covariance_), nonlinearFunction(nonlinearFunction_), seed{ 0 } {}
    SimpleGaussianField(void * cov, void * nonlin):
	SimpleGaussianField(reinterpret_cast<double(*)(double*)>(cov), reinterpret_cast<double(*)(double)>(nonlin)){}
    SimpleGaussianField(double(*cov)(double*), double(*nonlin)(double)):
	covariance([cov](Point<DIM> x){return cov(x.data());}), nonlinearFunction(nonlin), seed{ 0 } {}
    SimpleGaussianField(Interf_FuncPointer i1, Interf_FuncPointer i2): SimpleGaussianField(reinterpret_cast<void*>(i1.address), reinterpret_cast<void*>(i2.address)) {}
    //! covariance function
    std::function<double(Point<DIM>)> covariance;
    //! nonlinear function
    std::function<double(double)> nonlinearFunction;
    //! seed for the standard white noise
    int seed;

public:
    template<bool use_openmp_for_fft>
    static SimpleGaussianField create(std::function<double(Point<DIM>)> covariance_,
        std::function<double(double)> nonlinearFunction_){
	if constexpr (use_openmp_for_fft){
	    throw runtime_error("You shall compile with another option for using this function. See USE_OPENMP_FOR_FFT");
	} else {
	    return SimpleGaussianField(covariance_, nonlinearFunction_);
	}
    }



public:
    //! setter (need a copy of the gaussian::RV, which is of abstract type)
    template<class GaussianRV>
    void setNonlinearTransform(const GaussianRV& nonlinearTransform_);
};

//! @brief object for voxellizing the covariance field
template<unsigned short DIM>
class NumericalCovariance {
public:
    NumericalCovariance(std::function<double(Point<DIM>)> cov) : covariance(cov) {};
    //! covariance function
    std::function<double(Point<DIM>)> covariance;
};

// template version
template<unsigned short DIM, class FFTScal, class T>
vector<double> createField_T(const T& gaussianField, const Grid& grid, size_t size_);
//! create a field
vector<double> createField(const NumericalCovariance<3>& covariance, const Grid& grid, size_t size_);
vector<double> createField(const NumericalCovariance<2>& covariance, const Grid& grid, size_t size_);
vector<double> createField(const SimpleGaussianField<3>& gaussianField, const Grid& grid, size_t size_);
vector<double> createField(const SimpleGaussianField<2>& gaussianField, const Grid& grid, size_t size_);


} // namespace gaussianField
} // namespace merope

#include "../Field/GaussianField.ixx"

#endif /* GAUSSIANFIELD_HXX_ */

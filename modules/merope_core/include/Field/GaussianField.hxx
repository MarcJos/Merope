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
    //! covariance function
    std::function<double(Point<DIM>)> covariance;
    //! nonlinear function
    std::function<double(double)> nonlinearFunction;
    //! seed for the standard white noise
    int seed;


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

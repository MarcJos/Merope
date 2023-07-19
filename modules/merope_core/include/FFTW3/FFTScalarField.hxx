//! Copyright : see license.txt
//!
//! \brief Scalar FFT fields with FFTW
//!

#ifndef _FFTSCALARFIELD_HXX
#define _FFTSCALARFIELD_HXX 1

#include<functional>
#include<array>

#include "../FFTW3/FFTField.hxx"



#include "../MeropeNamespace.hxx"


namespace merope {

namespace gaussianField {
class CovSum;
} // gaussianField

class VTKstream;
//! Scalars field based on FFT fields
class FFTScalarField : public FFTField {
public:
    //! Constructor
    //! \param grid Grid
    //! \param IP True for in place memory allocation
    //! \param flags plane construction method
    FFTScalarField(const Grid& grid, bool IP, unsigned flags = FFTW_ESTIMATE);
    //! Constructor with covariance
    //! \param grid Grid
    //! \param cs Covariance function
    //! \param IP True for in place memory allocation
    //! \param flags plane construction method
    FFTScalarField(const Grid& grid, const gaussianField::CovSum& cs, bool IP, int seed, unsigned flags =
        FFTW_ESTIMATE);
    //! Constructor with covariance
    //! \param grid : Grid
    //! \param cs : Covariance function
    //! \param IP : True for in place memory allocation
    //! \param showCovariance : only show the covariance function, not the gaussian field
    //! \param seed : seed for random generator
    //! \param flags : plane construction method
    FFTScalarField(const Grid& grid, const std::function<double(std::array<double, 3>)>& cs, bool IP,
        int seed, bool showCovariance = false,
        unsigned flags = FFTW_ESTIMATE);
    FFTScalarField(const Grid& grid, const std::function<double(std::array<double, 2>)>& cs, bool IP,
        int seed, bool showCovariance = false,
        unsigned flags = FFTW_ESTIMATE);
    //! Writes a VTK CELL file
    //! \param fvtk VTK file stream
    //! \param cname vector coordinate name
    void toVTKCELL(VTKstream& fvtk, const char* cname) const;
private:
    //! Fill the field with a covariance function
    //! \param Lx,Ly,Lz Grid sizes
    //! \param cs Covariance function
    template<class C>
    void setCov(double Lx, double Ly, double Lz, const C& cs);
    //! Transform a covariance into a Random Function
    void randFunc(int seed);
    //! @see : constructor
    template<class COVARIANCE_TYPE>
    void build(const Grid& grid, const COVARIANCE_TYPE& cs, bool IP, int seed,
        unsigned flags, bool showCovariance);
    //! put the Covariance in Fourier (only retaining the positive part)
    void prepareCovarianceInFourier();
    //! apply the function on each fourier mode
    template<class FUNCTION>
    void loopOnFrequencies(const FUNCTION& function);
};

} // namespace merope

#include "../FFTW3/FFTScalarField.ixx"

#endif /* _FFTSCALARFIELD_HXX */


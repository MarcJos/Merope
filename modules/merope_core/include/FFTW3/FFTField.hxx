//! Copyright : see license.txt
//!
//! \brief FFT fields based on FFTW library

#ifndef _FFT_FIELD_HXX
#define _FFT_FIELD_HXX 1

#include<vector>
#include<iostream>

#if defined(_WIN32) || defined(WIN32)  // ugly for Eclipse
#include "../../../FFTW_Win/fftw3.h"
#else
#include<fftw/fftw3.h>
#endif

#include<assert.h>

#include"../FFTW3/NumericTypes.h"
#include "../localMFront.h"



#include "../MeropeNamespace.hxx"


namespace merope {

class Grid;

//----------------------------------------------
//! FFT field mother class
class FFTField {
protected:
    unsigned char d;  //!< Space dimension
    unsigned char nv; //!< Number of variables
    size_t fSize;  //!< Number of spatial nodes
    size_t FSize;  //!< Number of spectral nodes
    bool isSpatial; //!< Current representation (spatial or spectral)

    int n[3]; //!< Grid dimensions
    unsigned nz, nzb; //!< Last real grid dimension (for in-place FFT)

    rfloat* f; //!< Spatial field
    cfloat* F; //!< Spectral field

private:
    FFTW_PREF(plan) forwardPlan;  //!< FFT forward  plane
    FFTW_PREF(plan) backwardPlan; //!< FFT backward plane
    unsigned flags; //!< FFT plane type

public:
    void forward();  //!< Direct transformation (with normalization)
    void backward(); //!< Reverse transformation

    //! Print the field values
    //! \param os Output flux
    void print(std::ostream& os) const;
    //! Initialize a uniform field in spatial space
    void setSpatialField();
    //! Initialize a uniform field in spectral space
    void setSpectralField();
    //! Initialize a uniform field in spectral space
    //! \param values values to set at null frequency
    void setSpectralField(const CastemReal* values);
    //! Cancel the field mean in spectral representation
    void SpectralZeroMean();
    //! Reset spectral field
    //! \param values values to set at null frequency
    //! \param k coeff to apply to all other frequencies values
    void resetSpectralField(const CastemReal* values, const double k);

    cfloat* getSpectralField(); //!< getter on spectral field
    const cfloat* getSpectralField() const; //!< const getter on spectral field
    //! Spatial vector at node i
    //! \param i Node number
    //! \return Spatial vector at node i
    rfloat* getSpatialField(const unsigned i);
    const rfloat* getSpatialField(const unsigned i) const;
    //! \return number of variables
    unsigned char getNv() const;
    //! \return Number of spatial nodes
    size_t getfSize() const;

    //! affectation operator that allocate plans if not done yet
    //! This affectation is to be used only
    //!    for fields that have the sames sizes
    FFTField& operator=(const FFTField& field);
    //! Subtraction of two fields
    //! \param ff1,ff2 FFTFields
    void minus(const FFTField& ff1, const FFTField& ff2);
    //! Subtraction to the current field
    //! \param ff2 FFTField to substract
    FFTField& operator-=(const FFTField& ff2);
    //! Addition of two fields
    //! \param ff1,ff2 FFTFields
    void plus(const FFTField& ff1, const FFTField& ff2);
    //! Addition to the current field
    //! \param ff2 FFTField to add
    FFTField& operator+=(const FFTField& ff2);
    //! Field comparison
    //! \param ff2 Second fft field to be compared
    //! \return The mean square difference
    long double compare(const FFTField& ff2) const;
    //! Scalar product
    //! \param ff2 Other FFTfield
    //! \return the scalar product
    long double scalarProduct(const FFTField& ff2) const;
    //! Field mean square
    //! \return The mean square
    long double norm2() const;
    //! Linear combination of FFTFields
    //! \param V Table of FFT Fields
    //! \param l Weights
    //! \param N Combination length
    void linComb(const FFTField* const* V, const long double* l,
        unsigned char N);
    //! Destructor
    virtual ~FFTField();
protected:
    //! Contructor
    //! \param grid Grid
    //! \param flags plane construction method
    FFTField(const Grid& grid, unsigned flags);
    //! Memory allocation
    //! \param IP True for in place memory allocation
    void alloc(bool IP);

    //! Get the grid dimensions
    //! \param nx,ny First two grid dimensions
    void getDim(size_t& nx, size_t& ny) const;

    //!< Spatial  field allocation
    void spatialAlloc();
    //!< Spectral field allocation
    void spectralAlloc();

    //!< check that field is in spatial  representation, and exists
    void checkSpatial(const std::string& mName) const;
    //!< check that field is in spectral representation, and exists
    void checkSpectral(const std::string& mName) const;

    //!< assert that field is in spatial  representation, and exists
    inline void assertSpatial() const {
        assert(isSpatial && std::string("FILE: ") + std::string(__FILE__) == "");
        assert(f);
    }

    //!< assert that field is in spectral representation, and exists
    inline void assertSpectral() const {
        assert(!isSpatial);
        assert(F);
    }

private:
    FFTField(const FFTField&); //!< Forbid copy constructor

    void buildForwardPlan();  //!< Direct  transformation plane construction
    void buildBackwardPlan(); //!< Reverse transformation plane construction
    //! Spatial scalar product
    //! \param ff2 Other FFTfield
    //! \return the scalar product
    long double SpatialScalarProduct(const FFTField& ff2) const;
    //! Spectral scalar product
    //! \param ff2 Other FFTfield
    //! \return the scalar product
    long double SpectralScalarProduct(const FFTField& ff2) const;
};

} // namespace merope

#endif // _FFT_FIELD_HXX


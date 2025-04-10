//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

namespace merope {
namespace gaussianField {

//! Abstract class for a covariance function
class CovIso {
protected:
    double r;  //!< Range
    double s;  //!< Sill
public:
    //! Default Destructor
    virtual ~CovIso();
    //! 1D Covariance Function
    //! \param x Coordinate
    double cov(double x) const;
    //! 2D Covariance Function
    //! \param x,y Coordinates
    double cov(double x, double y) const;
    //! 3D Covariance Function
    //! \param x,y,z Coordinates
    double cov(double x, double y, double z) const;

protected:
    //! Default constructor
    //! \param r Range
    //! \param s Sill
    CovIso(double r, double s);
    //! Isotropic covariance function
    //! \param h Distance
    virtual double covIso(double h) const = 0;
};

enum class CovType {
    Spherical,
    Gaussian,
    Exponential
};

//! Exponential covariance function
class CovExpo : public CovIso {
public:
    //! Default constructor
    //! \param r Range
    //! \param s Sill
    CovExpo(double r, double s);
    //! Isotropic covariance function
    //! \param h Distance
    double covIso(double h) const;
};

//! Spherical covariance function
class CovSpherical : public CovIso {
public:
    //! Default constructor
    //! \param r Range
    //! \param s Sill
    CovSpherical(double r, double s);
    //! Isotropic covariance function
    //! \param h Distance
    double covIso(double h) const;
};

//! Gaussian covariance function
class CovGaussian : public CovIso {
public:
    //! Default constructor
    //! \param r Range
    //! \param s Sill
    CovGaussian(double r, double s);
    //! Isotropic covariance function
    //! \param h Distance
    double covIso(double h) const;
};

//! Covariance sum
class CovSum : public std::vector<CovIso*> {
public:
    //! destructor
    ~CovSum();
    //! Add a covariance function
    //! \param covType : Covariance function
    //! \param parameters : model parameters
    void add(CovType covType, vector<double> modelParameters);
    //! 1D Covariance Function
    //! \param x Coordinate
    double cov(double x) const;
    //! 2D Covariance Function
    //! \param x,y Coordinates
    double cov(double x, double y) const;
    //! 3D Covariance Function
    //! \param x,y,z Coordinates
    double cov(double x, double y, double z) const;

private:
    //! Add a covariance function
    //! \param c Covariance function
    void add(CovIso* c);
};

}  // namespace gaussianField
}  // namespace merope



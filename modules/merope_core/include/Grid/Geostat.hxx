//! Copyright : see license.txt
//!
//! \brief Geostatistic data needed to simulate a periodic medium
//!

#ifndef _GEOSTAT_HXX
#define _GEOSTAT_HXX 1


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {
namespace gaussianField {

//! virtual Random Variable
class RV {
public:
    //! Destructor
    virtual ~RV();
    //! Inverse anamorphosis
    //! \param x Value
    double inverseA(double x) const;
    //! Initialize the stepwise approximation
    virtual void setApprox() const = 0;
private:
    //! Inverse of the cumulative distribution function
    virtual double inverseF(double) const = 0;
};

//! Positive step Probability distribution
class StepDis: public RV {
    std::vector<double> xi; //!< Interval limits
    std::vector<double> Fi; //!< Cumulative distribution function
public:
    //! Default constructor
    StepDis();
    //! Constructor by reading a file
    //! \param nom File name
    StepDis(std::string nom);
    //! \return its size
    size_t size() const;
    //! Clear the distribution
    void clear();
    //! Stepwise approximation of a PDF
    //! M RV type
    //! \param mod Random variable
    //! \param x0 Discretizing Beginning
    //! \param dx Discretizing step
    //! \param eps Precision
    template<typename M>
    void FromModel(const M& mod, double x0, double dx, double eps = 1.e-5);
    //! Initialize the stepwise approximation
    void setApprox() const;
    //! \return the inverse of the cumulative distribution function
    //! \param x Value
    double inverseF(double x) const;
    //! get
    double getMin() const { return xi[0]; };
    double getMax() const { return xi[xi.size() - 1]; };
private:
    //! Read from file
    //! \param nom File name
    void read(std::string nom);
};

//! Gaussian Random Variable
class GaussRV {
protected:
    double m; //!< Mean
    double s; //!< standard deviation
    double w;  //!< Weight in a sum
public:
    //! Default constructor
    //! \param m Mean
    //! \param s standard deviation
    //! \param w Weight in a sum
    GaussRV(double m, double s, double w = 1);
    //! Destructor
    virtual ~GaussRV();
    //! Probability density function
    //! \param x Value
    virtual double PDF(double x) const;
    //! Cumulative distribution function
    //! \param x Value
    virtual double CDF(double x) const;
};

//! Truncated Gaussian
class TGauss: public GaussRV {
    double l; //!< Threshhold
    double A, Fl; //!< Norm
public:
    //! Default constructor of a truncated Gaussian random variable
    //! \param m Mean
    //! \param s standard deviation
    //! \param w Weight in a sum
    TGauss(double m, double s, double w = 1);
    //! Probability density function
    //! \param x Value
    double PDF(double x) const;
    //! Cumulative distribution function
    //! \param x Value
    double CDF(double x) const;
private:
    //! Auxiliary parameters set up
    void MAJ();
};

//! Lognormal Random Variable
class LogNorm: public GaussRV {
public:
    //! Default constructor of a Lognormal random variable
    //! \param m Mean
    //! \param s standard deviation
    //! \param w Weight in a sum
    LogNorm(double, double, double = 1);
    //! Probability density function
    //! \param x Value
    double PDF(double x) const;
    //! Cumulative distribution function
    //! \param x Value
    double CDF(double x) const;
};

//! Sum of Elementary Random Variables
class SofERV: public RV, public std::vector<GaussRV*> {
    mutable StepDis sd; //!< Step wise approximation
public:
    //! Destructor
    virtual ~SofERV();
    //! Add a new Elementary Random variable
    //! \param e Elementary Random variable
    void add(GaussRV* e);
    //! Probability density function
    //! \param x Value
    double PDF(double x) const;
    //! Cumulative distribution function
    //! \param x Value
    double CDF(double x) const;
    //! Initialize the stepwise approximation
    void setApprox() const;
private:
    //! \param Inverse of the cumulative distribution function
    //! \param x Value
    double inverseF(double x) const;
};

} // namespace gaussianField
} // namespace merope

#endif /* _GEOSTAT_HXX */


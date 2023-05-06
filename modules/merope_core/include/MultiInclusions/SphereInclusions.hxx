//! Copyright : see license.txt
//!
//! \brief

#ifndef SPHEREINCLUSIONS_HXX_
#define SPHEREINCLUSIONS_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/Interface.hxx"
#include "../Geometry/GeomTools.hxx"
#include "../MesoStructure/InsideTorus.hxx"
#include "../MultiInclusions/SphereSeeds.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM>
class SphereInclusions: public SphereSeeds<DIM> {
    //! \brief Class for modelling spherical inclusions inside a matrix
    //! 0 is the matrix phase
    //! 1->N are the inclusion phases
public:
    // default constructor
    SphereInclusions():
        SphereSeeds<DIM>() {}

    //! return the covariances of spheres
    void covX(const unsigned nx, std::ostream& fout) const; //!< Theoretical covariance (X direction)
    void covY(const unsigned ny, std::ostream& fout) const; //!< Theoretical covariance (Y direction)
    void covZ(const unsigned nz, std::ostream& fout) const; //!< Theoretical covariance (Z direction)
    //! Getter on number of phases in the geometry.
    inline unsigned long getNumberOfPhases() const {
        return 1 + this->theSpheres.getNbPhases();
    }
    // +1 corresponds to the matrix phase
private:
    //! Theoretical covariance in direction i
    //! \param nx Spatial discretisation
    //! \param fout Output stream
    void covDirection(size_t direction, const unsigned n,
        std::ostream& fout) const;

};

} // namespace merope


#include "../MultiInclusions/SphereInclusions.ixx"
#endif /* SPHEREINCLUSIONS_HXX_ */

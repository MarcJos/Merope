//! Copyright : see license.txt
//!
//! \brief

#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/GeomTools.hxx"

#include "../../../AlgoPacking/include/Interface.hxx"

#include "../MesoStructure/InsideTorus.hxx"
#include "../MultiInclusions/SphereSeeds.hxx"


namespace merope {

template<unsigned short DIM>
class SphereInclusions final : public SphereSeeds<DIM> {
    //! \briefClass for modelling spherical inclusions inside a matrix
    //! 0 is the matrix phase
    //! 1->N are the inclusion phases
public:
    // default constructor
    SphereInclusions() :
        SphereSeeds<DIM>() {}

    //! return the covariances of spheres
    void covX(const unsigned nx, std::ostream& fout) const;  //!< Theoretical covariance (X direction)
    void covY(const unsigned ny, std::ostream& fout) const;  //!< Theoretical covariance (Y direction)
    void covZ(const unsigned nz, std::ostream& fout) const;  //!< Theoretical covariance (Z direction)
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

}  // namespace merope


#include "../MultiInclusions/SphereInclusions.ixx"


//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../../AlgoPacking/src/StdHeaders.hxx"


#include "../localMFront.h"


#include "../MeropeNamespace.hxx"


namespace merope {

namespace auxiRotations {
constexpr int PRINTED_DOUBLE_PRECISION = 13;
}

//! 2D Rotation
class Rotation2D {
    //! Rotation Matrix in the column-major order
    //! (FORTRAN/CAST3M convention)
    //! First vector Second vector
    //!  a[0]=a[1,1]  a[2]=a[1,2]
    //!  a[1]=a[2,1]  a[3]=a[2,2]
    CastemReal a[4];
public:
    //! Default constructor
    Rotation2D();
    //! Constructor from one vector
    //! \param V 1rst base vector of the material frame
    //! \param verif true if a verification is needed
    Rotation2D(const CastemReal* V, bool verif = false);
    //! Constructor from an angle
    //! \param theta Angle [rad]
    Rotation2D(double theta);
    //! Constructor from the Euler angles
    //! \param eng Random number engine that generates pseudo-random numbers
    void randomEuler(std::default_random_engine& eng);
    //! From Euler angles to rotaion matrix
    //! \param theta Euler angle [rad]
    void EulerToAxis(double theta);
    //! Print the first axis
    //! \param f Output file
    void print(std::ofstream& f) const;
private:
    //! \return the rotation matrix
    const CastemReal* getMat() const;
    //! Frame change (direct)
    //! \param vm Vector in the material frame
    //! \param vg Vector in the global frame
    void direct(const CastemReal* vm, CastemReal* vg) const;
    //! Frame change (reverse)
    //! \param vg Vector in the global frame
    //! \param vm Vector in the material frame
    void reverse(const CastemReal* vg, CastemReal* vm) const;
    //! Vectorial product
    void VecProd();
};

//! 3D Rotation
class Rotation3D {
    //! Rotation Matrix in the column-major order
    //! (FORTRAN/CAST3M convention)
    //! First vector Second vector Third vector
    //!  a[0]=a[1,1]   a[3]=a[1,2]  a[6]=a[1,3]
    //!  a[1]=a[2,1]   a[4]=a[2,2]  a[7]=a[2,3]
    //!  a[2]=a[3,1]   a[5]=a[3,2]  a[8]=a[3,3]
    CastemReal a[9];
public:
    //! Default constructor
    Rotation3D();
    //! Constructor from 2 vectors
    //! \param V 1rst and 2nd base vectors of the material frame
    //! \param verif true if a verification is needed
    Rotation3D(const CastemReal* V, bool verif = false);
    //! Constructor from the Euler angles
    //! \param phi,theta,psi Euler angles [rad]
    Rotation3D(double phi, double theta, double psi);
    //! Constructor from the Euler angles
    //! \param eng Random number engine that generates pseudo-random numbers
    void randomEuler(std::default_random_engine& eng);
    //! From Euler angles to rotation matrix
    //! \param phi,theta,psi Euler angles [rad]
    void EulerToAxis(double phi, double theta, double psi);
    //! Print first 2 axes
    //! \param f Output file
    void print(std::ofstream& f) const;
private:
    //! \return the rotation matrix
    const CastemReal* getMat() const;
    //! Frame change (direct)
    //! \param vm Vector in the material frame
    //! \param vg Vector in the global frame
    void direct(const CastemReal* vm, CastemReal* vg) const;
    //! Frame change (reverse)
    //! \param vg Vector in the global frame
    //! \param vm Vector in the material frame
    void reverse(const CastemReal* vg, CastemReal* vm) const;
    //! Vectorial product
    void VecProd();
};

}  // namespace merope




//! Copyright : see license.txt
//!
//! \brief Gathers tools for the surrounding shape
//! to be refactored : no virtual method should be here!

#ifndef GLOBALSHAPE_HXX_
#define GLOBALSHAPE_HXX_

#include "StdHeaders.hxx"

#include "AmbiantSpace.hxx"
#include "Geometry/GeomTools_1.hxx"

namespace sac_de_billes {
using namespace std;

namespace AmbiantSpace {

enum class NameShape {
    Tore, Cube, Sphere, Cylinder
};

inline NameShape readShape(string s) {
    if (s == "Tore") {
        return NameShape::Tore;
    } else if (s == "Cube") {
        return NameShape::Cube;
    } else if (s == "Sphere") {
        return NameShape::Sphere;
    } else if (s == "Cylinder") {
        return NameShape::Cylinder;
    } else {
        throw invalid_argument("readShapes");
    }
}

template<short unsigned DIM>
class BigShape {
    //! For the volume in which the spheres are placed
    //! SHOULD BE CONVEX!
public:
    //! space dimension
    short d;
    //! the lengths of the sides of the cube
    Point<DIM> L;

    //! Default Constructor
    BigShape(Point<DIM> L_) :
        d{ DIM }, L{ L_ }, boundaryExclusionDistance{ 0 } {}

    virtual ~BigShape() {}
    //! \set BoundaryExclusionDistance
    void setBoundaryExclusionDistance(const double& distance) {
        boundaryExclusionDistance = distance;
    }

    //! \return the norm of a vector
    double normeCarre(const Point<DIM>& vec) const;
    double norme(const Point<DIM>& vec) const;

    //------------------------------------
    // Virtual methods
    //------------------------------------
    //! \return Volume of the BigShape
    virtual double volume() const = 0;
    //! \return the squared distance between two points
    //! \param x1, x2 two points
    virtual double distanceCarre(const Point<DIM>& x1,
        const Point<DIM>& x2) const = 0;
    //! \return whether the point is inside the BigShape (taking into account the minRadius)
    virtual bool isInside(const Point<DIM>& point,
        const double& minRadius) const = 0;
    //! \return whether a sphere is inside the BigShape
    //! coincides with isInside above BUT for periodic geometry
    virtual bool isInside(const Sphere<DIM>& sph) const = 0;
    //! \return the vector p1p2 (nontrivial for periodic geometries)
    virtual Point<DIM> geomVector(const Point<DIM>& p1,
        const Point<DIM>& p2) const = 0;
    //! \return the sphere after it has bounced on the boundary
    virtual void bounce(Sphere<DIM>& sph) const = 0;
    //! type
    virtual NameShape type() const = 0;

    //------------------------------------
    // Deduced from the virtual methods
    //------------------------------------
    //! \return whether spheres are intersected
    bool areSphereIntersected(const Sphere<DIM>& sph1,
        const Sphere<DIM>& sph2) const;

protected:
    double boundaryExclusionDistance;
};

template<short unsigned DIM>
class Cube : public BigShape<DIM> {
public:
    // constructor
    //! \param L_ : the lengths of the sides of the cube
    Cube(Point<DIM> L_);

    //! \return the volume of the cube
    double volume() const override;
    //! \return the squared PERIODIC distance between two points
    //! \param x1, x2 two points
    double distanceCarre(const Point<DIM>& x1,
        const Point<DIM>& x2) const override;
    //! \return whether the point is inside the cube
    bool isInside(const Point<DIM>& point, const double& minRadius) const override;
    //! a sphere may cross the boundary of the cube
    bool isInside(const Sphere<DIM>& sph) const override;
    //! \see BigShape::geomVector
    Point<DIM> geomVector(const Point<DIM>& p1,
        const Point<DIM>& p2) const override;
    //! \return the projection of a point inside the Cube
    virtual void projection(Point<DIM>& point) const;
    //! \see BigShape::bounce
    void bounce(Sphere<DIM>& sph) const override;
    //! \see BigShape::type
    NameShape type() const override {
        return NameShape::Cube;
    }
    //! return the area of a given face
    double faceArea(int direction);
};

template<short unsigned DIM>
class Tore : public Cube<DIM> {
public:
    // constructor
    Tore(Point<DIM> L_);

    //! \return the squared PERIODIC distance between two points
    //! \param x1, x2 two points
    double distanceCarre(const Point<DIM>& x1, const Point<DIM>& x2) const
        override;
    //! \return whether the point is inside the bigShape
    bool isInside(const Point<DIM>& point, const double& minRadius) const
        override;
    //! a sphere is always inside the bigShape
    bool isInside(const Sphere<DIM>& sph) const override;
    //! \see BigShape::geomVector
    Point<DIM> geomVector(const Point<DIM>& p1,
        const Point<DIM>& p2) const override;
    //! \see BigShape::type
    NameShape type() const override {
        return NameShape::Tore;
    }

    //! \return the projection of a point inside the Cube
    void projection(Point<DIM>& point) const override;
    //! \see BigShape::bounce
    void bounce(Sphere<DIM>& sph) const override;
};

//! warning : why does BigSphere inherit from Cube ?
template<short unsigned DIM>
class BigSphere : public Cube<DIM> {
private:
    Point<DIM> center;
    double radius;
public:
    // constructor
    //! \param L_: surrounding box.
    //! the sphere is centered of maximal radius in the box
    BigSphere(Point<DIM> L_);

    //! \return the volume of the sphere
    double volume() const override;
    //! \return whether the point is inside the bigShape
    bool isInside(const Point<DIM>& point, const double& minRadius) const override;
    //! a sphere is always inside the bigShape
    bool isInside(const Sphere<DIM>& sph) const override;
    //! \see BigShape::bounce
    void bounce(Sphere<DIM>& sph) const override;
    //! \see BigShape::type
    NameShape type() const override {
        return NameShape::Sphere;
    }
    //! not applicable
    double faceArea(int direction) = delete;
};

template<short unsigned DIM>
class BigCylinder : public Cube<DIM> {
private:
    double height;
    BigSphere<2> circle;
public:
    // constructor
    //! \param L_: surrounding box.
    //! the cylinder has its axis in the z direction, and is maximal in all directions wrt the box
    BigCylinder(Point<DIM> L_);

    //! \return the volume of the cylinder
    double volume() const override;
    //! \return whether the point is inside the bigShape
    bool isInside(const Point<DIM>& point, const double& minRadius) const
        override;
    //! a sphere is always inside the bigShape
    bool isInside(const Sphere<DIM>& sph) const override;
    //! \see BigShape::bounce
    void bounce(Sphere<DIM>& sph) const override;
    //! \see BigShape::type
    NameShape type() const override {
        return NameShape::Cylinder;
    }
};

template<unsigned short DIM>
inline unique_ptr<AmbiantSpace::BigShape<DIM>> createShape(
    AmbiantSpace::NameShape nameShape, Point<DIM> L) {
    unique_ptr<AmbiantSpace::BigShape<DIM>> ptr;
    if (nameShape == AmbiantSpace::NameShape::Tore) {
        ptr.reset(new AmbiantSpace::Tore<DIM>(L));
    } else if (nameShape == AmbiantSpace::NameShape::Cube) {
        ptr.reset(new AmbiantSpace::Cube<DIM>(L));
    } else if (nameShape == AmbiantSpace::NameShape::Sphere) {
        ptr.reset(new AmbiantSpace::BigSphere<DIM>(L));
    } else if (nameShape == AmbiantSpace::NameShape::Cylinder) {
        ptr.reset(new AmbiantSpace::BigCylinder<DIM>(L));
    } else {
        throw invalid_argument(__PRETTY_FUNCTION__);
    }
    return ptr;
}

} // namespace AmbiantSpace
} // namespace sac_de_billes

#include "GlobalShape.ixx"
#endif /* GLOBALSHAPE_HXX_ */

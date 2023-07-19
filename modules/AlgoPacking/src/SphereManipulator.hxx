//! Copyright : see license.txt
//!
//! \brief For manipulating spheres placed in a torus. (Applying isometries on it.)
//!

#ifndef SPHEREMANIPULATOR_HXX_
#define SPHEREMANIPULATOR_HXX_

#include "StdHeaders.hxx"

#include "AlgoPacking.hxx"
#include "AmbiantSpace.hxx"
#include "GlobalShape.hxx"
#include "SphereContainer.hxx"


namespace sac_de_billes {
using namespace std;

template<unsigned DIM>
class SphereManipulator : public SphereContainer<DIM> {
private:
    //! Global shape (only periodicity is relevant)
    unique_ptr<AmbiantSpace::Tore<DIM>> torus;
    //!
    vector<Sphere<DIM>> originalSpheres;

public:
    //! Puts a collection of spheres inside the SphereManipulator
    SphereManipulator(vector<Sphere<DIM>> theSpheres_,
        Point<DIM> length);
    //! Puts a collection of spheres inside the SphereManipulator
    SphereManipulator(vector<SimpleSphereFormat> theSpheres_,
        Point<DIM> length);
    //! Puts a collection of spheres inside the from a proceeded AlgoRSA<DIM>
    template<class T>
    SphereManipulator(const AlgoInterface<DIM, T>& myAlgo);
    //! translates (periodically) the sphere collection
    //! \param vector for translation
    void translate(const Point<DIM>& vector);
    void translate_back(const Point<DIM>& vector);
    //! \return the distance between the spheres and the artificial boundary of the cube
    //! \see SphereManipulator<DIM>::distanceToFaceS, SphereManipulator<DIM>::distanceToEdgeS , SphereManipulator<DIM>::distanceToCornerS
    double distance_cube_spheres() const;
    //! randomly tries to maximize the distance_cube_spheres
    //! \param nb_tries : number of random points picked
    double random_search(size_t nb_tries);
    //! return an upper bound on the minimal distance from the spheres to the cube that may be attained
    //! uses here the optimal solution for the distance to the faces
    double upper_bound_on_best_distmin();
    //! get the dimension of the surrounding box
    Point<DIM> getLength() const override;

private:
    //! stores the vector by which the collection of spheres has been translated
    Point<DIM> vector_of_translation;
    //! \return the distance of the parallel tangential plane to the face
    double distanceToFaceS(const Sphere<DIM>& sphere) const;
    //! \return the distance : radius - distance(center,edge)
    double distanceToEdgeS(const Sphere<DIM>& sphere) const;
    //! \return the distance : raidus - distance(center, corner)
    double distanceToCornerS(const Sphere<DIM>& sphere) const;
    //! \return true if the corner is within a distance 'radius +dist' of the sphere center
    bool cornerWithinSphere(const Sphere<DIM>& sphere, double dist) const;
};

} // namespace sac_de_billes

#include "SphereManipulator.ixx"

#endif /* SPHEREMANIPULATOR_HXX_ */

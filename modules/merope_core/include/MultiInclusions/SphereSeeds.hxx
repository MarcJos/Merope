//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/Interface.hxx"
#include "../MesoStructure/InsideTorus.hxx"
#include "../Geometry/GeomTools.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

//! Ecriture pour Cast3M d'une longue ligne qui risque d'avoir plus que 72 caracteres
void writeLongLine(const std::vector<std::string>&, std::ostream&);
// Read a keyword in the input file
// \param keyword Keyword
// \param fileStreamo
void readKeyword(const std::string keyword, ifstream& f);

template<unsigned short DIM>
class SphereSeeds : public SphereCollection<DIM>, public InsideTorus<DIM> {
    //! class holding, printing and manipulating collections of spheres
public:
    SphereSeeds();

    //! reads from File
    //! VER format
    //! fixme (quasi-obsolete + buggy)
    void fromFile(string f);
    // fixme
    void fromFile(string f1, string f2);
    //! set the spheres from a spherecollection
    void setSpheres(const vector<Sphere<DIM>> sph) {
        this->fill(sph);
    }
    //! set the spheres from sphereSeeds
    void setSpheres(const SphereSeeds<DIM>& sphereSeed);
    //! Spheres generation according to 2 given radius histograms
    //! \param seed : random seed
    //! \param dist : Spheres distribution method
    //! \param di : Repulsion distance (inactive for boolean schemes)
    //! \param desiredRPhi : desired (radius,volume_fraction)
    //! \param tabPhases_ (optional) = corresponding phases
    void fromHisto(unsigned seed, const algoSpheres::TypeAlgo& dist,
        const double di, vector<array<double, 2>> desiredRPhi,
        vector<PhaseType> tabPhases = { });

private:
    //! Verify that the spheres are not touching each others (including a repulsion distance)
    //! \param di Repulsion distance
    bool isClose(const double di);
    //! default_torus
    static constexpr array<double, DIM> DEFAULT_L = create_array<DIM>(1.);
};

}  // namespace merope


#include "../MultiInclusions/SphereSeeds.ixx"



//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "StdHeaders.hxx"

#include "AmbiantSpace.hxx"
#include "GlobalShape.hxx"

namespace sac_de_billes {
using namespace std;

namespace auxi_SphereCollection {
//! tuple for (phase, volume fraction)
template<class PHASE_TYPE>
class PhaseFrac {
public:
    PHASE_TYPE phase;
    double fracVol;  // volume fraction
    PhaseFrac(PHASE_TYPE p = 0, double f = 0) :
        phase{ p }, fracVol{ f } {}
    bool operator<(PhaseFrac<PHASE_TYPE>& phfv2) {
        return this->phase < phfv2.phase;
    }
};

//! tuple for (phase, volume_fraction,radius)
class PhaseFracRadN : public PhaseFrac<PhaseType> {
public:
    PhaseFracRadN(unsigned short p = 0, double f = 0, double r = 0) :
        PhaseFrac(p, f), radius{ r }, number{ 0 } {}
    double radius;  // radius (of spheres)
    size_t number;
};
}  // namespace auxi_SphereCollection

template<unsigned short DIM>
class SphereContainer {
    // Common class for outputs of a sphere collection inside a cube
protected:
    //! Sphere collection
    vector<Sphere<DIM>> theSpheres;

public:
    SphereContainer() :
        theSpheres{}, dictionary{ nullptr } {}
    virtual ~SphereContainer() {
        delete dictionary;
    }

    //! output for Ovito
    void printDump(string nameFile) const;
    //! output csv of spheres (for easy use in other applications, such as gmsh)
    void printCSV(string nameFile) const;
    void printCSV_space(string nameFile) const;
    //! output .pos for Combs
    void printPos(string nameFile) const;
    //! output in the old VER format
    void printVER(string nameFile) const;
    //! \return the center and radii of all placed spheres (first indices for the center, last for the radius)
    vector<SimpleSphereFormat> getPlacedSpheres_nonT() const;
    //! translate all the spheres
    //! \param translationVector : vector to translate
    void translate(const Point<DIM>& translationVector);
    //! \return the phases related to the placedSpheres
    vector<PhaseType> getPhases() const;
    //! \return the popsition and radii of all placed spheres
    const vector<Sphere<DIM>>& getSpheres() const { return theSpheres; }
    //! total volume
    double totalVolume() const;
    //! return the number of spheres
    size_t size() const;
    //! sets phases from n -> N (from possibily non-consecutive phases)
    void resetPhases(PhaseType n);
    //! get the dimension of the surrounding box
    virtual Point<DIM> getLength() const = 0;
    //! set the dictionary for names of phases
    void setNamePhase(map<PhaseType, string> dico);

private:
    //! see printCSV
    void printCSV_with_separator(string nameFile, string separator) const;
    //! print on ost the spheres in 4 double (center, radius) + phase
    //! \param ost : on which we write
    //! \param sphere : the sphere
    //! \param separator : string to separate the data
    void printSphere(ofstream& ost, const Sphere<DIM>& sphere,
        string separator) const;
    //! dictionary for names of phases
    map<PhaseType, string>* dictionary;
};

template<unsigned short DIM>
class SphereCollection : public SphereContainer<DIM> {
public:
    //! default constructor
    SphereCollection() :
        SphereContainer<DIM>(), bigShape{ nullptr }, nbPhases{ 0 } {
    }
    //! acquires a bigShape
    void setShape(AmbiantSpace::BigShape<DIM>* bs) {
        bigShape = bs;
    }
    //! initializes the collection of Spheres
    void fill(vector<Sphere<DIM>> vs);
    //! \see  SphereContainer<DIM>::getLength
    inline Point<DIM> getLength() const override {
        return bigShape->L;
    }
    //! subscripting
    inline const Sphere<DIM>& operator[](const size_t& i) const {
        return this->theSpheres[i];
    }
    //! get the number of phases
    size_t getNbPhases() const;
    //! return the array of (phase, radius, volume_fraction)
    vector<auxi_SphereCollection::PhaseFracRadN> tab_PhaseFracRad() const;
    //! print the volume fractions of spheres in smili-VER format
    void printFracVol(string filename) const;

protected:
    //! refers to the global shape
    AmbiantSpace::BigShape<DIM>* bigShape;

private:
    //! avoids to re-compute the nb of Phases at each time
    size_t nbPhases;
};

namespace auxi_SphereCollection {
//! reads the spheres from an imput file containing lines of x, y, z, radius
//! fixme
template<unsigned short DIM>
SphereCollection<DIM> fromFile(istream& fileStream, PhaseType phase,
    AmbiantSpace::BigShape<DIM>* bigShape_);
}
}  // namespace sac_de_billes

#include "SphereContainer.ixx"


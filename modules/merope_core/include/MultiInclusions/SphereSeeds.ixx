//! Copyright : see license.txt
//!
//! \brief
//

#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM>
inline SphereSeeds<DIM>::SphereSeeds() :
    InsideTorus<DIM>() {
    this->setShape(&this->tore);
}

template<unsigned short DIM>
inline void SphereSeeds<DIM>::setSpheres(const SphereSeeds<DIM>& sphereSeed) {
    this->setLength(sphereSeed.getLength());
    this->theSpheres = sphereSeed.getPlacedSpheres();
}

template<unsigned short DIM>
inline void SphereSeeds<DIM>::fromFile(string fname1, string fname2) {
    SphereSeeds <DIM> mi1{ };
    SphereSeeds <DIM> mi2{ };
    mi1.fromFile(fname1);
    mi2.fromFile(fname2);
    // set the length
    this->setLength(mi1.tore.L);
    // get the sphetes
    vector < Sphere <DIM>> vSph1 = mi1.getSpheres();  //phase = 1
    vector < Sphere <DIM>> vSph2 = mi2.getSpheres();  //phase = 1
    for (auto& sph : vSph2) {
        sph.phase = 2;
    }
    vSph1.insert(vSph1.end(), vSph2.begin(), vSph2.end());
    this->theSpheres = vSph1;
}

template<unsigned short DIM>
inline void SphereSeeds<DIM>::fromFile(string fname) {
    // fixme
    ifstream fileStream(fname);
    if (!fileStream) {
        cerr << "Can't open file" + string(fname) << endl;
        throw(runtime_error(__PRETTY_FUNCTION__));
    }
    string str;
    // Parse dimension
    readKeyword("'Dimension'", fileStream);
    unsigned short d;
    fileStream >> d;
    if (d != DIM) {
        throw runtime_error("Inconsistent dimensions");
    }
    // Parse nb of seeds
    readKeyword("'Points'", fileStream);
    size_t seedsNb;
    fileStream >> seedsNb;
    vector < Sphere <DIM>> listSpheres{ };

    // Parse lengths
    array<double, DIM> length;
    readKeyword("'Lx'", fileStream);
    fileStream >> length[0];
    if (d > 1) {  // d= 2 ou d=3
        readKeyword("'Ly'", fileStream);
        fileStream >> length[1];
    }
    if (d > 2) {  // d=3
        readKeyword("'Lz'", fileStream);
        fileStream >> length[2];
    }
    //
    this->setLength(length);
    auto sphereColl = auxi_SphereCollection::fromFile<DIM>(fileStream, 1, &(this->tore));
    this->theSpheres = sphereColl.getSpheres();
}

template<unsigned short DIM>
inline void SphereSeeds<DIM>::fromHisto(unsigned seed,
    const algoSpheres::TypeAlgo& dist, const double di,
    vector<array<double, 2>> desiredRPhi, vector<PhaseType> tabPhases) {
    if (tabPhases.size() != desiredRPhi.size()) {
        cerr << "The phase table is incorrectly set : hence, by default, all the sphere phases is 1" << endl;
        tabPhases = vector<PhaseType>(desiredRPhi.size(), 1);
    }
    const vector<Sphere<DIM>>& listSpheres =
        algoSpheres::throwSpheres < DIM
        >(dist, AmbiantSpace::NameShape::Tore, this->tore.L, seed, desiredRPhi, tabPhases, di);
    this->theSpheres = listSpheres;
}

template<unsigned short DIM>
inline bool SphereSeeds<DIM>::isClose(const double di) {
    for (const auto& sph1 : this->theSpheres) {
        for (const auto& sph2 : this->theSpheres) {
            if (this->tore.distanceCarrre(sph1.center, sph2.center)
                < auxi_function::puissance < 2
                >(sph1.radius = sph2.radius + di)) {
                return true;
            }
        }
    }
    return false;
}

}  // namespace merope




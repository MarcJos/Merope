//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../GenericMerope/StdHeaders.hxx"


namespace sac_de_billes {

/// SphereContainer<DIM>

template<unsigned short DIM>
void SphereContainer<DIM>::translate(const Point<DIM>& translationVector) {
    for (auto& sphere : theSpheres) {
        sphere.center = translationVector + sphere.center;
    }
}


template<unsigned short DIM>
inline vector<SimpleSphereFormat> SphereContainer<DIM>::getPlacedSpheres_nonT() const {
    return sphereTools::vector_fromSphere2Array(getSpheres());
}

template<unsigned short DIM>
inline vector<PhaseType> SphereContainer<DIM>::getPhases() const {
    vector<PhaseType> phases(getSpheres().size());
    for (size_t i = 0; i < getSpheres().size(); i++) {
        phases[i] = getSpheres()[i].phase;
    }
    return phases;
}

template<unsigned short DIM>
void SphereContainer<DIM>::setNamePhase(map<PhaseType, string> dico) {
    delete dictionary;
    dictionary = new map<PhaseType, string>(dico);
}

template<unsigned short DIM>
void SphereContainer<DIM>::printSphere(ofstream& ost, const Sphere<DIM>& sphere,
    string separator) const {
    for (size_t i = 0; i < 3; i++) {
        if (i < DIM)
            ost << sphere.center[i] << separator;
        else ost << 0 << separator;
    }
    ost << sphere.radius << separator;
    if (dictionary) {
        ost << dictionary->at(sphere.phase);
    } else {
        ost << sphere.phase;
    }
}


template<unsigned short DIM>
void SphereContainer<DIM>::printVTK(string nameFile) const {
    ofstream file(nameFile);
    if (!file.is_open()) return;
    file << "# vtk DataFile Version 2.0" << endl;
    file << "SphereCollection" << endl;
    file << "ASCII" << endl;
    file << "DATASET UNSTRUCTURED_GRID" << endl;
    int numSpheres = getSpheres().size();
    file << "POINTS " << numSpheres << " float" << endl;
    for (const auto& sphere : getSpheres()) {
        file << sphere.center[0] << " " << sphere.center[1] << " ";
        if constexpr (DIM == 3) {
            file << sphere.center[2];
        } else {
            file << 0.;
        }
        file << endl;
    }
    file << "POINT_DATA " << numSpheres << endl;
    file << "SCALARS radius float 1" << endl;
    file << "LOOKUP_TABLE default" << endl;
    for (const auto& sphere : getSpheres()) {
        file << sphere.radius << "\n";
    }
}

template<unsigned short DIM>
void SphereContainer<DIM>::printDump(string nameFile) const {
    Point<DIM> loc_L = getLength();
    ofstream ost{ nameFile };
    ost << "ITEM: TIMESTEP" << endl;
    ost << 0 << endl;
    ost << "ITEM: NUMBER OF ATOMS" << endl;
    ost << getSpheres().size() << endl;
    ost << "ITEM: BOX BOUNDS ff ff ff" << endl;
    ost << 0 << " " << loc_L[0] << endl;
    ost << loc_L[1] << " " << 0 << endl;
    if constexpr (DIM == 3) {
        ost << loc_L[2] << " " << 0 << endl;
    } else {
        ost << loc_L[1] / 10. << " " << 0 << endl;
    }
    ost << "ITEM: ATOMS id x y z Radius Phase" << endl;
    for (unsigned long i = 0; i < getSpheres().size(); i++) {
        ost << (i + 1) << " ";
        printSphere(ost, getSpheres()[i], " ");
        ost << endl;
    }
}

template<unsigned short DIM>
void SphereContainer<DIM>::printCSV_space(string nameFile) const {
    printCSV_with_separator(nameFile, " ");
}

template<unsigned short DIM>
void SphereContainer<DIM>::printCSV(string nameFile) const {
    printCSV_with_separator(nameFile, " , ");
}

template<unsigned short DIM>
void SphereContainer<DIM>::printCSV_with_separator(string nameFile,
    string separator) const {
    ofstream ost{ nameFile };
    Point<DIM> loc_L = getLength();
    ost << loc_L[0] << separator << loc_L[1];
    if (DIM == 3) {
        ost << separator << loc_L[2];
    }
    ost << endl;
    for (const auto& sphere : getSpheres()) {
        printSphere(ost, sphere, separator);
        ost << endl;
    }
}

template<unsigned short DIM>
void SphereContainer<DIM>::printPos(string nameFile) const {
    ofstream ost{ nameFile };
    // 1ere ligne
    ost << "Combs 1.2.0 " << endl;
    // 2eme ligne
    /// Dimensions du tore
    for (size_t i = 0; i < DIM; i++) {
        ost << getLength()[i] << " ";
    }
    if (DIM == 2) {
        ost << 0.1 * getLength()[0] << " ";
    }
    /// Incompris
    ost << "cas1 matrice pbc 0" << endl;
    // Liste des spheres
    for (const auto& sph : getSpheres()) {
        printSphere(ost, sph, " ");
        ost << endl;
    }
}

template<unsigned short DIM>
void SphereContainer<DIM>::printVER(string nameFile) const {
    ofstream os{ nameFile };
    os << "'Dimension' " << DIM << endl;
    os << "'Points' " << theSpheres.size() << endl;
    os << "'Lx' " << getLength()[0] << endl;
    if (DIM > 1) os << "'Ly' " << getLength()[1] << endl;
    if (DIM > 2) os << "'Lz' " << getLength()[2] << endl;
    os << std::setprecision(6);
    for (const auto& sph : theSpheres) {
        auto coord = sphereTools::fromSphere2Array(sph);
        for (size_t i = 0; i < 4; i++) {
            os << coord[i] << " ";
        }
        os << endl;
    }
}

template<unsigned short DIM>
double SphereContainer<DIM>::totalVolume() const {
    return geomTools::volume_all(theSpheres);
}

template<unsigned short DIM>
size_t SphereContainer<DIM>::size() const {
    return theSpheres.size();
}

template<unsigned short DIM>
void SphereContainer<DIM>::resetPhases(PhaseType n) {
    sphereTools::sort(theSpheres);
    PhaseType current_phase = theSpheres[0].phase;
    PhaseType desired_phase = n;
    for (auto& sph : theSpheres) {
        if (sph.phase != current_phase) {
            desired_phase++;
            current_phase = sph.phase;
        }
        sph.phase = desired_phase;
    }
}

/// SphereCollection<DIM>
template<unsigned short DIM>
inline void SphereCollection<DIM>::printFracVol(string filename) const {
    vector < auxi_SphereCollection::PhaseFracRadN > phase_frac_rad =
        tab_PhaseFracRad();
    ofstream fout(filename);
    fout << "# Phase / Radius / Volume Fraction / Number of inclusions\n";
    for (const auto& pfr : phase_frac_rad) {
        fout << pfr.phase << " " << pfr.radius << " " << pfr.fracVol << " "
            << " " << pfr.number << endl;
    }
}

template<unsigned short DIM>
void SphereCollection<DIM>::fill(vector<Sphere<DIM>> vs) {
    this->theSpheres = vs;
}

template<unsigned short DIM>
size_t SphereCollection<DIM>::getNbPhases() const {
    set<PhaseType> allPhases{ };
    for (const auto& s : this->theSpheres) {
        allPhases.insert(s.phase);
    }
    nbPhases = allPhases.size();
    return nbPhases;
}

template<unsigned short DIM>
vector<auxi_SphereCollection::PhaseFracRadN> SphereCollection<DIM>::tab_PhaseFracRad() const {
    vector < Sphere <DIM>> vecSph = this->theSpheres;
    sphereTools::sort(vecSph);
    vector < auxi_SphereCollection::PhaseFracRadN > result{ };
    PhaseType currentPhase = vecSph[0].phase;
    double currentRadius = 0.;
    bool firstTime = true;  // special case of entering the loop
    size_t i = 0;
    for (const auto& sph : vecSph) {
        if (firstTime) {
            firstTime = false;
            i = 0;
            result.push_back(
                auxi_SphereCollection::PhaseFracRadN(sph.phase, 0.,
                    sph.radius));
        } else if (sph.phase != currentPhase
            or abs(sph.radius - currentRadius) > 0.0001 * currentRadius
            or firstTime) {
            i++;
            result.push_back(
                auxi_SphereCollection::PhaseFracRadN(sph.phase, 0.,
                    sph.radius));
        }
        currentPhase = sph.phase;
        currentRadius = sph.radius;
        result[i].fracVol += geomTools::volume(sph) / bigShape->volume();
        result[i].number++;
    }
    return result;
}

// auxi_SphereCollection

template<unsigned short DIM>
SphereCollection<DIM> auxi_SphereCollection::fromFile(istream& fileStream, PhaseType phase, AmbiantSpace::BigShape<DIM>* bigShape_) {
    vector<Sphere<DIM>> vecSph{};
    Sphere<DIM> sphere{};
    while (sphereTools::fromLine<DIM>(fileStream, phase, sphere)) {
        vecSph.push_back(sphere);
    }
    SphereCollection <DIM> sc{ };
    sc.fill(vecSph);
    sc.setShape(bigShape_);
    return sc;
}

}  // namespace sac_de_billes



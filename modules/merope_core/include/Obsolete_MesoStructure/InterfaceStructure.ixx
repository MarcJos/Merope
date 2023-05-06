//! Copyright : see license.txt
//!
//! \brief 
//!

#ifndef INTERFACESTRUCTURE_IXX_
#define INTERFACESTRUCTURE_IXX_

#include "../MeropeNamespace.hxx"


namespace merope {

//! Pre_Interface
template<unsigned short DIM>
Pre_InterfaceMultiInclusions<DIM>::Pre_InterfaceMultiInclusions():
    MicroType(), SphereSeeds<DIM>(), WithAspratio<DIM>() {}

template<unsigned short DIM>
inline vector<PhaseType> Pre_InterfaceMultiInclusions<DIM>::getPhaseList() const {
    vector < PhaseType > phaseListCollSphe(this->getSpheres().size());
    for (size_t i = 0; i < this->getSpheres().size(); i++) {
        phaseListCollSphe[i] = this->getSpheres()[i].phase;
    }
    sort(phaseListCollSphe.begin(), phaseListCollSphe.end());
    auto new_end = unique(phaseListCollSphe.begin(), phaseListCollSphe.end());
    phaseListCollSphe.resize(distance(phaseListCollSphe.begin(), new_end));
    return phaseListCollSphe;
}

template<unsigned short DIM>
inline vector<Identifier> Pre_InterfaceMultiInclusions<DIM>::getAllIdentifiers() const {
    return auxi_MultiInclusions::getAllIdentifiers(this->getSpheres().size());
}

////////////////////////////
//! InterfaceMultiInclusions
////////////////////////////

template<unsigned short DIM>
inline InterfaceMultiInclusions<DIM>::InterfaceMultiInclusions():
    InterfaceMultiInclusions(Pre_InterfaceMultiInclusions<DIM>()) {}

template<unsigned short DIM>
inline InterfaceMultiInclusions<DIM>::InterfaceMultiInclusions(
    Pre_InterfaceMultiInclusions<DIM> pre_InterfaceMultiInclusions) :
    Pre_InterfaceMultiInclusions<DIM>(pre_InterfaceMultiInclusions), layerList{ } {
}

template<unsigned short DIM>
inline MultiInclusions<DIM> InterfaceMultiInclusions<DIM>::build() {
    auto multiInclusions = this->buildSimple();
    multiInclusions.addLayer(layerList);
    return multiInclusions;
}

template<unsigned short DIM>
void InterfaceMultiInclusions<DIM>::setLayerList(
    const vector<smallShape::LayerInstructions>& layerList_) {
    this->layerList = layerList_;
}

template<unsigned short DIM>
const vector<smallShape::LayerInstructions>& InterfaceMultiInclusions<DIM>::getLayerList() const {
    return layerList;
}

template<unsigned short DIM>
inline MultiInclusions<DIM> InterfaceMultiInclusions<DIM>::buildSimple() {
    MultiInclusions <DIM> multiInclusions{ };
    if (this->typeCrystal == TypeCrystal::Spheres) {
        SphereInclusions <DIM> sphInc{ };
        sphInc.setLength(this->getL());
        sphInc.setSpheres(this->getSpheres());
        multiInclusions.setInclusions(sphInc);
    }
    else if (this->typeCrystal == TypeCrystal::Laguerre
        or this->typeCrystal == TypeCrystal::Voronoi) {
        if (this->typeCrystal == TypeCrystal::Voronoi) {
            for (auto& s : this->theSpheres) { // all the spheres shoud have radius 1
                s.radius = 1.;
            }
        }
        LaguerreTess <DIM> laguerreTess(this->getL(), this->getSpheres());
        laguerreTess.setAspRatio(this->aspratio);
        multiInclusions.setInclusions(laguerreTess);
    }
    else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Undefined behavior");
    }
    return multiInclusions;
}

//! InterfacePolyCrystal<DIM>

template<unsigned short DIM>
inline InterfaceStructure<DIM>::InterfaceStructure():
    Colorize(),
    mainInclusions(), secdInclusions(), innerPhase{ DEFAULT_INNERPHASE },
    erosionPhase{ DEFAULT_EROSIONPHASE }, erosionInclusionsPhase{ DEFAULT_EROSIONINCLUSIONSPHASE }, erosionWidth{ 0. } {}

template<unsigned short DIM>
inline void InterfaceStructure<DIM>::setErosionPhase(
    PhaseType erosionPhase_) {
    if (colorization == ColorMaterialID::Erode) {
        cerr << "__PRETTY_FUNCTION__";
        throw runtime_error("Impossible to set the erosionPhase for a material of type Erode. "
            "It is equal to the number of spheres + 1."
            "Consider use Erode2Mat instead");
    }
    erosionPhase = erosionPhase_;
}

template<unsigned short DIM>
inline void InterfaceStructure<DIM>::setErosionInclusionsPhase(
    PhaseType erosioNinclusionsPhase_) {
    erosionInclusionsPhase = erosioNinclusionsPhase_;
}

template<unsigned short DIM>
inline void InterfaceStructure<DIM>::setErosionWidth(
    double erosionWidth_) {
    this->erosionWidth = erosionWidth_;
}

template<unsigned short DIM>
inline Structure<DIM> InterfaceStructure<DIM>::build() {
    preBuild();
    Structure<DIM>structure(buildMainInclusions());
    postBuild(structure);
    return structure;
}

template<unsigned short DIM>
inline void InterfaceStructure<DIM>::preBuild() {
    secdInclusions.setTypeCrystal(TypeCrystal::Spheres);
    switch (this->colorization) {
    case (ColorMaterialID::Poly):
    case (ColorMaterialID::Erode):
    {
        auto theSpheres = mainInclusions.getSpheres();
        for (size_t i = 0; i < theSpheres.size(); i++) {
            theSpheres[i].phase = i;
        }
        mainInclusions.setSpheres(theSpheres);
        break;
    }
    case (ColorMaterialID::Erode2Mat):
    {
        auto coll = mainInclusions.getSpheres();
        for (auto& sph : coll) {
            sph.phase = innerPhase;
        }
        mainInclusions.setSpheres(coll);
        break;
    }
    case (ColorMaterialID::Erode3Mat):
    {
        auto coll = mainInclusions.getSpheres();
        for (auto& sph : coll) {
            sph.phase = DEFAULT_INNERPHASE;
        }
        mainInclusions.setSpheres(coll);
        break;
    }
    default:
        throw runtime_error("Unknown ColorMaterialID");
    }
}

template<unsigned short DIM>
inline void InterfaceStructure<DIM>::postBuild(Structure<DIM>& structure) {
    switch (this->colorization) {
    case (ColorMaterialID::Poly):
    case (ColorMaterialID::Erode):
    case (ColorMaterialID::Erode2Mat):
        break;
    case (ColorMaterialID::Erode3Mat):
        structure = Structure<DIM>(structure, buildSecdInclusions(), this->functionPorositySpheres(), vector<PhaseType> {erosionInclusionsPhase, erosionPhase });
        break;
    default:
        throw runtime_error("Unknown ColorMaterialID");
    }
}


template<unsigned short DIM>
inline void InterfaceStructure<DIM>::setInnerPhase(PhaseType innerPhase_) {
    this->innerPhase = innerPhase_;
}

template<unsigned short DIM>
inline MultiInclusions<DIM> InterfaceStructure<DIM>::buildMainInclusions() {
    InterfaceMultiInclusions<DIM> multiInclusion(mainInclusions);
    switch (this->colorization) {
    case(ColorMaterialID::Poly):
        break;
    case(ColorMaterialID::Erode):
    case(ColorMaterialID::Erode2Mat):
    case(ColorMaterialID::Erode3Mat):
        multiInclusion.setLayerList(buildLayerInstructions(temporaryErosionPhase()));
        break;
    default:
        throw runtime_error("Unknown ColorMaterialID");
    }
    return multiInclusion.build();
}

template<unsigned short DIM>
inline MultiInclusions<DIM> InterfaceStructure<DIM>::buildSecdInclusions() {
    InterfaceMultiInclusions<DIM> multiInclusion(secdInclusions);
    return multiInclusion.build();
}

template<unsigned short DIM>
inline vector<smallShape::LayerInstructions> InterfaceStructure<DIM>::buildLayerInstructions(PhaseType erosionPhase_) const {
    vector<smallShape::LayerInstructions> instructions{};
    for (size_t i = 0; i < mainInclusions.getAllIdentifiers().size(); i++) {
        instructions.push_back(smallShape::LayerInstructions{ i, erosionPhase_, erosionWidth });
    }
    return instructions;
}

template<unsigned short DIM>
inline PhaseType InterfaceStructure<DIM>::temporaryErosionPhase() const {
    switch (this->colorization) {
    case(ColorMaterialID::Poly):
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Should not be used.");
    case(ColorMaterialID::Erode):
        return mainInclusions.getSpheres().size();
        break;
    case(ColorMaterialID::Erode2Mat):
        return erosionPhase;
        break;
    case(ColorMaterialID::Erode3Mat):
        return DEFAULT_EROSIONPHASE;
        break;
    default:
        throw runtime_error("Unknown ColorMaterialID");
    }
}

template<unsigned short DIM>
inline void InterfaceStructure<DIM>::setLength(array<double, DIM> L) {
    mainInclusions.setLength(L);
    secdInclusions.setLength(L);
}

template<unsigned short DIM>
inline std::function<PhaseType(PhaseType, PhaseType)> InterfaceStructure<DIM>::functionPorositySpheres() const {
    PhaseType defaultInnerPh = DEFAULT_INNERPHASE;
    PhaseType defaultErosionPh = DEFAULT_EROSIONPHASE;
    PhaseType erosionPh = erosionPhase;
    PhaseType erosionInclusionPh = erosionInclusionsPhase;
    PhaseType innerPh = innerPhase;
    std::function<PhaseType(PhaseType, PhaseType)> func = [defaultInnerPh, defaultErosionPh, innerPh, erosionPh, erosionInclusionPh](PhaseType phi1, PhaseType phi2) {
        if (phi1 == defaultInnerPh) {
            return innerPh;
        }
        else if (phi1 == defaultErosionPh and phi2 == 0) {
            return erosionPh;
        }
        else if (phi1 == defaultErosionPh and phi2 == 1) {
            return erosionInclusionPh;
        }
        else {
            cerr << "###############" << endl;
            cerr << __PRETTY_FUNCTION__ << endl;
            cerr << "First phase  : " << phi1 << endl;
            cerr << "Secnd phase  : " << phi2 << endl;
            cerr << "###############" << endl;
            throw runtime_error("Unexpected phases");
        }
    };
    return func;
}

template<unsigned short DIM>
void InterfaceStructure<DIM>::frac2erosionWidth(double phi) {
    setErosionWidth(0.5 * calculeErodeeJdGn(phi));
}

template<unsigned short DIM>
double InterfaceStructure<DIM>::calculeErodeeJdGn(double phi) const {
    switch (DIM) {
    case 2:
    {
        double eJdGn = phi * (0.185384 * phi + 0.502622);
        return eJdGn * sqrt(mainInclusions.tore.volume() / mainInclusions.getSpheres().size());
    }
    case 3:
    {
        double eJdGn = phi * (0.166028 * phi + 0.348287);
        return eJdGn * cbrt(mainInclusions.tore.volume() / mainInclusions.getSpheres().size());
    }
    default:
        throw logic_error("Seeds::calculeErodeeJdGn : dimension should be 2 or 3");
    }
    // Voir les fonctions e2D(x) et e3D(x) dans le fichier /home/mi234124/Travail/PolyErode/plot.gp
}

template<unsigned short DIM>
double InterfaceStructure<DIM>::erosionWidth2frac() {
    switch (DIM) {
    case 2:
    {
        double eJdGn = 2. * erosionWidth / sqrt(mainInclusions.tore.volume() / mainInclusions.getSpheres().size());
        return eJdGn * (1.96789 - 1.05251 * eJdGn);
    }
    case 3:
    {
        double eJdGn = 2. * erosionWidth / cbrt(mainInclusions.tore.volume() / mainInclusions.getSpheres().size());
        return eJdGn * (2.82648 - 2.62721 * eJdGn);
    }
    default:
        throw logic_error("VoronoiErode3Mat::eJdG2frac : only works in 2D and 3D");
    }
    // Voir les fonctions phi2D0(x) et phi3D0(x) dans le fichier /home/mi234124/Travail/PolyErode/plot.gp
}

} // namespace merope


#endif /* INTERFACESTRUCTURE_IXX_ */

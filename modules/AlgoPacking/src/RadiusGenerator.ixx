//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef RADIUSGENERATOR_IXX_
#define RADIUSGENERATOR_IXX_

namespace sac_de_billes {
using namespace std;

namespace algoRSA_aux {
// Sets default phases if possible
template<class T>
inline void checkPhase(const vector<T>& tabRadii_, vector<PhaseType>& phases) {
    if (phases.size() != tabRadii_.size()) {
        if (phases.size() == 0) {
            phases = vector<PhaseType>(tabRadii_.size(),
                algoRSA_aux::DEFAULT_PHASE);
        }
        else {
            throw invalid_argument(
                "Incoherency between list of phases and the list of radii (different size)");
        }
    }
}

inline vector<algoRSA_aux::RadiusPhase> createTabRadiiPhase(
    const vector<double>& tabRadii_, vector<PhaseType> phases_) {
    checkPhase(tabRadii_, phases_);
    vector < algoRSA_aux::RadiusPhase > tabRadiiPhase = { };
    tabRadiiPhase.reserve(tabRadii_.size());
    for (size_t i = 0; i < tabRadii_.size(); i++) {
        tabRadiiPhase.push_back(
            algoRSA_aux::RadiusPhase(tabRadii_[i], phases_[i]));
    }
    return tabRadiiPhase;
}

inline vector<algoRSA_aux::RPhiPhase> createTabRPhiPhase(
    const vector<array<double, 2>>& desiredRPhi, vector<PhaseType> phases_) {
    checkPhase(desiredRPhi, phases_);
    vector < algoRSA_aux::RPhiPhase > tabRPhiPhase = { };
    tabRPhiPhase.reserve(desiredRPhi.size());
    for (size_t i = 0; i < desiredRPhi.size(); i++) {
        tabRPhiPhase.push_back(
            algoRSA_aux::RPhiPhase(desiredRPhi[i][0], desiredRPhi[i][1],
                phases_[i]));
    }
    return tabRPhiPhase;
}

template<unsigned short DIM>
inline vector<array<double, 2>> reComputeRPhi(
    vector<array<double, 2>> desiredRPhi, double exclusionDistance) {
    vector < array<double, 2> > new_desiredRPhi { };
    double new_R = 0;
    double new_Phi = 0;
    for (auto phiR : desiredRPhi) {
        new_R = phiR[0] + 0.5 * exclusionDistance;
        new_Phi = phiR[1] * pow(new_R, DIM) / pow(phiR[0], DIM);
        new_desiredRPhi.push_back(array<double, 2> { new_R, new_Phi });
    }
    return new_desiredRPhi;
}

inline vector<algoRSA_aux::RadiusPhase> readRadii(string nameFile) {
    vector<double> tabRadii_ { };
    vector<PhaseType> tabPhases_ { };
    ifstream input(nameFile);
    string line;
    while (getline(input, line)) {
        istringstream iss(line);
        double radius;
        PhaseType phase;
        iss >> radius;
        tabRadii_.push_back(radius);
        if (iss >> phase) { // the 2nd identifier phase is not mandatory
            tabPhases_.push_back(phase);
        }
    }
    cout << "Number of radii  = " << tabRadii_.size() << endl;
    return createTabRadiiPhase(tabRadii_, tabPhases_);
}
}

template<unsigned short DIM>
inline bool algoRSA_aux::RadiusGenerator<DIM>::nextRadius() {
    indexRadius++;
    return indexRadius != tabRadii.size();
}

template<unsigned short DIM>
inline double algoRSA_aux::RadiusGenerator<DIM>::getRadius() const {
    return tabRadii.at(indexRadius).radius;
}

template<unsigned short DIM>
inline PhaseType algoRSA_aux::RadiusGenerator<DIM>::getPhase() const {
    return tabRadii.at(indexRadius).phase;
}

template<unsigned short DIM>
inline void algoRSA_aux::RadiusGenerator<DIM>::emptyTabRadii() {
    if (tabRadii.size() == 0) {
        throw runtime_error(
            "algoRSA_aux::RadiusGenerator<DIM>::RadiusGenerator : No spheres should be placed!");
    }
    tabRadii.shrink_to_fit();
}

template<unsigned short DIM>
inline algoRSA_aux::RadiusGenerator<DIM>::RadiusGenerator(
    AmbiantSpace::BigShape<DIM>* bigShape_,
    vector<array<double, 2>> desiredRPhi_, vector<PhaseType> phases,
    double exclusionDistance):
    RadiusGenerator(bigShape_) {
    auto normalized_desiredRPhi = reComputeRPhi < DIM
    >(desiredRPhi_, exclusionDistance);
    setDesiredRPhi(createTabRPhiPhase(normalized_desiredRPhi, phases));
    emptyTabRadii();
}

template<unsigned short DIM>
inline algoRSA_aux::RadiusGenerator<DIM>::RadiusGenerator(
    AmbiantSpace::BigShape<DIM>* bigShape_, vector<double> tabRadii_,
    vector<PhaseType> phases_, double exclusionDistance):
    RadiusGenerator(bigShape_, createTabRadiiPhase(tabRadii_, phases_),
        exclusionDistance) {}

template<unsigned short DIM>
inline algoRSA_aux::RadiusGenerator<DIM>::RadiusGenerator(
    AmbiantSpace::BigShape<DIM>* bigShape_, string nameFile,
    double exclusionDistance) :
    RadiusGenerator(bigShape_, readRadii(nameFile), exclusionDistance) {}

template<unsigned short DIM>
inline algoRSA_aux::RadiusGenerator<DIM>::RadiusGenerator(
    AmbiantSpace::BigShape<DIM>* bigShape_,
    vector<RadiusPhase> tabRadPhase_, double exclusionDistance) :
    RadiusGenerator(bigShape_) {
    setTabRadii(tabRadPhase_, exclusionDistance);
    emptyTabRadii();
}

template<unsigned short DIM>
inline void algoRSA_aux::RadiusGenerator<DIM>::shrink(double factor) {
    for (auto& rphi : tabRadii) {
        rphi.radius *= factor;
    }
}

template<unsigned short DIM>
inline void algoRSA_aux::RadiusGenerator<DIM>::setTabRadii(
    vector<algoRSA_aux::RadiusPhase> tabRadii_, double exclusionDistance) {
    tabRadii = tabRadii_;
    for (auto& radiusPhase : tabRadii) {
        radiusPhase.radius += 0.5 * exclusionDistance;
    }
    // sorts the vector tabRadii in decreasing order w.r.t. the radii
    sort(tabRadii.begin(), tabRadii.end(), [](auto tr1, auto tr2) {
        return tr1.radius > tr2.radius;
        });
}

template<unsigned short DIM>
inline void algoRSA_aux::RadiusGenerator<DIM>::setDesiredRPhi(
    vector<RPhiPhase> desiredRPhi) {
    sort(desiredRPhi.begin(), desiredRPhi.end(),
        [](const RPhiPhase& rp1, const RPhiPhase& rp2) {
            return rp1.radius > rp2.radius;
        });
    // find the final size of the vector tabRadii
    for (size_t i = 0; i < desiredRPhi.size(); i++) {
        long numberOfSpheres = (desiredRPhi[i].volFrac * bigShape->volume()
            / sphereTools::volumeSphere <DIM>(desiredRPhi[i].radius));
        tabRadii.reserve(tabRadii.size() + numberOfSpheres);
        for (long j = 0; j != numberOfSpheres; j++) {
            tabRadii.push_back(
                RadiusPhase(desiredRPhi[i].radius, desiredRPhi[i].phase));
        }
    }
}

template<unsigned short DIM>
inline double algoRSA_aux::RadiusGenerator<DIM>::maxRadius() const {
    return (tabRadii[0]).radius;
}

template<unsigned short DIM>
inline double algoRSA_aux::RadiusGenerator<DIM>::minRadius() const {
    return (tabRadii[tabRadii.size() - 1]).radius;
}

template<unsigned short DIM>
vector<array<double, 2>> algoRSA_aux::RadiusGenerator<DIM>::achievedRPhi() const {
    vector<array<double, 2>> res = { };
    if (indexRadius == 0) {
        return res;
    }
    double currentRadius = tabRadii[0];
    double currentPhi = 0;
    for (size_t i = 0; i < indexRadius; i++) {
        if (abs((tabRadii[i]).radius - currentRadius)
            < 0.00001 * currentRadius) { // consider it is the same radius
            currentPhi += sphereTools::volumeSphere < DIM
            >(currentRadius) / bigShape->volumeTore();
        }
        else {
            res.push_back(array<double, 2> { currentRadius, currentPhi });
            currentRadius = tabRadii[i];
            currentPhi = sphereTools::volumeSphere < DIM
            >(currentRadius) / bigShape->volumeTore();
        }
    }
    res.push_back(array<double, 2> { currentRadius, currentPhi });
    return res;
}

template<unsigned short DIM>
vector<array<double, 2>> createRPhi_Bool(vector<array<double, 2>> desiredRPhi) {
    // verify volume fraction not too large
    double totVolFrac = 0.;
    for (const auto& rphi : desiredRPhi) {
        totVolFrac += rphi[1];
    }
    if (totVolFrac >= 1.) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Impossible to have a final volume fraction larger than 1 with BOOL algorithm!");
    }
    //
    totVolFrac = 0.;
    for (auto& rphi : desiredRPhi) {
        //double radius = rphi[0];
        double phi = rphi[1];
        rphi[1] = log((1 - totVolFrac) / (1 - totVolFrac - phi));
        totVolFrac += phi;
    }
    return desiredRPhi;
}

} // namespace sac_de_billes

#endif /* RADIUSGENERATOR_IXX_ */

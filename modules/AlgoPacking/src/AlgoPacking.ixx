//! Copyright : see license.txt
//!
//! \brief
//!
#pragma once


namespace sac_de_billes {
using namespace std;

/// AlgoPacking<DIM>

template<unsigned short DIM>
inline Point<DIM> AlgoPacking<DIM>::getLength() const {
    if (not isSet_bigShape) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("No shape set!");
    }
    return bigShape.get()->L;
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::setBigShape(vector<double> L_,
    AmbiantSpace::NameShape nameShape_) {
    Point<DIM> L{ };
    if (L_.size() != DIM) {
        throw invalid_argument(__PRETTY_FUNCTION__);
    }
    for (size_t i = 0; i < L_.size(); i++) {
        L[i] = L_[i];
    }
    setBigShape(L, nameShape_);
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::setBigShape(Point<DIM> L_,
    AmbiantSpace::NameShape nameShape_) {
    verifyBigShape();
    bigShape = AmbiantSpace::createShape <DIM>(nameShape_, L_);
    bigShape.get()->setBoundaryExclusionDistance(boundaryExclusionDistance);
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::setBigShape(vector<double> L_,
    string nameShape_) {
    setBigShape(L_, AmbiantSpace::readShape(nameShape_));
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::setExclusionDistance(double exclusionDistance_) {
    verifyExclusionDistance();
    exclusionDistance = exclusionDistance_;
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::setBoundaryExclusionDistance(
    double boundaryExclusionDistance_) {
    if (isSet_bigShape) {
        throw runtime_error(
            "You cannot set the boundaryExclusionDistance after the bigShape");
    }
    isSet_exclusionDistance = true;
    boundaryExclusionDistance = boundaryExclusionDistance_
        - 0.5 * exclusionDistance;
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::verifyExclusionDistance() {
    if (isSet_radiusGen) {
        throw runtime_error(
            "You cannot set the exclusionDistance after the radiusGenerator");
    }
    if (not isSet_exclusionDistance) {
        isSet_exclusionDistance = true;
    } else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error(
            "The exclusionDistance should be set only once AND before the boundaryExclusionDistance!");
    }
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::verifyRadiusGenerator() {
    if (not isSet_radiusGen and isSet_bigShape) {
        isSet_radiusGen = true;
    } else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error(
            "The radiusGenerator should be set only once! It is not possible to set the radius generator before the geometry!");
    }
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::verifyBigShape() {
    if (not isSet_bigShape) {
        isSet_bigShape = true;
    } else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("The bigShape should be set only once!");
    }
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::setRadiusGenerator(
    vector<array<double, 2> > desiredRPhi_) {
    setRadiusGenerator(desiredRPhi_, vector<PhaseType> { });
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::setRadiusGenerator(
    vector<array<double, 2> > desiredRPhi_, vector<PhaseType> tabPhases_) {
    this->verifyRadiusGenerator();
    if (this->getTypeAlgo() == algoSpheres::TypeAlgo::BOOL) {
        desiredRPhi_ = createRPhi_Bool<DIM>(desiredRPhi_);
    }
    radiusGen.reset(
        new algoRSA_aux::RadiusGenerator<DIM>(this->bigShape.get(),
            desiredRPhi_, tabPhases_, this->exclusionDistance));
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::setRadiusGenerator(vector<double> tabRadii_) {
    setRadiusGenerator(tabRadii_, vector<PhaseType> { });
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::setRadiusGenerator(vector<double> tabRadii_,
    vector<PhaseType> tabPhases_) {
    this->verifyRadiusGenerator();
    radiusGen.reset(
        new algoRSA_aux::RadiusGenerator<DIM>(this->bigShape.get(),
            tabRadii_, tabPhases_, this->exclusionDistance));
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::setRadiusGenerator(string nameFile) {
    this->verifyRadiusGenerator();
    radiusGen.reset(
        new algoRSA_aux::RadiusGenerator<DIM>(this->bigShape.get(),
            nameFile, this->exclusionDistance));
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::setRadiusGenerator(
    const algoRSA_aux::RadiusGenerator<DIM>& radiusGenerator_) {
    this->verifyRadiusGenerator();
    radiusGen.reset(new algoRSA_aux::RadiusGenerator<DIM>(radiusGenerator_));
}

template<unsigned short DIM>
bool AlgoPacking<DIM>::verifySphere() const {
    auto placedSpheres = this->getSpheres();
    for (unsigned i = 0; i < placedSpheres.size(); i++) {
        for (unsigned j = 0; j < placedSpheres.size(); j++) {
            if (not (i == j)) {
                if (this->bigShape->areSphereIntersected(placedSpheres[i],
                    placedSpheres[j])) {
                    cout << "####################################" << endl;
                    //printConsole();
                    cout << i << " " << j << endl;
                    Sphere <DIM> sph1 = placedSpheres[i];
                    Sphere <DIM> sph2 = placedSpheres[j];
                    for (unsigned short k = 0; k < DIM; k++) {
                        cout << sph1.center[k] << " ";
                    }
                    cout << endl;
                    for (unsigned short k = 0; k < DIM; k++) {
                        cout << sph2.center[k] << " ";
                    }
                    cout << endl;
                    cout << "####################################" << endl;
                    cout << "Distance = "
                        << sqrt(
                            bigShape.get()->distanceCarre(sph1.center,
                                sph2.center)) << endl;
                    cout << "Radii = " << sph1.radius << " , " << sph2.radius
                        << endl;
                    cout << "####################################" << endl;
                    return false;
                }
            }
        }
    }
    return true;
}

template<unsigned short DIM>
double AlgoPacking<DIM>::volumeFraction() const {
    return this->totalVolume() / bigShape->volume();
}


template<unsigned short DIM>
inline map<string, string> AlgoPacking<DIM>::proceed(unsigned seed, unsigned method) {
    //
    map<string, string> outputMessage = {};
    if (not this->isSet_radiusGen or not this->isSet_bigShape) {
        cerr << __PRETTY_FUNCTION__ << " 1 " << endl;
        throw runtime_error("The radiusGenerator or the bigShape is not set!");
    } else if (this->has_proceeded) {
        cerr << __PRETTY_FUNCTION__ << " 2 " << endl;
        throw runtime_error("I can proceed only once!");
    } else {
        outputMessage = proceed_loc(seed, method);
        has_proceeded = true;
    }

    setSpheres();
    printVolumeFraction();
    return outputMessage;
}

template<unsigned short DIM>
inline void AlgoPacking<DIM>::printVolumeFraction() {
    cout << "Volume fraction : " << volumeFraction() << endl;
    cout << "----------------------------" << endl;
}

/// AlgoInterface<DIM,T>

template<unsigned short DIM, class T>
inline void AlgoInterface<DIM, T>::setSpheres() {
    auto corePlacedSpheres = *(algo->placedSpheres);  // beware, they take into account the exclusionDistance
    double exclusionRadius = 0.5 * this->exclusionDistance;
#pragma omp parallel for firstprivate(exclusionRadius)
    for (size_t i = 0; i < corePlacedSpheres.size(); i++) {
        corePlacedSpheres[i].radius -= exclusionRadius;
    }
    this->theSpheres = corePlacedSpheres;
}

template<unsigned short DIM, class T>
inline map<string, string> AlgoInterface<DIM, T>::proceed_loc(unsigned seed, unsigned method) {
    algo = new T(this->bigShape.get(), this->radiusGen.get(), seed);
    return algo->proceed(method);
}

}  // namespace sac_de_billes



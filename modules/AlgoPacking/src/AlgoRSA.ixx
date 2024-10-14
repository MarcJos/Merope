//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace sac_de_billes {
using namespace std;

template <unsigned short DIM>
template <unsigned short NUM_METHOD>
void algoRSA_aux::MotherGrid<DIM>::addSphere(const Sphere<DIM>& sphere,
    MotherVoxel<DIM>* motherVoxel) {
    placedSpheres.push_back(sphere);
    size_t iSphe = placedSpheres.size() - 1;
    motherVoxel->updateCovered(sphere, minRadius);
    if constexpr (NUM_METHOD == 1) {
        linkSphereToVoxel(motherVoxel, sphere, iSphe);
    } else if constexpr (NUM_METHOD == 2) {
        motherVoxel->spheresInside.push_back(iSphe);
    } else {
        throw invalid_argument(__PRETTY_FUNCTION__);
    }
}

template <unsigned short DIM>
void algoRSA_aux::MotherGrid<DIM>::linkSphereToVoxel(
    MotherVoxel<DIM>* motherVoxel, const Sphere<DIM>& sphere, size_t iSphe) {
    array<int, DIM> iMin{};
    array<int, DIM> iMax{};
    Point<DIM> posMother = motherVoxel->position();
    bool sphereOnBoundary = false;
    // find the candidate surrounding squares that may intersect the sphere
    for (unsigned short j = 0; j < DIM; j++) {
        iMin[j] = floor(
            (sphere.center[j] - sphere.radius - maxRadius - posMother[j]) * inverse_voxelLength);
        iMax[j] = floor(
            (sphere.center[j] + sphere.radius + maxRadius - posMother[j]) * inverse_voxelLength) +
            1;
    }
    if (motherVoxel->isClose2Boundary) {
        sphereOnBoundary = correctForBoundaries(motherVoxel, iMin, iMax);
    }
    // copies the sphere inside the square
    auto mapped_i_voxel = surroundingVoxels(motherVoxel->discreteCoordinates,
        iMin, iMax);
    for (const auto& [i, voxel] : mapped_i_voxel) {
        auxi_linkSphereToVoxel(voxel, iSphe, i, sphereOnBoundary);
    }
}

//\fixme : the fact that we are or not near the boundary could be used as a template parameter.
template <unsigned short DIM>
inline vector<tuple<array<int, DIM>, algoRSA_aux::MotherVoxel<DIM>*>> algoRSA_aux::MotherGrid<
    DIM>::surroundingVoxels(DiscPoint<DIM> centralVoxel,
        array<int, DIM> iMin, array<int, DIM> iMax) {
    vector<tuple<array<int, DIM>, algoRSA_aux::MotherVoxel<DIM>*>> result;
    array<int, DIM> i({ 0 });

    auto getIndex = [&centralVoxel, &i, this](size_t j) {
        return auxi_function::fast_modulo(centralVoxel[j] + i[j], this->sizes[j]);
        };  // to get the j-th coordinate of the voxel placed at position i wrt central voxel
    DiscPoint<DIM> index({ 0 });

    if constexpr (DIM == 3) {
        size_t size_0 = iMax[0] - iMin[0];
        size_t size_1 = iMax[1] - iMin[1];
        size_t size_2 = iMax[2] - iMin[2];
        result.resize(size_0 * size_1 * size_2);
        size_t idx = 0;
        for (i[0] = iMin[0]; i[0] < iMax[0]; i[0]++) {
            index[0] = getIndex(0);
            for (i[1] = iMin[1]; i[1] < iMax[1]; i[1]++) {
                index[1] = getIndex(1);
                for (i[2] = iMin[2]; i[2] < iMax[2]; i[2]++, idx++) {
                    index[2] = getIndex(2);
                    result[idx] = make_tuple(i, getVoxel(index));
                }
            }
        }
    } else if constexpr (DIM == 2) {
        size_t size_0 = iMax[0] - iMin[0];
        size_t size_1 = iMax[1] - iMin[1];
        result.resize(size_1 * size_0);
        size_t idx = 0;
        for (i[0] = iMin[0]; i[0] < iMax[0]; i[0]++) {
            index[0] = getIndex(0);
            for (i[1] = iMin[1]; i[1] < iMax[1]; i[1]++, idx++) {
                index[1] = getIndex(1);
                result[idx] = make_tuple(i, getVoxel(index));
            }
        }
    } else {
        throw invalid_argument(__PRETTY_FUNCTION__);
    }

    return result;
}

template <unsigned short DIM>
inline vector<tuple<array<int, DIM>, algoRSA_aux::MotherVoxel<DIM>*>> algoRSA_aux::MotherGrid<
    DIM>::surroundingVoxels(algoRSA_aux::MotherVoxel<DIM>* centralVoxel,
        int layerWidth) {
    array<int, DIM> iMin{};
    array<int, DIM> iMax{};
    for (size_t j = 0; j < DIM; j++) {
        iMin[j] = -layerWidth;
        iMax[j] = layerWidth + 1;
    }
    correctForBoundaries(centralVoxel, iMin, iMax);
    return surroundingVoxels(centralVoxel->discreteCoordinates, iMin, iMax);
}

template <unsigned short DIM>
DiscPoint<DIM> algoRSA_aux::MotherGrid<DIM>::discCoordinate_withoutModulo(
    const Point<DIM>& point) const {
    DiscPoint<DIM> res{};
    for (unsigned short i = 0; i < DIM; i++) {
        res[i] = floor(point[i] * inverse_voxelLength);
    }
    return res;
}

template <unsigned short DIM>
inline array<size_t, DIM> algoRSA_aux::MotherGrid<DIM>::discCoordinate(
    const Point<DIM>& doubleCoord) const {
    return discCoordinate(discCoordinate_withoutModulo(doubleCoord));
}

template <unsigned short DIM>
inline array<size_t, DIM> algoRSA_aux::MotherGrid<DIM>::discCoordinate(
    const DiscPoint<DIM>& discCoord) const {
    array<size_t, DIM> result{};
    for (size_t i = 0; i < DIM; i++) {
        result[i] = auxi_function::fast_modulo(discCoord[i], sizes[i]);
    }
    return result;
}

template <unsigned short DIM>
inline size_t algoRSA_aux::MotherGrid<DIM>::getNbVoxels() {
    size_t nbVox = 1;
    for (size_t i = 0; i < DIM; i++) {
        nbVox *= sizes[i];
    }
    return nbVox;
}

template <unsigned short DIM>
inline bool algoRSA_aux::MotherGrid<DIM>::correctForBoundaries(
    algoRSA_aux::MotherVoxel<DIM>* motherVoxel, array<int, DIM>& iMin,
    array<int, DIM>& iMax) {
    bool nearBoundary = false;
    for (size_t j = 0; j < DIM; j++) {
        const auto discCoord_j = motherVoxel->discreteCoordinates[j];
        if (discCoord_j + iMax[j] >= sizes[j] - 2) {
            iMax[j]++;
            nearBoundary = true;
        }
        if (discCoord_j + iMin[j] <= 1) {
            iMin[j]--;
            nearBoundary = true;
        }
    }
    return nearBoundary;
}

template <unsigned short DIM>
void algoRSA_aux::MotherGrid<DIM>::auxi_linkSphereToVoxel(
    algoRSA_aux::MotherVoxel<DIM>* currentVoxel, const size_t& iSphe,
    const array<int, DIM>& i, const bool& sphereOnBoundary) {
    if (not currentVoxel->covered) {  /// unecessary. Avoid to add spheres where they will not be tested anymore.
        currentVoxel->spheresInside.push_back(iSphe);
        array<int, DIM> relativePos{};
        // By default, we cannot infer easily the relative position of the voxel wrt the sphere in that case
        // thus, we take the worst case, considering the sphere is inside the voxel
        if (not sphereOnBoundary) {  //< in that case, we may compute more precisely the relative position
            for (size_t j = 0; j < DIM; j++) {
                relativePos[j] = -i[j];
            }
        }
        currentVoxel->spheresRelativePosition.push_back(
            Path::get().fromCorner2Index<DIM>(relativePos));
    }
}

template <unsigned short DIM>
vector<Sphere<DIM>> algoRSA_aux::MotherGrid<DIM>::neighborSpheres(
    const MotherVoxel<DIM>* voxel) const {
    size_t nbSph = voxel->spheresInside.size();
    vector<Sphere<DIM>> lisSpheres(nbSph);
    for (size_t i = 0; i < nbSph; i++) {
        lisSpheres[i] = voxel->getSphere(i);
    }
    return lisSpheres;
}

template <unsigned short DIM>
template <unsigned short NUM_METHOD>
bool algoRSA_aux::MotherGrid<DIM>::isPlaceable(const Sphere<DIM>& sph1,
    const MotherVoxel<DIM>* mother) {
    if constexpr (NUM_METHOD == 1) {
        return isPlaceable_1(sph1, mother);
    } else if constexpr (NUM_METHOD == 2) {
        return isPlaceable_2(sph1, mother);
    } else {
        throw invalid_argument(
            "algoRSA_aux::MotherGrid<DIM>::isPlaceable Incorrect NUM_METHOD");
    }
}

template <unsigned short DIM>
bool algoRSA_aux::MotherGrid<DIM>::isPlaceable_1(const Sphere<DIM>& sph1,
    const MotherVoxel<DIM>* mother) const {
    size_t nbSph = mother->spheresInside.size();
    for (size_t i = 0; i < nbSph; i++) {
        if (bigShape->areSphereIntersected(sph1, mother->getSphere(i))) {
            return false;
        }
    }

    return true;
}

template <unsigned short DIM>
bool algoRSA_aux::MotherGrid<DIM>::isPlaceable_2(const Sphere<DIM>& sphere,
    const MotherVoxel<DIM>* mother) {
    if (mother->isClose2Boundary) {
        return checksOtherSpheres_2<7>(mother->discreteCoordinates, sphere, Path::get().pathForNeighbors<DIM, 7>());
    } else {
        return checksOtherSpheres_2<5>(mother->discreteCoordinates, sphere, Path::get().pathForNeighbors<DIM, 5>());
    }
}

template <unsigned short DIM>
template <size_t taille>
bool algoRSA_aux::MotherGrid<DIM>::checksOtherSpheres_2(
    const DiscPoint<DIM>& discCoord, const Sphere<DIM>& sph1,
    const array<array<int, DIM>, nbSubcubes<DIM>(taille)>& tabPath) {
    DiscPoint<DIM> index{};
    Sphere<DIM> sph2;
    for (const auto& i : tabPath) {
        for (size_t j = 0; j < DIM; j++) {
            index[j] = auxi_function::fast_modulo(discCoord[j] + i[j], sizes[j]);
        }
        for (size_t iSphe : (getVoxel(index))->spheresInside) {
            sph2 = placedSpheres[iSphe];
            if (bigShape->areSphereIntersected(sph1, sph2)) {
                return false;
            }
        }
    }
    return true;
}

template <unsigned short DIM>
template <typename T>
algoRSA_aux::MotherVoxel<DIM>* algoRSA_aux::MotherGrid<DIM>::getVoxel(
    array<T, DIM> index) {
    return &(tabVoxels[index]);
}

template <unsigned short DIM>
algoRSA_aux::MotherGrid<DIM>::MotherGrid(std::array<size_t, DIM> sizes_,
    AmbiantSpace::BigShape<DIM>* bigShape_, double voxelLength_,
    double minRadius_, double maxRadius_)
    : sizes(sizes_), bigShape(bigShape_), tabVoxels(), voxelLength(voxelLength_), inverse_voxelLength(1. / voxelLength_), minRadius(minRadius_), maxRadius(maxRadius_) {
    // Inialize the voxel grid
    merope::vox::ArrayDimensions<DIM> arrayDim(sizes_);
    tabVoxels = TYPE_TABVOXEL(arrayDim, MotherVoxel<DIM>());
    // Define each voxel
    DiscPoint<DIM> i{};
    if constexpr (DIM == 3) {
        for (i[0] = 0; i[0] < sizes[0]; i[0]++) {
            for (i[1] = 0; i[1] < sizes[1]; i[1]++) {
                for (i[2] = 0; i[2] < sizes[2]; i[2]++) {
                    initializeVoxel(i);
                }
            }
        }
    } else if constexpr (DIM == 2) {
        for (i[0] = 0; i[0] < sizes[0]; i[0]++) {
            for (i[1] = 0; i[1] < sizes[1]; i[1]++) {
                initializeVoxel(i);
            }
        }
    } else {
        throw invalid_argument(
            "algoRSA_aux::MotherGrid<DIM>::MotherGrid : incorrect DIM");
    }
    placedSpheres.reserve(
        ceil(
            0.5 * bigShape->volume() / sphereTools::volumeSphere<DIM>(minRadius)));  //< Initialize with the maximal possible size.
}

template <unsigned short DIM>
void algoRSA_aux::MotherGrid<DIM>::initializeVoxel(const DiscPoint<DIM>& i) {
    auto currentVoxel = getVoxel(i);
    (*currentVoxel) = MotherVoxel<DIM>(i, voxelLength, this);
    currentVoxel->mother = currentVoxel;
    currentVoxel->spheresInside.reserve(NB_EXPECTED_SPHERES_IN_CUBE);
    currentVoxel->isClose2Boundary = currentVoxel->isItClose2Boundary(bigShape);
    if (currentVoxel->isClose2Boundary) {
        currentVoxel->covered = currentVoxel->isOutside(bigShape, minRadius);
    }
}

template <unsigned short DIM>
algoRSA_aux::CurrentGrid<DIM>::CurrentGrid(MotherGrid<DIM>* motherGrid_,
    mt19937* randGen_) : motherGrid(motherGrid_), randGen(randGen_), nbUncoveredVoxels(0) {
    size_t compteur = 0;
    if constexpr (DIM == 3) {
        size_t nbVoxels = motherGrid->sizes[0] * motherGrid->sizes[1] * motherGrid->sizes[2];
        // defines the grid of voxels
        vector<Voxel<DIM>> futureTabVoxels(nbVoxels,
            motherGrid->tabVoxels[0]);
        loop<true>(motherGrid->sizes[0], motherGrid->sizes[1], motherGrid->sizes[2],
            [this, &compteur, &futureTabVoxels](
                const array<size_t, 3u>& i) {
                    futureTabVoxels[compteur] = *(motherGrid->getVoxel(i));
                    compteur++;
            });
        setTabVoxels(std::move(futureTabVoxels));
    } else if constexpr (DIM == 2) {
        size_t nbVoxels = motherGrid->sizes[0] * motherGrid->sizes[1];
        // defines the grid of voxels
        vector<Voxel<DIM>> futureTabVoxels(nbVoxels, motherGrid->tabVoxels[0]);
        loop<true>(motherGrid->sizes[0], motherGrid->sizes[1],
            [this, &compteur, &futureTabVoxels](
                const array<size_t, 2u>& i) {
                    futureTabVoxels[compteur] = *(motherGrid->getVoxel(i));
                    compteur++;
            });
        setTabVoxels(std::move(futureTabVoxels));
    } else {
        throw invalid_argument("CurrentGrid<DIM>::CurrentGrid. Uncorrect DIM.");
    }
}

template <unsigned short DIM>
template <unsigned short NUM_METHOD>
void algoRSA_aux::CurrentGrid<DIM>::subdivide() {
    // estimates the size of the new tabVoxel
    vector<Voxel<DIM>> futureTabVoxel{};
    futureTabVoxel.reserve(ceil(pow(2, DIM - 1) * nbUncoveredVoxels));  // a reasonable new size
    // objects for the loop
    for (Voxel<DIM>& voxel : tabVoxels) {
        if (not voxel.covered and not voxel.mother->covered) {
            voxel.Voxel<DIM>::template subdivide<NUM_METHOD>(futureTabVoxel);
        }
    }
    setTabVoxels(std::move(futureTabVoxel));
}

template <unsigned short DIM>
algoRSA_aux::Voxel<DIM>* algoRSA_aux::CurrentGrid<DIM>::pickUncoveredVoxel() {
    bool oneMoreTime = true;
    size_t indexVoxel = 0;
    while (oneMoreTime) {
        indexVoxel = distribution(*randGen);
        oneMoreTime = tabVoxels[indexVoxel].covered;
    }
    return &(tabVoxels[indexVoxel]);
}

template <unsigned short DIM>
algoRSA_aux::AlgoRSA_Template<DIM>::AlgoRSA_Template(
    AmbiantSpace::BigShape<DIM>* T_, RadiusGenerator<DIM>* radiusGen_,
    unsigned seed) : bigShape{ T_ }, radiusGen{ radiusGen_ } {
    // radius and volume fraction
    minRadius = radiusGen->minRadius();
    maxRadius = radiusGen->maxRadius();
    // Intialize the random generator, cf https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution
    mt19937 gen(seed);  // Standard mersenne_twister_engine seeded with rd()
    randGenerator = gen;
    randomReal = uniform_real_distribution<>(0., 1.);
    // Define the sizes parameter
    std::array<size_t, DIM> sizes;
    double voxelLength = maxRadius;  // difficult to find the optimal choice...
    for (size_t i = 0; i < DIM; i++) {
        sizes[i] = floor(bigShape->L[i] / voxelLength) + 1;
    }
    motherGrid = new MotherGrid<DIM>(sizes, bigShape, voxelLength, minRadius,
        maxRadius);
    currentGrid = new CurrentGrid<DIM>(motherGrid, &randGenerator);
    placedSpheres = &(motherGrid->placedSpheres);
    ///
    nbMaxVoxels = MULT_MAX * currentGrid->tabVoxels.size();
}

template <unsigned short DIM>
algoRSA_aux::AlgoRSA_Template<DIM>::~AlgoRSA_Template() {
    delete motherGrid;
    delete currentGrid;
}

template <unsigned short DIM>
inline map<string, string> algoRSA_aux::AlgoRSA_Template<DIM>::proceed(unsigned short method) {
    bool packed = false;
    switch (method) {
    case 1:
        packed = proceed_T<1>();
        break;
    case 2:
        packed = proceed_T<2>();
        break;
    default:
        throw invalid_argument(__PRETTY_FUNCTION__);
    }
    return map<string, string>{ { "Packed", packed ? "True" : "False" }};
}

template <unsigned short DIM>
template <unsigned short NUM_METHOD>
bool algoRSA_aux::AlgoRSA_Template<DIM>::proceed_T() {
    size_t cptSubDiv = 0;
    // printStep(); -> for testing
    while (currentGrid->nbUncoveredVoxels > 0 and cptSubDiv < NB_MAX_SUB_DIV) {
        for (long i = 0; i < 1 + PARAM_A * currentGrid->nbUncoveredVoxels;
            i++) {
            if (throwDart<NUM_METHOD>()) {
                if (not radiusGen->nextRadius()) {  /// not more radius to be placed
                    printFinalMessage(true);
                    return false;
                }
            }
        }
        currentGrid->CurrentGrid<DIM>::template subdivide<NUM_METHOD>();
        cptSubDiv++;
        // printStep(); -> for testing
        if (currentGrid->nbUncoveredVoxels > nbMaxVoxels and currentGrid->nbUncoveredVoxels > NB_MAX_SUB_VOXELS) {
            string s =
                "algoRSA<DIM>::proceed_Template. The grid has been too much subdivided (more than NB_MAX_SUB_VOXELS voxels). \n "
                "It is likely that : \n"
                "1) You try to insert too much spheres, and fail to insert spheres that are not of minimal size.\n"
                "     To check 1), suppress the smallest spheres in your requirement and try once more.\n"
                "2) The constant algoRSA<DIM>::PARAM_A is too low.\n"
                "3) The gap between the maximal and minimal radii is too large. \n"
                "4) There is a bug.";
            throw runtime_error(s);
        }
    }
    printFinalMessage();
    return (currentGrid->nbUncoveredVoxels == 0);
}

template <unsigned short DIM>
template <unsigned short NUM_METHOD>
bool algoRSA_aux::AlgoRSA_Template<DIM>::throwDart() {
    double radius = radiusGen->getRadius();
    PhaseType phase = radiusGen->getPhase();
    bool isNewBallAdded = false;
    Voxel<DIM>* pickedVoxel = currentGrid->pickUncoveredVoxel();
    if (pickedVoxel->mother->covered) {
        pickedVoxel->covered = true;
        currentGrid->nbUncoveredVoxels--;
    } else {
        Point<DIM> point = pickedVoxel->pickPoint(randGenerator, randomReal);
        Sphere<DIM> sphere(point, radius, phase);
        if (not pickedVoxel->mother->isClose2Boundary or (bigShape->isInside(point, minRadius) and bigShape->isInside(sphere))) {
            if (motherGrid->MotherGrid<DIM>::template isPlaceable<NUM_METHOD>(sphere,
                pickedVoxel->mother)) {
                motherGrid->MotherGrid<DIM>::template addSphere<NUM_METHOD>(sphere,
                    pickedVoxel->mother);
                isNewBallAdded = true;
                pickedVoxel->updateCovered(sphere, minRadius);
                if (pickedVoxel->covered) {
                    currentGrid->nbUncoveredVoxels--;
                }
            }
        }
    }
    return isNewBallAdded;
}

template <unsigned short DIM>
void algoRSA_aux::CurrentGrid<DIM>::setTabVoxels(
    vector<Voxel<DIM>>&& tabVoxels_) {
    tabVoxels = std::move(tabVoxels_);
    uniform_int_distribution<> distrib(0, tabVoxels.size() - 1);
    distribution = distrib;
    nbUncoveredVoxels = tabVoxels.size();
}

//////////////////////////////////////////////////////////
//// Beginning of functions unimportant for performance.//
//////////////////////////////////////////////////////////

template <unsigned short DIM>
void algoRSA_aux::AlgoRSA_Template<DIM>::printFinalMessage(
    bool noMoreSpheres) const {
    cout << "----------------------------" << endl;
    cout << "AlgoRSA::proceed()" << endl;
    if (noMoreSpheres) {
        cout << "Stopped because no more spheres should be placed." << endl;
    }
    if (currentGrid->nbUncoveredVoxels > 0) {
        cout << "Not fully packed" << endl;
    } else {
        cout << "Fully packed" << endl;
    }
    cout << "Nb of spheres : " << motherGrid->placedSpheres.size() << endl;
}

template <unsigned short DIM>
void algoRSA_aux::AlgoRSA_Template<DIM>::printStep() const {
    cout << "Uncovered voxels  = " << currentGrid->nbUncoveredVoxels << endl;
    cout << "Nb of spheres = " << (*placedSpheres).size() << endl;
}
}  // namespace sac_de_billes



//! Copyright : see license.txt
//!
//! \brief Implementations for class Voxel and MotherVoxel
//
#ifndef VOXEL_IXX_
#define VOXEL_IXX_


namespace sac_de_billes {
using namespace std;

template<unsigned short DIM>
algoRSA_aux::Voxel<DIM>::Voxel(const DiscPoint<DIM>& discreteCoordinates_,
    double length_, MotherGrid<DIM>* motherGrid_):
    discreteCoordinates(discreteCoordinates_), length(length_), covered(
        false), motherGrid(motherGrid_), mother(nullptr) {}

template<unsigned short DIM>
Point<DIM> algoRSA_aux::Voxel<DIM>::pickPoint(mt19937& randGenerator,
    uniform_real_distribution<>& randomReal) {
    Point <DIM> point({ 0. });
    for (size_t i = 0; i < DIM; i++) {
        point[i] = (discreteCoordinates[i] + randomReal(randGenerator))
            * length;
    }
    return point;
}

template<unsigned short DIM>
void algoRSA_aux::Voxel<DIM>::updateCovered(const Sphere<DIM>& sphere,
    double minRadius) {
    covered = covered or isInSphere(sphere, minRadius);
}

template<unsigned short DIM>
Point<DIM> algoRSA_aux::Voxel<DIM>::corner(int indice) const {
    if constexpr (DIM == 3) {
        return (corner(Path_auxi::TABCORNER3D[indice]));
    }
    else if constexpr (DIM == 2) {
        return (corner(Path_auxi::TABCORNER2D[indice]));
    }
    else {
        throw invalid_argument(
            "algoRSA_aux::Voxel<DIM>::corner. Incorrect DIM");
    }
}

template<unsigned short DIM>
Point<DIM> algoRSA_aux::Voxel<DIM>::corner(
    const array<unsigned short, DIM>& discCoord) const {
    Point <DIM> coordCorner{ };
    for (size_t i = 0; i < DIM; i++) {
        coordCorner[i] = (discreteCoordinates[i] + discCoord[i]) * length;
    }
    return coordCorner;
}

template<unsigned short DIM>
bool algoRSA_aux::Voxel<DIM>::isInSphere(const Sphere<DIM>& sphere,
    const double& minRadius,
    const vector<array<unsigned short, DIM>>& pathForCorner) const {
    double sphereInfluenceRadiusSquared = (sphere.radius + minRadius)
        * (sphere.radius + minRadius);
    for (const auto& i : pathForCorner) {
        if (mother->motherGrid->bigShape->distanceCarre(sphere.center,
            corner(i)) > sphereInfluenceRadiusSquared) {
            return false;
        }
    }
    return true;
}

template<unsigned short DIM>
template<unsigned short NUM_METHOD>
void algoRSA_aux::Voxel<DIM>::subdivide(
    vector<algoRSA_aux::Voxel<DIM>>& futureTabVoxel) const {
    algoRSA_aux::VoxelSubdivision <DIM> subdivision(this);
    if constexpr (NUM_METHOD == 1) {
        subdivision.findSpheres();
    }
    if (not subdivision.isEmpty()) {
        if constexpr (NUM_METHOD == 1) {
            subdivision.checkCoveredVoxel_1();
        }
        else if constexpr (NUM_METHOD == 2) {
            subdivision.checkCoveredVoxel_2();
        }
        else {
            throw invalid_argument(__PRETTY_FUNCTION__);
        }
    }
    subdivision.fillNewVoxels(futureTabVoxel);
}

template<unsigned short DIM>
bool algoRSA_aux::Voxel<DIM>::isOutside(AmbiantSpace::BigShape<DIM>* bigShape,
    double minRadius) {
    Point <DIM> point{ };
    bool answer = true;
    if constexpr (DIM == 3) {
        for (auto i : getIndices<2, 2, 2>()) {
            for (size_t j = 0; j < DIM; j++) {
                point[j] = (discreteCoordinates[j] + i[j]) * length;
            }
            answer = answer and not bigShape->isInside(point, minRadius);
        }
    }
    else if constexpr (DIM == 2) {
        for (auto i : getIndices<2, 2>()) {
            for (size_t j = 0; j < DIM; j++) {
                point[j] = (discreteCoordinates[j] + i[j]) * length;
            }
            answer = answer and not bigShape->isInside(point, minRadius);
        }
    }
    return answer;
}

template<unsigned short DIM>
Point<DIM> algoRSA_aux::Voxel<DIM>::center() const {
    Point <DIM> centerPoint{ };
    for (size_t i = 0; i < DIM; i++) {
        centerPoint[i] = (discreteCoordinates[i] + 0.5) * length;
    }
    return centerPoint;
}

template<unsigned short DIM>
bool algoRSA_aux::Voxel<DIM>::isCovered(const vector<Sphere<DIM>>* sphereList,
    double minRadius) const {
    for (const auto& sphere : *sphereList) {
        if (isInSphere(sphere, minRadius)) {
            return true;
        }
    }
    return false;
}

template<unsigned short DIM>
Point<DIM> algoRSA_aux::Voxel<DIM>::position() const {
    Point <DIM> pos;
    for (size_t j = 0; j < DIM; j++) {
        pos[j] = discreteCoordinates[j] * length;
    }
    return pos;
}

// MotherVoxel

template<unsigned short DIM>
algoRSA_aux::MotherVoxel<DIM>::MotherVoxel(
    const DiscPoint<DIM>& discreteCoordinates_, double length_,
    MotherGrid<DIM>* motherGrid_):
    algoRSA_aux::Voxel<DIM>(discreteCoordinates_, length_, motherGrid_), spheresInside(
        { }), spheresRelativePosition({ }), isClose2Boundary(false) {}

template<unsigned short DIM>
inline Sphere<DIM> algoRSA_aux::MotherVoxel<DIM>::getSphere(int i) const {
    return this->motherGrid->placedSpheres[spheresInside[i]];
}

template<unsigned short DIM>
inline bool algoRSA_aux::MotherVoxel<DIM>::isItClose2Boundary(
    AmbiantSpace::BigShape<DIM>* bigShape) {
    Point <DIM> p{ }; // this will be the coordinates of the surrounding cube
    bool answer = false;
    if constexpr (DIM == 3) {
        loop(2, 2, 2,
            [&](const array<size_t, 3u>& i) {
                for (size_t j = 0; j < DIM; j++) {
                    p[j] = (this->discreteCoordinates[j] + 3 * i[j]
                        - 2 * (1 - i[j])) * this->length;
                }
        answer = answer or not (bigShape->isInside(p, 0)); // tests if the corner of the big cube is inside the shape
            });
    }
    else if constexpr (DIM == 2) {
        loop(2, 2,
            [&](const array<size_t, 2u>& i) {
                for (size_t j = 0; j < DIM; j++) {
                    p[j] = (this->discreteCoordinates[j] + 3 * i[j]
                        - 2 * (1 - i[j])) * this->length;
                }
        answer = answer or not (bigShape->isInside(p, 0)); // tests if the corner of the big cube is inside the shape
            });
    }
    else {
        throw invalid_argument(__PRETTY_FUNCTION__);
    }
    return answer;
}

// VoxelSubdivision

template<unsigned short DIM>
void algoRSA_aux::VoxelSubdivision<DIM>::initialize() {
    NbUncoveredVoxels = SIZE;
    algoRSA_aux::Voxel <DIM> newVoxel = *fatherVoxel;
    newVoxel.length = 0.5 * fatherVoxel->length;
    size_t compteur = 0;
    if constexpr (DIM == 3) {
        loop<0, 2, 0, 2, 0, 2>(
            [&](const array<long, 3u>& i) {
                for (size_t j = 0; j < DIM; j++) {
                    newVoxel.discreteCoordinates[j] = 2
                        * fatherVoxel->discreteCoordinates[j] + i[j];
                }
        tabNewVoxels[compteur] = newVoxel;
        tabUnCovered[compteur] = true;
        compteur++;
            });
    }
    else if constexpr (DIM == 2) {
        loop<0, 2, 0, 2>(
            [&](const array<long, 2u>& i) {
                for (size_t j = 0; j < DIM; j++) {
                    newVoxel.discreteCoordinates[j] = 2
                        * fatherVoxel->discreteCoordinates[j] + i[j];
                }
        tabNewVoxels[compteur] = newVoxel;
        tabUnCovered[compteur] = true;
        compteur++;
            });
    }
    else {
        throw invalid_argument(
            "VoxelSubdivision<DIM>::initialize. DIM not 2 or 3.");
    }

    for (size_t j = 0; j < SIZE; j++) {
        if (fatherVoxel->mother->isClose2Boundary
            and tabNewVoxels[j].isOutside(fatherVoxel->motherGrid->bigShape,
                minRadius)) {
            setFalse(j);
        }
    }
}

template<unsigned short DIM>
void algoRSA_aux::VoxelSubdivision<DIM>::findSpheres() {
    closeSpheres = fatherVoxel->motherGrid->neighborSpheres(
        fatherVoxel->mother);
    closeSpheresRelativePos = fatherVoxel->mother->spheresRelativePosition;
}

template<unsigned short DIM>
bool algoRSA_aux::VoxelSubdivision<DIM>::isEmpty() {
    return (NbUncoveredVoxels == 0);
}

template<unsigned short DIM>
void algoRSA_aux::VoxelSubdivision<DIM>::setFalse(int j) {
    tabUnCovered[j] = false;
    NbUncoveredVoxels--;
}

template<unsigned short DIM>
void algoRSA_aux::VoxelSubdivision<DIM>::checkCoveredVoxel_1() {
    Sphere <DIM> sphere{ };
    int relativePos{ };
    for (size_t i = 0; i < closeSpheres.size(); i++) {
        sphere = closeSpheres[i];
        relativePos = closeSpheresRelativePos[i];
        if (hitsMidPoint(sphere)) {
            // the sphere has to intersect the middle point of the Voxel to contain any subvoxel
            for (size_t j = 0; j < SIZE; j++) {
                if (testSphere < 1 >(sphere, j, relativePos)) {
                    return;
                }
            }
        }
    }
}

template<unsigned short DIM>
void algoRSA_aux::VoxelSubdivision<DIM>::checkCoveredVoxel_2() {
    if (not fatherVoxel->mother->isClose2Boundary) {
        testSubVoxels<5>();
    }
    else {
        testSubVoxels<7>();
    }
}

template<unsigned short DIM>
void algoRSA_aux::VoxelSubdivision<DIM>::fillNewVoxels(
    vector<algoRSA_aux::Voxel<DIM>>& futureTabVoxel) {
    for (size_t j = 0; j < SIZE; j++) {
        if (tabUnCovered[j]) {
            futureTabVoxel.push_back(tabNewVoxels[j]);
        }
    }
}

template<unsigned short DIM>
bool algoRSA_aux::VoxelSubdivision<DIM>::checkSphere(const int& j,
    const Sphere<DIM>& sphere,
    const vector<array<short unsigned, DIM>>& pathForCorner) {
    if (tabUnCovered[j]) {
        if (tabNewVoxels[j].isInSphere(sphere, minRadius, pathForCorner)) {
            setFalse(j);
            return isEmpty();
        }
    }
    return false;
}

template<unsigned short DIM>
bool algoRSA_aux::VoxelSubdivision<DIM>::hitsMidPoint(Sphere<DIM>& sphere) {
    return (fatherVoxel->motherGrid->bigShape->distanceCarre(sphere.center,
        midPoint)
        < (sphere.radius + minRadius) * (sphere.radius + minRadius));
}

template<unsigned short DIM>
template<unsigned short NB_NGHB>
void algoRSA_aux::VoxelSubdivision<DIM>::testSubVoxels() {
    const array<array<int, DIM>, nbSubcubes<DIM>(NB_NGHB)>& pathForNeighb =
        Path::get().pathForNeighbors<DIM, NB_NGHB>();
    Sphere <DIM> sphere;
    DiscPoint<DIM> index({ 0 });

    for (size_t relPos = 0; relPos < pathForNeighb.size(); relPos++) {
        for (size_t j = 0; j < DIM; j++) {
            index[j] = auxi_function::fast_modulo(
                fatherVoxel->mother->discreteCoordinates[j]
                + pathForNeighb[relPos][j],
                fatherVoxel->motherGrid->sizes[j]);
        }
        for (size_t iSphe : fatherVoxel->motherGrid->getVoxel(index)->spheresInside) {
            sphere = fatherVoxel->motherGrid->placedSpheres[iSphe];
            if (hitsMidPoint(sphere)) {
                // the sphere has to intersect the middle point of the Voxel to contain any subvoxel
                for (size_t j = 0; j < SIZE; j++) {
                    if (testSphere<2, NB_NGHB>(sphere, j, relPos)) {
                        return;
                    }
                }
            }
        }
    }
}

template<unsigned short DIM>
template<unsigned short NUM_METHOD, unsigned short NB_NGHB>
bool algoRSA_aux::VoxelSubdivision<DIM>::testSphere(const Sphere<DIM>& sphere,
    const size_t j, const int relativePosition) {
    if constexpr (DIM == 3) {
        if constexpr (NUM_METHOD == 1) {
            return checkSphere(j, sphere,
                Path::get().pathFor27Corners[relativePosition]);
        }
        else if constexpr (NUM_METHOD == 2) {
            if constexpr (NB_NGHB == 5) {
                return (checkSphere(j, sphere,
                    Path::get().pathForCorners[relativePosition]));
            }
            else {
                return (checkSphere(j, sphere, Path::get().pathForCorners[0]));
            }
        }
    }
    else if constexpr (DIM == 2) { //\fixme not optimized for 2D.
        return checkSphere(j, sphere, Path_auxi::get().TABCORNER<DIM>());
    }
    else {
        throw invalid_argument(
            "VoxelSubdivision<DIM>::testSubVoxels. Incorrect DIM");
    }
}
} // namespace sac_de_billes

#endif /* VOXEL_IXX_ */

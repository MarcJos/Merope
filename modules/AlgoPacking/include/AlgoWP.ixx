//! Copyright : see license.txt
//!
//! \brief
//

/// IntersectedSpheres<DIM>

// algoWP_aux::Pair

#pragma once

namespace sac_de_billes {
inline bool algoWP_aux::Pair::operator<(const algoWP_aux::Pair& pair2) const {
    return std::lexicographical_compare(this->begin(), this->end(), pair2.begin(), pair2.end());
}

/// IntersectedSpheres

template<unsigned short DIM>
inline double algoWP_aux::IntersectedSpheres<DIM>::energy_J() const {
    double J = 0., deltaij_;
    for (const auto& pair : listOfPairs) {
        const Sphere<DIM>& sph1 = motherGrid->placedSpheres[listOf[pair[0]]];
        const Sphere<DIM>& sph2 = motherGrid->placedSpheres[listOf[pair[1]]];
        deltaij_ = algoWP_aux::deltaij(bigShape, sph1, sph2);
        J += deltaij_ * deltaij_;
    }
    return J;
}

template<unsigned short DIM>
inline void algoWP_aux::IntersectedSpheres<DIM>::buildLists() {
    __gnu_parallel::sort(listOfPairs.begin(), listOfPairs.end());
    listOfPairs.erase(std::unique(listOfPairs.begin(), listOfPairs.end()), listOfPairs.end());
    listOf.resize(listOfPairs.size() * 2);

#pragma omp parallel for
    for (size_t idx = 0; idx < listOfPairs.size(); idx++) {
        listOf[2 * idx] = listOfPairs[idx][0];
        listOf[2 * idx + 1] = listOfPairs[idx][1];
    }

    __gnu_parallel::sort(listOf.begin(), listOf.end());
    listOf.erase(std::unique(listOf.begin(), listOf.end()), listOf.end());
}

template<unsigned short DIM>
inline void algoWP_aux::IntersectedSpheres<DIM>::updatePositions() {
    //! \fixme : not really equivalent to the other version!
    vector<Point<DIM>> displacements(listOfPairs.size());
    // compute displacements
#pragma omp parallel for
    for (size_t l = 0; l < listOfPairs.size(); l++) {
        const auto& pair = listOfPairs[l];
        size_t i = pair[0];
        size_t j = pair[1];
        if (i < j) {
            const Sphere <DIM>& sph1 = motherGrid->placedSpheres[i];
            const Sphere <DIM>& sph2 = motherGrid->placedSpheres[j];
            const auto& geomVector = bigShape->geomVector(sph1.center,
                sph2.center);
            double distance = sqrt(bigShape->normeCarre(geomVector));
            double factor = h_multiplier * 2 * algoWP_aux::deltaij(sph1, sph2, distance)
                / distance;
            displacements[l] = factor * geomVector;
        }
    }
    // apply displacements
#pragma omp parallel for
    for (size_t l = 0; l < listOfPairs.size(); l++) {
        const auto& pair = listOfPairs[l];
        size_t i = pair[0];
        size_t j = pair[1];
        if (i < j) {
            Sphere <DIM>& sph1 = motherGrid->placedSpheres[i];
            Sphere <DIM>& sph2 = motherGrid->placedSpheres[j];
            for (size_t k = 0; k < DIM; k++) {
#pragma omp atomic
                sph1.center[k] -= displacements[l][k];
#pragma omp atomic
                sph2.center[k] += displacements[l][k];
            }
        }
    }
    for (auto i : listOf) {
        bigShape->bounce(motherGrid->placedSpheres[i]);
    }
}

template<unsigned short DIM>
inline void algoWP_aux::IntersectedSpheres<DIM>::setListOfPairs(vector<algoWP_aux::Pair>&& a_listOfPairs) {
    this->listOfPairs = std::move(a_listOfPairs);
    buildLists();
}

template<unsigned short DIM>
inline vector<Point<DIM> > algoWP_aux::IntersectedSpheres<DIM>::energy_NablaJ() const {
    map < size_t, size_t > fromIndexToListOf({ });
    for (size_t i = 0; i < listOf.size(); i++) {
        fromIndexToListOf[listOf[i]] = i;
    }
    vector<Point<DIM>> nablaJ(listOf.size(), Point<DIM> { });
    for (const auto& pair : listOfPairs) {
        size_t i = pair[0];
        size_t j = pair[1];
        if (i < j) {
            const Sphere<DIM>& sph1 = motherGrid->placedSpheres[i];
            const Sphere<DIM>& sph2 = motherGrid->placedSpheres[j];
            const auto& geomVector = bigShape->geomVector(sph1.center,
                sph2.center);
            double distance = sqrt(bigShape->normeCarre(geomVector));
            double factor = 2 * algoWP_aux::deltaij(sph1, sph2, distance)
                / distance;
            for (auto k = 0; k < DIM; k++) {
                double nablaJ_ij_k = factor * geomVector[k];
                nablaJ[fromIndexToListOf[i]][k] += nablaJ_ij_k;
                nablaJ[fromIndexToListOf[j]][k] -= nablaJ_ij_k;
            }
        }
    }
    return nablaJ;
}

template<unsigned short DIM>
inline vector<Sphere<DIM>> algoWP_aux::IntersectedSpheres<DIM>::newPositions() const {
    vector < Sphere <DIM>> newSpheres(listOf.size());
    auto nablaJ = energy_NablaJ();

    // preparing the parallel loop
    double h_multiplier_loc = h_multiplier;
    //
#pragma omp for firstprivate(h_multiplier_loc)
    for (size_t i = 0; i < listOf.size(); i++) {
        const Sphere<DIM>& oldsph = motherGrid->placedSpheres[listOf[i]];
        Sphere <DIM>& newsph = newSpheres[i];
        newsph.radius = oldsph.radius;
        for (auto k = 0; k < DIM; k++) {
            newsph.center[k] = oldsph.center[k] - h_multiplier * nablaJ[i][k];
        }
        bigShape->bounce(newsph);
    }
    return newSpheres;
}

/// WPGrid<DIM>

template<unsigned short DIM>
inline void WPGrid<DIM>::shrink(double shrinkFactor) {
    Sphere <DIM> oldSphere{ };
    size_t imax = this->placedSpheres.size();
#pragma omp parallel for firstprivate(oldSphere, shrinkFactor)
    for (size_t i = 0; i < imax; i++) {
        Sphere <DIM>& sph = this->placedSpheres[i];
        oldSphere = sph;
        sph.radius *= shrinkFactor;
        this->bigShape->bounce(sph);
        updateSphereIndex(i, oldSphere, sph);
    }
}

template<unsigned short DIM>
inline vector<size_t> WPGrid<DIM>::getNearSpheres(
    algoRSA_aux::MotherVoxel<DIM>* centralvoxel) {
    auto mapped_i_voxel = this->surroundingVoxels(centralvoxel, 1);
    vector < size_t > theSpheres{ };
    for (const auto& [i, voxel] : mapped_i_voxel) {
        for (size_t sphe : (voxel->spheresInside)) {
            theSpheres.push_back(sphe);
        }
    }
    return theSpheres;
}

template<unsigned short DIM>
inline WPGrid<DIM>::WPGrid(DiscPoint<DIM> sizes_,
    AmbiantSpace::BigShape<DIM>* bigShape_, double voxelLength_,
    double minRadius_, double maxRadius_) :
    algoRSA_aux::MotherGrid<DIM>(auxi_function::convertArray<DIM, long, size_t>(sizes_), bigShape_, voxelLength_,
        minRadius_, maxRadius_), activatedVoxels{ } {
    //\fixme : maybe shrink should decrease when reaching the desired radius?
}

template<unsigned short DIM>
inline void WPGrid<DIM>::resetActivatedVoxels() {
    activatedVoxels.resize(this->getNbVoxels());
    for (size_t index = 0; index < this->tabVoxels.size(); index++) {
        activatedVoxels[index] = &(this->tabVoxels[index]);
    }
}

template<unsigned short DIM>
inline algoWP_aux::IntersectedSpheres<DIM> WPGrid<DIM>::findIntersectedSpheres() {
    algoWP_aux::IntersectedSpheres<DIM> intersectedSpheres(this, this->bigShape);
    intersectedSpheres.setListOfPairs(findIntersectedSpheres_listOfPairs());
    return intersectedSpheres;
}

template<unsigned short DIM>
vector<algoWP_aux::Pair> WPGrid<DIM>::findIntersectedSpheres_listOfPairs() {
    vector<algoWP_aux::Pair> listOfPairs = {};
    {
#pragma omp declare reduction (merge : vector<algoWP_aux::Pair> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp parallel for reduction(merge : listOfPairs)
        for (const auto& voxel : activatedVoxels) {
            vector < size_t > nearSpheresIndex = this->getNearSpheres(voxel);
            for (size_t sphe0Index : voxel->spheresInside) {
                const Sphere<DIM>& sphe0 = this->placedSpheres[sphe0Index];
                for (size_t sphe1Index : nearSpheresIndex) {
                    const Sphere<DIM>& sphe1 = this->placedSpheres[sphe1Index];
                    if ((sphe0Index != sphe1Index) && this->bigShape->areSphereIntersected(sphe0, sphe1)) {
                        listOfPairs.emplace_back(algoWP_aux::Pair(min(sphe0Index, sphe1Index), max(sphe0Index, sphe1Index)));
                    }
                }
            }
        }
    }
    return listOfPairs;
}

template<unsigned short DIM>
inline void WPGrid<DIM>::removeOverlap() {
    resetActivatedVoxels();
    size_t compteur = 0;
    while (compteur < MAX_ITERATIONS and activatedVoxels.size() != 0) {
        moveOverlappingSpheres<NUM_METHOD_OVERLAP>();
        compteur++;
    }
    /// Display
    cout << "WPGrid<DIM>::removeOverlap : Number of iterations = " << compteur << endl;
    if (compteur == MAX_ITERATIONS) {
        cerr << "Number of spheres = " << this->placedSpheres.size() << endl;
        cerr << "Radius = " << this->placedSpheres[0].radius << endl;
        throw runtime_error("I can't pack!"
            " The desired volume fraction is probably too high!");
    }
}

template<unsigned short DIM>
inline void WPGrid<DIM>::updateSphereIndex(size_t indexSphere,
    const Sphere<DIM>& oldSphere, const Sphere<DIM>& newSphere) {
    auto departureVoxel = this->getVoxel(
        this->discCoordinate_withoutModulo(oldSphere.center));
    auto arrivalVoxel = this->getVoxel(
        this->discCoordinate_withoutModulo(newSphere.center));
    updateSphereIndex(indexSphere, departureVoxel, arrivalVoxel);
}

template<unsigned short DIM>
inline void WPGrid<DIM>::updateSphereIndex(size_t indexSphere,
    algoRSA_aux::MotherVoxel<DIM>* departureVoxel,
    algoRSA_aux::MotherVoxel<DIM>* arrivalVoxel) {
    if (arrivalVoxel != departureVoxel) {
#pragma omp critical
        {
            bool indexFound = false;
            for (size_t i = 0; i < departureVoxel->spheresInside.size(); i++) {
                if (departureVoxel->spheresInside[i] == indexSphere) {
                    indexFound = true;
                    departureVoxel->spheresInside.erase(
                        departureVoxel->spheresInside.begin() + i);
                    break;
                }
            }
            if (not indexFound) {
                throw runtime_error(__PRETTY_FUNCTION__);
            }
            arrivalVoxel->spheresInside.push_back(indexSphere);
        }
    }
}

template<unsigned short DIM>
template<unsigned short NUM_METHOD_OVERLAP_>
inline void WPGrid<DIM>::moveOverlappingSpheres() {
    if constexpr (NUM_METHOD_OVERLAP_ == 1) {
        return moveOverlappingSpheres_1();
    } else if constexpr (NUM_METHOD_OVERLAP_ == 2) {
        return moveOverlappingSpheres_2();
    } else {
        throw invalid_argument(__PRETTY_FUNCTION__);
    }
}

template<unsigned short DIM>
inline void WPGrid<DIM>::moveOverlappingSpheres_1() {
    auto intersectedSpheres = findIntersectedSpheres();
    vector < Sphere <DIM>> oldSpheres(intersectedSpheres.listOf.size());

#pragma omp parallel for
    for (size_t i = 0; i < oldSpheres.size(); i++) {
        oldSpheres[i] = this->placedSpheres[intersectedSpheres.listOf[i]];
    }

    intersectedSpheres.updatePositions();
    activatedVoxels.resize(intersectedSpheres.listOf.size());

#pragma omp parallel for
    for (size_t i = 0; i < intersectedSpheres.listOf.size(); i++) {
        size_t indexSphere = intersectedSpheres.listOf[i];
        const Sphere<DIM>& oldSphere = oldSpheres[i];
        const Sphere<DIM>& newSphere = this->placedSpheres[indexSphere];
        auto departureVoxel = this->getVoxel(
            this->discCoordinate_withoutModulo(oldSphere.center));
        auto arrivalVoxel = this->getVoxel(
            this->discCoordinate_withoutModulo(newSphere.center));
        activatedVoxels[i] = arrivalVoxel;
        updateSphereIndex(indexSphere, departureVoxel, arrivalVoxel);
    }
}

template<unsigned short DIM>
inline void WPGrid<DIM>::moveOverlappingSpheres_2() {
    auto intersectedSpheres = findIntersectedSpheres();
    auto newSpheres = intersectedSpheres.newPositions();
    activatedVoxels.resize(intersectedSpheres.listOf.size());
    for (size_t i = 0; i < intersectedSpheres.listOf.size(); i++) {
        size_t indexSphere = intersectedSpheres.listOf[i];
        const Sphere<DIM>& oldSphere = this->placedSpheres[indexSphere];
        const Sphere<DIM>& newSphere = newSpheres[i];
        auto departureVoxel = this->getVoxel(
            this->discCoordinate_withoutModulo(oldSphere.center));
        auto arrivalVoxel = this->getVoxel(
            this->discCoordinate_withoutModulo(newSphere.center));
        this->placedSpheres[indexSphere] = newSpheres[i];
        activatedVoxels[i] = arrivalVoxel;
        updateSphereIndex(indexSphere, departureVoxel, arrivalVoxel);
    }
}

template<unsigned short DIM>
inline map<string, string> algoWP_aux::AlgoWP_Template<DIM>::proceed(unsigned short method) {
    initialize(method);
    double Delta_RenormVolumeFraction = (1 - pow(RADIUS_MULT, DIM))
        / NUMBER_OF_STEPS;
    double RenormVolumeFraction = pow(RADIUS_MULT, DIM);
    double RenormVolumeFraction_next = pow(RADIUS_MULT, DIM);
    for (size_t i = 0; i < NUMBER_OF_STEPS; i++) {
        RenormVolumeFraction_next += Delta_RenormVolumeFraction;
        double shrinkFactor = pow(
            RenormVolumeFraction_next / RenormVolumeFraction, 1. / DIM);
        RenormVolumeFraction = RenormVolumeFraction_next;
        grid.get()->shrink(shrinkFactor);
        grid.get()->removeOverlap();
    }
    printFinalMessage();
    return map<string, string>{ { "Packed", "False" }};
}

template<unsigned short DIM>
inline algoWP_aux::AlgoWP_Template<DIM>::AlgoWP_Template(
    AmbiantSpace::BigShape<DIM>* bigShape_,
    algoRSA_aux::RadiusGenerator<DIM>* radiusGen_, unsigned seed_) :
    placedSpheres{ nullptr }, bigShape{ bigShape_ }, radiusGen{
            radiusGen_ }, seed{ seed_ } {
    DiscPoint <DIM> sizes;
    double maxRadius = radiusGen->maxRadius();
    double minRadius = radiusGen->minRadius();
    double voxelLength = 2.0001 * maxRadius;  // difficult to find the optimal choice...
    if (voxelLength < 2 * maxRadius) {
        throw runtime_error(
            "Algorithm badly parametrized"
            "It is important that potentially intersecting spheres are only in the closest voxels");
    }
    for (size_t i = 0; i < DIM; i++) {
        sizes[i] = ceil(bigShape->L[i] / voxelLength);
    }
    grid.reset(
        new WPGrid<DIM>(sizes, bigShape, voxelLength, minRadius,
            maxRadius));
    placedSpheres = &(grid.get()->placedSpheres);
}

template<unsigned short DIM>
inline void algoWP_aux::AlgoWP_Template<DIM>::printFinalMessage() {
    cout << "----------------------------" << endl;
    cout << "AlgoWP::proceed()" << endl;
}

template<unsigned short DIM>
inline void algoWP_aux::AlgoWP_Template<DIM>::initialize(
    unsigned short method) {
    auto placedSpheresRSA = initPlacedSpheres(method);
    *placedSpheres = placedSpheresRSA;
    for (size_t i = 0; i < placedSpheresRSA.size(); i++) {
        const auto& sph = placedSpheresRSA[i];
        auto voxel = grid.get()->getVoxel(
            grid.get()->discCoordinate(sph.center));
        voxel->spheresInside.push_back(i);
    }
}

template<unsigned short DIM>
vector<Sphere<DIM>> algoWP_aux::AlgoWP_Template<DIM>::initPlacedSpheres(
    unsigned short method) {
    AlgoRSA <DIM> algoRSABase{ };
    Point<DIM> L_loc = bigShape->L;
    algoRSABase.setBigShape(vector<double>(L_loc.begin(), L_loc.end()),
        bigShape->type());
    auto radiusGenCopy = *radiusGen;
    radiusGenCopy.shrink(RADIUS_MULT);
    algoRSABase.setRadiusGenerator(radiusGenCopy);
    algoRSABase.proceed(seed, method);
    return algoRSABase.getSpheres();
}
}  // namespace sac_de_billes



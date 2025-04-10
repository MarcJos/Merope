//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace sac_de_billes {

template<unsigned short DIM>
unique_ptr<AlgoPacking<DIM>> algoSpheres::newAlgo(
    const algoSpheres::TypeAlgo& typeAlgo) {
    unique_ptr < AlgoPacking <DIM>> algo;
    switch (typeAlgo) {
    case algoSpheres::TypeAlgo::RSA:
        algo = unique_ptr < AlgoPacking <DIM>>(new AlgoRSA<DIM>());
        break;
    case algoSpheres::TypeAlgo::WP:
        algo = unique_ptr < AlgoPacking <DIM>>(new AlgoWP<DIM>());
        break;
    case algoSpheres::TypeAlgo::BOOL:
        algo = unique_ptr < AlgoPacking <DIM>>(new AlgoBool<DIM>());
        break;
    default:
        throw runtime_error(__PRETTY_FUNCTION__);
    }
    return algo;
}

template<unsigned short DIM>
vector<Sphere<DIM>> algoSpheres::throwSpheres(
    const algoSpheres::TypeAlgo& typeAlgo, AmbiantSpace::NameShape nameShape,
    Point<DIM> L, unsigned seed,
    vector<array<double, 2>> desiredRPhi, vector<PhaseType> tabPhases,
    double mindist) {
    auto algo = newAlgo <DIM>(typeAlgo);
    algo->setExclusionDistance(mindist);
    algo->setBigShape(L, nameShape);
    algo->setRadiusGenerator(desiredRPhi, tabPhases);
    algo->proceed(seed, 1);
    auto theSpheres = algo->getSpheres();
    return theSpheres;
}

template<unsigned short DIM>
vector<Sphere<DIM>> algoSpheres::throwSpheres(
    const algoSpheres::TypeAlgo& typeAlgo, AmbiantSpace::NameShape nameShape,
    Point<DIM> L, unsigned seed,
    vector<double> tabRadii, vector<PhaseType> tabPhases,
    double mindist) {
    auto algo = newAlgo <DIM>(typeAlgo);
    auto algo_ptr = algo.get();
    algo_ptr->setExclusionDistance(mindist);
    algo_ptr->setBigShape(L, nameShape);
    algo_ptr->setRadiusGenerator(tabRadii, tabPhases);
    algo_ptr->proceed(seed, 1);
    auto theSpheres = algo_ptr->getSpheres();
    return theSpheres;
}

template<unsigned short DIM>
double algoSpheres::monoDispReferenceRadius(Point<DIM> L,
    size_t seedsNb, double phi, algoSpheres::TypeAlgo typeAlgo) {
    AmbiantSpace::Tore <DIM> torus(L);
    switch (typeAlgo) {
    case algoSpheres::TypeAlgo::RSA:
    case algoSpheres::TypeAlgo::WP:
    {
        double volumeRef = phi * torus.volume()
            / (sphereTools::volumeSphere <DIM>(1.) * seedsNb);
        return pow(volumeRef, 1. / DIM);
    }
    case algoSpheres::TypeAlgo::BOOL:
    {
        double volumeRef = -log(1. - phi) * torus.volume()
            / (sphereTools::volumeSphere <DIM>(1.) * seedsNb);
        return pow(volumeRef, 1. / DIM);
    }
    default:
        throw runtime_error("Unknown TypeAlgo");
    }
}

template<unsigned short DIM>
inline size_t algoSpheres::monoDispNbSpheres(Point<DIM> L,
    double radius, double phi, algoSpheres::TypeAlgo typeAlgo) {
    AmbiantSpace::Tore <DIM> torus(L);
    switch (typeAlgo) {
    case algoSpheres::TypeAlgo::RSA:
    case algoSpheres::TypeAlgo::WP:
    {
        return static_cast<size_t>(phi * torus.volume()
            / sphereTools::volumeSphere <DIM>(radius));
    }
    case algoSpheres::TypeAlgo::BOOL:
    {
        return static_cast<size_t>(-log(1 - phi) * torus.volume()
            / sphereTools::volumeSphere <DIM>(radius));
    }
    default:
        throw runtime_error("Unknown TypeAlgo");
    }
}

template<unsigned short DIM>
inline vector<Sphere<DIM>> algoSpheres::fillMaxRSA(
    AmbiantSpace::NameShape nameShape, Point<DIM> L,
    size_t NbSpheres, unsigned seed, double mindist) {
    double phi = fractionVolThMaxRSA(DIM);
    double radius = monoDispReferenceRadius < DIM
    >(L, NbSpheres, phi, algoSpheres::TypeAlgo::RSA) - 0.5 * mindist;
    if (radius < 0) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Radius cannot be negative!");
    }
    vector < Sphere <DIM>> spheres;
    vector < PhaseType > tabPhases(NbSpheres, 1);
    size_t counter = 0;
    constexpr size_t COUNTER_MAX = 10;
    while (counter < COUNTER_MAX) {
        vector<double> tabRadii(NbSpheres, radius);
        spheres =
            throwSpheres < DIM
            >(algoSpheres::TypeAlgo::RSA, nameShape, L, seed, tabRadii, tabPhases, mindist);
        if (spheres.size() == NbSpheres) {
            break;
        }
        radius *= 0.95;
        counter++;
    }
    if (counter == COUNTER_MAX) {
        throw runtime_error("I cannot achieve the number of spheres.");
    }
    return spheres;
}

}  // namespace sac_de_billes



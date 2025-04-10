//! Copyright : see license.txt
//!
//! \brief
//

namespace sac_de_billes {

template<unsigned DIM>
inline SphereManipulator<DIM>::SphereManipulator(
    vector<Sphere<DIM> > theSpheres_, Point<DIM> length) :
    originalSpheres(theSpheres_), vector_of_translation{
    Point<DIM> { } } {
    this->theSpheres = theSpheres_;
    torus.reset(new AmbiantSpace::Tore<DIM>(length));
}

template<unsigned DIM>
template<class T>
inline SphereManipulator<DIM>::SphereManipulator(
    const AlgoInterface<DIM, T>& myAlgo) :
    SphereManipulator<DIM>(myAlgo.getSpheres(), myAlgo.getLength()) {}

template<unsigned DIM>
inline SphereManipulator<DIM>::SphereManipulator(
    vector<SimpleSphereFormat> theSpheres_, Point<DIM> length) :
    SphereManipulator<DIM>(
        sphereTools::vector_fromArray2Sphere <DIM>(theSpheres_),
        length) {}

template<unsigned DIM>
inline void SphereManipulator<DIM>::translate_back(const Point<DIM>& vector) {
    Point <DIM> vector_inverse = vector;
    for (size_t i = 0; i < DIM; i++) {
        vector_inverse[i] *= -1;
    }
    translate(vector_inverse);
}

template<unsigned DIM>
inline void SphereManipulator<DIM>::translate(const Point<DIM>& vector) {
    Sphere <DIM> mySphere;
    for (size_t j = 0; j < originalSpheres.size(); j++) {
        mySphere = originalSpheres[j];
        for (size_t i = 0; i < DIM; i++) {
            mySphere.center[i] += vector[i];
        }
        torus.get()->projection(mySphere.center);
        this->theSpheres[j] = mySphere;
    }
    for (size_t i = 0; i < DIM; i++) {
        vector_of_translation[i] = vector[i];
    }
}

template<unsigned DIM>
inline double SphereManipulator<DIM>::distance_cube_spheres() const {
    double dist = numeric_limits<double>::max();
    for (const auto& sphere : this->theSpheres) {
        dist = min(dist, distanceToFaceS(sphere));
        dist = min(dist, distanceToEdgeS(sphere));
        dist = min(dist, distanceToCornerS(sphere));
        // avoid the case where the corner is within a sphere
        dist = (cornerWithinSphere(sphere, sphere.radius * 0.5) ? 0. : dist);
    }
    return dist;
}

template<unsigned DIM>
inline double SphereManipulator<DIM>::distanceToFaceS(
    const Sphere<DIM>& sphere) const {
    double dist = numeric_limits<double>::max();
    for (size_t i = 0; i < DIM; i++) {
        dist = min(dist, sphere.center[i]);						// lower face
        dist = min(dist, abs(torus.get()->L[i] - sphere.center[i]));  // upper face
    }
    return dist;
}

template<unsigned DIM>
inline double SphereManipulator<DIM>::distanceToEdgeS(
    const Sphere<DIM>& sphere) const {
    double dist = numeric_limits<double>::max();			// result
    Point < 2 > pt_loc;  // 2D representation of the vector from edge to the center
    double distance_loc;			// desired distance
    for (size_t i = 0; i < DIM; i++) {					// choose the 1st face
        for (size_t j = i + 1; j < DIM; j++) {			// choose the 2nd face
            for (size_t i1 = 0; i1 < 2; i1++) {
                for (size_t j1 = 0; j1 < 2; j1++) {
                    pt_loc[0] = sphere.center[i] - i1 * torus.get()->L[i];
                    pt_loc[1] = sphere.center[j] - j1 * torus.get()->L[j];
                    distance_loc = abs(
                        sphere.radius
                        - sqrt(
                            pt_loc[0] * pt_loc[0]
                            + pt_loc[1] * pt_loc[1]));
                    dist = min(dist, distance_loc);
                }
            }
        }
    }
    return dist;
}

template<unsigned DIM>
inline double SphereManipulator<DIM>::random_search(size_t nb_tries) {
    upper_bound_on_best_distmin();
    Point<DIM> vector_translate{ };
    Point<DIM> best_vector{ };
    double distance = 0;
    double distmin = 0;
    for (size_t i = 0; i < nb_tries; i++) {
        distance = distance_cube_spheres();
        if (distmin < distance) {
            best_vector = vector_of_translation;
            distmin = distance;
        }
        for (size_t j = 0; j < DIM; j++) {
            vector_translate[j] = getLength()[j] * rand() * 1. / RAND_MAX;
        }
        translate(vector_translate);
    }
    // sets to best_vector
    translate(best_vector);
    // verifies coherency
    distance = distance_cube_spheres();
    if (distance < distmin * (0.999)) {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
    return distance;
}

template<unsigned DIM>
inline double SphereManipulator<DIM>::upper_bound_on_best_distmin() {
    Point<DIM> best_vector;  // candidate best_vector
    double distance_to_faces = numeric_limits<double>::max();
    for (size_t coord = 0; coord < DIM; coord++) {
        vector<double> table_extrema{ };
        vector<double> table_steps{ };
        vector < size_t > table_indices{ };
        for (size_t i = 0; i < this->theSpheres.size(); i++) {
            auto sphere = this->theSpheres[i];
            table_extrema.push_back(sphere.center[coord] - sphere.radius);  // beware, not necessarily in [0,L]
            table_extrema.push_back(sphere.center[coord] + sphere.radius);  // beware, not necessarily in [0,L]
        }
        for (size_t i = 0; i < table_extrema.size(); i++) {
            geomTools::projection_periodic_1D(table_extrema[i],
                getLength()[coord]);  // get back to [0,L]
            table_indices.push_back(i);
        }
        sort(table_extrema.begin(), table_extrema.end());  // sorts the extrema
        for (size_t i = 0; i < table_extrema.size() - 1; i++) {
            table_steps.push_back(table_extrema[i + 1] - table_extrema[i]);
        }
        table_steps[table_extrema.size() - 1] = getLength()[coord]
            + table_extrema[0] - table_extrema[table_extrema.size() - 1];  // periodicity
        sort(table_indices.begin(), table_indices.end(),
            [&](size_t k, size_t l) {
                return table_steps[k] > table_steps[l];
            });
        // finds the larger steps between spheres
        size_t i = 0;
        if (table_indices[i] != table_extrema.size() - 1) {
            best_vector[coord] = 0.5
                * (table_extrema[table_indices[i]]
                    + table_extrema[table_indices[i] + 1]);  // in the middle!
        } else {
            best_vector[coord] = 0.5
                * (table_extrema[0]
                    + table_extrema[table_extrema.size() - 1]
                    + getLength()[coord]);  // in the middle!
        }
        translate(best_vector);
        distance_to_faces = min(0.5 * table_steps[table_indices[0]],
            distance_to_faces);
    }
    return distance_to_faces;
}

template<unsigned DIM>
inline Point<DIM> SphereManipulator<DIM>::getLength() const {
    return torus.get()->L;
}

template<unsigned DIM>
inline double SphereManipulator<DIM>::distanceToCornerS(
    const Sphere<DIM>& sphere) const {
    double dist = numeric_limits<double>::max();
    if constexpr (DIM == 3) {
        Point <DIM> corner = { 0, 0, 0 };  // there is only 1 corner with periodicity!
        dist = abs(
            sphere.radius
            - sqrt(
                torus.get()->distanceCarre(sphere.center,
                    corner)));
    }
    return dist;
}

template<unsigned DIM>
inline bool SphereManipulator<DIM>::cornerWithinSphere(
    const Sphere<DIM>& sphere, double margin) const {
    if constexpr (DIM == 3) {
        ;  // there is only 1 corner with periodicity!
        return (sphere.radius + margin
            >= sqrt(
                torus.get()->distanceCarre(sphere.center,
                    Point <DIM> {0, 0, 0})));
    } else {
        return false;
    }
}
}  // namespace sac_de_billes

//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef SPHERE_IXX_
#define SPHERE_IXX_

#include "../AmbiantSpace.hxx"

namespace sac_de_billes {

template<unsigned short DIM>
inline double Sphere<DIM>::volume() const {
    return sphereTools::volumeSphere<DIM>(radius);
}

template<unsigned short DIM>
inline bool Sphere<DIM>::isInside(const Point<DIM>& x) const {
    return geomTools::normeCarre<DIM>(x - center) < radius * radius;
}

template<unsigned short DIM>
bool areIntersected(const Sphere<DIM>& sph1, const Sphere<DIM>& sph2) {
    return geomTools::distanceCarre<DIM>(sph1.center, sph2.center)
        < (sph1.radius + sph2.radius) * (sph1.radius + sph2.radius);
}

template<unsigned short DIM>
inline double sphereTools::volumeSphere(const double& R) {
    if constexpr (DIM == 1) {
        return 2 * R;
    } else if constexpr (DIM == 2) {
        return m_PI * R * R;
    } else if constexpr (DIM == 3) {
        return 4. / 3. * m_PI * R * R * R;
    } else {
        throw invalid_argument(
            "BigShape<DIM>::volumeBoule: La dimension doit etre 1, 2 ou 3");
    }
}

template<unsigned short DIM>
inline vector<SimpleSphereFormat> sphereTools::vector_fromSphere2Array(
    const vector<Sphere<DIM>>& vec_sphere) {
    vector<SimpleSphereFormat> vec_array{ };
    for (const auto& sphere : vec_sphere) {
        vec_array.push_back(fromSphere2Array(sphere));
    }
    return vec_array;
}

template<unsigned short DIM>
inline vector<Sphere<DIM>> sphereTools::vector_fromArray2Sphere(
    const vector<SimpleSphereFormat>& vec_array) {
    vector<Sphere<DIM>> vec_sphere{ };
    for (const auto& array_ : vec_array) {
        vec_sphere.push_back(fromArray2Sphere<DIM>(array_));
    }
    return vec_sphere;

}

template<unsigned short DIM>
inline SimpleSphereFormat sphereTools::fromSphere2Array(
    const Sphere<DIM>& sphere) {
    SimpleSphereFormat result{ };
    for (size_t i = 0; i < DIM; i++) {
        result[i] = sphere.center[i];
    }
    result[3] = sphere.radius;
    result[4] = sphere.phase;
    return result;
}

template<unsigned short DIM>
inline Sphere<DIM> sphereTools::fromArray2Sphere(
    const SimpleSphereFormat& array_) {
    Point<DIM> center(create_array<DIM>(0.));
    for (size_t i = 0; i < DIM; i++) {
        center[i] = array_[i];
    }
    double radius = array_[3];
    PhaseType phase = array_[4];
    return Sphere<DIM>(center, radius, phase);
}

template<unsigned short DIM>
double sphereTools::volInter(const double& R, const double& r,
    const double& d) {
    if (d >= R + r) {
        return 0;
    }
    if (d <= 1e-7 * (R + r)) {
        return sphereTools::volumeSphere<DIM>(min(r, R));
    }
    if constexpr (DIM == 1) {
        return R + r - d;
    } else if constexpr (DIM == 2) {
        double theta1 = (d * d + r * r - R * R) / (2. * d * r);
        double theta2 = (d * d + R * R - r * r) / (2. * d * R);
        double d_tilde = sqrt(
            (-d + r + R) * (d - r + R) * (d + r - R) * (d + r + R));
        return r * r * acos(theta1) + R * R * acos(theta2) - 0.5 * d_tilde;
    } else if constexpr (DIM == 3) {
        return m_PI * auxi_function::puissance<2>(R + r - d)
            * (d * (d + 2 * r + 2 * R)
                - 3 * auxi_function::puissance<2>(r - R)) / (12. * d);
    } else {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
}

template<unsigned short DIM>
void sphereTools::sort(vector<Sphere<DIM>>& theSpheres) {
    std::sort(theSpheres.begin(), theSpheres.end(),
        [](const auto& s1, const auto& s2) {
            if (s1.phase != s2.phase) {
                return (s1.phase < s2.phase);
            } else {
                return (s1.radius < s2.radius);
            }
        });
}

template<unsigned short DIM>
double sphereTools::volInter(const Sphere<DIM>& s1, const Sphere<DIM>& s2,
    const AmbiantSpace::BigShape<DIM>* const bigShape) {
    // see https://mathworld.wolfram.com/Sphere-SphereIntersection.html
    // see https://mathworld.wolfram.com/Circle-CircleIntersection.html
    double d = bigShape->norme(bigShape->geomVector(s1.center, s2.center));
    return volInter<DIM>(s1.radius, s2.radius, d);
}

template<unsigned short DIM>
bool sphereTools::fromLine(istream& fileStream, size_t phase, Sphere<DIM>& sphere) {
    string line = "";
    if (std::getline(fileStream, line)) {
        std::stringstream ss;
        ss.str(line);
        vector<double> data{};
        double oneData;
        while (ss >> oneData) {
            data.push_back(oneData);
        }
        if (data.size() == 0) { // there should be a void line
            return fromLine(fileStream, phase, sphere);
        }
        for (size_t i = 0; i < DIM; i++) {
            sphere.center[i] = data[i];
        }
        sphere.radius = data[data.size() - 1];
        sphere.phase = phase;
        return true;
    } else {
        return false;
    }
}

} // namespace sac_de_billes

#endif /* SPHERE_IXX_ */

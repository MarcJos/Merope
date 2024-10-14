//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace sac_de_billes {

template<unsigned short DIM, typename T>
inline bool Cylinder<DIM, T>::isInside(const Point<DIM>& point) const {
    Point<DIM> proj_point = geomTools::projection<DIM>(this->axis, point);
    if (geomTools::prodScal<DIM>(proj_point - this->axis[0], this->axis[1] - this->axis[0]) < 0 or
        geomTools::prodScal<DIM>(proj_point - this->axis[1], this->axis[0] - this->axis[1]) < 0
        ) {
        return false;
    }
    if (geomTools::normeCarre<DIM>(point - proj_point) > this->radius * this->radius) {
        return false;
    }
    return true;
}

template<unsigned short DIM, typename T>
inline void Cylinder<DIM, T>::print(std::ostream& f) const {
    f << "Cylinder [";
    f << "Axis : ";
    auxi_function::writeVectorToString(this->axis[0], f); f << " , ";
    auxi_function::writeVectorToString(this->axis[1], f); f << " ;";
    f << "Radius : " << this->radius << ";";
    f << "Phase  :" << this->phase << "] ";
}

}  // namespace sac_de_billes



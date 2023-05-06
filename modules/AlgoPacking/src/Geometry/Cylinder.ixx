//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef ALGOPACKING_SRC_GEOMETRY_CYLINDER_IXX_
#define ALGOPACKING_SRC_GEOMETRY_CYLINDER_IXX_

namespace sac_de_billes {

template<unsigned short DIM>
inline bool Cylinder<DIM>::isInside(const Point<DIM>& point) const {
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

template<unsigned short DIM>
inline void Cylinder<DIM>::print(std::ostream& f) const {
    f << "Cylinder [";
    f << "Axis : ";
    auxi_function::writeVectorToString(this->axis[0], f); f << " , ";
    auxi_function::writeVectorToString(this->axis[1], f); f << " ;";
    f << "Radius : " << this->radius << "] ";
}

} // namespace sac_de_billes

#endif /* ALGOPACKING_SRC_GEOMETRY_CYLINDER_IXX_ */

//! Copyright : see license.txt
//!
//! \brief
#pragma once


#include "../Grid/ConvertGrix.hxx"


namespace merope {
namespace vox {
namespace voxellizer {

template<unsigned short DIM>
template<class COMPOSITE_CELL>
const CartesianGrid<DIM, COMPOSITE_CELL>& GridRepresentation<DIM>::get() const {
    try {
        return std::get<CartesianGrid<DIM, COMPOSITE_CELL>>(internal_field);
    }
    catch (const std::bad_variant_access& ex) {
        cerr << ex.what() << endl;
        raise_error_internal_type_message(__PRETTY_FUNCTION__);
    }
    // for type correctness
    return std::get<CartesianGrid<DIM, COMPOSITE_CELL>>(internal_field);
}

template<unsigned short DIM>
template<class ERROR_MESSAGE_TYPE>
void GridRepresentation<DIM>::raise_error_internal_type_message(ERROR_MESSAGE_TYPE error_message) const {
    cerr << error_message << endl << endl;
    cerr << "Type is : "; print_type(cerr); cerr << endl << endl;
    throw runtime_error("GridRepresentation is not in a correct state to return the desired field");
}

template<unsigned short DIM, class INTERNAL_FIELD, class COMPOSITE>
void local_apply_coefficients(const CartesianGrid<DIM, COMPOSITE>& field,
    INTERNAL_FIELD& internal_field, const vector<double>& coefficients) {
    if constexpr (composite::is_composite<COMPOSITE>) {
        if constexpr (is_same_v<composite::Basic_Phase_Type<COMPOSITE>, PhaseType>) {
            internal_field = vox::voxellizer::apply_coefficients<DIM>(field, coefficients);
        }
    }
}

template<unsigned short DIM>
void GridRepresentation<DIM>::apply_coefficients(vector<double> coefficients) {
    auto transform = [this, &coefficients](const auto& field) {
        local_apply_coefficients(field, this->internal_field, coefficients);
        };
    try_on_all(transform);
}

template<unsigned short DIM>
void GridRepresentation<DIM>::apply_homogRule(homogenization::Rule homogRule) {
    auto field_iso = this->get<composite::Iso<double>>();
    if (homogRule == homogenization::Rule::Reuss) {
        internal_field = vox::voxellizer::applyHomogRule_T<homogenization::Rule::Reuss>(field_iso);
        return;
    }
    if (homogRule == homogenization::Rule::Voigt) {
        internal_field = vox::voxellizer::applyHomogRule_T<homogenization::Rule::Voigt>(field_iso);
        return;
    }
    if (homogRule == homogenization::Rule::Smallest) {
        internal_field = vox::voxellizer::applyHomogRule_T<homogenization::Rule::Smallest>(field_iso);
        return;
    }
    if (homogRule == homogenization::Rule::Largest) {
        internal_field = vox::voxellizer::applyHomogRule_T<homogenization::Rule::Largest>(field_iso);
        return;
    }
    Merope_error_impossible();
}

template<unsigned short DIM>
void GridRepresentation<DIM>::apply_homogRule(homogenization::Rule homogRule, vector<double> coefficients) {
    auto field_iso = this->get<composite::Iso<PhaseType>>();
    if (homogRule == homogenization::Rule::Reuss) {
        internal_field = vox::voxellizer::applyHomogRule_T<homogenization::Rule::Reuss>(field_iso, coefficients);
        return;
    }
    if (homogRule == homogenization::Rule::Voigt) {
        internal_field = vox::voxellizer::applyHomogRule_T<homogenization::Rule::Voigt>(field_iso, coefficients);
        return;
    }
    if (homogRule == homogenization::Rule::Smallest) {
        internal_field = vox::voxellizer::applyHomogRule_T<homogenization::Rule::Smallest>(field_iso, coefficients);
        return;
    }
    if (homogRule == homogenization::Rule::Largest) {
        internal_field = vox::voxellizer::applyHomogRule_T<homogenization::Rule::Largest>(field_iso, coefficients);
        return;
    }
    Merope_error_impossible();
}

template<unsigned short DIM>
template<class INPUT, class OUTPUT>
void GridRepresentation<DIM>::apply_texture(const std::function<OUTPUT(Point<DIM>, INPUT)>& texturer) {
    auto apply_texture_single_grid = [&](const auto& grid) {
        using COMPOSITE = remove_cvref_t<decltype(grid)>::value_type;
        if constexpr (composite::is_composite<COMPOSITE>) {
            if constexpr (is_same_v<composite::Basic_Phase_Type<COMPOSITE>, INPUT>) {
                internal_field = voxellizer::apply_texture<DIM>(grid, texturer);
            }
        }
        };
    try_on_all(apply_texture_single_grid);
}

template<unsigned short DIM>
template<class FUNCTION>
void GridRepresentation<DIM>::try_on_all(FUNCTION function) const {
    try_on_single<vox::composite::Pure<PhaseType>>(function);
    try_on_single<vox::composite::Pure<double>>(function);
    try_on_single<vox::composite::Iso<PhaseType>>(function);
    try_on_single<vox::composite::Iso<double>>(function);
    try_on_single<vox::composite::AnIso<DIM, PhaseType>>(function);
    try_on_single<vox::composite::AnIso<DIM, double>>(function);
    try_on_single<vox::composite::PolyGeom<DIM, PhaseType>>(function);
    try_on_single<vox::composite::PolyGeom<DIM, double>>(function);

    try_on_single<vox::composite::stl_format_Iso<PhaseType>>(function);
    try_on_single<vox::composite::stl_format_Iso<double>>(function);
    try_on_single<vox::composite::stl_format_AnIso<DIM, PhaseType>>(function);
    try_on_single<vox::composite::stl_format_AnIso<DIM, double>>(function);
}

template<unsigned short DIM>
template<class FUNCTION>
void GridRepresentation<DIM>::try_on_all(FUNCTION function) {
    try_on_single<vox::composite::Pure<PhaseType>>(function);
    try_on_single<vox::composite::Pure<double>>(function);
    try_on_single<vox::composite::Iso<PhaseType>>(function);
    try_on_single<vox::composite::Iso<double>>(function);
    try_on_single<vox::composite::AnIso<DIM, PhaseType>>(function);
    try_on_single<vox::composite::AnIso<DIM, double>>(function);
    try_on_single<vox::composite::PolyGeom<DIM, PhaseType>>(function);
    try_on_single<vox::composite::PolyGeom<DIM, double>>(function);

    try_on_single<vox::composite::stl_format_Iso<PhaseType>>(function);
    try_on_single<vox::composite::stl_format_Iso<double>>(function);
    try_on_single<vox::composite::stl_format_AnIso<DIM, PhaseType>>(function);
    try_on_single<vox::composite::stl_format_AnIso<DIM, double>>(function);
}

template<unsigned short DIM>
template<VoxelRule VOXEL_RULE, bool Assume_no_Intersection, class PHASE_TYPE, class STRUCTURE>
void GridRepresentation<DIM>::construct(const STRUCTURE& structure, GridParameters<DIM> preSubGrid, VoxelPolicy<VOXEL_RULE, Assume_no_Intersection, PHASE_TYPE> voxelPolicy) {
    internal_field = voxellizer::transformStructIntoGrid<DIM>(structure, preSubGrid, voxelPolicy);
}

template<unsigned short DIM>
template<class COMPOSITE_TYPE, class FUNCTION>
void GridRepresentation<DIM>::try_on_single(FUNCTION function) const {
    auto pointer = get_if<COMPOSITE_TYPE>();
    if (pointer) {
        function(*pointer);
    }
}

template<unsigned short DIM>
template<class COMPOSITE_TYPE, class FUNCTION>
void GridRepresentation<DIM>::try_on_single(FUNCTION function) {
    auto pointer = get_if<COMPOSITE_TYPE>();
    if (pointer) {
        function(*pointer);
    }
}

template<unsigned short DIM>
void GridRepresentation<DIM>::print_type(std::ostream& f) const {
    f << to_string_type();
}

template<unsigned short DIM>
string GridRepresentation<DIM>::to_string_type() const {
    string res = "";
    auto string_loc = [&res](const auto& elem) {
        res = name_type_of(elem);
        };
    try_on_all(string_loc);
    return res;
}

template<unsigned short DIM>
void GridRepresentation<DIM>::removeUnusedPhase() {
    auto pointer = get_if<PhaseType>();
    if (pointer) {
        vector<PhaseType> coefficients = {};
        vox::convertGrid::removeUnusedPhase(*pointer, coefficients);
    } else {
        raise_error_internal_type_message(__PRETTY_FUNCTION__);
    }
}

template<unsigned short DIM, class INTERNAL_FIELD, class COMPOSITE>
void local_convert_to_stl_format(const CartesianGrid<DIM, COMPOSITE>& field,
    INTERNAL_FIELD& internal_field) {
    if constexpr (vox::composite::is_Pure<COMPOSITE>
        or vox::composite::is_Iso<COMPOSITE>
        or vox::composite::is_AnIso<COMPOSITE>) {
        internal_field = vox::convertGrid::convert_to_stl_format<DIM>(field);
    }
}

template<unsigned short DIM>
void GridRepresentation<DIM>::convert_to_stl_format() {
    auto conversion = [this](const auto& field) {
        local_convert_to_stl_format(field, this->internal_field);
        };
    try_on_all(conversion);
}

template<unsigned short DIM, class COMPOSITE_OUT, class COMPOSITE_IN,
    class COMPOSITE_FIELD, class INTERNAL_FIELD>
void local_convert_to(const CartesianGrid<DIM, COMPOSITE_FIELD>& field,
    INTERNAL_FIELD& internal_field) {
    if constexpr (is_same_v<COMPOSITE_IN, COMPOSITE_FIELD>) {
        internal_field = vox::convertGrid::convert_to<DIM, COMPOSITE_OUT, COMPOSITE_IN>(field);
    }

}

template<unsigned short DIM>
template<class COMPOSITE_OUT, class COMPOSITE_IN>
void GridRepresentation<DIM>::convert_to() {
    auto conversion = [this](const auto& field) {
        local_convert_to<DIM, COMPOSITE_OUT, COMPOSITE_IN>(field, this->internal_field);
        };
    try_on_all(conversion);
}

template<unsigned short DIM>
void GridRepresentation<DIM>::convert_to_Iso_format() {
    convert_to<vox::composite::Iso<PhaseType>, vox::composite::AnIso<DIM, PhaseType>>();
    convert_to<vox::composite::Iso<double>, vox::composite::AnIso<DIM, double>>();
    //
    convert_to<vox::composite::Iso<PhaseType>, vox::composite::PolyGeom<DIM, PhaseType>>();
    convert_to<vox::composite::Iso<double>, vox::composite::PolyGeom<DIM, double>>();
    //
    convert_to<vox::composite::Iso<PhaseType>, vox::composite::Pure<PhaseType>>();
    convert_to<vox::composite::Iso<double>, vox::composite::Pure<double>>();
}

template<unsigned short DIM>
void GridRepresentation<DIM>::convert_to_AnIso_format() {
    convert_to<vox::composite::AnIso<DIM, PhaseType>, vox::composite::PolyGeom<DIM, PhaseType>>();
    convert_to<vox::composite::AnIso<DIM, double>, vox::composite::PolyGeom<DIM, double>>();
}

template<unsigned short DIM, class COMPOSITE>
string name_type_of(CartesianGrid<DIM, COMPOSITE>) {
    if constexpr (is_same_v<COMPOSITE, vox::composite::Pure<PhaseType>>) {
        return "vox::composite::Pure<PhaseType>";
    } else if constexpr (is_same_v<COMPOSITE, vox::composite::Pure<double>>) {
        return "vox::composite::Pure<double>";
    } else if constexpr (is_same_v<COMPOSITE, vox::composite::Iso<PhaseType>>) {
        return "vox::composite::Iso<PhaseType>";
    } else if constexpr (is_same_v<COMPOSITE, vox::composite::Iso<double>>) {
        return "vox::composite::Iso<double>";
    } else if constexpr (is_same_v<COMPOSITE, vox::composite::AnIso<DIM, PhaseType>>) {
        return "vox::composite::AnIso<DIM, PhaseType>";
    } else if constexpr (is_same_v<COMPOSITE, vox::composite::AnIso<DIM, double>>) {
        return "vox::composite::AnIso<DIM, double>";
    } else if constexpr (is_same_v<COMPOSITE, vox::composite::PolyGeom<DIM, PhaseType>>) {
        return "vox::composite::PolyGeom<DIM, PhaseType>";
    } else if constexpr (is_same_v<COMPOSITE, vox::composite::PolyGeom<DIM, double>>) {
        return "vox::composite::PolyGeom<DIM, double>";
    } else if constexpr (is_same_v<COMPOSITE, vox::composite::stl_format_Iso<PhaseType>>) {
        return "vox::composite::stl_format_Iso<PhaseType>";
    } else if constexpr (is_same_v<COMPOSITE, vox::composite::stl_format_Iso<double>>) {
        return "vox::composite::stl_format_Iso<double>";
    } else if constexpr (is_same_v<COMPOSITE, vox::composite::stl_format_AnIso<DIM, PhaseType>>) {
        return "vox::composite::stl_format_AnIso<DIM, PhaseType>";
    } else if constexpr (is_same_v<COMPOSITE, vox::composite::stl_format_AnIso<DIM, double>>) {
        return "vox::composite::stl_format_AnIso<DIM, double>";
    } else Merope_static_error(COMPOSITE, "Unknown type");
}

}  // namespace  voxellizer
} /* namespace vox */
}  // namespace merope
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
    cerr << error_message << endl;
    print_type(cerr); cerr << endl;
    throw runtime_error("GridRepresentation is not in a correct state to return the deisred field");
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
    cerr << __PRETTY_FUNCTION__ << endl;
    throw runtime_error("Impossible");
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
    cerr << __PRETTY_FUNCTION__ << endl;
    throw runtime_error("Impossible");
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

    try_on_single<vox::composite::stl_format_Iso<PhaseType>>(function);
    try_on_single<vox::composite::stl_format_Iso<double>>(function);
    try_on_single<vox::composite::stl_format_AnIso<DIM, PhaseType>>(function);
    try_on_single<vox::composite::stl_format_AnIso<DIM, double>>(function);
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
    if constexpr (composite::is_composite<COMPOSITE>) {
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

template<unsigned short DIM, class COMPOSITE>
string name_type_of(CartesianGrid<DIM, COMPOSITE>) {
    if constexpr (is_same_v<COMPOSITE, vox::composite::Pure<PhaseType>>) {
        return "vox::composite::Pure<PhaseType>";
    }
    if constexpr (is_same_v<COMPOSITE, vox::composite::Pure<double>>) {
        return "vox::composite::Pure<double>";
    }
    if constexpr (is_same_v<COMPOSITE, vox::composite::Iso<PhaseType>>) {
        return "vox::composite::Iso<PhaseType>";
    }
    if constexpr (is_same_v<COMPOSITE, vox::composite::Iso<double>>) {
        return "vox::composite::Iso<double>";
    }
    if constexpr (is_same_v<COMPOSITE, vox::composite::AnIso<DIM, PhaseType>>) {
        return "vox::composite::AnIso<DIM, PhaseType>";
    }
    if constexpr (is_same_v<COMPOSITE, vox::composite::AnIso<DIM, double>>) {
        return "vox::composite::AnIso<DIM, double>";
    }
    if constexpr (is_same_v<COMPOSITE, vox::composite::stl_format_Iso<PhaseType>>) {
        return "vox::composite::stl_format_Iso<PhaseType>";
    }
    if constexpr (is_same_v<COMPOSITE, vox::composite::stl_format_Iso<double>>) {
        return "vox::composite::stl_format_Iso<double>";
    }
    if constexpr (is_same_v<COMPOSITE, vox::composite::stl_format_AnIso<DIM, PhaseType>>) {
        return "vox::composite::stl_format_AnIso<DIM, PhaseType>";
    }
    if constexpr (is_same_v<COMPOSITE, vox::composite::stl_format_AnIso<DIM, double>>) {
        return "vox::composite::stl_format_AnIso<DIM, double>";
    }
    throw runtime_error("Unknown type");
}

}  // namespace  voxellizer
} /* namespace vox */
}  // namespace merope
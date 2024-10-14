//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../../AlgoPacking/src/StdHeaders.hxx"
#include "../Voxellation/Voxellizer.hxx"

namespace merope {
namespace vox {
namespace voxellizer {

template<unsigned short DIM>
class GridRepresentation {
public:
    //! @brief constructor
    template<class STRUCTURE>
    GridRepresentation(const STRUCTURE& structure, GridParameters<DIM> preSubGrid, vox::VoxelRule voxelRule) : internal_field(nullptr) {
        if (voxelRule == VoxelRule::Center) {
            internal_field = voxellizer::transformStructIntoGrid<DIM, VoxelRule::Center>(structure, preSubGrid);
        }
        if (voxelRule == VoxelRule::Average) {
            internal_field = voxellizer::transformStructIntoGrid<DIM, VoxelRule::Average>(structure, preSubGrid);
        }
        if (voxelRule == VoxelRule::Laminate) {
            internal_field = voxellizer::transformStructIntoGrid<DIM, VoxelRule::Laminate>(structure, preSubGrid);
        }
    }

    //! @brief transform the composite<int> grid into a composite<double> grid
    //! according to the rule i -> coefficients[i]
    void apply_coefficients(vector<double> coefficients);
    //! @brief apply the homogenization rule to the grid
    void apply_homogRule(homogenization::Rule homogRule);
    //! @brief equivalent to
    //! this->apply_coefficients(coefficients);
    //! this->apply_homogRule(homogRule);
    void apply_homogRule(homogenization::Rule homogRule, vector<double> coefficients);
    //! @brief apply a texturer to the structure
    //! @param texturer : function of (x = spatial position, phase)
    template<class INPUT, class OUTPUT>
    void apply_texture(const std::function<OUTPUT(Point<DIM>, INPUT)>& texturer);

    //! @brief : convert to the adequate stl format
    void convert_to_stl_format();
    //! @brief : change the phase id so that all phases are inside [0, N] with each phase present at least in one voxel.
    //! preserves the phase order
    void removeUnusedPhase();
    //! @return the internal field
    //! throws error if type incoherent with internal state
    //! @tparam COMPOSITE_CELL : type of cell
    template<class COMPOSITE_CELL>
    const CartesianGrid<DIM, COMPOSITE_CELL>& get() const;
    //! @return a pointer to the internal field
    //! return null pointer if type incoherent with internal state
    //! @tparam COMPOSITE_CELL : type of cell
    template<class COMPOSITE_CELL>
    const CartesianGrid<DIM, COMPOSITE_CELL>* get_if() const { return std::get_if<CartesianGrid<DIM, COMPOSITE_CELL>>(&internal_field); }
    template<class COMPOSITE_CELL>
    CartesianGrid<DIM, COMPOSITE_CELL>* get_if() { return std::get_if<CartesianGrid<DIM, COMPOSITE_CELL>>(&internal_field); }
    //! @return : a string describing the internal type
    string to_string_type() const;
    //! @brief : print the internal type on f
    void print_type(std::ostream& f) const;
    //! @brief : throws an error, with type of the internal state displayed
    template<class ERROR_MESSAGE_TYPE>
    void raise_error_internal_type_message(ERROR_MESSAGE_TYPE error_message) const;

    //! @brief given a lambda function, try it on all variant types
    //! @param function : given lambda function
    template<class FUNCTION>
    void try_on_all(FUNCTION function);
    template<class FUNCTION>
    void try_on_all(FUNCTION function) const;

private:
    template<class COMPOSITE_TYPE, class FUNCTION>
    void try_on_single(FUNCTION function);
    template<class COMPOSITE_TYPE, class FUNCTION>
    void try_on_single(FUNCTION function) const;

    std::variant<
        CartesianGrid<DIM, vox::composite::Pure<PhaseType>>,
        CartesianGrid<DIM, vox::composite::Pure<double>>,
        CartesianGrid<DIM, vox::composite::Iso<PhaseType>>,
        CartesianGrid<DIM, vox::composite::Iso<double>>,
        CartesianGrid<DIM, vox::composite::AnIso<DIM, PhaseType>>,
        CartesianGrid<DIM, vox::composite::AnIso<DIM, double>>,

        //CartesianGrid<DIM, vox::composite::stl_format_Pure<PhaseType>>, SAME AS ::Pure
        // CartesianGrid<DIM, vox::composite::stl_format_Pure<double>>, SAME AS ::Pure
        CartesianGrid<DIM, vox::composite::stl_format_Iso<PhaseType>>,
        CartesianGrid<DIM, vox::composite::stl_format_Iso<double>>,
        CartesianGrid<DIM, vox::composite::stl_format_AnIso<DIM, PhaseType>>,
        CartesianGrid<DIM, vox::composite::stl_format_AnIso<DIM, double>>,
        void*
    > internal_field;
};

template<unsigned short DIM, class COMPOSITE>
string name_type_of(CartesianGrid<DIM, COMPOSITE>);

}  // namespace  voxellizer
} /* namespace vox */
}  // namespace merope

#include "../Voxellation/DynamicVoxellizer.ixx"
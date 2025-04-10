//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../../Geometry/include/VolumeInCube.hxx"

namespace merope {
namespace vox {
namespace composite {

template<unsigned short DIM, class PHASE_OUT, class COMPOSITE,
    class CONVERSION, class TEST>
Change_Type_Composite<DIM, COMPOSITE, PHASE_OUT> transform_phase(const COMPOSITE& composite,
    const CONVERSION& conversion) {
    if constexpr (is_Pure<COMPOSITE>) {
        return conversion(composite);
    }
    //

    //
    else if constexpr (is_Iso<COMPOSITE>) {
        Change_Type_Composite<DIM, COMPOSITE, PHASE_OUT> result{};
        for (const auto& pf : composite) {
            result.push_back(PhaseFrac<PHASE_OUT>(conversion(pf.phase), pf.fracVol));
        }
        return result;
    }
    //

    //
    else if constexpr (is_AnIso<COMPOSITE>) {
        Change_Type_Composite<DIM, COMPOSITE, PHASE_OUT> result{};
        for (const auto& pf : composite) {
            result.push_back(PhaseFracNormal<DIM, PHASE_OUT>(conversion(pf.phase), pf.fracVol, pf.normal));
        }
        return result;
    }
    //

    //
    else if constexpr (is_PolyGeom<COMPOSITE>) {
        Change_Type_Composite<DIM, COMPOSITE, PHASE_OUT> result{};
        result.halfSpaces = composite.halfSpaces;
        result.phases.resize(composite.phases.size());
        for (size_t i = 0; i < composite.phases.size(); i++) {
            result.phases[i] = conversion(composite.phases[i]);
        }
        if (composite.is_there_matrix()) {
            result.setMatrixPhase(conversion(composite.getMatrixPhase()));
        }
        return result;
    }
    //

    //
    else {
        Merope_static_error(COMPOSITE, "Impossible");
    }
}

template<unsigned short DIM, class COMPOSITE, class TEXTURER,
    class PHASE_OUT, class TEST>
Change_Type_Composite<DIM, COMPOSITE, PHASE_OUT> apply_texture_loc(const COMPOSITE& composite,
    const TEXTURER& texturer, const Point<DIM>& x) {
    auto conversion = [&](const auto& phase) {
        return texturer(x, phase);
        };
    return transform_phase<DIM, PHASE_OUT>(composite, conversion);
}


template<unsigned short DIM, class COMPOSITE_OUT, class COMPOSITE_IN>
COMPOSITE_OUT convert_to(const COMPOSITE_IN& composite_in) {
    return auxi::Helper_convert_to<DIM, COMPOSITE_OUT, COMPOSITE_IN>::eval(composite_in);
}

template<unsigned short DIM, class VOXEL_POLICY, class PHASE_TYPE>
void pure_phase_inserter(OutputFormat<VOXEL_POLICY::voxelRule, DIM, PHASE_TYPE>& vox, PHASE_TYPE purePhase, const VOXEL_POLICY& voxelPolicy) {
    constexpr vox::VoxelRule VOXEL_RULE = VOXEL_POLICY::voxelRule;
    if constexpr (VOXEL_POLICY::Assume_no_Intersection) {
        if constexpr (VOXEL_RULE == VoxelRule::Center) {
            vox = purePhase;
        }
        //
        else if constexpr (VOXEL_RULE == VoxelRule::Average
            or VOXEL_RULE == VoxelRule::Laminate
            or VOXEL_RULE == VoxelRule::PolyGeom) {
            vox.set_pure(purePhase);
        }
        //
        else {
            Merope_static_error(PHASE_TYPE, "Impossible");
        }
    }
    // case not(VOXEL_POLICY::Assume_no_Intersection)
    else {
        auto newVoxel = OutputFormat<VOXEL_RULE, DIM, PHASE_TYPE>(purePhase);
        voxel_filler<DIM>(vox, newVoxel, voxelPolicy);
    }
}

template<unsigned short DIM, class VOXEL_TYPE, class VOXEL_POLICY>
void voxel_filler(VOXEL_TYPE& voxel_to_be_filled, const VOXEL_TYPE& voxel_to_include, const VOXEL_POLICY& voxelPolicy) {
    if constexpr (vox::composite::is_Pure<VOXEL_TYPE>) {
        if constexpr (VOXEL_POLICY::Assume_no_Intersection) {
            voxel_to_be_filled = voxel_to_include;
        } else {
            voxel_to_be_filled = (is_nan_like(voxel_to_be_filled)) ? voxel_to_include : voxelPolicy.getIntersectionRule()(voxel_to_be_filled, voxel_to_include);
        }
    } else if constexpr (vox::composite::is_Iso<VOXEL_TYPE> or vox::composite::is_AnIso<VOXEL_TYPE>) {
        if constexpr (VOXEL_POLICY::Assume_no_Intersection) {
            for (const auto& phfv : voxel_to_include) {
                voxel_to_be_filled.push_back(phfv);
            }
        } else {
            if (voxel_to_be_filled.size() == 0) {
                voxel_to_be_filled = voxel_to_include;
            } else {
                Merope_assert(false, "If intersection is assumed, impossible to decided how to mix 2 voxels with this rule.");
            }
        }
    } else if constexpr (vox::composite::is_PolyGeom<VOXEL_TYPE>) {
        if constexpr (VOXEL_POLICY::Assume_no_Intersection) {
            voxel_to_be_filled.merge_with(voxel_to_include);
        } else {
            voxel_to_be_filled.merge_with(voxel_to_include, voxelPolicy.getIntersectionRule());
        }
    } else {
        Merope_static_error(VOXEL_TYPE, "VOXEL_TYPE incorrect");
    }
}

template<unsigned short DIM, class VOXEL_TYPE, class VOXEL_POLICY>
void postProcess(VOXEL_TYPE& voxel, const VOXEL_POLICY&, const auto& matrixPhaseHolder) {
    if constexpr (VOXEL_POLICY::voxelRule == VoxelRule::Average
        or VOXEL_POLICY::voxelRule == VoxelRule::Laminate
        or VOXEL_POLICY::voxelRule == VoxelRule::PolyGeom) {
        voxel.postProcess(matrixPhaseHolder);
    }
    if constexpr (VOXEL_POLICY::voxelRule == VoxelRule::Center) {
        voxel = (is_nan_like(voxel)) ? matrixPhaseHolder.getMatrixPhase() : voxel;
    }
}




template<class VOXEL_TYPE, class VOXEL_POLICY>
VOXEL_TYPE make_default_voxel(const VOXEL_POLICY& voxelPolicy, const auto& matrixPhaseHolder) {
    if constexpr (VOXEL_POLICY::voxelRule == VoxelRule::Center) {
        if constexpr (VOXEL_POLICY::Assume_no_Intersection) {
            return matrixPhaseHolder.getMatrixPhase();
        } else {
            return make_nan_like<composite::Basic_Phase_Type<VOXEL_TYPE>>();
        }
    }
    //
    else if constexpr (VOXEL_POLICY::voxelRule == VoxelRule::PolyGeom) {
        VOXEL_TYPE voxel{};
        if (matrixPhaseHolder.is_there_matrix()) {
            voxel.setMatrixPhase(matrixPhaseHolder.getMatrixPhase());
        }
        return voxel;
    }
    //
    else {
        return VOXEL_TYPE{};
    }
}

namespace auxi {

template<unsigned short DIM, class COMPOSITE_OUT, class COMPOSITE_IN>
COMPOSITE_OUT Helper_convert_to<DIM, COMPOSITE_OUT, COMPOSITE_IN>::eval(const COMPOSITE_IN& composite_in) {
    using BasicPhaseType = composite::Basic_Phase_Type<COMPOSITE_IN>;

    //
    if constexpr (Helper_convert_to::areTypes<AnIso<DIM, BasicPhaseType>, PolyGeom<DIM, BasicPhaseType>>) {
        AnIso<DIM, BasicPhaseType> result{};
        //
        for (size_t i = 0; i < composite_in.phases.size(); i++) {
            auto p = composite_in.phases[i];
            auto f = geomTools::volume_in_cube(composite_in.halfSpaces[i]);
            Point<DIM> normal = create_array<DIM>(1.);
            if (composite_in.halfSpaces[i].size() > 0) {
                normal = composite_in.halfSpaces[i][0].vec();
            }
            geomTools::renormalize<DIM>(normal);
            result.add(PhaseFracNormal<DIM, BasicPhaseType>(p, f, normal));
        }
        result.postProcess(composite_in);
        //
        return result;
    }
    //

    //
    else if constexpr (Helper_convert_to::areTypes<Iso<BasicPhaseType>, PolyGeom<DIM, BasicPhaseType>>) {
        Iso<BasicPhaseType> result{};
        //
        for (size_t i = 0; i < composite_in.phases.size(); i++) {
            auto p = composite_in.phases[i];
            auto f = geomTools::volume_in_cube(composite_in.halfSpaces[i]);
            result.add(PhaseFrac<BasicPhaseType>(p, f));
        }
        result.postProcess(composite_in);
        //
        return result;
    }
    //

    //
    else if constexpr (Helper_convert_to::areTypes<Iso<BasicPhaseType>, AnIso<DIM, BasicPhaseType>>) {
        Iso<BasicPhaseType> result{};
        for (const auto& pfn : composite_in) {
            result.add(PhaseFrac<BasicPhaseType>(pfn.phase, pfn.fracVol));
        }
        return result;
    }
    //

    //
    else if constexpr (composite::is_Pure<COMPOSITE_IN>) {
        return COMPOSITE_OUT(composite_in);
    }
    //

    //
    else if constexpr (Helper_convert_to::areTypes<to_stl_format<COMPOSITE_IN>, COMPOSITE_IN>) {
        return auxi::convert_to_stl_format(composite_in);
    }
    //

    //
    else {
        Merope_static_error(COMPOSITE_IN, "Impossible");
    }
}


template<class COMPOSITE, unsigned short DIM>
to_stl_format<COMPOSITE> convert_to_stl_format(const COMPOSITE& composite_) {
    static_assert(is_Pure<COMPOSITE> or is_Iso<COMPOSITE> or is_AnIso<COMPOSITE>);
    using PHASE_TYPE = Basic_Phase_Type<COMPOSITE>;
    if constexpr (is_Pure<COMPOSITE>) {
        return composite_;
    }
    if constexpr (is_Iso<COMPOSITE>) {
        composite::stl_format_Iso<PHASE_TYPE> result = {};
        for (const auto& phaseFrac : composite_) {
            result.push_back(make_tuple(phaseFrac.phase, phaseFrac.fracVol));
        }
        return result;
    }
    if constexpr (is_AnIso<COMPOSITE>) {
        composite::stl_format_Iso<PHASE_TYPE> result_phaseFrac = {};
        for (const auto& phaseFrac : composite_) {
            result_phaseFrac.push_back(make_tuple(phaseFrac.phase, phaseFrac.fracVol));
        }
        Point<DIM> normal = vox::composite::compute_global_normal(composite_);
        composite::stl_format_AnIso<DIM, PHASE_TYPE> result = { result_phaseFrac, normal };
        return result;
    }
}
}  // namespace  auxi


}  // namespace  composite
}  // namespace vox
}  // namespace merope



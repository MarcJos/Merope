//! Copyright : see license.txt
//!
//! \briefClass for voxellization of inclusions
//
#pragma once


#include "../../GenericMerope/StdHeaders.hxx"
#include "../../GenericMerope/MeropeNamespace.hxx"

#include "../../Geometry/include/GeomTools.hxx"
#include "../../Grid/CartesianGrid.hxx"
#include "../../Grid/ConvertGrix.hxx"
#include "../../Grid/Grid_VER.hxx"
#include "../../Grid/GridManipulations.hxx"
#include "../../Grid/GridTypes.hxx"
#include "../../Grid/VtkHistogramm.hxx"
#include "../../Voxellation/VoxRecurStructure.hxx"
#include "../../VTKinout/VTKStream.hxx"
#include "../../Voxellation/Voxellizer.hxx"


namespace merope {
namespace vox {

//! data structure for Voxellation
template<unsigned short DIM>
class Voxellation_data : public GridParameters<DIM> {
public:
    //! main constructor
    Voxellation_data(const Structure<DIM>& structure);
    Voxellation_data(const FieldStructure<DIM>& fieldStructure);
    //! destructor
    virtual ~Voxellation_data() {}
    //! forbid copy
    Voxellation_data(const Voxellation_data& other) = delete;
    Voxellation_data& operator=(const Voxellation_data& other) = delete;
    Voxellation_data(Voxellation_data&& other) = delete;
    Voxellation_data& operator=(Voxellation_data&& other) = delete;

    // setter
    //! sets the material coefficients
    //! \param pureCoeffs:  correspond to pure phase 0 (=matrix), 1->N spheres
    //! \warning : also sets the coefficients (which may be modifier afterwards due to the homogenization rule)
    void setPureCoeffs(vector<double> pureCoeffs_);
    //! chooses the homogenization rule
    void setHomogRule(homogenization::Rule rule) { homogRule = rule; }
    //! \param nbNodes : size of the large grid
    //! \param nMin : index of the lower left corner of the slice
    //! \param nMax : index of the upper right corner of the slice
    void setDiscretization(array<size_t, DIM> nbNodes_, array<size_t, DIM> nMin_, array<size_t, DIM> nMax_);
    void setDiscretization(array<size_t, DIM> nbNodes_) { this->setDiscretization(nbNodes_, create_array<DIM, size_t>(0), nbNodes_); }
    //! set the rule for choosing which phase is in the voxel
    void setVoxelRule(VoxelRule voxelRule_) { voxelRule = voxelRule_; }

protected:
    //! determines which rule to fill the voxel
    VoxelRule voxelRule;
    //! homogenization rule for building voxellization
    homogenization::Rule homogRule;
    //! inner inclusions ref
    unique_ptr<const Structure<DIM>> structure;
    //! inner field ref
    unique_ptr<const FieldStructure<DIM>> fieldStructure;
    //! stores the coefficient of pure materials
    vector<double> pureCoeffs;
    //! says whether pureCoeffs is set
    bool pureCoeffsIsSet;
    //! auxiliary constructor
    Voxellation_data(array<double, DIM> L);
    //! warning : does not seem well-placed!
    //! checks if the pureCoeffs are set. If not, set them.
    void checkPureCoeffs();
    //! @brief check whether the voxel rule is consistent
    template<VoxelRule VOXEL_RULE>
    void check_voxel_rule(const char* pretty_function) const;
    //! @brief check whether the homogenization rule is consistent
    template<homogenization::Rule HOMOG_RULE>
    void check_homog_rule(const char* pretty_function) const;
    //! @brief check whether the structure is consistent
    template<class STRUCTURE = bool>
    void check_structure(const char* pretty_function) const;
};

//! class for turning geometrical structures to voxellation
template<unsigned short DIM>
class Voxellation final : public Voxellation_data<DIM>, private auxi_function::Deprecated {
public:
    //! main constructor
    Voxellation(const Structure<DIM>& structure_) :Voxellation_data<DIM>(structure_) {}
    Voxellation(const MultiInclusions<DIM>& mIncl) :Voxellation_data<DIM>(Structure<DIM>(mIncl)) {}
    Voxellation(const FieldStructure<DIM>& fieldStructure_) : Voxellation_data<DIM>(fieldStructure_) {}
    //! destructor
    virtual ~Voxellation() {}


    //! \param nbNodes : size of the large grid
    //! \param nMin : index of the lower left corner of the slice
    //! \param nMax : index of the upper right corner of the slice
    void proceed(array<size_t, DIM> nbNodes_, array<size_t, DIM> nMin_, array<size_t, DIM> nMax_) {
        this->setDiscretization(nbNodes_, nMin_, nMax_);
    }
    //! \param nbNodes : sizes of the voxellation grid
    void proceed(array<size_t, DIM> nbNodes_) { this->proceed(nbNodes_, create_array<DIM, size_t>(0), nbNodes_); }


    //! computes a field
    GridField<DIM> getField();
    //! compute a phaseGrid
    vox::CartesianGrid<DIM, PhaseType> computePurePhaseGrid(vector<double>& coefficients_);
    //! for outputs in python Format
    //! \return a vector of vectors of (phase, volume fraction)
    vector<composite::stl_format_Iso<PhaseType>> computePhaseGrid();
    //! \return a vector of vectors of composite voxels with volume fraction and normal
    vector<composite::stl_format_AnIso<DIM, PhaseType>> computeCompositeGrid();
    //! @tparam VOXEL_RULE 
    //! @param nbNodes_ : size of the grid
    //! @return the cartesian Grid containing the composite voxels
    template<vox::VoxelRule VOXEL_RULE>
    CartesianGrid<DIM, composite::OutputFormat<VOXEL_RULE, DIM, PhaseType>> computeCartesianGrid();

    // output
    //! print the output with composite voxels
    //! \param fileVTK : output file for the field table
    void printFieldFile(string fileVTK);
    //! print the output with composite voxels
    //! \param fileVTK : output file for the phase table (\see Grid)
    //! \param fileCoeff : output storing the coefficient associated with each phase (compatible with AMITEX format)
    void printFile(string fileVTK, string fileCoeff);
    //! @return print the inner phases of the grid
    string print();

private:
    //! sets all the voxels in the grid from phaseFracVol, with homogenization rule
    template<class VOXEL_TYPE>
    vox::GridField<DIM> applyHomogRule(const vox::CartesianGrid<DIM, VOXEL_TYPE>& phaseFracVol);
    //! sets all the voxels in the grid from phaseFracVol, with homogenization rule
    template<homogenization::Rule, class VOXEL_TYPE>
    vox::GridField<DIM> applyHomogRule_T(const vox::CartesianGrid<DIM, VOXEL_TYPE>& phaseFracVol);
    //! \see computeField
    template<class STRUCTURE>
    GridField<DIM> computeField(const STRUCTURE& structure);

    //! forbid copy
    Voxellation(const Voxellation& other) = delete;
    Voxellation& operator=(const Voxellation& other) = delete;
    Voxellation(Voxellation&& other) = delete;
    Voxellation& operator=(Voxellation&& other) = delete;
};

} /* namespace vox */
}  // namespace merope


#include "../../Voxellation/Deprecated/Voxellation.ixx"


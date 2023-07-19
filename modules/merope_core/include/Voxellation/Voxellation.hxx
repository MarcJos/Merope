//! Copyright : see license.txt
//!
//! \brief Class for voxellization of inclusions
//
#ifndef VOXELLATION_HXX_
#define VOXELLATION_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Geometry/GeomTools.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "../Grid/ConvertGrix.hxx"
#include "../Grid/Grid_VER.hxx"
#include "../Grid/GridManipulations.hxx"
#include "../Grid/GridTypes.hxx"
#include "../Grid/VtkHistogramm.hxx"
#include "../Voxellation/VoxRecurStructure.hxx"
#include "../Voxellation/PreVoxellation.hxx"
#include "../VTKinout/VTKStream.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

//! data structure for Voxellation
template<unsigned short DIM>
class Voxellation_data : public  PreVoxellation<DIM> {
public:
    //! main constructor
    Voxellation_data(const Structure<DIM>&);
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

    // output
    //! print the output with composite voxels
    //! \param fileVTK : output file for the phase table (\see Grid)
    //! \param fileCoeff : output storing the coefficient associated with each phase (compatible with AMITEX format)
    void printFile(string fileVTK, string fileCoeff) const;

protected:
    //! homogenization rule for building voxellization
    homogenization::Rule homogRule;
    //! inner inclusions ref
    unique_ptr<const Structure<DIM>> structure;
    //! inner field ref
    unique_ptr<const FieldStructure<DIM>> fieldStructure;
    //! stores the coefficient for composite voxels
    vector<double> coefficients;
    //! stores the coefficient of pure materials
    vector<double> pureCoeffs;
    //! says whether pureCoeffs is set
    bool pureCoeffsIsSet;
    //! auxiliary constructor
    Voxellation_data(array<double, DIM> L);
    // verifications
        //! checks if the structure exists. If not, throws an error
    bool checkStructure(const char* pretty_function);
    //! warning : does not seem well-placed!
    //! checks if the pureCoeffs are set. If not, set them.
    void checkPureCoeffs();
private:
    // output
        //! print the coefficients
    void printCoeffs(string fileCoeff) const;
};

//! class for turning geometrical structures to voxellation
template<unsigned short DIM>
class Voxellation final : public Voxellation_data<DIM> {
public:
    //! main constructor
    Voxellation(const Structure<DIM>& structure_) :Voxellation_data<DIM>(structure_) {}
    Voxellation(const MultiInclusions<DIM>& mIncl) :Voxellation_data<DIM>(Structure<DIM>(mIncl)) {}
    Voxellation(const FieldStructure<DIM>& fieldStructure_) : Voxellation_data<DIM>(fieldStructure_) {}
    //! destructor
    virtual ~Voxellation() {}
    //! forbid copy
    Voxellation(const Voxellation& other) = delete;
    Voxellation& operator=(const Voxellation& other) = delete;
    Voxellation(Voxellation&& other) = delete;
    Voxellation& operator=(Voxellation&& other) = delete;

    //! computes a slice of the grid
    //! \param nbNodes : size of the large grid
    //! \param nMin : index of the lower left corner of the slice
    //! \param nMax : index of the upper right corner of the slice
    void proceed(array<size_t, DIM> nbNodes, array<size_t, DIM> nMin, array<size_t, DIM> nMax);
    //! computes the grid
    //! \param nbNodes : sizes of the voxellation grid
    void proceed(array<size_t, DIM> nbNodes);

    //! builds the array of phase concentrations
    template<VoxelRule VOXEL_RULE>
    vox::CartesianGrid<DIM, OutputFormat<VOXEL_RULE>> getPhaseFracField();
    // output
    //! print the output with composite voxels
    //! \param fileVTK : output file for the phase table (\see Grid)
    //! \param fileCoeff : output storing the coefficient associated with each phase (compatible with AMITEX format)
    void printFieldFile(string fileVTK);
    //! computes a field
    GridField<DIM> getField();
    //! compute a phaseGrid
    vox::CartesianGrid<DIM, vox::VTK_PHASE> getPhaseGrid();
    //! for outputs in python Format
    //! \return a vector of vectors of (phase, volume fraction)
    vector<vector<tuple<vox::VTK_PHASE, double>>> computePhaseGrid(array<size_t, DIM> nbNodes_);

private:
    //! sets all the voxels in the grid from phaseFracVol, with homogenization rule
    vox::GridField<DIM> applyHomogRule(const vox::CartesianGrid<DIM, VoxelPhaseFrac>& phaseFracVol);
    //! sets all the voxels in the grid from phaseFracVol, with homogenization rule
    template<homogenization::Rule>
    vox::GridField<DIM> applyHomogRule_T(const vox::CartesianGrid<DIM, VoxelPhaseFrac>& phaseFracVol);
    //! \see getField
    GridField<DIM> getField_Structure();
    //! \see getField
    GridField<DIM> getField_FieldStructure();
};

} /* namespace vox */
} // namespace merope


#include "../Voxellation/Voxellation.ixx"
#endif /* VOXELLATION_HXX_ */

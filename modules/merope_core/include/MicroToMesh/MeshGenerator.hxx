//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MICROTOMESH_MESHGENERATOR_HXX_
#define MICROTOMESH_MESHGENERATOR_HXX_

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Mesh/MeshStructure.hxx"
#include "../MultiInclusions/MultiInclusions.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {
namespace mesh {
namespace generator {

class MeshGenerator {
public:
    //! constructor
    MeshGenerator() : meshSize{ 0.05 }, meshOrder{ 1 },
        adimensionnalMergeDistance_0{ 1.e-5 },
        adimensionnalMergeDistance_1{ 1.e-5 },
        binaryOutput{ false },
        multiInclusions{ nullptr } {}
    //! getter
    size_t getMeshOrder() const { return meshOrder; }
    //! setter
    void setMeshOrder(size_t meshOrder_) { this->meshOrder = meshOrder_; }
    //! write the gmsh input into a stream
    //! \param f : outputstream
    void write(std::ostream& f) const;
    //! write the gmsh input into a file
    //! \param nameFile : the name of the file
    void write(string nameFile) const;
    //! setter
    void setMultiInclusions(const MultiInclusions<3>& multiInclusions_) { this->multiInclusions = &multiInclusions_; }
    //! getter
    double getMeshSize() const { return meshSize; }
    //! setter
    void setMeshSize(double meshSize_) { this->meshSize = meshSize_; }
    //! setter
    void setAdimMergeDistance0(double adimensionnalMergeDistance0_) { adimensionnalMergeDistance_0 = adimensionnalMergeDistance0_; }
    //! setter
    void setAdimMergeDistance1(double adimensionnalMergeDistance1_) { adimensionnalMergeDistance_1 = adimensionnalMergeDistance1_; }
    //! setter
    void setBinaryOutput(bool binaryOutput_) { binaryOutput = binaryOutput_; }

private:
    //! Options
    //! size (=step) of the mesh elements
    double meshSize;
    //! order of the mesh elements
    size_t meshOrder;
    //! triggers the merge distance
    double adimensionnalMergeDistance_0;
    //! triggers the merge distance
    double adimensionnalMergeDistance_1;
    //! decides whether the output format is binary
    bool binaryOutput;

    //! inner copy of a multiInclusions
    const MultiInclusions<3>* multiInclusions;
private:
    //! write spherical inclusions
    //! \param f : outputstream
    void writeSphericalInclusions(std::ostream& f) const;
    //! write laguerre tessellation
    //! \param f : outputstream
    void writeLaguerreTess(std::ostream& f) const;
};

} // namespace generator
} // namespace mesh
} // namespace merope

#endif /* MICROTOMESH_MESHGENERATOR_HXX_ */

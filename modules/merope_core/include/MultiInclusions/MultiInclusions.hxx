//! Copyright : see license.txt
//!
//! \brief
//!

#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../GenericTools/SOA.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"

#include "../SingleVoxel/SingleVoxel_Headers.hxx"

#include "../MultiInclusions/LaguerreTess.hxx"
#include "../MultiInclusions/Rectangle.hxx"
#include "../MultiInclusions/SphereInclusions.hxx"
#include "../MultiInclusions/ObjectInclusions.hxx"

#include "../MicroInclusion/MicroInclusion.hxx"



namespace merope {

template<unsigned short DIM>
using SOA_for_MultiI = std::enable_if_t<DIM == 2 or DIM == 3,
        std::conditional_t<DIM == 2,
        SOA <smallShape::SphereInc<DIM>, smallShape::ConvexPolyhedronInc<DIM>, smallShape::EllipseInc<DIM>, smallShape::SpheroPolyhedronInc<DIM>>,
        SOA <smallShape::SphereInc<DIM>, smallShape::ConvexPolyhedronInc<DIM>, smallShape::EllipseInc<DIM>, smallShape::SpheroPolyhedronInc<DIM>, smallShape::CylinderInc<3>>
        >>;

template<unsigned short DIM>
class MultiInclusions final :
        public InsideTorus<DIM>,
        private  SOA_for_MultiI<DIM>,
        public MatrixPhaseHolder<PhaseType> {
        //! class for NON-INTERSECTING inclusions inside a matrix
public:
        using SOA_type = SOA_for_MultiI<DIM>;
        //! constructor
        MultiInclusions();
        //! from a PolyInclusions
        void setInclusions(const PolyInclusions<DIM>& polyInc);
        //! from a tessellation
        void setInclusions(LaguerreTess<DIM> polyX);
        //! from spherical inclusions
        void setInclusions(const SphereInclusions<DIM>& sphereI);
        //! from ObjectInclusions
        template<class ObjectInc>
        void setInclusions(const ObjectInclusions<DIM, ObjectInc>& objectInclusion);
        //! from Object vectors
        template<class ObjectInc>
        void setInclusions(const vector<ObjectInc>& vectInclusions, Point<DIM> L);
        //! from a rectangle
        void setInclusions(const Rectangle<DIM>& rect);
        //! enlarge each inclusion
        //! \param identifiers : concerned identifiers of inclusions
        //! \param width : width of the enlargement
        void enlarge(const vector<Identifier>& identifiers, const vector<double>& width);
        void enlarge(const vector<Identifier>& identifiers, double width) {
                this->enlarge(identifiers, vector<double>(identifiers.size(), width));
        }
        //! add a layer
        //! \param layersToAdd : instructions for adding a layer
        void addLayer(vector<smallShape::LayerInstructions> layersToAdd);
        //! add a layer
        //! \param identifiers : concerned identifiers of inclusions
        //! \param newPhases : phases identifiers to be given to the new layer
        //! \param width : widths of the layer to be added
        void addLayer(const vector<Identifier>& identifiers,
                const vector<PhaseType>& newPhase, const vector<double>& width);
        void addLayer(const vector<Identifier>& identifiers, PhaseType newPhase,
                double width);
        //! change the phase of the inclusions
        void changePhase(const vector<Identifier>& identifiers,
                const vector<PhaseType>& newPhase);
        //! sets inner shapes
        //! \param shapes : vector of inclusions
        template<class C>
        void setInnerShapes(const C& inclusions);
        //! \return the identifiers corresponding to a list of phases
        vector<Identifier> getIdentifiers(vector<PhaseType> phases) const;
        //! \return the identifiers for all inclusions
        vector<Identifier> getAllIdentifiers() const;
        //! vector of all phases
        vector<PhaseType> getAllPhases() const;
        //! vector of all phases at a specific layer
        //! \param : layer_index, is the number of the layer (0=close to boundary, +1 = farther from the boundary ...)
        vector<PhaseType> getAllPhases(size_t layer_index) const;
        //! vector of all phases verifying a specific test
        template<class TEST_AND_FILL_FUNCTION, class TEST_FUNCTION>
        vector<PhaseType> getAllPhases(TEST_AND_FILL_FUNCTION test_and_fill_function,
                TEST_FUNCTION insert_matrix_phase) const;
        //! get all the centers of the inclusions
        vector<Point<DIM>> getAllCenters() const;
        //! get inclusions
        using SOA_type::get;
        using SOA_type::apply_on_all;
        //! @return : the inner laguerre tesselation if applicable
        const LaguerreTess<DIM>& getLaguerreTess() const;
        //! @return : is a Laguerre tessellation
        bool isLaguerreTess() const { return laguerreTess != nullptr; }
private:
        //! check whether the identifier is actually in the structure
        //! if not, throws an error
        //! \param layersToAdd : layer that we wish to add
        bool checkAddLayer(vector<smallShape::LayerInstructions> layersToAdd) const;
        //! inner Laguerre tessellation if applicable
        shared_ptr<const LaguerreTess<DIM>> laguerreTess;
};

namespace auxi_MultiInclusions {
template<class C>
using Instruction = std::function<void(C*, smallShape::LayerInstructions)>;

//! get all the identifiers of the inclusions
//! \param NbOfSeeds : number of seeds in the Structure
static vector<Identifier> getAllIdentifiers(size_t NbOfSeeds);

//! \see MultiInclusions<DIM>::enlarge
template<class INCLUSIONVECTOR>
void enlarge_T(vector<smallShape::LayerInstructions> layerInstructions,
        INCLUSIONVECTOR& inclusions);
//! \see MultiInclusions<DIM>::changePhase
template<class INCLUSIONVECTOR>
void changePhase_T(vector<smallShape::LayerInstructions> layerInstructions,
        INCLUSIONVECTOR& inclusions);
//! \see MultiInclusions<DIM>::addLayer
template<class INCLUSIONVECTOR>
void addLayer_T(vector<smallShape::LayerInstructions> instructions,
        INCLUSIONVECTOR& inclusions);
//! adds layers to inclusions accorded to the LayerInstructions
//! both inputs are sorted wrt identifier property
//! \param sortedLayerInstructions : layersInstruction
//! \param sortedPointerInclusionList : vector of pointer to inclusions
template<class C>
void applyLayerInstruction_T(const vector<smallShape::LayerInstructions>& sortedLayerInstructions,
        vector<C*>& sortedPointerInclusionList, Instruction<C> applyInstruction);
//! \return all the identifiers of inclusions that display a certain phase
//! \param phase : reseached phases
//! \param inclusionList : vector of inclusion
template<class INCLUSIONVECTOR>
vector<Identifier> getIdentifiers(vector<PhaseType> phases, const INCLUSIONVECTOR& inclusionList);
//!
template<class C>
vector<C*> sortInclusionAndInstructions(vector<C>& inclusions,
        vector<smallShape::LayerInstructions>& layerInstructions);
}  // namespace auxi_MultiInclusions

}  // namespace merope

#include "../MultiInclusions/MultiInclusions.ixx"



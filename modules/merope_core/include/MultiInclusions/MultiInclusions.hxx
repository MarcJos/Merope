//! Copyright : see license.txt
//!
//! \brief 
//!

#ifndef MULTIINCLUSIONS_HXX_
#define MULTIINCLUSIONS_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"
#include "../MeropeNamespace.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../MultiInclusions/EllipseInclusions.hxx"
#include "../MultiInclusions/LaguerreTess.hxx"
#include "../MultiInclusions/Rectangle.hxx"
#include "../MultiInclusions/SphereInclusions.hxx"
#include "../MultiInclusions/SpheroPolyhedronInclusions.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"



namespace merope {

template<unsigned short DIM>
class MultiInclusions: public InsideTorus<DIM> {
        //! class for NON-INTERSECTING inclusions inside a matrix
public:
        //! constructor
        MultiInclusions();
        //! from a PolyInclusions
        void setInclusions(const PolyInclusions<DIM>& polyInc);
        //! from a tessellation
        void setInclusions(LaguerreTess<DIM> polyX);
        //! from spherical inclusions
        void setInclusions(const SphereInclusions<DIM>& sphereI);
        //! from elliptical inclusions
        void setInclusions(const EllipseInclusions<DIM>& ellipseI);
        //! from a rectangle
        void setInclusions(const Rectangle<DIM>& rect);
        //! from spheroPolyhedrons
        void setInclusions(const vector<smallShape::SpheroPolyhedronInc<DIM>>& spheroPolyhedrons, Point<DIM> L);
        //! from spheroPolyhedrons
        void setInclusions(const SpheroPolyhedronInclusions<DIM>& spheroPolyhedrons);
        //! add a layer
        //! \param layersToAdd : instructions for adding a layer
        void addLayer(vector<smallShape::LayerInstructions> layersToAdd);
        //! add a layer
        //! \param identifiers : concerned phases of inclusions
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
        //! get all the centers of the inclusions
        vector<Point<DIM>> getAllCenters() const;
        //! set the matrix phase
        void setMatrixPhase(PhaseType matrixPhase_) { this->matrixPhase = matrixPhase_; }
        //! get inclusions
        template<class C>
        vector<C>& getInclusions();

        //! getter
        const vector<smallShape::EllipseInc<DIM> >& getEllipseInc() const { return ellipseInc; }
        //! getter
        const vector<smallShape::ConvexPolyhedronInc<DIM> >& getPolyhedrons() const { return polyhedrons; }
        //! getter
        const vector<smallShape::SphereInc<DIM> >& getSphereInc() const { return sphereInc; }
        //! getter
        const vector<smallShape::SpheroPolyhedronInc<DIM> >& getSpheroPolyhedrons() const { return spheroPolyhedrons; }
        //! getter
        PhaseType getMatrixPhase() const { return matrixPhase; }

private:
        //! stores the polyhedrons inclusions
        vector<smallShape::ConvexPolyhedronInc<DIM>> polyhedrons;
        //! stores the spherical inclusions
        vector<smallShape::SphereInc<DIM>> sphereInc;
        //! stores the ellipsoid inclusions
        vector<smallShape::EllipseInc<DIM>> ellipseInc;
        //! stores the ellipsoid inclusions
        vector<smallShape::SpheroPolyhedronInc<DIM>> spheroPolyhedrons;
        //! Phase of the matrix
        PhaseType matrixPhase;
private:
        //! check whether the identifier is actually in the structure
        //! if not, throws an error
        //! \param layersToAdd : layer that we wish to add
        bool checkAddLayer(vector<smallShape::LayerInstructions> layersToAdd) const;
        //! todo
        template<class C>
        void setInclusions_T(const C& vectorOfInclusions);
        //! todo
        template<class LAMBDA_FUNCTION>
        void applyOnAllInclusions(LAMBDA_FUNCTION f);
        //! todo
        template<class LAMBDA_FUNCTION>
        void applyOnAllInclusions(LAMBDA_FUNCTION f) const;
};

namespace auxi_MultiInclusions {
template<class C>
using Instruction = std::function<void(C*, smallShape::LayerInstructions)>;

//! get all the identifiers of the inclusions
//! \param NbOfSeeds : number of seeds in the Structure
static vector<Identifier> getAllIdentifiers(size_t NbOfSeeds);

//! \see MultiInclusions<DIM>::changePhase
template<class INCLUSIONVECTOR>
void changePhase_T(vector<smallShape::LayerInstructions> layerInstructions,
        INCLUSIONVECTOR& inclusions);
//! \see MultiInclusions<DIM>::addLayer
template<class INCLUSIONVECTOR>
void addLayer_T(vector<smallShape::LayerInstructions> layerInstructions,
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
}

} // namespace merope

#include "../MultiInclusions/MultiInclusions.ixx"

#endif /* MULTIINCLUSIONS_HXX_ */

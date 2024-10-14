//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../MultiInclusions/LaguerreTess.hxx"
#include "../MultiInclusions/MultiInclusions.hxx"
#include "../MultiInclusions/SphereSeeds.hxx"
#include "../MesoStructure/Structure.hxx"
#include "../Obsolete_MesoStructure/MicroType.hxx"
#include "../../../AlgoPacking/src/AuxiFunctions.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM>
class Pre_InterfaceMultiInclusions :
    public MicroType,
    public SphereSeeds<DIM>,
    public WithAspratio<DIM>,
    private auxi_function::Deprecated {
public:
    //! constructor
    Pre_InterfaceMultiInclusions();
    //! get the list of the phases of spheres in the sphere collection
    vector<PhaseType> getPhaseList() const;
    //! return the list of Identifiers to the microstructure
    vector<Identifier> getAllIdentifiers() const;
};

template<unsigned short DIM>
class InterfaceMultiInclusions final : public Pre_InterfaceMultiInclusions<DIM> {
public:
    //! constructor
    InterfaceMultiInclusions();
    //! constructor
    InterfaceMultiInclusions(Pre_InterfaceMultiInclusions<DIM>);
    //! setter
    void setLayerList(const vector<smallShape::LayerInstructions>& layerList);
    //! getter
    const vector<smallShape::LayerInstructions>& getLayerList() const;
    //! builds the MultiInclusions
    MultiInclusions<DIM> build();

protected:
    //! layers to be added to the inclusions
    vector<smallShape::LayerInstructions> layerList;
private:
    //! build the MultiInclusions
    MultiInclusions<DIM> buildSimple();
};


template<unsigned short DIM>
class InterfaceStructure : public Colorize, private auxi_function::Deprecated {
public:
    //! constructor
    InterfaceStructure();
    //! main MultiInclusions
    Pre_InterfaceMultiInclusions<DIM> mainInclusions;
    //! secondary MultiInclusions
    Pre_InterfaceMultiInclusions<DIM> secdInclusions;
    //! \return a structure parametrized
    Structure<DIM> build();
    //! setter
    void setErosionPhase(PhaseType erosionPhase);
    void setErosionInclusionsPhase(PhaseType erosioNinclusionsPhase);
    void setErosionWidth(double widthOfErosion_);
    void setInnerPhase(PhaseType innerPhase_);
    void setLength(array<double, DIM> L);

    //! fixme : too much empirical?
    //! Compute the thickness of grain boundaries from its volumic or areal fraction
    //! using fitted formulas on Voronoi microstructures
    //! \param phi Volumic or areal fraction of grain boundaries
    void frac2erosionWidth(double phi);
    //! fixme : too much empirical?
    //! Compute the volumic or areal fraction from the thickness of grain boundaries
    //! using fitted formulas on Voronoi microstructures
    double erosionWidth2frac();

private:
    //! PhaseTypes
    constexpr static PhaseType DEFAULT_INNERPHASE = 0;
    PhaseType innerPhase;
    constexpr static PhaseType DEFAULT_EROSIONPHASE = 1;
    PhaseType erosionPhase;
    constexpr static PhaseType DEFAULT_EROSIONINCLUSIONSPHASE = 2;
    PhaseType erosionInclusionsPhase;
    //! width of the layer on the inclusions
    double erosionWidth;
    //! layer instructions for the main phase
    //! erosionPhase_ : erosion phase to be imposed
    vector<smallShape::LayerInstructions> buildLayerInstructions(PhaseType erosionPhase_) const;
    //! build the main inclusion
    MultiInclusions<DIM> buildMainInclusions();
    //! build the secondary inclusion
    MultiInclusions<DIM> buildSecdInclusions();
    //! \return a function that helps for Erode3Mat, deciding which pahse to change
    std::function<PhaseType(PhaseType, PhaseType)> functionPorositySpheres() const;
    //! fixme : too much empirical?
    //! \return the thickness of grain boundaries from its volumic or areal fraction
    //! using fitted formulas on Voronoi microstructures (! not on Laguerre or Johnson-Mehl!)
    //! \param phi : Volumic or areal fraction of grain boundaries
    double calculeErodeeJdGn(double phi) const;
    //! performs some operations on the InterfaceStructure to make things coherent
    //! forces the definition of mainInclusion to take innerPhase into account
    void preBuild();
    //! performs some operations on the outputStructure
    void postBuild(Structure<DIM>& structure);
    //! return the erosionPhase that is necessary to compute the structure
    PhaseType temporaryErosionPhase() const;
};

}  // namespace merope


#include "../Obsolete_MesoStructure/InterfaceStructure.ixx"



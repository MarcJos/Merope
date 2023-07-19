//! Copyright : see license.txt
//!
//! \brief A simple grid at 1,2 or 3 dimensions
//!

#ifndef _GRID_HXX
#define _GRID_HXX 1


#include "../../../AlgoPacking/src/StdHeaders.hxx"
#include "../localMFront.h"
#include "../MeropeNamespace.hxx"


namespace merope {

//! type definition for a list of node indexes.
typedef std::vector<size_t> NodesList;

class VTKRead;
class VTKstream;
class TIFFRead;

//! Class representing a grid, possibly in 1D, 2D or 3D.
class Grid {

protected:
    unsigned char d; //!< Grid dimension (1D, 2D or 3D).

    size_t ng; //!< Number of nodes in the grid \f$n_T=n_x \times n_y \times n_z\f$.
    size_t nx; //!< Number of nodes in x direction (spatial discretization).
    size_t ny; //!< Number of nodes in y direction (spatial discretization).
    size_t nz; //!< Number of nodes in z direction (spatial discretization).

    double lx;   //!< Length of the grid in x direction.
    double ly;   //!< Length of the grid in y direction.
    double lz;   //!< Length of the grid in z direction.

public:
    //! table to retrieve grid index \f$i_g\f$ from medium position \f${\alpha,i_\alpha}: i_g=mat2Grid[\alpha][i_\alpha]\f$
    std::vector<NodesList> phases;

public:
    //! Default constructor
    Grid();
    //! 1D void constructor
    //! \param lx Length
    //! \param nx Grid size
    Grid(double lx, size_t nx);
    //! 2D void constructor
    //! \param lx,ly Length
    //! \param nx,ny Grid size
    Grid(double lx, double ly, size_t nx, size_t ny);
    //! 3D void constructor
    //! \param lx,ly,lz Length
    //! \param nx,ny,nz Grid size
    Grid(double lx, double ly, double lz, size_t nx, size_t ny, size_t nz);
    //! Grid constructor by reading a periodical continuous geometry
    //! \param geometry the periodical discrete geometry
    Grid(VTKRead& geometry);
    //! Grid constructor by reading a TIFF File
    //! \param tf: TIFF File
    //! \param dx: Grid step
    Grid(TIFFRead& tf, double);
    //! Print the internal variables
    //! \param os: Output stream
    void print(std::ostream& os = std::cout) const;

    unsigned char getDim() const; //!< getter for grid nb of nodes.
    size_t getNx() const; //!< getter for nb of node in x direction.
    size_t getNy() const; //!< getter for nb of node in y direction.
    size_t getNz() const; //!< getter for nb of node in z direction.
    size_t getNg() const; //!< getter for nb of node in grid (nx x ny x nz)
    double getLx() const; //!< getter for x length.
    double getLy() const; //!< getter for y length.
    double getLz() const; //!< getter for z length.

    unsigned short getNbOfPhases() const; //!< getter for number of phases \f$\alpha\f$.

    //! getter to retrieve a list of points for a material (ie a phase).
    //! \param alpha the phase number \f$\alpha\f$
    const NodesList* getPhase(unsigned short alpha) const;
    //! getter to retrieve all points for all phases
    const std::vector<NodesList>& getPhases() const;
    //! Return phase fraction
    //! \param i Phase index starting by 0
    double phaseFraction(unsigned long i) const;
    //! Grid to VTK order conversion
    //! \param Ig Grid order = k+nz*(j+i*ny)
    //! \return Iv VTK order = i+nx*(j+k*ny)
    unsigned Go2VTKo(const size_t Ig) const;

    //! Correction of Divergence Norm
    //! \return Correction factor sqrt(N)/L
    double DivNormFactor() const;

    //! VTK CELL header file
    //! \param fvtk VTK output stream
    void VTKheaderCELL(VTKstream& fvtk) const;
    //! Writes a VTK CELL file.
    //! \param fvtk VTK output stream
    void toVTKCELL(VTKstream& fvtk) const;

    //! Variogramm in the x direction
    //! \param v Variable field
    //! \param fout Output stream
    void varioX(const CastemReal* v, std::ostream& fout) const;
    //! Variogramm in the y direction (2 & 3 D)
    //! \param v Variable field
    //! \param fout Output stream
    void varioY(const CastemReal* v, std::ostream& fout) const;
    //! Variogramm in the z direction (3D only)
    //! \param v Variable field
    //! \param fout Output stream
    void varioZ(const CastemReal* v, std::ostream& fout) const;
    //! Covariance in the x direction
    //! \param a,b Phases
    //! \param fout Output stream
    void covX(unsigned short a, unsigned short b, std::ostream& fout) const;
    //! Covariance in the y direction
    //! \param a,b Phases
    //! \param fout Output stream
    void covY(unsigned short a, unsigned short b, std::ostream& fout) const;
    //! Covariance in the z direction
    //! \param a,b Phases
    //! \param fout Output stream
    void covZ(unsigned short a, unsigned short b, std::ostream& fout) const;

    Grid(const Grid&);            //!< Forbid use of copy constructor.
    Grid& operator=(const Grid&); //!< Forbid use of affectation.

protected:
    //! Common initialisation
    void initCommon();
    //! VTK File Grid reader
    //! \param vf: VTK File Grid
    //! \param T: Type of the material ID
    template<typename T>
    void readVTK(VTKRead& vf);
    //! One material
    void oneMat();

protected:
    //! Fill a 1/2/3D table with phaseIdx saved in each point of the grid
    //! \param ids Filled table defining geometry phases on the grid
    void vtkReorderMaterialIdx(unsigned short* const ids) const;
    //! VTK CELL header DIM-D file
    //! \param fvtk VTK file stream
    template<unsigned short DIM>
    void VTKheaderCELL_T(VTKstream& fvtk) const;
    //! Writes a VTK CELL DIM-D file, for VTK generation purpose.
    //! \param fvtk VTK file stream
    template<unsigned short DIM>
    void toVTKCELL_T(VTKstream& fvtk) const;
};

} // namespace merope


#include "../Grid/Grid.ixx"

#endif // _GRID_HXX


//! Copyright : see license.txt
//!
//! \briefReader for a periodic discrete medium
//! given by a VTK File (DATASET STRUCTURED_POINTS)
//

#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {

//! Component
struct Comp {
    string name;
    vector<streampos> p;
    //! Default constructor
    Comp() {}
    //! Complete constructor
    //! \param name Field name
    Comp(const char* name);
    //! Print for test
    //! \param fout Output stream
    void print(ostream& fout = cout) const;
};

//! VTK reading structure
class VTKRead : public ifstream {
public:
    //! Material Id data type
    enum MType {
        UC, US, ST, CH
    };
protected:
    size_t nx, ny, nz;  //!< Grid dimensions
    size_t ng;
    unsigned char d;  //!< Space dimension
    double dx, dy, dz;  //!< Spacing
    vector<Comp> cmps;  //!< Variable components
    Comp vel;  //!< Velocity
    streampos pMId;  //!< Material Id position
    MType mt;  //!< Material Id data type
    bool isList;  //!< True after the component parsing
    bool isFor;  //!< True for a formated VTK
public:
    //! VTK reader constructor
    //! \param fname File name
    VTKRead(const char* fname);
    // Print for test
    // fout: Output stream
    void print(ostream& fout = cout) const;
    //! Read the material data
    //! \param tab Table value, in BigEndian
    //! \param T: Type of the material ID
    template<typename T>
    void readMat(vector<T>& tab);
    //! Transfers parameters to a Grid object
    //! \param d Grid dimension
    //! \param nx,ny,nz Grid sizes
    //! \param lx,ly,lz Periodical geometry dimensions
    //! \param mt Material Id data type
    void getParamG(unsigned char& d, size_t& nx, size_t& ny, size_t& nz,
        double& lx, double& ly, double& lz, MType& mt) const;
    //! Get the component list
    //! \return the component list
    const vector<Comp>& getCmps() const;

    //! Read one component
    //! \param c Component
    //! \param tabv Table value, in BigEndian
    void readComp(const Comp& c, vector<float>& tabv);
    //! Read one component
    //! \param name Component name
    //! \param tabv Table value, in BigEndian
    //! \return the component number
    unsigned char readComp(const char* name, vector<float>& tabv);

    //! Parse the variable list
    //! \param strainN Strain tensor field name
    //! \param stressN Stress tensor field name
    void variables(const char* strainN, const char* stressN);
protected:
    //! CELL DATA field reading
    //! \param p Position in the stream
    //! \param tv Table value, in BigEndian
    void fromCELL(const streampos p, float* tv);
    //! Gather the components
    //! \param strainN Strain tensor field name
    //! \param stressN Stress tensor field name
    void gatherL(const char* strainN, const char* stressN);
};

}  // namespace merope


#include "../VTKinout/VTKRead.ixx"




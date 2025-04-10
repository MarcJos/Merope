//! Copyright : see license.txt
//!
//! \briefReader for a periodic discrete medium
//! given by a VTK File (DATASET STRUCTURED_POINTS)
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"


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
};

}  // namespace merope


#include "../VTKinout/VTKRead.ixx"




//! Copyright : see license.txt
//!
//! \brief Defines a VTK stream format
//

#ifndef _VTKSTREAM_HXX
#define _VTKSTREAM_HXX 1


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

//! VTK file managment
class VTKstream: public std::ofstream {
private:
        bool modeIsSet;   //!< first call of setCELL or setPOINT set the mode
        bool modeIsCELL;  //!< writing mode is CELL_DATA or POINT_DATA

public:
        //! Constructor
        //! \param filename Name of the file
        //! \param comment optionnal comment
        //! \param mode Flags describing the requested i/o mode for the file
        VTKstream(const char* filename, const char* comment = NULL,
                ios_base::openmode mode = ios_base::out);
        //! Opening a file
        //! \param filename Name of the file
        //! \param comment optionnal comment
        //! \param mode Flags describing the requested i/o mode for the file
        void open(const char* filename, const char* comment = NULL,
                ios_base::openmode mode = ios_base::out);

        //! Generate start of a 2D file of type STRUCTURED_POINTS
        //! \param nx,ny grid sizes
        //! \param dx,dy grid steps
        void STRUCTURED_POINTS(const unsigned nx, const unsigned ny,
                const double dx, const double dy);
        //! Generate start of a 3D file of type STRUCTURED_POINTS
        //! \param nx,ny,nz grid sizes
        //! \param dx,dy,dz grid steps
        void STRUCTURED_POINTS(const unsigned nx, const unsigned ny,
                const unsigned nz, const double dx, const double dy,
                const double dz);

        //! Set writing mode as "CELL_DATA n"
        //! \param n : writing mode
        void setCELL(const unsigned n);
        //! Set writing mode as "POINT_DATA n"
        //! \param n : writing mode
        void setPOINT(const unsigned n);

        //! Binary writing (Big Endian)
        //! \param data Data to write
        template<typename DATA>
        void write(DATA data);
private:
        //! Binary VTK header
        //! \param comment optionnal comment
        void headerBIN(const char* comment = NULL);
};

} // namespace merope


#include "../VTKinout/VTKStream.ixx"
#endif // _VTKSTREAM_HXX


//! Copyright : see license.txt
//!
//! \briefDefines a VTK stream format
//

#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

//! VTK file managment
class VTKstream : public std::ofstream {
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

        template<unsigned short DIM>
        void STRUCTURED_POINTS(const array<size_t, DIM>& n, const array<double, DIM>& dx) {
                if constexpr (DIM == 3) {
                        STRUCTURED_POINTS(n[0], n[1], n[2], dx[0], dx[1], dx[2]);
                } else if constexpr (DIM == 2) {
                        STRUCTURED_POINTS(n[0], n[1], dx[0], dx[1]);
                } else {
                        cerr << __PRETTY_FUNCTION__ << endl; throw std::invalid_argument("DIM");
                }
        }

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

        //! @brief : write vector data
        //! @param vectorData : vector of data
        //! @param n : vector containing number of voxels in each 3 dimensions
        template<class PHASE_OUT, class VECTOR_DATA, class ARRAY_DIM>
        void writeVector(const VECTOR_DATA& vectorData, ARRAY_DIM n);

private:
        //! Binary VTK header
        //! \param comment optionnal comment
        void headerBIN(const char* comment = NULL);
};

}  // namespace merope


#include "../VTKinout/VTKStream.ixx"



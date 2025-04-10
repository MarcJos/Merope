//! Copyright : see license.txt
//!
//! \briefdefines a VTK stream format
//!


#include "VTKinout/VTKStream.hxx"


namespace merope {

//----------------------------------------------
inline void VTKstream::headerBIN(const char* const comment) {
    *this << "# vtk DataFile Version 4.0" << endl;
    if (comment) {
        *this << comment << endl;
    } else {
        *this << "Periodic GridGenerator" << endl;
    }
    *this << "BINARY" << endl;
}

VTKstream::VTKstream(const char* const filename, const char* const comment,
    ios_base::openmode mode) :
    ofstream(filename, mode) {
    if (mode != ofstream::app) {
        modeIsSet = false;
        headerBIN(comment);
    } else {
        modeIsSet = true;
        modeIsCELL = true;
    }
}

void VTKstream::open(const char* const filename, const char* const comment,
    ios_base::openmode mode) {
    ofstream::open(filename, mode);
    if (mode != ofstream::app) {
        modeIsSet = false;
        headerBIN(comment);
    } else {
        modeIsSet = true;
        modeIsCELL = true;
    }
}

void VTKstream::STRUCTURED_POINTS(const unsigned nx, const unsigned ny,
    const double dx, const double dy) {
    STRUCTURED_POINTS(nx, ny, 0, dx, dy, 1.);
}

void VTKstream::STRUCTURED_POINTS(const unsigned nx, const unsigned ny,
    const unsigned nz, const double dx, const double dy, const double dz) {
    *this << "DATASET STRUCTURED_POINTS" << endl;
    *this << "DIMENSIONS " << nx + 1 << " " << ny + 1 << " " << nz + 1 << endl;
    *this << "ORIGIN  0. 0. 0." << endl;
    *this << "SPACING " << dx << " " << dy << " " << dz << endl;
}

void VTKstream::setCELL(const unsigned n) {
    if (!modeIsSet) {
        modeIsSet = true;
        modeIsCELL = false;
    }
    if (!modeIsCELL) {
        *this << "CELL_DATA " << n << endl;
        modeIsCELL = true;
    }
}

void VTKstream::setPOINT(const unsigned n) {
    if (!modeIsSet) {
        modeIsSet = true;
        modeIsCELL = true;
    }
    if (modeIsCELL) {
        *this << "POINT_DATA " << n << endl;
        modeIsCELL = false;
    }
}

}  // namespace merope


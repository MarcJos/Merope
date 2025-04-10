//! Copyright : see license.txt
//!
//! \briefReader for a periodic discrete medium
//! given by a VTK File (DATASET STRUCTURED_POINTS)
//!
//!


#include "VTKinout/VTKRead.hxx"


namespace merope {
// Define a constant for the maximum length of a VTK line
constexpr size_t VTK_LINE_MAX = 256;

Comp::Comp(const char* const name_i) :
    name(name_i), p(1, 0) {}

void Comp::print(std::ostream& fout) const {
    fout << name << " ";
    for (vector<streampos>::const_iterator i = p.begin(); i != p.end(); ++i) {
        fout << *i << " ";
    }
    fout << endl;
}

VTKRead::VTKRead(const char* const fname) :
    vel(""), isList(false) {
    open(fname);
    if (!(*this)) {
        throw(logic_error(
            "VTKread::VTKread: Can't open file '" + string(fname) + "'"));
    }

    // VTK header
    char line[VTK_LINE_MAX];
    getline(line, VTK_LINE_MAX);
    getline(line, VTK_LINE_MAX);

    // Data type
    getline(line, VTK_LINE_MAX);
    if ("BINARY" != string(line) and "ASCII" != string(line)) {
        throw(logic_error(
            "VTKread::VTKread: A VTK file is either BINARY or ASCII"));
    }
    if ("ASCII" == string(line))
        isFor = true;
    else
        isFor = false;

    // VTK dataset Type
    getline(line, VTK_LINE_MAX);
    if ("DATASET STRUCTURED_POINTS" != string(line)) {
        throw(logic_error(
            "VTKread::VTKread: The VTK File must be a DATASET STRUCTURED_POINTS"));
    }

    // VTK dimensions
    getline(line, VTK_LINE_MAX);
    if (sscanf(line, "DIMENSIONS%zu%zu%zu", &nx, &ny, &nz) != 3) {
        throw(logic_error(
            "VTKread::VTKread: In your VTK File '" + string(fname)
            + "' DIMENSIONS are missing"));
    }
    // Points -> Cells correction
    --nx;
    --ny;
    if (nz > 1) {
        d = 3;
        --nz;
    } else {
        d = 2;
    }

    ng = nx * ny * nz;

    // VTK Origin and spacing
    getline(line, VTK_LINE_MAX);
    getline(line, VTK_LINE_MAX);
    if (sscanf(line, "SPACING%lf%lf%lf", &dx, &dy, &dz) != 3) {
        string msg("VTKread::VTKread: ");
        msg += "In your VTK File '" + string(fname) + "' SPACING is missing";
        throw(logic_error(msg));
    }

    // CELL Data line
    getline(line, VTK_LINE_MAX);
    size_t nT;
    if (sscanf(line, "CELL_DATA%zu", &nT) != 1) {
        string msg("VTKRead::VTKRead: ");
        msg += "In VTK File '" + string(fname) + "' CELL_DATA expected";
        throw(logic_error(msg));
    }
    if (nT != ng) {
        string msg("VTKRead::VTKRead: ");
        msg += "In VTK File '" + string(fname) + "' the cell number is wrong";
        throw(logic_error(msg));
    }

    // Scalars
    getline(line, VTK_LINE_MAX);
    char dataType[VTK_LINE_MAX];
    if (sscanf(line, "SCALARS%*s%s", dataType) != 1) {
        string msg("VTKRead::VTKRead: ");
        msg += "In VTK File '" + string(fname) + "' SCALARS expected";
        throw(logic_error(msg));
    }
    // Table reading
    string dT(dataType);
    if ("unsigned_char" == dT) {
        mt = UC;
    } else if ("unsigned_short" == dT) {
        mt = US;
    } else if ("short" == dT) {
        mt = ST;
    } else if ("char" == dT) {
        mt = CH;
    } else {
        string msg("VTKRead::VTKRead: ");
        msg +=
            "In VTK File '" + string(fname)
            + "' SCALARS type must be char, unsigned_char, short or unsigned_short";
        throw(logic_error(msg));
    }

    // Material Id data position
    pMId = tellg();
    // Pass an optional line
    getline(line, VTK_LINE_MAX);
    if (!strncmp(line, "LOOKUP_TABLE", 12)) {
        pMId = tellg();
    }
}

void VTKRead::print(std::ostream& fout) const {
    fout << nx << " " << ny << " " << nz << endl;
    fout << dx << " " << dy << " " << dz << endl;
    fout << pMId << endl;
    cout << "Components :" << endl;
    for (vector<Comp>::const_iterator i = cmps.begin(); i != cmps.end(); ++i) {
        i->print(fout);
    }
    if (vel.p[0]) {
        cout << "Velocity :" << endl;
        vel.print(fout);
    }
}

void VTKRead::getParamG(unsigned char& d_o, size_t& nx_o, size_t& ny_o,
    size_t& nz_o, double& lx, double& ly, double& lz, MType& mt_o) const {
    d_o = d;
    nx_o = nx;
    ny_o = ny;
    nz_o = nz;
    lx = nx * dx;
    ly = ny * dy;
    lz = nz * dz;
    mt_o = mt;
}
}  // namespace merope


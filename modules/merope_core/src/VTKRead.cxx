//! Copyright : see license.txt
//!
//! \briefReader for a periodic discrete medium
//! given by a VTK File (DATASET STRUCTURED_POINTS)
//!
//!


#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "VTKinout/VTKRead.hxx"


#include "MeropeNamespace.hxx"


namespace merope {

#define VTK_LINE_MAX 256

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

void VTKRead::variables(const char* const strainN, const char* const stressN) {
    if (isList) return;
    char line[VTK_LINE_MAX];
    // Cell and Point field dimension
    size_t NTP;
    switch (d) {
    case 1:
        NTP = nx + 1;
        break;
    case 2:
        NTP = (nx + 1) * (ny + 1);
        break;
    default:
        NTP = (nx + 1) * (ny + 1) * (nz + 1);
    }

    // Position after Material data
    size_t N = ng;
    seekg(pMId);
    switch (mt) {
    case UC:
        seekg(N * sizeof(unsigned char), ios_base::cur);
        break;
    case US:
        seekg(N * sizeof(unsigned short), ios_base::cur);
        break;
    case ST:
        seekg(N * sizeof(short), ios_base::cur);
        break;
    case CH:
        seekg(N * sizeof(char), ios_base::cur);
        break;
    }
    getline(line, VTK_LINE_MAX);

    istringstream is;
    string w;
    vel.p[0] = 0;
    do {
        getline(line, VTK_LINE_MAX);
        if (eof()) break;
        is.clear();
        is.str(line);
        is >> w;
        if ("CELL_DATA" == w) {
            N = ng;
            getline(line, VTK_LINE_MAX);
            is.clear();
            is.str(line);
            is >> w;
        } else if ("POINT_DATA" == w) {
            N = NTP;
            getline(line, VTK_LINE_MAX);
            is.clear();
            is.str(line);
            is >> w;
        }

        Comp* cmp = NULL;
        if ("SCALARS" == w) {
            is >> w;
            cmps.push_back(Comp(w.c_str()));
            cmp = &cmps.back();
            is >> w;
            if ("float" != w) {
                throw(logic_error("VTKRead::variables: float expected"));
            }
            // LOOKUP_TABLE
            getline(line, VTK_LINE_MAX);
            cmp->p[0] = tellg();
            seekg(N * sizeof(float), ios_base::cur);
        } else if ("VECTORS" == w) {
            is >> w;
            vel.name = w;
            is >> w;
            if ("float" != w) {
                throw(logic_error("VTKRead::variables: float expected"));
            }

            vel.p[0] = tellg();
            seekg(3 * N * sizeof(float), ios_base::cur);
        } else {
            cout << "Line = " << line << endl;
            throw(logic_error(
                "VTKRead::variables: Type must be SCALARS or VECTORS"));
        }
        getline(line, VTK_LINE_MAX);
    } while (true);
    clear();
    isList = true;

    // Gather the components
    gatherL(strainN, stressN);
}

void VTKRead::gatherL(const char* const strainN, const char* const stressN) {
    vector<Comp> cm2;
    Comp c2("");

    for (vector<Comp>::const_iterator c = cmps.begin(); c != cmps.end();) {
        c2.p.clear();
        // Stress
        if ("SXX" == c->name) {
            c2.name = stressN;
            c2.p.push_back(c->p[0]);
            ++c;
            if ("SYY" != c->name) {
                throw(logic_error(
                    "VTKRead::Lana: unknown component :" + c->name));
            }
            c2.p.push_back(c->p[0]);
            ++c;
            if ("SZZ" != c->name) {
                throw(logic_error(
                    "VTKRead::Lana: unknown component :" + c->name));
            }
            c2.p.push_back(c->p[0]);
            ++c;
            if ("SXY" != c->name) {
                throw(logic_error(
                    "VTKRead::Lana: unknown component :" + c->name));
            }
            c2.p.push_back(c->p[0]);
            if (3 == d) {
                ++c;
                if ("SXZ" != c->name) {
                    throw(logic_error(
                        "VTKRead::Lana: unknown component :" + c->name));
                }
                c2.p.push_back(c->p[0]);
                ++c;
                if ("SYZ" != c->name) {
                    throw(logic_error(
                        "VTKRead::Lana: unknown component :" + c->name));
                }
                c2.p.push_back(c->p[0]);
            }
            ++c;
        }
        // Strain
        else if ("EXX" == c->name) {
            c2.name = strainN;
            c2.p.push_back(c->p[0]);
            ++c;
            if ("EYY" != c->name) {
                throw(logic_error(
                    "VTKRead::Lana: unknown component :" + c->name));
            }
            c2.p.push_back(c->p[0]);
            ++c;
            if ("EZZ" != c->name) {
                throw(logic_error(
                    "VTKRead::Lana: unknown component :" + c->name));
            }
            c2.p.push_back(c->p[0]);
            ++c;
            if ("GXY" != c->name) {
                throw(logic_error(
                    "VTKRead::Lana: unknown component :" + c->name));
            }
            c2.p.push_back(c->p[0]);
            if (3 == d) {
                ++c;
                if ("GXZ" != c->name) {
                    throw(logic_error(
                        "VTKRead::Lana: unknown component :" + c->name));
                }
                c2.p.push_back(c->p[0]);
                ++c;
                if ("GYZ" != c->name) {
                    throw(logic_error(
                        "VTKRead::Lana: unknown component :" + c->name));
                }
                c2.p.push_back(c->p[0]);
            }
            ++c;
        }
        // Vectorial component
        else {
            size_t n = c->name.find_first_of('[');
            if (n != string::npos) {
                string deb = c->name.substr(0, n);
                c2.name = deb;
                while (c != cmps.end() and deb == c->name.substr(0, n)) {
                    c2.p.push_back(c->p[0]);
                    ++c;
                }
            }
            // Simple component
            else {
                c2.name = c->name;
                c2.p.push_back(c->p[0]);
                ++c;
            }
        }

        cm2.push_back(c2);
    }
    cm2.swap(cmps);
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

void VTKRead::fromCELL(const streampos p, float* const tv) {
    seekg(p);
    read((char*)tv, ng * sizeof(float));
}

void VTKRead::readComp(const Comp& c, vector<float>& tabv) {
    if (isFor)
        throw(logic_error(
            "VTKRead::readComp: Components can only be read in a binary file"));

    tabv.resize(ng * c.p.size());
    float* v = &tabv[0];
    unsigned char i = 0;
    for (; i < c.p.size(); ++i, v += ng) {
        fromCELL(c.p[i], v);
    }
}

unsigned char VTKRead::readComp(const char* const name, vector<float>& tabv) {
    if (isFor)
        throw(logic_error(
            "VTKRead::readComp: Components can only be read in a binary file"));

    for (vector<Comp>::const_iterator c = cmps.begin(); c != cmps.end(); ++c) {
        if (c->name == name) {
            readComp(*c, tabv);
            return static_cast<unsigned char>(c->p.size());
        }
    }
    throw(logic_error(
        "VTKRead::readComp: " + string(name) + ": unknown component"));
}

const vector<Comp>& VTKRead::getCmps() const {
    return cmps;
}

}  // namespace merope

//! Copyright : see license.txt
//!
//! \briefMedium from a TIFF File
//!


#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "VTKinout/TIFFRead.hxx"
#include "VTKinout/VTKStream.hxx"


#include "MeropeNamespace.hxx"


namespace merope {

TIFFRead::TIFFRead(const char* const nom) :
    ifstream(nom) {
    Header();
}

inline void TIFFRead::read16(unsigned short& i16) {
    read((char*)&i16, 2);
    if (isMM) vtkByteSwapper<2>::Swap(&i16);
}

inline void TIFFRead::read32(unsigned& i32) {
    read((char*)&i32, 4);
    if (isMM) vtkByteSwapper<4>::Swap(&i32);
}

void TIFFRead::Header() {
    char c[2];
    read(c, 2);
    if (('M' == c[0]) and ('M' == c[1]))
        isMM = true;
    else if (('I' == c[0]) and ('I' == c[1]))
        isMM = false;
    else
        throw(logic_error("TIFFRead::Header: Début difficile"));

    // Lecture de l'offset de l'Image File Directory (IFD)
    unsigned i32;
    seekg(4, ios_base::beg);
    read32(i32);

    // Dimensionnement du tableau des Directory Entries
    unsigned short i16;
    seekg(i32, ios_base::beg);
    read16(i16);
    DEs.resize(i16);

    // Lecture des Directory Entries
    vector<DE>::iterator i;
    for (i = DEs.begin(); i != DEs.end(); ++i) {
        read16(i->tag);
        read16(i->type);
        read32(i->count);
        read32(i->value);
    }
}

unsigned TIFFRead::TIFFGetField(const unsigned short tag) const {
    vector<DE>::const_iterator i;
    for (i = DEs.begin(); i != DEs.end(); ++i) {
        if (i->tag == tag) return i->value;
    }
    throw logic_error("TIFFRead::TIFFGetField : Le champ n'existe pas");
}

void TIFFRead::readValues(unsigned short* v) {
    // Positionnement sur le tableau de valeurs
    seekg(TIFFGetField(TIFFTAG_STRIPOFFSETS), ios_base::beg);
    // Lecture
    read((char*)v, TIFFGetField(TIFFTAG_STRIPBYTECOUNTS));
}

void TIFFRead::print(ostream& fout) const {
    fout << "Nombre de Directory Entries : " << DEs.size() << endl;
    vector<DE>::const_iterator i;
    for (i = DEs.begin(); i != DEs.end(); ++i) {
        fout << i->tag << " ";
        fout << i->type << " ";
        fout << i->count << " ";
        fout << i->value << endl;
    }
}

}  // namespace merope


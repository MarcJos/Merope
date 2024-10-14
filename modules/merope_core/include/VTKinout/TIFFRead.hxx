//! Copyright : see license.txt
//!
//! \briefMedium from a TIFF File
//

#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {

#define TIFFTAG_IMAGEWIDTH      256 // Image width in pixels
#define TIFFTAG_IMAGELENGTH     257 // Image height in pixels
#define TIFFTAG_STRIPOFFSETS    273 // Offsets to data strips
#define TIFFTAG_STRIPBYTECOUNTS 279 // Bytes counts for strips

//! TIFF file 8 Bit
class TIFFRead : public std::ifstream {
    bool isMM;  // Détection du Big Endian
    // Directory Entry
    struct DE {
        unsigned short tag, type;
        unsigned count, value;
    };

    // Directory Entries
    std::vector<DE> DEs;
public:
    //! Constructeur
    //! \param nom : Nom du fichier
    TIFFRead(const char* const nom);
    //! Lecture des valeurs
    //! \param v : Tableau (16bits) pour la récupération des valeurs
    void readValues(unsigned short* v);
    //! Recherche d'un champ
    //! \param tag : Label du champ
    unsigned TIFFGetField(const unsigned short tag) const;
    //! Pour test : impression
    //! \param fout : Flux de sortie
    void print(std::ostream& fout = std::cout) const;
private:
    //! Lecture d'un entier 16 bits
    //! i16 : Entier à écrire
    void read16(unsigned short& i16);
    //! Lecture d'un entier 32 bits
    //! \param i32 : Entier à écrire
    void read32(unsigned& i32);
    //! Lecture de l'entête
    void Header();
};

}  // namespace merope





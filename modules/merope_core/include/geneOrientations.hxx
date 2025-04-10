//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../GenericMerope/StdHeaders.hxx"


namespace merope {
namespace geneOrientation {

void readOrient0(const char* const, std::ofstream&, std::ofstream&);
void printRandomOrient_3D(const unsigned short N, const unsigned short seed_i, std::string nameFile);
void randomOrient0(const unsigned short N, const unsigned short seed_i, std::ofstream&);
//! Difference entre 3 methodes de generation de nombres reels aleatoires
void test();
void EulerToAxis(const double* const ang, double* const vec);
//! \return une erreur si vec ne contient pas une base orthonormale
//! \param vec constitutu� de 3 vecteurs concat�n�s {v0,v1,v2,u0,u1,u2,w0,w1,w2}
void checkOrth(const double* const vec);
//! \return un point dans le premier cadran de la sphere unite
//! \param v: Coordonnees d'un vecteur dans le premier cadran
//!     0 <= vy <= vx <= vz
void PremCad(const double* const vec, double* const v);
//! \return un r�sultat imprime
void printPremCad(const double* const ang, const double* const vec, ofstream& f);
void printAllCad(const double* const ang, const double* const vec, ofstream& f);
void printLicos(const unsigned short n, const double* const vec, ofstream& f);
void printTmfft(const unsigned short n, const double* const ang, ofstream& f);
void printNoFormatting(const unsigned short n, const double* const vec, ofstream& f);

}  // geneOrientation
}  // namespace merope





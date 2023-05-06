//! Copyright : see license.txt
//!
//! localMFront.h
//! Format des lois de comportement et materiau
//!

#ifndef _LOCALMFRONT_H
#define _LOCALMFRONT_H 1

/* Integer type */
#ifdef UNIX32
typedef int CastemInt;
#else
#ifdef UNIX64
typedef long CastemInt;
#else
typedef long CastemInt;
#endif
#endif

/* Real type */
#ifdef MFRONT_USE_LONGDOUBLE
typedef long double CastemReal;
#else
typedef double CastemReal;
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Prototype for material law */
typedef CastemReal (*MatFctPtr)(const CastemReal* const);

/* Prototype for mfront fonctions */
typedef void (*MFRONTFctPtr)(CastemReal* const, /* stress                   */
CastemReal* const, /* internal state variables */
CastemReal* const, /* tangent operator         */
CastemReal* const, CastemReal* const, CastemReal* const, CastemReal* const,
		CastemReal* const, CastemReal* const, CastemReal* const,
		const CastemReal* const, /* strain tensor    */
		const CastemReal* const, /* strain increment */
		const CastemReal* const, const CastemReal* const, /* time increment   */
		const CastemReal* const, /* temperature      */
		const CastemReal* const, /* temperature increment    */
		const CastemReal* const, /* external state variables */
		const CastemReal* const, /* external state variables increments   */
		const char* const, const CastemInt* const, /* modelling hypothesis                  */
		const CastemInt* const, const CastemInt* const, /* number of components of tensors       */
		const CastemInt* const, /* number of internal state variables    */
		const CastemReal* const, /* material properties               */
		const CastemInt* const, /* number of material properties         */
		const CastemReal* const, const CastemReal* const, /* rotation matrix                       */
		CastemReal* const, /* estimation of the next time increment */
		const CastemReal* const, const CastemReal* const,
		const CastemReal* const, const CastemInt* const, const CastemInt* const,
		const CastemInt* const, const CastemInt* const, const CastemInt* const,
		CastemInt* const, const int /* hidden fortran parameter */);

#ifdef __cplusplus
}
#endif

#endif /* _LOCALMFRONT_H */


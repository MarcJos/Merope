//! Copyright : see license.txt
//!
//! \brief
//
#pragma once



#include "../MeropeNamespace.hxx"


namespace merope {

//! Hyperbole
class Hyper {
    //! Point milieu
    double xc, yc, zc;
    //! Coordonnees du vecteur u dans l'ancien repere (x,y)
    double ux, uy, uz;
    //! Pente de l'asymptote montante et ses puissances 2 et 4
    double m, m2, m4;
    //! Abscisse du sommet DR/2
    double DRs2;
    //! Precision
    double EPS2;
public:
    //! \return Constructeur 2D
    //! \param xa,ya: Point A
    //! \param  xb,yb: Point B
    //! \param  Ra,Rb: Rayons
    //! \param  EPS2: Precision
    Hyper(const double xa, const double ya, const double xb, const double yb,
        const double Ra, const double Rb, const double = 1e-20);
    // \return Constructeur d'hyperbole (3D)
    // \param xa,ya,za: Point A
    // \param xb,yb,zb: Point B
    // \param Ra,Rb: Rayons
    // \param EPS2: Precision
    Hyper(const double xa, const double ya, const double za, const double xb,
        const double yb, const double zb, const double Ra, const double Rb,
        const double = 1e-20);
    //! Distance a l'hyperbole dans le repere d'origine (2D)
    //! \param xm,ym: Point a projeter sur l'hyperbole
    double distanceRO(double xm, double ym) const;
    //! Distance a l'hyperbole dans le repere d'origine (3D)
    //! \param xm,ym,zm: Point a projeter sur l'hyperbole
    double distanceRO(double xm, double ym, double zm) const;
private:
    //! \return Projection perpendiculaire a l'asymptote montante
    //! \param um,vm: Point a projeter
    //! \param up,vp: Point projete
    void projM(const double um, const double vm, double& up, double& vp) const;
    //! \return Projection perpendiculaire a l'asymptote descendante
    //! \param um,vm: Point a projeter
    //! \param up,vp: Point projete
    void projD(const double um, const double vm, double& up, double& vp) const;
    //! \return Initialisation de la recherche
    //! \param um,vm: Point a projeter
    //! \param up,vp: Point de depart de la recherche
    void init(const double um, const double vm, double& up, double& vp) const;
    //! Une iteration de Newton
    //! \param um,vm: Point a projeter
    //! \param u,v: Ancienne et nouvelle iteration
    //! \param retourne la distance de la derniere iteration
    double oneIter(const double um, const double vm, double& u,
        double& v) const;
    //! Distance a l'hyperbole
    //! \param um,vm: Point a projeter
    double distance(const double um, const double vm) const;
};

}  // namespace merope




//! Copyright : see license.txt
//!
//! \brief

#pragma once


namespace merope {

template<unsigned short DIM>
inline void SphereInclusions<DIM>::covX(const unsigned nx,
    std::ostream& fout) const {
    covDirection(0, nx, fout);
}

template<unsigned short DIM>
inline void SphereInclusions<DIM>::covY(const unsigned ny,
    std::ostream& fout) const {
    covDirection(1, ny, fout);
}

template<unsigned short DIM>
inline void SphereInclusions<DIM>::covZ(const unsigned nz,
    std::ostream& fout) const {
    covDirection(2, nz, fout);
}

template<unsigned short DIM>
inline void SphereInclusions<DIM>::covDirection(size_t direction,
    const unsigned n_xi, std::ostream& fout) const {
    // fixme
    double V = this->tore.volume();
    double volFracCarre = auxi_function::puissance<2>(geomTools::volume_all(this->theSpheres) / V);
    double dxi = this->tore.L(direction) / n_xi;
    array<double, DIM> x = { 0, 0, 0 };
    unsigned n = n_xi / 2 + 1, i;
    for (x[direction] = 0, i = 0; i < n; ++i, x[direction] += dxi) {
        /*fout << x[direction] << " " << cS.volInter(x, getL()) / V - fv2 << endl;
         */
        Merope_error_not_done();
    }
}

}  // namespace merope




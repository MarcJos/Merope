//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef PATH_IXX_
#define PATH_IXX_

namespace sac_de_billes {
using namespace std;

template<unsigned short DIM>
inline const vector<array<unsigned short, DIM> >& Path_auxi::TABCORNER() const {
    static_assert(DIM == 2 or DIM == 3);
    if constexpr (DIM == 3) {
        return TABCORNER3D;
    } else if constexpr (DIM == 2) {
        return TABCORNER2D;
    }
}

template<unsigned short DIM, size_t NB_NGHB>
inline array<array<int, DIM>, nbSubcubes<DIM>(NB_NGHB)> Path::createParcours(){
   array<array<int, DIM>, nbSubcubes<DIM>(NB_NGHB)> path_ = { };
   int compteur = 0;
   constexpr auto ix = IX;
   if constexpr (DIM == 3) {
       loop<0, NB_NGHB, 0, NB_NGHB, 0, NB_NGHB>([&](const array<long, 3u> &i) {
           path_[compteur] = array<int, DIM> { ix[i[0]], ix[i[1]], ix[i[2]] };
           compteur++;
       });
   } else if constexpr (DIM == 2) {
       loop<0, NB_NGHB, 0, NB_NGHB>([&](const array<long, 2u> &i) {
           path_[compteur] = array<int, DIM> { ix[i[0]], ix[i[1]] };
           compteur++;
       });
   } else {
       throw invalid_argument("path::createParcours(). Incorrect DIM.");
   }

   auto norme2 = [](auto x) {
       double res = 0;
       for (int j = 0; j < DIM; j++) {
           res += x[j] * x[j];
       }
       return res;
   };

   auto criterionSort = [norme2](auto x, auto y) {
       return norme2(x) < norme2(y);
   };
   sort(path_.begin(), path_.end(), criterionSort);
   return path_;
}

template<unsigned DIM, size_t NB_NGHB>
inline const array<array<int, DIM>, nbSubcubes<DIM>(NB_NGHB)>& Path::pathForNeighbors() const {
    if constexpr (DIM == 3) {
        if constexpr (NB_NGHB == 5) {
            return pathNghb_3D_5;
        } else if constexpr (NB_NGHB == 7) {
            return pathNghb_3D_7;
        }
    } else if constexpr (DIM == 2) {
        if constexpr (NB_NGHB == 5) {
            return pathNghb_2D_5;
        } else if constexpr (NB_NGHB == 7) {
            return pathNghb_2D_7;
        }
    }
    throw invalid_argument("Path::pathForNeighbors Incorrect DIM");
}

template<unsigned short DIM>
inline int Path::fromCorner2Index(array<int, DIM> posRelative) const {
    int result = 0;
    array<int, 3> exponents( { 9, 3, 1 });
    for (int j = 0; j < DIM; j++) {
        if (posRelative[j] > 0) {
            result += 2 * exponents[j];
        } else if (posRelative[j] == 0) {
            result += 1 * exponents[j];
        }
        // else if(posRelative[j]<0){
        //  ;
        // }
    }
    return result;
}

} // namespace sac_de_billes

#endif /* PATH_IXX_ */

//! Copyright : see license.txt
//!
//! \brief
//
#include "StdHeaders.hxx"

#include "Path.hxx"

namespace sac_de_billes {
using namespace std;

Path_auxi::Path_auxi():
        TABCORNER3D {
                array<unsigned short, 3> { 0, 0, 0 },
                array<unsigned short, 3> { 1, 1, 1 },
                array<unsigned short, 3> { 0, 0, 1 },
                array<unsigned short, 3> { 1, 1, 0 },
                array<unsigned short, 3> { 0, 1, 0 },
                array<unsigned short, 3> { 1, 0, 1 },
                array<unsigned short, 3> { 0, 1, 1 },
                array<unsigned short, 3> { 1, 0, 0 } },
        TABCORNER2D {
                array<unsigned short, 2> { 0, 0 },
                array<unsigned short, 2> { 1, 1 },
                array<unsigned short, 2> { 0, 1 },
                array<unsigned short, 2> { 1, 0 } }{
}

Path::Path() :
        pathForCorners{buildPathForCoin()},
        pathFor27Corners{build27Corners()},
        CORNER_27_DIRS{create27Directions()},
        pathNghb_3D_5{createParcours<3, 5>()},
        pathNghb_3D_7{createParcours<3, 7>()},
        pathNghb_2D_5{createParcours<2, 5>()},
        pathNghb_2D_7{createParcours<2, 7>()}{

}

inline const vector<array<unsigned short, 3> > & Path::myCornersEfficient(
        int indexPosRelative)  const{
    return pathFor27Corners[indexPosRelative];
}

/// private methods

vector<array<unsigned short, 3> > Path::myCorners(array<int, 3> posRelative) const {
    constexpr unsigned short DIM = 3;
    array<unsigned short, DIM> i = { 0, 0, 0 };
    array<int, DIM> iMin = { 1, 1, 1 };
    array<int, 3> iMax = { 1, 1, 1 };
    for (int j = 0; j < 3; j++) {
        if (-posRelative[j] <= 0) {
            iMin[j] = 0;
        }
        if (-posRelative[j] >= 0) {
            iMax[j] = 2;
        }
    }
    vector<array<unsigned short, 3>> answer = { };
    for (i[0] = iMin[0]; i[0] < iMax[0]; i[0]++) {
        for (i[1] = iMin[1]; i[1] < iMax[1]; i[1]++) {
            for (i[2] = iMin[2]; i[2] < iMax[2]; i[2]++) {
                answer.push_back(i);
            }
        }
    }

    // reordering the corners in a more efficient way (in the spirit of TABCORNER)
    if (answer.size() == 8) {
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 3; k++) {
                answer[j][k] = Path_auxi::get().TABCORNER<DIM>()[j][k];
            }
        }
    }
    if (answer.size() == 4) {
        swap(answer[1], answer[3]);
    }
    return answer;
}

vector<vector<array<unsigned short, 3>>> Path::buildPathForCoin()  {
    vector<vector<array<unsigned short, 3>>> pathForCorner = { };
    for (int i = 0; i < 5 * 5 * 5; i++) {
        pathForCorner.push_back(myCorners(pathNghb_3D_5[i]));
    }
    return pathForCorner;
}

array<array<int, 3>, 27> Path::create27Directions()  {
    array<array<int, 3>, 27> res { };
    array<int, 3> index { };
    int compteur = 0;
    loop<0, 3, 0, 3, 0, 3>([&](const array<long, 3u> &i) {
        for (int j = 0; j < 3; j++) {
            index[j] = static_cast<int>(i[j] - 1);
        }
        res[compteur] = index;
        compteur++;
    });
    return res;
}

vector<vector<array<unsigned short, 3>>> Path::build27Corners()  {
    vector<vector<array<unsigned short, 3>>> pathForCorner = { };
    for (auto relPos : CORNER_27_DIRS) {
        pathForCorner.push_back(myCorners(relPos));
    }
    return pathForCorner;
}

const Path& Path::get() {
    static Path path_{};
    return path_;
}

const Path_auxi& Path_auxi::get() {
    static Path_auxi path_{};
    return path_;
}

} // namespace sac_de_billes

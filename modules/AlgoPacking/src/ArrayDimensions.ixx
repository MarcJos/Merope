#include "ArrayDimensions.hxx"
namespace merope {

// SubArrayDimensions<DIM>
template<unsigned short DIM>
vox::SubArrayDimensions<DIM>::SubArrayDimensions(array<size_t, DIM> nbNodes_, array<size_t, DIM> nMin_, array<size_t, DIM> nMax_):
    nbNodes(nbNodes_), coverTorus{ true }, nMin{}, nMax{}, nbNodeSubgrid{} {
    this->setSubGridIndices(nMin_, nMax_);
}

template<unsigned short DIM>
vox::SubArrayDimensions<DIM>::SubArrayDimensions(array<size_t, DIM> nbNodes_):
    SubArrayDimensions(nbNodes_, create_array<DIM, size_t>(0), nbNodes_) {}

template<unsigned short DIM>
inline bool vox::SubArrayDimensions<DIM>::testCoherent() const {
    for (size_t i = 0; i < DIM; i++) {
        if (not (0 <= this->nMin[i] and this->nMin[i] < this->nMax[i] and this->nMax[i] <= this->getNbNodeBigGrid()[i])) {
            return false;
        }
    }
    return true;
}

template<unsigned short DIM>
void vox::SubArrayDimensions<DIM>::setNbNodes(array<size_t, DIM> nbNodes_) {
    nbNodes = nbNodes_;
    this->setSubGridIndices(create_array<DIM, size_t>(0), nbNodes_);
}

template<unsigned short DIM>
inline void vox::SubArrayDimensions<DIM>::setSubGridIndices(array<size_t, DIM> nMin_,
    array<size_t, DIM> nMax_) {
    nMin = nMin_;
    nMax = nMax_;
    this->setNodeSubGrid();
    this->coverTorus = this->computeCoverTorus();
    if (not testCoherent()) {
        throw invalid_argument("Incompatible grid extraction!");
    }
}

template<unsigned short DIM>
inline void vox::SubArrayDimensions<DIM>::setNodeSubGrid() {
    for (size_t i = 0; i < DIM; i++) {
        this->nbNodeSubgrid[i] = this->getNMax()[i] - this->getNMin()[i];
    }
}

template<unsigned short DIM>
inline bool vox::SubArrayDimensions<DIM>::computeCoverTorus() const {
    for (size_t i = 0; i < DIM; i++) {
        if (this->getNMin()[i] > 0 or this->getNMax()[i] < this->getNbNodeBigGrid()[i]) {
            return false;
        }
    }
    return true;
}

}// namespace merope
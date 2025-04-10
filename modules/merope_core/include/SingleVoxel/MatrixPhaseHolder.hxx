//! Copyright : see license.txt
//!
//! \brief

#pragma once

namespace merope {

//! @brief : class containing a matrix phase, which can be defined or not.
template<class MatrixPhase>
class MatrixPhaseHolder {
public:
    MatrixPhaseHolder(bool matrixPresence_, MatrixPhase matrixPhase_) : matrixPresence{ matrixPresence_ }, matrixPhase(matrixPhase_) {}
    MatrixPhaseHolder(MatrixPhase matrixPhase_) :matrixPresence{ true }, matrixPhase(matrixPhase_) {}
    MatrixPhaseHolder() :matrixPresence{ false }, matrixPhase{} {}

    ~MatrixPhaseHolder() {}
    //! set the matrix phase
    void setMatrixPhase(MatrixPhase matrixPhase_) {
        matrixPresence = true;
        this->matrixPhase = matrixPhase_;
    }
    void setMatrixPhase_if_not_present(MatrixPhase matrixPhase_) {
        if (not matrixPresence) setMatrixPhase(matrixPhase_);
    }
    //! getter
    MatrixPhase getMatrixPhase() const { return matrixPhase; }
    //! is there a matrix ?
    bool is_there_matrix() const { return matrixPresence; }
private:
    //! is there a matrix ?
    bool matrixPresence;
    //! Phase of the matrix
    MatrixPhase matrixPhase;
};

}  // namespace  merope

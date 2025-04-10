//! Copyright : see license.txt
//!
//! \briefDiscrete periodical Geometries
//
#pragma once


#include "../VTKinout/VTKStream.hxx"


namespace merope {

template<typename T>
void VTKRead::readMat(std::vector<T>& tab) {
    seekg(pMId);
    if (isFor) {
        for (typename std::vector<T>::iterator i = tab.begin(); i != tab.end();
            ++i) {
            *this >> *i;
        }
    } else {
        for (typename std::vector<T>::iterator i = tab.begin(); i != tab.end();
            ++i) {
            read((char*)&(*i), sizeof(T));
            vtkByteSwap(*i);
        }
    }
}

}  // namespace merope




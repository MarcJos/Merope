//! Copyright : see license.txt
//!
//! \brief Defines a VTK stream format
//

#ifndef _VTKSTREAM_IXX
#define _VTKSTREAM_IXX 1


#include "../MeropeNamespace.hxx"


namespace merope {

#if defined(BIGENDIAN)
//
// vtkWrite definition for BIG ENDIAN case
//

template <class DATA>
inline void vtkWrite(const DATA data, std::ostream* os) {
    os->write((char*)(&data), sizeof(DATA));
}

template <class DATA>
inline void vtkRead(DATA* data, std::istream& is) {
    is.read((char*)(data), sizeof(DATA));
}

#else
//
// vtkWrite definition for LITTLE ENDIAN case
//

// for little endian system, byte swap is needed!
template<size_t s>
struct vtkByteSwapper;

template<>
struct vtkByteSwapper<1> {
    static inline void Swap(void*) {}
};

template<>
struct vtkByteSwapper<2> {
    static inline void Swap(void* p) {
        char one_byte;
        char* data = static_cast<char*>(p);
        one_byte = data[0];
        data[0] = data[1];
        data[1] = one_byte;
    }
};

template<>
struct vtkByteSwapper<4> {
    static inline void Swap(void* p) {
        char one_byte;
        char* data = static_cast<char*>(p);
        one_byte = data[0];
        data[0] = data[3];
        data[3] = one_byte;
        one_byte = data[1];
        data[1] = data[2];
        data[2] = one_byte;
    }
};

template<>
struct vtkByteSwapper<8> {
    static inline void Swap(void* p) {
        char one_byte;
        char* data = static_cast<char*>(p);
        one_byte = data[0];
        data[0] = data[7];
        data[7] = one_byte;
        one_byte = data[1];
        data[1] = data[6];
        data[6] = one_byte;
        one_byte = data[2];
        data[2] = data[5];
        data[5] = one_byte;
        one_byte = data[3];
        data[3] = data[4];
        data[4] = one_byte;
    }
};

template<class DATA>
inline void vtkWrite(DATA data, std::ostream* os) {
    vtkByteSwapper<sizeof(DATA)>::Swap(&data);
    os->write((char*)&data, sizeof(DATA));
}

template<class DATA>
inline void vtkRead(DATA* data, std::istream& is) {
    is.read((char*)data, sizeof(DATA));
    vtkByteSwapper<sizeof(DATA)>::Swap(data);
}
#endif


template<class VECTOR_DATA, class ARRAY_DIM>
void VTKstream::writeVector(const VECTOR_DATA& vectorData, ARRAY_DIM n) {
    size_t n2 = 1;
    if (n.size() == 3) {
        n2 = n[2];
    }
    size_t i, j, k;
    //
    for (k = 0; k < n2; ++k) {
        for (j = 0; j < n[1]; ++j) {
            for (i = 0; i < n[0]; ++i) {
                write(vectorData[k + n2 * (j + i * n[1])]);
            }
        }
    }
}

// Binary writing (Big Endian)
// \param data: Data to write
template<typename DATA>
void VTKstream::write(DATA data) {
    vtkWrite(data, this);
}

// Binary reading (Big Endian)
// \param data: Data to write
template<typename DATA>
void vtkByteSwap(DATA& data) {
#if !defined(BIGENDIAN)
    vtkByteSwapper<sizeof(DATA)>::Swap(&data);
#endif
}

} // namespace merope


#endif // _VTKSTREAM_IXX

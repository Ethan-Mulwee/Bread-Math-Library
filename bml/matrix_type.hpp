#ifndef BML_MATRIX
#define BML_MATRIX

#include "vector_type.hpp"
#include <cstddef>

namespace bml {

    struct matrix2x2 {
        union {
            vector2 data[2];
            struct {
                vector2 i, j;
            };
            struct {
                vector2 x, y;
            };
        };

        inline vector2& operator[](const size_t i) {
            return data[i];
        }

        inline const vector2& operator[](const size_t i) const {
            return data[i];
        }
    };
    struct matrix3x3 {
        union {
            vector3 data[3];
            struct {
                vector3 i, j, k;
            };
            struct {
                vector3 x, y, z;
            };
        };
        
        inline vector3& operator[](const size_t i) {
            return data[i];
        }

        inline const vector3& operator[](const size_t i) const {
            return data[i];
        }
    };

    struct matrix4x4 {
        union {
            vector4 data[4];
            struct {
                vector4 i, j, k, l;
            };
            struct {
                vector4 x, y, z, w;
            };
        };

        inline vector4& operator[](const size_t i) {
            return data[i];
        }

        inline const vector4& operator[](const size_t i) const {
            return data[i];
        }
    };

}
#endif

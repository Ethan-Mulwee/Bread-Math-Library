#ifndef BML_MATRIX
#define BML_MATRIX

#include "vector_type.hpp"
#include <cstddef>

namespace bml {

    struct Matrix2x2 {
        union {
            Vector2 data[2];
            struct {
                Vector2 i, j;
            };
            struct {
                Vector2 x, y;
            };
        };

        inline Vector2& operator[](const size_t i) {
            return data[i];
        }

        inline const Vector2& operator[](const size_t i) const {
            return data[i];
        }
    };
    struct Matrix3x3 {
        union {
            Vector3 data[3];
            struct {
                Vector3 i, j, k;
            };
            struct {
                Vector3 x, y, z;
            };
        };
        
        inline Vector3& operator[](const size_t i) {
            return data[i];
        }

        inline const Vector3& operator[](const size_t i) const {
            return data[i];
        }
    };

    struct Matrix4x4 {
        union {
            Vector4 data[4];
            struct {
                Vector4 i, j, k, l;
            };
            struct {
                Vector4 x, y, z, w;
            };
        };

        inline Vector4& operator[](const size_t i) {
            return data[i];
        }

        inline const Vector4& operator[](const size_t i) const {
            return data[i];
        }
    };

}
#endif

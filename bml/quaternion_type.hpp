#ifndef BML_QUATERNION_TYPE
#define BML_QUATERNION_TYPE

#include <cmath>

#include "vector_type.hpp"

namespace bml {
    struct Quaternion {
        float x, y, z, w;

        float angle() const {
            return acos(w)*2.0f;
        }

        Vector3 axis() const {
            float angle = (*this).angle();
            float s = sin(angle/2);
            return Vector3{x/s,y/s,z/s};
        }

        float length() const {
            return sqrt(x*x+y*y+z*z+w*w);
        }

        void normalize() {
            float factor = 1.0f/(*this).length();
            x *= factor; y *= factor; z *= factor; w *= factor;
        }
    };
}

#endif
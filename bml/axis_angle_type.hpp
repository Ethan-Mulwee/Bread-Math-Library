#ifndef BML_AXIS_ANGLE_TYPE
#define BML_AXIS_ANGLE_TYPE

#include "vector_type.hpp"

namespace bml {
    struct axis_angle {
        float angle;
        vector3 axis;
    };
}

#endif
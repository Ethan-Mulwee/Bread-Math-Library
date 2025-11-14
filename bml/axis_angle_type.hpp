#ifndef BML_AXIS_ANGLE_TYPE
#define BML_AXIS_ANGLE_TYPE

#include "vector_type.hpp"

namespace bml {
    struct AxisAngle {
        float angle;
        Vector3 axis;
    };
}

#endif
#ifndef BML_EULER_TYPE
#define BML_EULER_TYPE

namespace bml {
    // TODO: figure out what to do with this
    enum euler_type{
        XYZ, XZY, YXZ, YZX, ZXY, ZYX
    };

    struct euler_xyz {
        float x, y, z;
    }; 
}

#endif
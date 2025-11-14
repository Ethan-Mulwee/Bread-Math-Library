#ifndef BML_EULER_TYPE
#define BML_EULER_TYPE

namespace bml {
    // TODO: figure out what to do with this
    enum EulerType{
        XYZ, XZY, YXZ, YZX, ZXY, ZYX
    };

    struct EulerXYZ {
        float x, y, z;
    }; 
}

#endif
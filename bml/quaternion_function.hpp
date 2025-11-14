#ifndef BML_QUATERNION_FUNCTION
#define BML_QUATERNION_FUNCTION

#include <cmath>

#include "quaternion_type.hpp"
#include "vector_function.hpp"
#include "bml_types.hpp"

namespace bml {
    // Multiply quaternions
    inline Quaternion operator*(const Quaternion &q, const Quaternion &p) {
        return Quaternion{
            q.w*p.x+q.x*p.w+q.y*p.z-q.z*p.y, // x
            q.w*p.y-q.x*p.z+q.y*p.w+q.z*p.x, // y
            q.w*p.z+q.x*p.y-q.y*p.x+q.z*p.w, // z
            q.w*p.w-q.x*p.x-q.y*p.y-q.z*p.z  // w
        };
    }

    // Add quaternions
    inline Quaternion operator+(const Quaternion &q, const Quaternion &p) {
        return Quaternion{q.x + p.x, q.y + p.y, q.z + p.z, q.w + p.w};
    }

    // Multiply quaternion by scalar
    inline Quaternion operator*(const Quaternion &q, const float s) {
        return Quaternion{q.x*s, q.y*s, q.z*s, q.w*s};
    }

    inline Quaternion normalize(const Quaternion &q) {
        float factor = 1.0f/q.length();
        return Quaternion{q.x*factor, q.y*factor, q.z*factor, q.w*factor};
    }

    // NOTE: this only works for quaternions of length 1 that represent rotations as this is the conjugate
    inline Quaternion inverse(const Quaternion &q) {
        return Quaternion{-q.x, -q.y, -q.z, q.w};
    }

    // TODO: Add rotate towards vector orientation? This one is weird 
    inline Quaternion quaternion_add_vector(const Quaternion &q, const Vector3 &v) {
        return normalize(Quaternion{
            q.x + 0.5f * (v.x * q.w + v.y * q.z - v.z * q.y),
            q.y + 0.5f * (v.y * q.w + v.z * q.x - v.x * q.z),
            q.z * 0.5f * (v.z * q.w + v.x * q.y - v.y * q.x),
            q.w + 0.5f * (-v.x * q.x - v.y * q.y - v.z * q.z),
        });
    }

    inline Vector3 quaternion_transform_vector3(const Quaternion &q, const Vector3 &v) {
        return Vector3{
            v.x*(q.x*q.x-q.y*q.y-q.z*q.z+q.w*q.w)+v.y*(2*q.x*q.y-2*q.w*q.z)+v.z*(2*q.x*q.z+2*q.w*q.y),
            v.x*(2*q.w*q.z+2*q.x*q.y)+v.y*(q.w*q.w-q.x*q.x+q.y*q.y-q.z*q.z)+v.z*(2*q.y*q.z-2*q.w*q.x),
            v.x*(2*q.x*q.z-2*q.w*q.y)+v.y*(2*q.w*q.x+2*q.y*q.z)+v.z*(q.w*q.w-q.x*q.x-q.y*q.y+q.z*q.z),
        };
    }

    // reinterpet a vector3 as a quaternion TODO: intution
    inline Quaternion quaternion_from_vector3(const Vector3 &v) {
        return Quaternion{v.x, v.y, v.z, 0};
    }

    inline Quaternion quaternion_from_axis_angle(const Vector3 &axis, const float angle) {
        Quaternion result;
        
        Vector3 axisNormalized = normalize(axis); 
        result.w = cosf(angle/2.0f);
        float s = sinf(angle/2.0f);
        result.x = s*axisNormalized.x;
        result.y = s*axisNormalized.y;
        result.z = s*axisNormalized.z;

        return result;
    }

    // Quaternion from (XYZ) euler angle
    inline Quaternion quaternion_from_euler_angles_XYZ(const float x, const float y, const float z) {
        Quaternion xRotation = quaternion_from_axis_angle(Vector3{1,0,0}, x);
        Quaternion yRotation = quaternion_from_axis_angle(Vector3{0,1,0}, y);
        Quaternion zRotation = quaternion_from_axis_angle(Vector3{0,0,1}, z);
        return xRotation * yRotation * zRotation;
    }

    // Quaternion from (ZYX) yaw, pitch, roll Euler Angles
    inline Quaternion quaternion_from_euler_angles_ZYX(const float yaw, const float pitch, const float roll) {
        Quaternion qYaw{0.0f, 0.0f, sinf(yaw/2.0f), cosf(yaw/2.0f)};
        Quaternion qPitch{0.0f, sinf(pitch/2.0f), 0.0f, cosf(pitch/2.0f)};
        Quaternion qRoll{sinf(roll/2.0f), 0.0f, 0.0f, cosf(roll/2.0f)};

        return qYaw * qPitch * qRoll;
    }

    // Adapted from https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
    inline Quaternion quaternion_from_matrix3x3(const Matrix3x3 &m) {
        Quaternion q;
        float trace = m[0][0] + m[1][1] + m[2][2];
        if( trace > 0 ) {
            float s = 0.5f / sqrtf(trace+ 1.0f);
            q.w = 0.25f / s;
            q.x = ( m[1][2] - m[2][1] ) * s;
            q.y = ( m[2][0] - m[0][2] ) * s;
            q.z = ( m[0][1] - m[1][0] ) * s;
        } 
        else {
            if ( m[0][0] > m[1][1] && m[0][0] > m[2][2] ) {
                float s = 2.0f * sqrtf( 1.0f + m[0][0] - m[1][1] - m[2][2]);
                q.w = (m[1][2] - m[2][1] ) / s;
                q.x = 0.25f * s;
                q.y = (m[1][0] + m[0][1] ) / s;
                q.z = (m[2][0] + m[0][2] ) / s;
            } 
            else if (m[1][1] > m[2][2]) {
                float s = 2.0f * sqrtf( 1.0f + m[1][1] - m[0][0] - m[2][2]);
                q.w = (m[2][0] - m[0][2] ) / s;
                q.x = (m[1][0] + m[0][1] ) / s;
                q.y = 0.25f * s;
                q.z = (m[2][1] + m[1][2] ) / s;
            } 
            else {
                float s = 2.0f * sqrtf( 1.0f + m[2][2] - m[0][0] - m[1][1] );
                q.w = (m[0][1] - m[1][0] ) / s;
                q.x = (m[2][0] + m[0][2] ) / s;
                q.y = (m[2][1] + m[1][2] ) / s;
                q.z = 0.25f * s;
            }
        }
        return q;
    }
}

#endif
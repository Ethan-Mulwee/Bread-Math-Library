#ifndef BML_VECTOR_FUNCTION
#define BML_VECTOR_FUNCTION

#include "vector_type.hpp"

namespace bml {

    /* -------------------------------------------------------------------------- */
    /*                                   vector2                                  */
    /* -------------------------------------------------------------------------- */

    inline vector2 operator+(const vector2 &a, const vector2 &b) {
        return vector2{a.x+b.x, a.y+b.y};
    }

    inline void operator+=(vector2 &a, const vector2 &b) {
        a.x += b.x;
        a.y += b.y;
    }

    inline vector2 operator-(const vector2 &a, const vector2 &b) {
        return vector2{a.x-b.x, a.y-b.y};
    }

    // Component-wise multipication of two vectors
    inline vector2 operator*(const vector2 &a, const vector2 &b) {
        return vector2{a.x*b.x, a.y*b.y};
    }

    // Scale vector by a scalar
    inline vector2 operator*(const vector2 &v, const float s) {
        return vector2{v.x*s, v.y*s};
    }

    // Scale vector by a scalar
    inline vector2 operator*(const float s, const vector2 &v) {
        return vector2{v.x*s, v.y*s};
    }


    inline vector2 operator/(const vector2 &a, const vector2 &b) {
        return vector2{a.x/b.x, a.y/b.y};
    }

    inline float dot(const vector2 &a, const vector2 &b) {
        return a.x*b.x+a.y*b.y;
    }

    /* -------------------------------------------------------------------------- */
    /*                                   vector3                                  */
    /* -------------------------------------------------------------------------- */

    inline vector3 operator+(const vector3 &a, const vector3 &b) {
        return vector3{a.x+b.x, a.y+b.y, a.z+b.z};
    }

    inline void operator+=(vector3 &a, const vector3 &b) {
        a.x += b.x;
        a.y += b.y;
        a.z += b.z;
    }

    inline vector3 operator-(const vector3 &a, const vector3 &b) {
        return vector3{a.x-b.x, a.y-b.y, a.z-b.z};
    }

    inline void operator-=(vector3 &a, const vector3 &b) {
        a.x -= b.x;
        a.y -= b.y;
        a.z -= b.z;
    }

    // Component-wise multipication of two vectors
    inline vector3 operator*(const vector3 &a, const vector3 &b) {
        return vector3{a.x*b.x, a.y*b.y, a.z*b.z};
    }

    // Scale vector by a scalar
    inline vector3 operator*(const vector3 &v, const float s) {
        return vector3{v.x*s, v.y*s, v.z*s};
    }

    inline void operator*=(vector3 &v, const float s) {
        v.x *= s; v.y *= s; v.z *= s;
    }

    // Scale vector by a scalar
    inline vector3 operator*(const float s, const vector3 &v) {
        return vector3{v.x*s, v.y*s, v.z*s};
    }

    inline vector3 operator/(const vector3 &v, const float s) {
        return vector3{v.x/s, v.y/s, v.z/s};
    }

    inline vector3 operator/(const float s, const vector3 &v) {
        return vector3{s/v.x, s/v.y, s/v.z};
    }

    inline vector3 operator/(const vector3 &a, const vector3 &b) {
        return vector3{a.x/b.x, a.y/b.y, a.z/b.z};
    }

    inline float dot(const vector3 &a, const vector3 &b) {
        return a.x*b.x+a.y*b.y+a.z*b.z;
    }

    inline vector3 cross(const vector3 &a, const vector3 &b) {
        return vector3{
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z, 
            a.x * b.y - a.y * b.x
        };
    }

    inline vector3 normalize(const vector3 &v) {
        float f = 1.0f/v.length();
        return vector3{v.x*f, v.y*f, v.z*f};
    }

    inline vector3 project(const vector3 &a, const vector3 &b) {
        return (dot(a,b)/dot(b,b))*b;
    }

    inline vector3 vector3_from_vector4(const vector4 &v) {
        return vector3{v.x, v.y, v.z};
    } 

    inline vector3 vector3_from_homogeneous_coordinate(const vector4 &v) {
        return vector3{v.x/v.w, v.y/v.w, v.z/v.w};
    }

    /* -------------------------------------------------------------------------- */
    /*                                   vector4                                  */
    /* -------------------------------------------------------------------------- */

    inline vector4 operator+(const vector4 &a, const vector4 &b) {
        return vector4{a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w};
    }

    inline void operator+=(vector4 &a, const vector4 &b) {
        a.x += b.x;
        a.y += b.y;
        a.z += b.z;
        a.w += b.w;
    }

    inline vector4 operator-(const vector4 &a, const vector4 &b) {
        return vector4{a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w};
    }

    // Component-wise multipication of two vectors
    inline vector4 operator*(const vector4 &a, const vector4 &b) {
        return vector4{a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w};
    }

    // Scale vector by a scalar
    inline vector4 operator*(const vector4 &v, const float s) {
        return vector4{v.x*s, v.y*s, v.z*s, v.w*s};
    }

    inline void operator*=(vector4 &v, const float s) {
        v.x *= s; v.y *= s; v.z *= s; v.w *= s;
    }

    inline void operator/=(vector4 &v, const float s) {
        v.x /= s; v.y /= s; v.z /= s; v.w /= s;
    }

    // Scale vector by a scalar
    inline vector4 operator*(const float s, const vector4 &v) {
        return vector4{v.x*s, v.y*s, v.z*s, v.w*s};
    }


    inline vector4 operator/(const vector4 &a, const vector4 &b) {
        return vector4{a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w};
    }

    inline float dot(const vector4 &a, const vector4 &b) {
        return a.x*b.x+a.y*b.y+a.z*b.z+a.w*b.w;
    }
    
}

#endif
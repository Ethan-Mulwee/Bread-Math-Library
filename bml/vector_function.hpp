#ifndef BML_VECTOR_FUNCTION
#define BML_VECTOR_FUNCTION

#include "vector_type.hpp"

namespace bml {

    /* -------------------------------------------------------------------------- */
    /*                                   vector2                                  */
    /* -------------------------------------------------------------------------- */

    inline Vector2 operator+(const Vector2 &a, const Vector2 &b) {
        return Vector2{a.x+b.x, a.y+b.y};
    }

    inline void operator+=(Vector2 &a, const Vector2 &b) {
        a.x += b.x;
        a.y += b.y;
    }

    inline Vector2 operator-(const Vector2 &a, const Vector2 &b) {
        return Vector2{a.x-b.x, a.y-b.y};
    }

    // Component-wise multipication of two vectors
    inline Vector2 operator*(const Vector2 &a, const Vector2 &b) {
        return Vector2{a.x*b.x, a.y*b.y};
    }

    // Scale vector by a scalar
    inline Vector2 operator*(const Vector2 &v, const float s) {
        return Vector2{v.x*s, v.y*s};
    }

    // Scale vector by a scalar
    inline Vector2 operator*(const float s, const Vector2 &v) {
        return Vector2{v.x*s, v.y*s};
    }


    inline Vector2 operator/(const Vector2 &a, const Vector2 &b) {
        return Vector2{a.x/b.x, a.y/b.y};
    }

    inline float dot(const Vector2 &a, const Vector2 &b) {
        return a.x*b.x+a.y*b.y;
    }

    /* -------------------------------------------------------------------------- */
    /*                                   vector3                                  */
    /* -------------------------------------------------------------------------- */

    inline Vector3 operator+(const Vector3 &a, const Vector3 &b) {
        return Vector3{a.x+b.x, a.y+b.y, a.z+b.z};
    }

    inline void operator+=(Vector3 &a, const Vector3 &b) {
        a.x += b.x;
        a.y += b.y;
        a.z += b.z;
    }

    inline Vector3 operator-(const Vector3 &a, const Vector3 &b) {
        return Vector3{a.x-b.x, a.y-b.y, a.z-b.z};
    }

    inline void operator-=(Vector3 &a, const Vector3 &b) {
        a.x -= b.x;
        a.y -= b.y;
        a.z -= b.z;
    }

    // Component-wise multipication of two vectors
    inline Vector3 operator*(const Vector3 &a, const Vector3 &b) {
        return Vector3{a.x*b.x, a.y*b.y, a.z*b.z};
    }

    // Scale vector by a scalar
    inline Vector3 operator*(const Vector3 &v, const float s) {
        return Vector3{v.x*s, v.y*s, v.z*s};
    }

    inline void operator*=(Vector3 &v, const float s) {
        v.x *= s; v.y *= s; v.z *= s;
    }

    // Scale vector by a scalar
    inline Vector3 operator*(const float s, const Vector3 &v) {
        return Vector3{v.x*s, v.y*s, v.z*s};
    }

    inline Vector3 operator/(const Vector3 &v, const float s) {
        return Vector3{v.x/s, v.y/s, v.z/s};
    }

    inline Vector3 operator/(const float s, const Vector3 &v) {
        return Vector3{s/v.x, s/v.y, s/v.z};
    }

    inline Vector3 operator/(const Vector3 &a, const Vector3 &b) {
        return Vector3{a.x/b.x, a.y/b.y, a.z/b.z};
    }

    inline float dot(const Vector3 &a, const Vector3 &b) {
        return a.x*b.x+a.y*b.y+a.z*b.z;
    }

    inline Vector3 cross(const Vector3 &a, const Vector3 &b) {
        return Vector3{
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z, 
            a.x * b.y - a.y * b.x
        };
    }

    inline Vector3 normalize(const Vector3 &v) {
        float f = 1.0f/v.length();
        return Vector3{v.x*f, v.y*f, v.z*f};
    }

    inline Vector3 project(const Vector3 &a, const Vector3 &b) {
        return (dot(a,b)/dot(b,b))*b;
    }

    inline Vector3 vector3_from_vector4(const Vector4 &v) {
        return Vector3{v.x, v.y, v.z};
    } 

    inline Vector3 vector3_from_homogeneous_coordinate(const Vector4 &v) {
        return Vector3{v.x/v.w, v.y/v.w, v.z/v.w};
    }

    inline Vector3 lerp(const Vector3 &a, const Vector3 &b, const float t) {
        return (1.0f-t) * a + t * b;
    }

    /* -------------------------------------------------------------------------- */
    /*                                   vector4                                  */
    /* -------------------------------------------------------------------------- */

    inline Vector4 operator+(const Vector4 &a, const Vector4 &b) {
        return Vector4{a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w};
    }

    inline void operator+=(Vector4 &a, const Vector4 &b) {
        a.x += b.x;
        a.y += b.y;
        a.z += b.z;
        a.w += b.w;
    }

    inline Vector4 operator-(const Vector4 &a, const Vector4 &b) {
        return Vector4{a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w};
    }

    // Component-wise multipication of two vectors
    inline Vector4 operator*(const Vector4 &a, const Vector4 &b) {
        return Vector4{a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w};
    }

    // Scale vector by a scalar
    inline Vector4 operator*(const Vector4 &v, const float s) {
        return Vector4{v.x*s, v.y*s, v.z*s, v.w*s};
    }

    inline void operator*=(Vector4 &v, const float s) {
        v.x *= s; v.y *= s; v.z *= s; v.w *= s;
    }

    inline void operator/=(Vector4 &v, const float s) {
        v.x /= s; v.y /= s; v.z /= s; v.w /= s;
    }

    // Scale vector by a scalar
    inline Vector4 operator*(const float s, const Vector4 &v) {
        return Vector4{v.x*s, v.y*s, v.z*s, v.w*s};
    }


    inline Vector4 operator/(const Vector4 &a, const Vector4 &b) {
        return Vector4{a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w};
    }

    inline float dot(const Vector4 &a, const Vector4 &b) {
        return a.x*b.x+a.y*b.y+a.z*b.z+a.w*b.w;
    }
    
}

#endif
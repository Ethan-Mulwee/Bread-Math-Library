#ifndef BML_MATRIX_FUNCTION_IMPL
#define BML_MATRIX_FUNCTION_IMPL

#include <cmath>

#include "smath_types.hpp"
#include "matrix_function.hpp"
#include "vector_function.hpp"

namespace bml {

    /* -------------------------------------------------------------------------- */
    /*                                  matrix2x2                                 */
    /* -------------------------------------------------------------------------- */

    inline void operator*=(Matrix2x2 &m, const float s) {
        m[0][0] *= s; m[1][0] *= s;
        m[0][1] *= s; m[1][1] *= s;
    }

    inline Matrix2x2 operator*(const Matrix2x2 &m, const float s) {
        return Matrix2x2{
            Vector2{m[0][0] * s, m[0][1] * s},
            Vector2{m[1][0] * s, m[1][1] * s},
        };
    }

    inline Matrix2x2 operator*(const float s, const Matrix2x2 &m) {
        return Matrix2x2{
            Vector2{m[0][0] * s, m[0][1] * s},
            Vector2{m[1][0] * s, m[1][1] * s},
        };
    }

    inline Matrix2x2 operator*(const Matrix2x2 &a, const Matrix2x2 &b) {
        Vector2 ihat = b.i;
        Vector2 jhat = b.j;

        ihat = a * ihat;
        jhat = a * jhat;

        return Matrix2x2{ihat, jhat};
    }

    inline Matrix2x2 inverse(const Matrix2x2 &m) {
        float factor = 1.0f/(m[0][0]*m[1][1]-m[1][0]*m[0][1]);

        return factor * Matrix2x2{
            m[1][1], -m[1][0],
            -m[0][1], m[0][0]
        };
    }

    inline Vector2 matrix2x2_transform_vector2(const Matrix2x2 &m, const Vector2 &v) {
        Vector2 x = v.x * m.i;
        Vector2 y = v.y * m.j;

        return x + y;
    }

    inline Vector2 operator*(const Matrix2x2 &m, const Vector2 &v) {
        Vector2 x = v.x * m.i;
        Vector2 y = v.y * m.j;

        return x + y;
    }

    inline Matrix2x2 matrix2x2_change_basis(const Matrix2x2 &m, const Matrix2x2 &b) {
        return inverse(b)*m*b;
    }

    inline Matrix2x2 matrix2x2_from_identity() {
        return Matrix2x2 {
            1, 0,
            0 ,1
        };
    }

    inline bool matrix2x2_is_inverse(const Matrix2x2 inverse, const Matrix2x2 matrix, float epsilon) {
        Matrix2x2 product = matrix*inverse;
        Matrix2x2 identity = matrix2x2_from_identity();

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                bool upperBound = product[i][j] < identity[i][j] + epsilon;
                bool lowerBound = product[i][j] > identity[i][j] - epsilon;

                if (!(upperBound && lowerBound))
                    return false;
            }
        }

        return true;
    }

    /* -------------------------------------------------------------------------- */
    /*                                 matrix3x3                                  */
    /* -------------------------------------------------------------------------- */
    
    inline void operator*=(Matrix3x3 &m, const float s) {
        m.i *= s; m.j *= s; m.k *= s;
    }
    
    inline Matrix3x3 operator*(const Matrix3x3 &m, const float s) {
        return Matrix3x3{m.i*s, m.j*s, m.k*s,};
    }

    inline Matrix3x3 operator*(const float s, const Matrix3x3 &m) {
        return Matrix3x3{m.i*s, m.j*s, m.k*s};
    }

    inline Matrix3x3 operator*(const Matrix3x3 &a, const Matrix3x3 &b) {

        Vector3 ihat = a * b.i;
        Vector3 jhat = a * b.j;
        Vector3 khat = a * b.k;

        return Matrix3x3{ihat, jhat, khat};
    }
    
    inline float determinant(const Matrix3x3 &m) {
        // a(ei-fh)
        float t1 = m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2]);
        // b(di-fg)
        float t2 = m[1][0]*(m[0][1]*m[2][2]-m[2][1]*m[0][2]);
        // c(dh-eg)
        float t3 = m[2][0]*(m[0][1]*m[1][2]-m[1][1]*m[0][2]);
        
        return t1-t2+t3;
    }

    inline Matrix3x3 transpose(const Matrix3x3 &m) {
        return {
            m[0][0], m[1][0], m[2][0],
            m[0][1], m[1][1], m[2][1],
            m[0][2], m[1][2], m[2][2],
        };
    }
    
    inline Matrix3x3 inverse(const Matrix3x3 &m) {
        // took this from glm
        float OneOverDeterminant = 1.0f / (
            +m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        );
            
        Matrix3x3 inverse;
        inverse[0][0] = +(m[1][1] * m[2][2] - m[1][2] * m[2][1]);
        inverse[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]);
        inverse[0][2] = +(m[0][1] * m[1][2] - m[0][2] * m[1][1]);
        inverse[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]);
        inverse[1][1] = +(m[0][0] * m[2][2] - m[0][2] * m[2][0]);
        inverse[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]);
        inverse[2][0] = +(m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        inverse[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]);
        inverse[2][2] = +(m[0][0] * m[1][1] - m[0][1] * m[1][0]);
        
        inverse *= OneOverDeterminant;
        return inverse;
    }

    inline Vector3 matrix3x3_transform_vector3(const Matrix3x3 &m, const Vector3 &v) {
        Vector3 x = v.x*m.i;
        Vector3 y = v.y*m.j;
        Vector3 z = v.z*m.k;
        return x+y+z;
    }

    inline Vector3 operator*(const Matrix3x3 &m, const Vector3 &v) {
        Vector3 x = v.x*m.i;
        Vector3 y = v.y*m.j;
        Vector3 z = v.z*m.k;
        return x+y+z;
    }

    inline Matrix3x3 matrix3x3_from_quaternion(const Quaternion &q) {
        Matrix3x3 m;

        // i hat
        m[0][0] = q.x*q.x-q.y*q.y-q.z*q.z+q.w*q.w;
        m[0][1] = 2*q.w*q.z+2*q.x*q.y;
        m[0][2] = 2*q.x*q.z-2*q.w*q.y;
        // j hat
        m[1][0] = 2*q.x*q.y-2*q.w*q.z;
        m[1][1] = q.w*q.w-q.x*q.x+q.y*q.y-q.z*q.z;
        m[1][2] = 2*q.w*q.x+2*q.y*q.z;
        // k hat
        m[2][0] = 2*q.x*q.z+2*q.w*q.y;
        m[2][1] = 2*q.y*q.z-2*q.w*q.x;
        m[2][2] = q.w*q.w-q.x*q.x-q.y*q.y+q.z*q.z;
    
        return m;
    }

    inline Matrix3x3 matrix3x3_from_ihat(const Vector3 &v) {
        Vector3 jhat, khat;

        Vector3 ihat = normalize(v);

        if (std::abs(ihat.y) < std::abs(ihat.z)) {
            khat = Vector3{-ihat.z, 0, ihat.x};
            jhat = cross(ihat, khat);
        }
        else {
            jhat = Vector3{-ihat.y, ihat.x, 0};
            khat = cross(ihat, jhat);
        }

        return Matrix3x3{ihat, jhat, khat};
    }

    inline Matrix3x3 matrix3x3_from_jhat(const Vector3 &v) {
        Vector3 ihat, khat;

        Vector3 jhat = normalize(v);

        if (std::abs(jhat.y) < std::abs(jhat.z)) {
            khat = Vector3{-jhat.z, 0, jhat.x};
            ihat = cross(khat, jhat);
        }
        else {
            ihat = Vector3{-jhat.y, jhat.x, 0};
            khat = cross(ihat, jhat);
        }

        return Matrix3x3{ihat, jhat, khat};
    }

    inline Matrix3x3 matrix3x3_from_khat(const Vector3 &v) {
        Vector3 ihat, jhat;

        Vector3 khat = normalize(v);

        if (std::abs(khat.y) < std::abs(khat.z)) {
            jhat = Vector3{-khat.z, 0, khat.x};
            ihat = cross(khat, jhat);
        }
        else {
            ihat = Vector3{-khat.y, khat.x, 0};
            jhat = cross(khat, ihat);
        }

        return Matrix3x3{ihat, jhat, khat};
    }

    // TODO:
    // matrix3x3 matrix3x3_from_euler(const euler_xyz &e) {

    // }

    // returns matrix in a different basis
    inline Matrix3x3 matrix3x3_change_basis(const Matrix3x3 &matrix, const Matrix3x3 &changeOfBasisMatrix) {
        return inverse(changeOfBasisMatrix)*matrix*changeOfBasisMatrix;
    }

    inline Matrix3x3 matrix3x3_change_basis_rotation(const Matrix3x3 &m, const Matrix3x3 &b) {
        return transpose(b)*m*b;
    }

    inline Matrix3x3 matrix3x3_from_matrix4x4(const Matrix4x4 &m) {
        return Matrix3x3{
            m[0][0], m[0][1], m[0][2],
            m[1][0], m[1][1], m[1][2],
            m[2][0], m[2][1], m[2][2]
        };
    }

    inline Matrix3x3 matrix3x3_from_identity() {
        return Matrix3x3{
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f
        };
    }

    inline Matrix3x3 matrix3x3_from_diagonal(const float s) {
        return Matrix3x3{
            s, 0.0f, 0.0f,
            0.0f, s, 0.0f,
            0.0f, 0.0f, s
        };
    }

    inline Matrix3x3 matrix3x3_normalize_basis(const Matrix3x3 &m) {
        Vector3 ihat = normalize(m.i);
        Vector3 jhat = normalize(m.j);
        Vector3 khat = normalize(m.k);
        return Matrix3x3{ihat, jhat, khat};
    }

    inline bool matrix3x3_is_inverse(const Matrix3x3 inverse, const Matrix3x3 matrix, float epsilon) {
        Matrix3x3 product = matrix*inverse;
        Matrix3x3 identity = matrix3x3_from_identity();

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                bool upperBound = product[i][j] < identity[i][j] + epsilon;
                bool lowerBound = product[i][j] > identity[i][j] - epsilon;

                if (!(upperBound && lowerBound))
                    return false;
            }
        }

        return true;
    }

    inline bool matrix3x3_is_orthonormal(const Matrix3x3 &m, float epsilon) {
        // float epsilon = 0.001f;

        // for (int i = 0; i < 3; i++) {
        //     vector3 basisVector = vector3_from_matrix3x3(m, i);
        //     bool upperBound = basisVector.length() < 1.0f + epsilon;
        //     bool lowerBound = basisVector.length() > 1.0f - epsilon;

        //     if (!(upperBound && lowerBound))
        //         return false;
        // }

        // return true;
        return matrix3x3_is_inverse(transpose(m), m);
    }

    inline bool matrix3x3_is_orthogonal(const Matrix3x3 &m, float epsilon) {

        float ij = dot(m.i,m.j);
        float jk = dot(m.j,m.k);
        float ik = dot(m.i,m.k);

        bool one = (ij < 0.0f + epsilon && ij > 0.0 - epsilon);
        bool two = (jk < 0.0f + epsilon && jk > 0.0 - epsilon);
        bool three = (ik < 0.0f + epsilon && ik > 0.0 - epsilon);

        return one && two && three;
    }

    /* -------------------------------------------------------------------------- */
    /*                                 matrix4x4                                  */
    /* -------------------------------------------------------------------------- */

    inline void operator*=(Matrix4x4 &m, const float s) {
        m.i *= s; m.j *= s; m.k *= s; m.l *= s;
    }

    inline void operator*=(Matrix4x4 &a, const Matrix4x4 &b) {
        a = a * b;
    }

    inline Matrix4x4 operator*(const float s, const Matrix4x4 &m) {
        return Matrix4x4{m.i*s, m.j*s, m.k*s, m.l*s};
    }

    inline Matrix4x4 operator*(const Matrix4x4 &m, const float s) {
        return Matrix4x4{m.i*s, m.j*s, m.k*s, m.l*s};
    }

    inline Matrix4x4 operator*(const Matrix4x4 &a, const Matrix4x4 &b)
    {

        Vector4 ihat = matrix4x4_transform_vector4(a, b.i);
        Vector4 jhat = matrix4x4_transform_vector4(a, b.j);
        Vector4 khat = matrix4x4_transform_vector4(a, b.k);
        Vector4 lhat = matrix4x4_transform_vector4(a, b.l);

        return Matrix4x4{ihat, jhat, khat, lhat};
    }

    inline float determinant(const Matrix4x4 &m) {
        return 
            -m[0][2]*m[1][1]*m[2][0]+
            m[0][1]*m[1][2]*m[2][0]+
            m[0][2]*m[1][0]*m[2][1]-
            m[0][0]*m[1][2]*m[2][1]-
            m[0][1]*m[1][0]*m[2][2]+
            m[0][0]*m[1][1]*m[2][2];
    }    

    inline Matrix4x4 transpose(const Matrix4x4 &m) {
        return {
            m[0][0], m[1][0], m[2][0], m[3][0],
            m[0][0], m[1][1], m[2][1], m[3][1],
            m[0][0], m[1][2], m[2][2], m[3][2],
            m[0][0], m[1][3], m[2][3], m[3][3],
        };
    }

    // adapted from https://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix#1148405
    inline Matrix4x4 inverse(const Matrix4x4 &m) {
        Matrix4x4 inv; 
        Matrix4x4 invOut; 
        float det;
        int i;

        inv[0][0] = m[1][1]  * m[2][2] * m[3][3] - 
                m[1][1]  * m[3][2] * m[2][3] - 
                m[1][2]  * m[2][1]  * m[3][3] + 
                m[1][2]  * m[3][1]  * m[2][3] +
                m[1][3] * m[2][1]  * m[3][2] - 
                m[1][3] * m[3][1]  * m[2][2];

        inv[0][1] = -m[0][1]  * m[2][2] * m[3][3] + 
                m[0][1]  * m[3][2] * m[2][3] + 
                m[0][2]  * m[2][1]  * m[3][3] - 
                m[0][2]  * m[3][1]  * m[2][3] - 
                m[0][3] * m[2][1]  * m[3][2] + 
                m[0][3] * m[3][1]  * m[2][2];

        inv[0][2]  = m[0][1]  * m[1][2] * m[3][3] - 
                m[0][1]  * m[3][2] * m[1][3] - 
                m[0][2]  * m[1][1] * m[3][3] + 
                m[0][2]  * m[3][1] * m[1][3] + 
                m[0][3] * m[1][1] * m[3][2] - 
                m[0][3] * m[3][1] * m[1][2] ;

        inv[0][3]  = -m[0][1]  * m[1][2] * m[2][3] + 
                m[0][1]  * m[2][2] * m[1][3] +
                m[0][2]  * m[1][1] * m[2][3] - 
                m[0][2]  * m[2][1] * m[1][3] - 
                m[0][3] * m[1][1] * m[2][2] + 
                m[0][3] * m[2][1] * m[1][2] ;

        inv[1][0] = -m[1][0]  * m[2][2] * m[3][3] + 
                m[1][0]  * m[3][2] * m[2][3] + 
                m[1][2]  * m[2][0] * m[3][3] - 
                m[1][2]  * m[3][0] * m[2][3] - 
                m[1][3] * m[2][0] * m[3][2] + 
                m[1][3] * m[3][0] * m[2][2] ;

        inv[1][1]  = m[0][0]  * m[2][2] * m[3][3] - 
                m[0][0]  * m[3][2] * m[2][3] - 
                m[0][2]  * m[2][0] * m[3][3] + 
                m[0][2]  * m[3][0] * m[2][3] + 
                m[0][3] * m[2][0] * m[3][2] - 
                m[0][3] * m[3][0] * m[2][2] ;

        inv[1][2]  = -m[0][0]  * m[1][2] * m[3][3] + 
                m[0][0]  * m[3][2] * m[1][3] + 
                m[0][2]  * m[1][0] * m[3][3] - 
                m[0][2]  * m[3][0] * m[1][3] - 
                m[0][3] * m[1][0] * m[3][2] + 
                m[0][3] * m[3][0] * m[1][2] ;

        inv[1][3]  = m[0][0]  * m[1][2] * m[2][3] - 
                m[0][0]  * m[2][2] * m[1][3] - 
                m[0][2]  * m[1][0] * m[2][3] + 
                m[0][2]  * m[2][0] * m[1][3] + 
                m[0][3] * m[1][0] * m[2][2] - 
                m[0][3] * m[2][0] * m[1][2] ;

        inv[2][0] = m[1][0]  * m[2][1] * m[3][3] - 
                m[1][0]  * m[3][1] * m[2][3] - 
                m[1][1]  * m[2][0] * m[3][3] + 
                m[1][1]  * m[3][0] * m[2][3] + 
                m[1][3] * m[2][0] * m[3][1] - 
                m[1][3] * m[3][0] * m[2][1] ;

        inv[2][1]  = -m[0][0]  * m[2][1] * m[3][3] + 
                m[0][0]  * m[3][1] * m[2][3] + 
                m[0][1]  * m[2][0] * m[3][3] - 
                m[0][1]  * m[3][0] * m[2][3] - 
                m[0][3] * m[2][0] * m[3][1] + 
                m[0][3] * m[3][0] * m[2][1] ;

        inv[2][2]  = m[0][0]  * m[1][1] * m[3][3] - 
                m[0][0]  * m[3][1] * m[1][3] - 
                m[0][1]  * m[1][0] * m[3][3] + 
                m[0][1]  * m[3][0] * m[1][3] + 
                m[0][3] * m[1][0] * m[3][1] - 
                m[0][3] * m[3][0] * m[1][1] ;

        inv[2][3]  = -m[0][0]  * m[1][1] * m[2][3] + 
                m[0][0]  * m[2][1] * m[1][3] + 
                m[0][1]  * m[1][0] * m[2][3] - 
                m[0][1]  * m[2][0] * m[1][3] - 
                m[0][3] * m[1][0] * m[2][1] + 
                m[0][3] * m[2][0] * m[1][1] ;

        inv[3][0] = -m[1][0] * m[2][1] * m[3][2] + 
                m[1][0] * m[3][1] * m[2][2] + 
                m[1][1] * m[2][0] * m[3][2] - 
                m[1][1] * m[3][0] * m[2][2] - 
                m[1][2] * m[2][0] * m[3][1] + 
                m[1][2] * m[3][0] * m[2][1] ;

        inv[3][1]  = m[0][0] * m[2][1] * m[3][2] - 
                m[0][0] * m[3][1] * m[2][2] - 
                m[0][1] * m[2][0] * m[3][2] + 
                m[0][1] * m[3][0] * m[2][2] + 
                m[0][2] * m[2][0] * m[3][1] - 
                m[0][2] * m[3][0] * m[2][1] ;

        inv[3][2]  = -m[0][0] * m[1][1] * m[3][2] + 
                m[0][0] * m[3][1] * m[1][2] + 
                m[0][1] * m[1][0] * m[3][2] - 
                m[0][1] * m[3][0] * m[1][2] - 
                m[0][2] * m[1][0] * m[3][1] + 
                m[0][2] * m[3][0] * m[1][1] ;

        inv[3][3]  = m[0][0] * m[1][1] * m[2][2] - 
                m[0][0] * m[2][1] * m[1][2] - 
                m[0][1] * m[1][0] * m[2][2] + 
                m[0][1] * m[2][0] * m[1][2] + 
                m[0][2] * m[1][0] * m[2][1] - 
                m[0][2] * m[2][0] * m[1][1] ;

        det = m[0][0] * inv[0][0] + m[1][0] * inv[0][1]  + m[2][0] * inv[0][2]  + m[3][0] * inv[0][3] ;

        det = 1.0 / det;

        for (i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                invOut[i][j] = inv[i][j] * det;

        return invOut;
    }

    inline Matrix4x4 invert_transform(const Matrix4x4 &m) {
        Matrix3x3 rotation = matrix3x3_from_matrix4x4(m);
        Matrix3x3 rotationInverse = transpose(rotation);

        Vector3 translation = m.l.xyz;
        Vector3 translationInverse = matrix3x3_transform_vector3(-1.0f*rotationInverse, translation);
    
        return matrix4x4_from_transformation(translationInverse, rotationInverse);
    }

    inline Vector3 matrix4x4_transform_vector3(const Matrix4x4 &m, const Vector3 &v) {
        Vector3 x = (v.x*m.i).xyz;
        Vector3 y = (v.y*m.j).xyz;
        Vector3 z = (v.z*m.k).xyz;
        Vector3 w = m.l.xyz;

        return x+y+z+w;
    }

    inline Vector4 matrix4x4_transform_vector4(const Matrix4x4 &m, const Vector4 &v) {
        Vector4 x = v.x*m.i;
        Vector4 y = v.y*m.j;
        Vector4 z = v.z*m.k;
        Vector4 w = v.w*m.l;

        return x+y+z+w;
    }

    inline Vector4 operator*(const Matrix4x4 &m, const Vector4 &v) {
        Vector4 x = v.x*m.i;
        Vector4 y = v.y*m.j;
        Vector4 z = v.z*m.k;
        Vector4 w = v.w*m.l;

        return x+y+z+w;
    }

    inline Matrix4x4 matrix4x4_from_perspective(float fovy, float aspect, float zNear, float zFar) {
        const float tanHalfFovy = tan(fovy/2.0f);

        Matrix4x4 m = {0};

        m[0][0] = 1.0f / (aspect * tanHalfFovy);
        m[1][1] = 1.0f / tanHalfFovy;
        m[2][2] = - (zFar + zNear) / (zFar - zNear);
        m[2][3] = -1.0f;
        m[3][2] = (-2.0f * zFar * zNear) / (zFar - zNear);
        return m;
    }

    inline Matrix4x4 matrix4x4_from_orthographic(float left, float right, float bottom, float top, float zNear, float zFar) {
        Matrix4x4 result{0};
        result[0][0] = 2.0f / (right - left);
        result[1][1] = 2.0f / (top - bottom);
        result[2][2] = - 2.0f / (zFar - zNear);
        result[3][0] = - (right + left) / (right - left);
        result[3][1] = - (top + bottom) / (top - bottom);
        result[3][2] = - (zFar + zNear) / (zFar - zNear);
        result[3][3] = 1.0f;

        return result;
    }

    inline Matrix4x4 matrix4x4_from_translation(const Vector3 &v) {
        return Matrix4x4{
            .i = {1, 0, 0, 0,},
            .j = {0, 1, 0, 0,},
            .k = {0, 0, 1, 0,},
            .l = {v.x, v.y, v.z, 1},
        };
    }

    inline Matrix4x4 matrix4x4_from_look(const Vector3 &direction, const Vector3 &center, const Vector3 &up) {
		Vector3 const f(normalize(center - direction));
		Vector3 const s(normalize(cross(f, up)));
		Vector3 const u(cross(s, f));

		Matrix4x4 result{0};
		result[0][0] = s.x;
		result[1][0] = s.y;
		result[2][0] = s.z;
		result[0][1] = u.x;
		result[1][1] = u.y;
		result[2][1] = u.z;
		result[0][2] =-f.x;
		result[1][2] =-f.y;
		result[2][2] =-f.z;
		result[3][0] =-dot(s, direction);
		result[3][1] =-dot(u, direction);
		result[3][2] = dot(f, direction);
        result[3][3] = 1.0f;
		return result;
    }

    inline Matrix4x4 matrix4x4_from_scale(const float s) {
        return Matrix4x4{
            s,    0.0f, 0.0f, 0.0f,
            0.0f, s,    0.0f, 0.0f,
            0.0f, 0.0f, s,    0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
        };
    }

    inline Matrix4x4 matrix4x4_from_scale(const Vector3 &s)
    {
        return Matrix4x4{
            s.x,  0.0f, 0.0f, 0.0f,
            0.0f, s.y,  0.0f, 0.0f,
            0.0f, 0.0f, s.z,  0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
        };
    }

    inline Matrix4x4 matrix4x4_from_rotation(const Quaternion &q) {
        return matrix4x4_from_matrix3x3(matrix3x3_from_quaternion(q));
    }

    inline Matrix4x4 matrix4x4_from_transformation(const Vector3 &translation, const Matrix3x3 &rotation) {
        return Matrix4x4 {
            rotation[0][0], rotation[0][1], rotation[0][2], 0,
            rotation[1][0], rotation[1][1], rotation[1][2], 0,
            rotation[2][0], rotation[2][1], rotation[2][2], 0,
            translation.x,   translation.y,  translation.z,  1
        };
    }

    // TODO:
    // matrix4x4 matrix4x4_from_transformation(const vector3 &translation, const matrix3x3 &rotation, const vector3 &scale) {
    //     return matrix4x4 {
    //         rotation[0][0]*scale.x,  
    //     }
    // }

    inline Matrix4x4 matrix4x4_from_matrix3x3(const Matrix3x3 &m) {
        return Matrix4x4{
            .i = {m.i.x, m.i.y, m.i.z, 0.0f},
            .j = {m.j.x, m.j.y, m.j.z, 0.0f},
            .k = {m.k.x, m.k.y, m.k.z, 0.0f},
            .l = { 0.0f,  0.0f,  0.0f, 1.0f},
        };
    }

    inline Matrix4x4 matrix4x4_from_identity() {
        return Matrix4x4{
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        };
    }

    inline Matrix4x4 matrix4x4_from_diagonal(const float s) {
        return Matrix4x4{
            s, 0, 0, 0,
            0, s, 0, 0,
            0, 0, s, 0,
            0, 0, 0, s
        };
    }

    inline bool matrix4x4_is_inverse(const Matrix4x4 inverse, const Matrix4x4 matrix, float epsilon) {
        Matrix4x4 product = matrix*inverse;
        Matrix4x4 identity = matrix4x4_from_identity();

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                bool upperBound = product[i][j] < identity[i][j] + epsilon;
                bool lowerBound = product[i][j] > identity[i][j] - epsilon;

                if (!(upperBound && lowerBound))
                    return false;
            }
        }

        return true;
    }
}

#endif
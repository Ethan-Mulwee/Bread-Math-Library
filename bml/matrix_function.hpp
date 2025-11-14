#ifndef BML_MATRIX_FUNCTION
#define BML_MATRIX_FUNCTION

#include <cmath>
#include "smath_types.hpp"
#include "vector_function.hpp"

namespace bml {

    /* -------------------------------------------------------------------------- */
    /*                                  matrix2x2                                 */
    /* -------------------------------------------------------------------------- */

    inline void operator*=(Matrix2x2 &m, const float s);

    inline Matrix2x2 operator*(const Matrix2x2 &m, const float s);

    inline Matrix2x2 operator*(const float s, const Matrix2x2 &m);

    inline Matrix2x2 operator*(const Matrix2x2 &a, const Matrix2x2 &b);
    
    inline float determinant(const Matrix2x2 &m);

    inline Matrix2x2 transpose(const Matrix2x2 &m);

    inline Matrix2x2 inverse(const Matrix2x2 &m);

    inline Vector2 matrix2x2_transform_vector2(const Matrix2x2 &m, const Vector2 &v);

    inline Vector2 operator*(const Matrix2x2 &m, const Vector2 &v);

    inline Matrix2x2 matrix2x2_change_basis(const Matrix2x2 &m, const Matrix2x2 &b);

    inline Matrix2x2 matrix2x2_from_identity();

    inline bool matrix2x2_is_inverse(const Matrix2x2 inverse, const Matrix2x2 matrix, float epsilon = 0.001f);

    /* -------------------------------------------------------------------------- */
    /*                                 matrix3x3                                  */
    /* -------------------------------------------------------------------------- */
    
    inline void operator*=(Matrix3x3 &m, const float s);
    
    inline Matrix3x3 operator*(const Matrix3x3 &m, const float s);

    inline Matrix3x3 operator*(const float s, const Matrix3x3 &m);

    inline Matrix3x3 operator*(const Matrix3x3 &a, const Matrix3x3 &b);
    
    inline float determinant(const Matrix3x3 &m);

    inline Matrix3x3 transpose(const Matrix3x3 &m);
    
    inline Matrix3x3 inverse(const Matrix3x3 &m);

    inline Vector3 matrix3x3_transform_vector3(const Matrix3x3 &m, const Vector3 &v);

    inline Vector3 operator*(const Matrix3x3 &m, const Vector3 &v);

    inline Matrix3x3 matrix3x3_from_quaternion(const Quaternion &q);

    inline Matrix3x3 matrix3x3_from_euler(const EulerXYZ &e);

    // Creates a 3x3 orthonormal matrix from ihat vector
    inline Matrix3x3 matrix3x3_from_ihat(const Vector3 &v); 

    // Creates a 3x3 orthonormal matrix from jhat vector
    inline Matrix3x3 matrix3x3_from_jhat(const Vector3 &v); 

    // Creates a 3x3 orthonormal matrix from khat vector
    inline Matrix3x3 matrix3x3_from_khat(const Vector3 &v); 

    // Simple change of basis function, NOTE: this uses a inverse which could be optimzied away as a tranpose in the case of a rotation matrix
    inline Matrix3x3 matrix3x3_change_basis(const Matrix3x3 &matrix, const Matrix3x3 &changeOfBasisMatrix);

    // for change of basis when the basis is a rotation matrix and so the inverse can be optimized to isntead be a simple tranpose
    inline Matrix3x3 matrix3x3_change_basis_rotation(const Matrix3x3 &m, const Matrix3x3 &b);

    inline Matrix3x3 matrix3x3_from_matrix4x4(const Matrix4x4 &m);

    inline Matrix3x3 matrix3x3_from_identity();

    inline Matrix3x3 matrix3x3_from_diagonal(const float s);

    inline Matrix3x3 matrix3x3_normalize_basis(const Matrix3x3 &m);

    // Checks if i is an inverse of m
    inline bool matrix3x3_is_inverse(const Matrix3x3 inverse, const Matrix3x3 matrix, float epsilon = 0.001f);

    inline bool matrix3x3_is_orthonormal(const Matrix3x3 &m, float epsilon = 0.001f); 

    inline bool matrix3x3_is_orthogonal(const Matrix3x3 &m, float epsilon = 0.001f);

    /* -------------------------------------------------------------------------- */
    /*                                 matrix4x4                                  */
    /* -------------------------------------------------------------------------- */

    inline void operator*=(Matrix4x4 &m, const float s);

    inline void operator*=(Matrix4x4 &a, const Matrix4x4 &b);
    
    inline Matrix4x4 operator*(const float s, const Matrix4x4 &m);

    inline Matrix4x4 operator*(const Matrix4x4 &m, const float s);

    inline Matrix4x4 operator*(const Matrix4x4 &a, const Matrix4x4 &b);

    inline float determinant(const Matrix4x4 &m);   

    inline Matrix4x4 transpose(const Matrix4x4 &m);

    inline Matrix4x4 inverse(const Matrix4x4 &m);
    
    // this assumes the bottom row is 0, 0, 0, 1
    // check if in the case that you have a 3x3 rot mat and a position that inverting is not much faster via transpose?
    inline Matrix4x4 invert_transform(const Matrix4x4 &m);

    inline Vector3 matrix4x4_transform_vector3(const Matrix4x4 &m, const Vector3 &v);

    inline Vector4 matrix4x4_transform_vector4(const Matrix4x4 &m, const Vector4 &v);

    inline Vector4 operator*(const Matrix4x4 &m, const Vector4 &v);

    inline Matrix4x4 matrix4x4_from_perspective(float fovy, float aspect, float zNear, float zFar);

    inline Matrix4x4 matrix4x4_from_orthographic(float left, float right, float bottom, float top, float zNear, float zFar);

    inline Matrix4x4 matrix4x4_from_translation(const Vector3 &v);

    inline Matrix4x4 matrix4x4_from_look(const Vector3 &direction, const Vector3 &center, const Vector3 &up);

    inline Matrix4x4 matrix4x4_from_scale(const float s);

    inline Matrix4x4 matrix4x4_from_scale(const Vector3 &s);

    inline Matrix4x4 matrix4x4_from_rotation(const Quaternion &q);

    inline Matrix4x4 matrix4x4_from_transformation(const Vector3 &translation, const Matrix3x3 &rotation);

    inline Matrix4x4 matrix4x4_from_transformation(const Vector3 &translation, const Matrix3x3 &rotation, const Vector3 &scale);

    inline Matrix4x4 matrix4x4_from_matrix3x3(const Matrix3x3 &m);

    inline Matrix4x4 matrix4x4_from_identity();

    inline Matrix4x4 matrix4x4_from_diagonal(const float s);

    inline bool matrix4x4_is_inverse(const Matrix4x4 inverse, const Matrix4x4 matrix, float epsilon = 0.001f);

    // TODO:
    // inline matrix4x4 matrix4x4_from_quaternion(const quaternion &q) {

    // }
}

#endif
#ifndef BML_MATRIX_FUNCTION
#define BML_MATRIX_FUNCTION

#include <cmath>
#include "smath_types.hpp"
#include "vector_function.hpp"

namespace bml {

    /* -------------------------------------------------------------------------- */
    /*                                  matrix2x2                                 */
    /* -------------------------------------------------------------------------- */

    inline void operator*=(matrix2x2 &m, const float s);

    inline matrix2x2 operator*(const matrix2x2 &m, const float s);

    inline matrix2x2 operator*(const float s, const matrix2x2 &m);

    inline matrix2x2 operator*(const matrix2x2 &a, const matrix2x2 &b);
    
    inline float determinant(const matrix2x2 &m);

    inline matrix2x2 transpose(const matrix2x2 &m);

    inline matrix2x2 inverse(const matrix2x2 &m);

    inline vector2 matrix2x2_transform_vector2(const matrix2x2 &m, const vector2 &v);

    inline vector2 operator*(const matrix2x2 &m, const vector2 &v);

    inline matrix2x2 matrix2x2_change_basis(const matrix2x2 &m, const matrix2x2 &b);

    inline matrix2x2 matrix2x2_from_identity();

    inline bool matrix2x2_is_inverse(const matrix2x2 inverse, const matrix2x2 matrix, float epsilon = 0.001f);

    /* -------------------------------------------------------------------------- */
    /*                                 matrix3x3                                  */
    /* -------------------------------------------------------------------------- */
    
    inline void operator*=(matrix3x3 &m, const float s);
    
    inline matrix3x3 operator*(const matrix3x3 &m, const float s);

    inline matrix3x3 operator*(const float s, const matrix3x3 &m);

    inline matrix3x3 operator*(const matrix3x3 &a, const matrix3x3 &b);
    
    inline float determinant(const matrix3x3 &m);

    inline matrix3x3 transpose(const matrix3x3 &m);
    
    inline matrix3x3 inverse(const matrix3x3 &m);

    inline vector3 matrix3x3_transform_vector3(const matrix3x3 &m, const vector3 &v);

    inline vector3 operator*(const matrix3x3 &m, const vector3 &v);

    inline matrix3x3 matrix3x3_from_quaternion(const quaternion &q);

    inline matrix3x3 matrix3x3_from_euler(const euler_xyz &e);

    // Creates a 3x3 orthonormal matrix from ihat vector
    inline matrix3x3 matrix3x3_from_ihat(const vector3 &v); 

    // Creates a 3x3 orthonormal matrix from jhat vector
    inline matrix3x3 matrix3x3_from_jhat(const vector3 &v); 

    // Creates a 3x3 orthonormal matrix from khat vector
    inline matrix3x3 matrix3x3_from_khat(const vector3 &v); 

    // Simple change of basis function, NOTE: this uses a inverse which could be optimzied away as a tranpose in the case of a rotation matrix
    inline matrix3x3 matrix3x3_change_basis(const matrix3x3 &matrix, const matrix3x3 &changeOfBasisMatrix);

    // for change of basis when the basis is a rotation matrix and so the inverse can be optimized to isntead be a simple tranpose
    inline matrix3x3 matrix3x3_change_basis_rotation(const matrix3x3 &m, const matrix3x3 &b);

    inline matrix3x3 matrix3x3_from_matrix4x4(const matrix4x4 &m);

    inline matrix3x3 matrix3x3_from_identity();

    inline matrix3x3 matrix3x3_from_diagonal(const float s);

    inline matrix3x3 matrix3x3_normalize_basis(const matrix3x3 &m);

    // Checks if i is an inverse of m
    inline bool matrix3x3_is_inverse(const matrix3x3 inverse, const matrix3x3 matrix, float epsilon = 0.001f);

    inline bool matrix3x3_is_orthonormal(const matrix3x3 &m, float epsilon = 0.001f); 

    inline bool matrix3x3_is_orthogonal(const matrix3x3 &m, float epsilon = 0.001f);

    /* -------------------------------------------------------------------------- */
    /*                                 matrix4x4                                  */
    /* -------------------------------------------------------------------------- */

    inline void operator*=(matrix4x4 &m, const float s);

    inline void operator*=(matrix4x4 &a, const matrix4x4 &b);
    
    inline matrix4x4 operator*(const float s, const matrix4x4 &m);

    inline matrix4x4 operator*(const matrix4x4 &m, const float s);

    inline matrix4x4 operator*(const matrix4x4 &a, const matrix4x4 &b);

    inline float determinant(const matrix4x4 &m);   

    inline matrix4x4 transpose(const matrix4x4 &m);

    inline matrix4x4 inverse(const matrix4x4 &m);
    
    // this assumes the bottom row is 0, 0, 0, 1
    // check if in the case that you have a 3x3 rot mat and a position that inverting is not much faster via transpose?
    inline matrix4x4 invert_transform(const matrix4x4 &m);

    inline vector3 matrix4x4_transform_vector3(const matrix4x4 &m, const vector3 &v);

    inline vector4 matrix4x4_transform_vector4(const matrix4x4 &m, const vector4 &v);

    inline vector4 operator*(const matrix4x4 &m, const vector4 &v);

    inline matrix4x4 matrix4x4_from_perspective(float fovy, float aspect, float zNear, float zFar);

    inline matrix4x4 matrix4x4_from_orthographic(float left, float right, float bottom, float top, float zNear, float zFar);

    inline matrix4x4 matrix4x4_from_translation(const vector3 &v);

    inline matrix4x4 matrix4x4_from_look(const vector3 &direction, const vector3 &center, const vector3 &up);

    inline matrix4x4 matrix4x4_from_scale(const float s);

    inline matrix4x4 matrix4x4_from_scale(const vector3 &s);

    inline matrix4x4 matrix4x4_from_rotation(const quaternion &q);

    inline matrix4x4 matrix4x4_from_transformation(const vector3 &translation, const matrix3x3 &rotation);

    inline matrix4x4 matrix4x4_from_transformation(const vector3 &translation, const matrix3x3 &rotation, const vector3 &scale);

    inline matrix4x4 matrix4x4_from_matrix3x3(const matrix3x3 &m);

    inline matrix4x4 matrix4x4_from_identity();

    inline matrix4x4 matrix4x4_from_diagonal(const float s);

    inline bool matrix4x4_is_inverse(const matrix4x4 inverse, const matrix4x4 matrix, float epsilon = 0.001f);

    // TODO:
    // inline matrix4x4 matrix4x4_from_quaternion(const quaternion &q) {

    // }
}

#endif
#ifndef BML_TRANSFORM_FUNCTION
#define BML_TRANSFORM_FUNCTION

#include "transform_type.hpp"
#include "vector_function.hpp"
#include "quaternion_function.hpp"

namespace bml {

    // // TODO: correctly invert scale
    // actually nah this is a nightmare and expensive
    // inline transform inverse(const transform &t) {
    //     transform inverseTransform;
    //     inverseTransform.scale = 1.0f / t.scale;
    //     inverseTransform.rotation = inverse(t.rotation);
    //     // inverseTransform.translation = -1.0f * quaternion_transform_vector(inverseTransform.rotation, inverseTransform.scale * t.translation);
    //     inverseTransform.translation = -1.0f * quaternion_transform_vector(inverseTransform.rotation, t.translation);

    //     return inverseTransform;
    // }

    inline Vector3 inverse_transform_vector(const Transform &t, const Vector3 &v) {
        Vector3 transformed_vector = v;
        transformed_vector = transformed_vector - t.translation;
        transformed_vector = quaternion_transform_vector3(inverse(t.rotation), transformed_vector);
        transformed_vector = transformed_vector / t.scale;

        return transformed_vector;
    }

    inline Vector3 transform_vector(const Transform &t, const Vector3 &v) {
        Vector3 transformed_vector = v;
        transformed_vector = transformed_vector * t.scale;
        transformed_vector = quaternion_transform_vector3(t.rotation, transformed_vector);
        transformed_vector = transformed_vector + t.translation;

        return transformed_vector;
    }

    inline Matrix4x4 matrix4x4_from_transform(const Transform &t) {
        Matrix3x3 rotation_matrix = matrix3x3_from_quaternion(t.rotation);
        Matrix3x3 scale_matrix = {
            t.scale.x, 0.0f, 0.0f,
            0.0f, t.scale.y, 0.0f,
            0.0f, 0.0f, t.scale.z
        };

        Matrix4x4 transform_matrix = matrix4x4_from_matrix3x3(rotation_matrix*scale_matrix);
        transform_matrix[3][0] = t.translation.x;
        transform_matrix[3][1] = t.translation.y;
        transform_matrix[3][2] = t.translation.z;

        return transform_matrix;
    }

    // Decompose matrix4x4 note this is very expensive and doesn't always work
    // inline transform transform_from_matrix4x4(const matrix4x4 &m) {
    //     const vector3 i = vector3_from_matrix4x4(m, 0); 
    //     const vector3 j = vector3_from_matrix4x4(m, 1); 
    //     const vector3 k = vector3_from_matrix4x4(m, 2); 

    //     vector3 scale = {i.length(), j.length(), k.length()};
    //     vector3 translation = vector3_from_matrix4x4(m,3);
    //     matrix3x3 rotationMatrix = matrix3x3_from_columns(
    //         normalized(i),
    //         normalized(j),
    //         normalized(k)
    //     );
    //     quaternion rotation = quaternion_from_matrix3x3(rotationMatrix);

    //     return transform{
    //         .translation = translation,
    //         .rotation = rotation,
    //         .scale = scale
    //     };
    // }

}

#endif
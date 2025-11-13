#include "smath.hpp"
#include "smath_iostream.hpp"

int main() {
    using namespace smath;

    std::cout << "\n-------------------------------------------------------------------------- \n"
                 "                              Transform Testing                              \n"
                 "--------------------------------------------------------------------------\n\n";

    transform t = {
        .translation = vector3{1.0f,1.0f,2.0f},
        .rotation = quaternion_from_euler_angles_XYZ(M_PI/2.0f,M_PI/3.0f,0.0f),
        .scale = vector3{1.0f,2.0f,3.0f}
    };

    // transform inverseTransform = inverse(t);

    vector3 v = {0.0f,1.0f,0.0f};

    vector3 transformedVector = transform_transform_vector3(t,v);

    std::cout << "Transform: " << t << "\n";
    std::cout << "Vector: " << v << "\n";

    std::cout << "Transformed Vector: " << transformedVector << "\n\n";

    matrix4x4 transformationMatrix = matrix4x4_from_transform(t);
    
    std::cout << "Transform to Matrix: " << transformationMatrix << "\n";
    std::cout << "Othrogonal Test: " << matrix3x3_is_orthogonal(matrix3x3_from_matrix4x4(transformationMatrix)) << "\n\n";
    std::cout << "Transformed Vector from Matrix: " << matrix4x4_transform_vector3(transformationMatrix, v) << "\n\n";


    // std::cout << "Matrix to Transform: " << transform_from_matrix4x4(transformationMatrix) << "\n";

    matrix4x4 inverseTransformationMatrix = inverse(transformationMatrix);

    
    std::cout << "Inverse Matrix: " << inverseTransformationMatrix << "\n";

    std::cout << "Inverse Test: " << matrix4x4_is_inverse(inverseTransformationMatrix, transformationMatrix) << "\n";

    std::cout << "Othrogonal Test: " << matrix3x3_is_orthogonal(matrix3x3_from_matrix4x4(inverseTransformationMatrix)) << "\n\n";
    
    // std::cout << "Inverse Matrix to Transform:" << transform_from_matrix4x4(inverseTransformationMatrix) << "\n";

    // std::cout << "Inverse Matrix Transfrom to Matrix" << matrix4x4_from_transform(transform_from_matrix4x4(inverseTransformationMatrix)) << "\n";

    // std::cout << "Inverse Transform: " << inverseTransform << "\n";

    // std::cout << "Untransformed Vector: " << transform_transform_vector3(inverseTransform1, transformedVector) << "\n";
}
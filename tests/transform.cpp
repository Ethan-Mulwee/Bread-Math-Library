#include "bml.hpp"
#include "bml_iostream.hpp"

int main() {
    using namespace bml;

    std::cout << "\n-------------------------------------------------------------------------- \n"
                 "                              Transform Testing                              \n"
                 "--------------------------------------------------------------------------\n\n";

    Transform t = {
        .translation = Vector3{1.0f,1.0f,2.0f},
        .rotation = quaternion_from_euler_angles_XYZ(M_PI/2.0f,M_PI/3.0f,0.0f),
        .scale = Vector3{1.0f,2.0f,3.0f}
    };

    // transform inverseTransform = inverse(t);

    Vector3 v = {0.0f,1.0f,0.0f};

    Vector3 transformed_vector = transform_vector(t,v);

    std::cout << "Transform: " << t << "\n";
    std::cout << "Vector: " << v << "\n";

    std::cout << "Transformed Vector: " << transformed_vector << "\n\n";

    Matrix4x4 transformation_matrix = matrix4x4_from_transform(t);
    
    std::cout << "Transform to Matrix: " << transformation_matrix << "\n";
    std::cout << "Othrogonal Test: " << matrix3x3_is_orthogonal(matrix3x3_from_matrix4x4(transformation_matrix)) << "\n\n";
    std::cout << "Transformed Vector from Matrix: " << matrix4x4_transform_vector3(transformation_matrix, v) << "\n\n";


    // std::cout << "Matrix to Transform: " << transform_from_matrix4x4(transformationMatrix) << "\n";

    Matrix4x4 inverse_transformation_matrix = inverse(transformation_matrix);

    
    std::cout << "Inverse Matrix: " << inverse_transformation_matrix << "\n";

    std::cout << "Inverse Test: " << matrix4x4_is_inverse(inverse_transformation_matrix, transformation_matrix) << "\n";

    std::cout << "Othrogonal Test: " << matrix3x3_is_orthogonal(matrix3x3_from_matrix4x4(inverse_transformation_matrix)) << "\n\n";
    
    // std::cout << "Inverse Matrix to Transform:" << transform_from_matrix4x4(inverseTransformationMatrix) << "\n";

    // std::cout << "Inverse Matrix Transfrom to Matrix" << matrix4x4_from_transform(transform_from_matrix4x4(inverseTransformationMatrix)) << "\n";

    // std::cout << "Inverse Transform: " << inverseTransform << "\n";

    // std::cout << "Untransformed Vector: " << transform_transform_vector3(inverseTransform1, transformedVector) << "\n";
}
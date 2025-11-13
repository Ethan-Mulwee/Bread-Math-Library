#include "smath.hpp"
#include "to_string.hpp"
#include "smath_iostream.hpp"
#include <iostream>

int main() {

    using namespace bml;

    std::cout << "\n--------------------------------------------------------------------------\n"
                 "                             \x1b[33m Matrix Testing \x1b[0m                             \n"
                 "--------------------------------------------------------------------------\n\n";

    std::cout << "-------------------Transformation Matrix Creation----------------------- \n \n";
    
    std::cout << "Perspective Matrix: FOVY=45 ASPECT=1 NEAR=0.1 FAR=100.0 \n";
    matrix4x4 perspectiveMatrix = matrix4x4_from_perspective(45.0f, 1.0f, 0.1f, 100.0f);
    std::cout << perspectiveMatrix << "\n \n";

    std::cout << "Translation Matrix: TRANSLATION=(1,2,3) \n";
    matrix4x4 translationMatrix = matrix4x4_from_translation(vector3{1,2,3});
    std::cout << translationMatrix << "\n \n";

    std::cout << "Othrographic Matrix: LEFT=-10 RIGHT=10 BOTTOM=-10 TOP=10 NEAR=-10 FAR=20 \n";
    matrix4x4 orthoMatrix = matrix4x4_from_orthographic(-10,10,-10,10,-10,20);
    std::cout << orthoMatrix << "\n \n";

    std::cout << "Look Matrix: DIRECTION=(normalize(-1,-1,-1)) CENTER = (0,0,0), UP = (0,1,0)"; 
    matrix4x4 lookMatrix = matrix4x4_from_look(normalize(vector3{-1.0f, -1.0f, -1.0f}), {0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f});
    std::cout << lookMatrix << "\n \n";

    std::cout << "-------------------Matrix3x3 Multiplication----------------------- \n \n";

    std::cout << "Matrix3x3 A: \n";
    matrix3x3 matrix3x3A = {
        1, 4, 3,
        4, 5, 6,
        7, 2, 9
    };
    std::cout << matrix3x3A << "\n \n";

    std::cout << "Matrix3x3 B: \n";
    matrix3x3 matrix3x3B = {
        1, 4, 7,
        2, 5, 8,
        3, 9, 9
    };
    std::cout << matrix3x3B << "\n \n";

    
    std::cout << "Matrix3x3 A * Matrix3x3 B: \n";
    std::cout << matrix3x3A*matrix3x3B << "\n \n";
    
    std::cout << "-------------------Matrix2x2 Multiplication----------------------- \n \n";

    std::cout << "Matrix2x2 A: \n";
    matrix2x2 matrix2x2A = {
        1, 2,
        7, 4
    };
    std::cout << matrix2x2A << "\n \n";

    std::cout << "Matrix2x2 B: \n";
    matrix2x2 matrix2x2B = {
        1, 1,
        2, 7
    };
    std::cout << matrix2x2B << "\n \n";

    
    std::cout << "Matrix2x2 A * Matrix2x2 B: \n";
    std::cout << matrix2x2A*matrix2x2B << "\n \n";

    std::cout << "-------------------Matrix2x2 Inverse----------------------- \n \n";

    {
        std::cout << "Matrix2x2 A: \n";
        matrix2x2 matrix2x2A = {
            1, 2,
            7, 4
        };
        std::cout << matrix2x2A << "\n \n";
        
        std::cout << "Matrix2x2 A inverse \n";
        std::cout << inverse(matrix2x2A) << "\n \n";
        std::cout << "Inverse Test: " << matrix2x2_is_inverse(inverse(matrix2x2A), matrix2x2A) << "\n\n";
    }

    std::cout << "-------------------Matrix2x2 Change of basis----------------------- \n \n";

    std::cout << "Change of basis matrix: \n";
    matrix2x2 changeOfBasisMatrix = {
        2, -1,
        1, 1
    };
    std::cout << changeOfBasisMatrix << "\n \n";

    std::cout << "Transformation matrix: \n";
    matrix2x2 transformationMatrix {
        0, -1,
        1, 0
    };
    std::cout << transformationMatrix << "\n \n";

    // TODO: this seems to be wrong at the moment
    std::cout << "Transformation matrix with changed basis: \n";
    std::cout << matrix2x2_change_basis(transformationMatrix, changeOfBasisMatrix) << "\n";

    std::cout << "-------------------Matrix3x3 Inverse----------------------- \n \n";

    {
        matrix3x3 m {
            0.6f, 0.2f, 1.0f,
            0.2f, 1.6f, 2.4f,
            1.3f, 0.2f, 0.1f
        };

        matrix3x3 i = inverse(m);

        std::cout << "Matrix: " << m << "\n";
        std::cout << "Inverse(Matrix): " << i << "\n";
        std::cout << "Inverse Test: " << matrix3x3_is_inverse(i, m) << "\n\n";
    }

    std::cout << "-------------------Matrix4x4 Inverse----------------------- \n \n";

    {
        matrix4x4 matrix = {
            0.1f, 5.2f, 1.2f, 0.6f,
            0.5f, 0.1f, 2.5f, 1.2f,
            3.4f, 1.1f, 0.6f, 0.2f,
            0.1f, 0.6f, 0.8f, 2.4f
        };

        matrix4x4 inverseMatrix = inverse(matrix);

        std::cout << "Matrix: " << matrix << "\n";
        std::cout << "Inverse(Matrix): " << inverseMatrix << "\n";
        std::cout << "Inverse Test: " << matrix4x4_is_inverse(inverseMatrix, matrix) << "\n\n";
    }

    std::cout << "-------------------Matrix4x4 Multiplication----------------------- \n \n";

    std::cout << "Matrix4x4 A:";
    matrix4x4 matrix4x4A = {
        1,2,3,4,
        5,6,7,8,
        9,10,11,12,
        13,14,15,16
    };
    std::cout << matrix4x4A << "\n \n";

    std::cout << "Matrix4x4 B:";
    matrix4x4 matrix4x4B = {
        1,5,9,13,
        2,6,10,14,
        3,7,11,15,
        4,8,12,16
    };
    std::cout << matrix4x4B << "\n \n";

    std::cout << "Matrix4x4 A * Matrix4x4 B:";
    std::cout << matrix4x4A*matrix4x4B << "\n \n";

    std::cout << "-------------------Matrix3x3 from Quaternion----------------------- \n \n";

    {
        quaternion q = {0.692537f, 0.0f, 0.0f, 0.721383f};
        std::cout << "Quaternion: " << q << "\n";
        std::cout << "Matrix3x3 from Quaternion:" << matrix3x3_from_quaternion(q) << "\n";
        std::cout << "Determinant: " << determinant(matrix3x3_from_quaternion(q)) << "\n";
        std::cout << "Orthonormal: " << matrix3x3_is_orthonormal(matrix3x3_from_quaternion(q)) << "\n\n";
    }

    std::cout << "-------------------Matrix3x3 Determinant----------------------- \n \n";

    {
        matrix3x3 m = {
            0.2f, 0.1f, 0.6f,
            0.5f, 0.6f, 0.0f,
            1.2f, 2.6f, 0.5f
        };

        std::cout << "Matrix:" << m << "\n";
        std::cout << "Determiant: " << determinant(m) << "\n\n";
    }

    std::cout << "-------------------Matrix3x3 Inverse Test----------------------- \n \n";

    {
        matrix3x3 m = {
            0.2f, 0.1f, 0.6f,
            0.5f, 0.6f, 0.0f,
            1.2f, 2.6f, 0.5f
        };

        matrix3x3 i = inverse(m);

        std::cout << "Matrix:" << m << "\n";
        std::cout << "Inverse:" << i << "\n";
        std::cout << "Product" << m*i << "\n";
        std::cout << "Inverse Test: " << matrix3x3_is_inverse(i, m) << "\n\n";

    }

    std::cout << "-------------------Matrix3x3 Orthonormal Test----------------------- \n \n";

    {
        matrix3x3 m = {
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f
        };

        std::cout << "Matrix:" << m << "\n";
        std::cout << "Orthonormal: " << matrix3x3_is_orthonormal(m) << "\n\n";
    }

    std::cout << "-------------------Matrix3x3 Orthogonal Test----------------------- \n \n";

    {
        matrix3x3 m = {
            1.0f, 0.0f, 0.0f,
            0.0f, 2.0f, 0.0f,
            0.0f, 0.0f, 1.0f
        };

        std::cout << "Matrix:" << m << "\n";
        std::cout << "Orthogonal: " << matrix3x3_is_orthogonal(m) << "\n";  
        std::cout << "Orthonormal: " << matrix3x3_is_orthonormal(m) << "\n\n";  
    }

    std::cout << "-------------------Transformation Matrix Inverse Test----------------------- \n \n";

    {
        quaternion rotation = quaternion_from_euler_angles_XYZ(M_PI/2.0f,0.0f,0.0f);
        matrix3x3 rotationMatrix = matrix3x3_from_quaternion(rotation);
        matrix3x3 scaleMatrix = matrix3x3_from_diagonal(2.0f);  
        vector3 translation = {1.0f,0.0f,2.0f};

        matrix4x4 transform = matrix4x4_from_matrix3x3(rotationMatrix*scaleMatrix);
        transform[0][3] = translation.x;
        transform[1][3] = translation.y;
        transform[2][3] = translation.z;

        vector3 testVector = {0,0,1};

        std::cout << "Vector: " << matrix4x4_transform_vector3(transform, testVector) << "\n";

        matrix4x4 inverseTransfrom = inverse(transform);

        std::cout << "Matrix: " << transform << "\n";
        std::cout << "Inverse(Matrix): " << inverseTransfrom << "\n";
        std::cout << "Inverse Test: " << matrix4x4_is_inverse(inverseTransfrom, transform) << "\n\n";
    }

    std::cout << "-------------------Memory Layout and Access----------------------- \n \n";

    {
        matrix3x3 a = {
            vector3{2.4f, 1.2f, 0.5f},
            vector3{0.1f, 0.0f, 4.6f},
            vector3{4.0f, 0.5f, 0.5f}
        };
        std::cout << a << "\n";
        std::cout << "ihat: " << a.i << "\n";
        std::cout << "jhat: " << a.j << "\n";
        std::cout << "khat: " << a.k << "\n";

        std::cout << "transpose: \n" << transpose(a) << "\n";

        matrix3x3 b = {0};
        std::cout << b << "\n";
    }

}
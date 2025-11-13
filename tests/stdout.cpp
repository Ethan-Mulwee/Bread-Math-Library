#include "smath.hpp"
#include "smath_iostream.hpp"

int main() {
    using namespace smath;

    std::cout << "\n-------------------------------------------------------------------------- \n"
                 "                            Standard Out Testing                            \n"
                 "--------------------------------------------------------------------------\n\n";

    vector2 a = {1,2};
    vector3 b = {1,2,3};
    vector3 c = {-1,2,-3};
    vector4 d = {1,2,3,4};

    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << "\n";
    std::cout << "c = " << c << "\n";
    std::cout << "d = " << d << "\n\n";

    matrix2x2 m2x2 = {
        1, 2,
        3, 4
    };

    matrix3x3 m3x3 = {
        1, 2, 3,
        4, 5, 6,
        7, 8, 9
    };

    matrix4x4 m4x4 = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16
    };

    std::cout << "m2x2: " << m2x2 << "\n";
    std::cout << "m3x3 = " << m3x3 << "\n";
    std::cout << "m4x4 = " << m4x4 << "\n";

    quaternion q = {1,2,3,4};

    std::cout << "q = " << q << "\n";

    std::cout << "\n";

}
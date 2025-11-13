#include "smath.hpp"
#include "to_string.hpp"
#include "smath_iostream.hpp"
#include <iostream>

int main() {
    using namespace bml;

    std::cout << "\n -------------------------------------------------------------------------- \n "
                 "                          \x1b[33m Quaternion Testing \x1b[0m                           \n"
                 " -------------------------------------------------------------------------- \n \n";

    {
        quaternion q = quaternion_from_axis_angle(vector3{1,0,0}, M_PI/2.0f);
        std::cout << "q: " << q << "\n";
        std::cout << "Axis: "<< q.axis() << "\n";
        std::cout << "Angle: " << q.angle() << "\n";
    }

    std::cout << "\n ---------------------------Vector transformation--------------------------- \n \n";

    {
        quaternion q = quaternion_from_axis_angle(vector3{1,0,0}, M_PI/2.0f);
        vector3 v{0.0f, 1.0f, 0.0f};
        std::cout << "Quaternion: " << q << "\n";
        std::cout << "Vector: " << v << "\n \n";

        std::cout << "Vector transformed by quaternion: " << quaternion_transform_vector(q, v) << "\n";
    }

    std::cout << "\n ---------------------------Quaterion from Euler--------------------------- \n \n";

    {
        vector3 euler = vector3{M_PI/2.0f, M_PI/2.0f, 0.0f};
        quaternion q = quaternion_from_euler_angles_XYZ(euler.x, euler.y, euler.z);
        std::cout << "Euler (XYZ): " << euler << "\n";
        std::cout << "Quaternion from Euler: " << q << "\n"; 
        std::cout << "Axis: " << q.axis() << "\n";
        std::cout << "Angle: " << q.angle() << "\n";
    }

    std::cout << "\n ---------------------------Quaternion Multiplcation--------------------------- \n \n";

    {
        quaternion q = quaternion_from_axis_angle(vector3{1.0f,0.0f,0.0f}, M_PI/2.0f);
        quaternion p = quaternion_from_axis_angle(vector3{0.0f,1.0f,0.5f}, M_PI/2.0f);
        std::cout << "q: " << q << "\n";
        std::cout << "p: " << p << "\n\n";

        std::cout << "q*p: " << q*p << "\n";
    }

    std::cout << "\n ---------------------------Quaternion Inverse--------------------------- \n \n";

    {
        quaternion q = quaternion_from_axis_angle(vector3{1.0f,0.0f,0.0f}, M_PI/2.0f);
        quaternion qInverse = inverse(q);
        std::cout << "q: " << q << "\n";
        std::cout << "inverse(q): " << qInverse << "\n";
        std::cout << "q*inverse(q): " << q*qInverse << "\n";
    }

    std::cout << "\n";
}
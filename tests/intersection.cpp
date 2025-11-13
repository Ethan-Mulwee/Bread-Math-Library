#include "smath.hpp"
#include "intersection.hpp"
#include "smath_iostream.hpp"

int main() {
    using namespace smath;

    std::cout << "\n-------------------------------------------------------------------------- \n"
                 "                            Intersection Testing                            \n"
                 "--------------------------------------------------------------------------\n\n";

    {
        vector3 ray_origin{2.6525f, 0.69536f, 1.4344f};
        vector3 ray_direction{-0.9545, -0.0597, -0.2923};
        vector3 aabb_min = {-1.0f,-1.0f,-1.0f};
        vector3 aabb_max = { 1.0f, 1.0f, 1.0f};
    
        matrix4x4 obb_matrix = matrix4x4_from_transform({
            .translation = vector3{0.0f, 0.0f, 0.0f},
            .rotation = quaternion{0.122f, 0.176f, 0.064f, 0.975f},
            .scale = vector3{1.0f, 1.0f, 1.0f}
        });
    
        float distance = 0.0f;
    
        if (ray_intersection_obb(ray_origin, ray_direction, aabb_min, aabb_max, obb_matrix, distance)) {
            std::cout << "hit \n";
            std::cout << "distance: " << distance << "\n";
            std::cout << "point: " << ray_origin + ray_direction * distance << "\n";
        } else {
          std::cout << "missed \n";  
        }
    }
}
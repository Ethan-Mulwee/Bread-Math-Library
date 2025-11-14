#include "bml.hpp"
#include "intersection.hpp"
#include "bml_iostream.hpp"

int main() {
    using namespace bml;

    std::cout << "\n-------------------------------------------------------------------------- \n"
                 "                            Intersection Testing                            \n"
                 "--------------------------------------------------------------------------\n\n";

    {
        Vector3 ray_origin{2.6525f, 0.69536f, 1.4344f};
        Vector3 ray_direction{-0.9545, -0.0597, -0.2923};
        Vector3 aabb_min = {-1.0f,-1.0f,-1.0f};
        Vector3 aabb_max = { 1.0f, 1.0f, 1.0f};
    
        Matrix4x4 obb_matrix = matrix4x4_from_transform({
            .translation = Vector3{0.0f, 0.0f, 0.0f},
            .rotation = Quaternion{0.122f, 0.176f, 0.064f, 0.975f},
            .scale = Vector3{1.0f, 1.0f, 1.0f}
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
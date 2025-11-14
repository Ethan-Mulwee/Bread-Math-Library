#ifndef BML_INTERSECTION
#define BML_INTERSECTION

#include "matrix_type.hpp"
#include "vector_function.hpp"
#include "matrix_function_impl.hpp"

// TODO: contains complex ray primitive intersection functions see: https://www.realtimerendering.com/intersections.html
namespace bml {
    // adapted from https://www.opengl-tutorial.org/miscellaneous/clicking-on-objects/picking-with-custom-ray-obb-function/
    inline bool ray_intersection_obb(Vector3 ray_origin, Vector3 ray_direction, Vector3 aabb_min, Vector3 aabb_max, Matrix4x4 transform, float& intersection_distance) {
        float t_min = 0.0f;
        float t_max = 100000.0f;

        Vector3 obb_position = transform.l.xyz;

        Vector3 delta = obb_position - ray_origin;

        {
            Vector3 x_axis = transform.i.xyz;
            float e = dot(x_axis, delta);
            float f = dot(ray_direction, x_axis);

            if (fabs(f) > 0.001f) {
                float t1 = (e+aabb_min.x)/f;
                float t2 = (e+aabb_max.x)/f;

                if (t1>t2){
                    float w=t1;t1=t2;t2=w;
                }

                if (t2 < t_max)
                    t_max = t2;
                if (t1 > t_min)
                    t_min = t1;

                if (t_max < t_min)
                    return false;
            } else {
                if(-e+aabb_min.x > 0.0f || -e+aabb_max.x < 0.0f)
                    return false;
            }
        }

        {
            Vector3 y_axis = transform.j.xyz;
            float e = dot(y_axis, delta);
            float f = dot(ray_direction, y_axis);

            if (fabs(f) > 0.001f) {
                float t1 = (e+aabb_min.y)/f;
                float t2 = (e+aabb_max.y)/f;

                if (t1>t2) {float w=t1;t1=t2;t2=w;}

                if (t2 < t_max)
                    t_max = t2;
                if (t1 > t_min)
                    t_min = t1;
                if (t_max < t_min)
                    return false;
            } else {
                if(-e+aabb_min.y > 0.0f || -e+aabb_max.y < 0.0f)
                    return false;
            }  
        }

        {
            Vector3 z_axis = transform.k.xyz;
            float e = dot(z_axis, delta);
            float f = dot(ray_direction, z_axis);

            if (fabs(f) > 0.001f) {
                float t1 = (e+aabb_min.z)/f;
                float t2 = (e+aabb_max.z)/f;

                if (t1>t2) {float w=t1;t1=t2;t2=w;}

                if (t2 < t_max)
                    t_max = t2;
                if (t1 > t_min)
                    t_min = t1;
                if (t_max < t_min)
                    return false;
            } else {
                if(-e+aabb_min.z > 0.0f || -e+aabb_max.z < 0.0f)
                    return false;
            }    
        }

        intersection_distance = t_min;
        return true;

    }
}

#endif
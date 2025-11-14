#ifndef BML_VECTOR
#define BML_VECTOR

#include <cmath>


namespace bml {

    struct Vector2 {
        union {
            struct {
                float x,y;            
            };
            float data[2];
        };
    
        inline float length() const {
            return sqrt(x*x+y*y);
        }
    
        inline float square_length() const {
            return x*x+y*y;
        }
    
        void normalize() {
            float factor = 1.0f/(*this).length();
            x *= factor; y *= factor;
        }

        inline float& operator[](const size_t i) {
            return data[i];
        }

        inline const float& operator[](const size_t i) const {
            return data[i];
        }
    };

    struct Vector3 {
        union {
            struct {
                float x,y,z;
            };
            float data[3];
        };
    
        inline float length() const {
            return sqrt(x*x+y*y+z*z);
        }
    
        inline float square_length() const {
            return x*x+y*y+z*z;
        }
    
        void normalize() {
            float factor = 1.0f/(*this).length();
            x *= factor; y *= factor; z *= factor;
        }

        inline float& operator[](const size_t i) {
            return data[i];
        }

        inline const float& operator[](const size_t i) const {
            return data[i];
        }
    };

    struct Vector4 {
        union {
            struct {
                float x, y, z, w;
            };
            float data[4];
            Vector3 xyz;
        };

        inline float length() const {
            return sqrt(x*x+y*y+z*z+w*w);
        }

        inline float square_length() const {
            return x*x+y*y+z*z+w*w;
        }

        void normalize() {
            float factor = 1.0f/(*this).length();
            x *= factor; y *= factor; z *= factor; w *= factor;
        }

        inline float& operator[](const size_t i) {
            return data[i];
        }

        inline const float& operator[](const size_t i) const {
            return data[i];
        }
    };
    
}



#endif
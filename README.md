# Bread Math Library
 A C++ header only math library. Originally existed as Simple-Math or `smath` but has been rebranded and refactored to better fit my other projects.

# Design 
Complexity is generally limited. And new features are written as I need them so some things are not yet implemented
- All types are aggregate types mean to be intialized with intalizer lists.
- Overloads are generally not used for functions with more than one parameter. 
- There is no usage of templates so far. 
- It is column major for easy use with Graphics APIs. 
- Methods are generally not used aside for a few basic operations that have no parameters
- Functions are snake case and names are somewhat verbose for the sake of code readability

Example code
```cpp
using namespace bml;

Vector3 a = {1,2,3};
Vector3 b = {2,3,4};

Vector3 c = a + b;

Matrix4x4 m = {
    .i = {1,0,0,0},
    .j = {0,1,0,0},
    .k = {0,0,1,0},
    .l = {0,0,0,1},
};

Matrix4x4 transpose = transpose(m);
Matrix4x4 inverse = inverse(m);

Vector4 i_hat = m.i;
Vector3 matrix_translation = m.l.xyz;

Vector3 transformed_vector = c * m;

Transform transform = {
    .translation = Vector3{1,2,3},
    .rotation = Quaternion{0,0,0,1},
    .scale = Vector3{1,2,1}
};
```

# Types
- Vector2, Vector3, Vector4
- Vector2i, Vector3i, Vector4i (planned)
- Vector2d, Vector3d, Vector4d (planned)
- Matrix2x2, Matrix3x3, Matrix4x4
- Quaternion
- Transform

types are float types unless otherwise denoted with a suffix such as i for integer and d for double.

# Printing & Debugging
Rich terminal printing functions are included in "bml_iostream.hpp" for easy debugging.
```cpp
#include "bml.hpp"
#include "bml_iostream.hpp"

bml::Vector3 v = {1,2,3};

int main() {
    std::cout << "vector: " << v << "\n";
}

// outputs: "vector: (1.000000, 2.000000, 3.000000)" with RGB coloring
```


# Testing
basic tests and examples for library functionally are included in `tests/` along with a makefile for easy use.


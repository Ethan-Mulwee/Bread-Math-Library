#include "bml.hpp"
#include "bml_iostream.hpp"
#include <iostream>

int main() {
    using namespace bml;

    std::cout << "\n-------------------------------------------------------------------------- \n"
                 "                              Vector Testing                              \n"
                 "--------------------------------------------------------------------------\n\n";

    Vector3 a = {1,2,3};
    Vector3 b = {3,2,1};

    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << "\n";
    std::cout << "a * b = " << a*b << "\n";

    std::cout << "a.y: " << a.y << "\n";
    std::cout << "a.data[1]: " << a.data[1] << "\n";

    std::cout << "\n";

}
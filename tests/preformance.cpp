#include <iostream>
#include <chrono>
#include <random>

#include "smath.hpp"
#include "to_string.hpp"

int main() {
    using namespace bml;

    vector3 a = {1,2,3};
    vector3 b = {4,5,6};

    int iterations = std::rand();

    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < iterations; i++) {
        a += b;
    }

    std::cout << to_string(a) << "\n";
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end -start;
    std::cout << "Elapsed time: " << elapsed.count() << "s" << "\n";
    std::cout << "Iterations: " << iterations << "\n";
}
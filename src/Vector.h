#pragma once
#include <vector>
#include <fstream>
namespace anita{
    template <typename T, unsigned int dimensions>
    class Vector{        
        T components[];
    public:
        Vector(T values[]);
    };
}

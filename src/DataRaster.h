#pragma once
#define _USE_MATH_DEFINES
#include <string>
#include <cmath>
#include <vector>
// #include <iostream>
// #include <fstream>
// #include <limits>
#include "PREM.h"

namespace anita{

    class DataRaster{
        static constexpr int DATA_ROWS = 6667;
        static constexpr int DATA_COLUMNS = DATA_ROWS;
        static constexpr int DATA_INTERVAL = 1000; // meters
        static constexpr int PROJECTION_PLANE_LAT = 71; // degrees
        
        static constexpr double Z_PLANE = sin(PROJECTION_PLANE_LAT*(M_PI/180.)) * anita::POLAR_EARTH_RADIUS;
        
        static constexpr double EPSILON = 0.001; // 1mm for numerical gradient
        static constexpr double RAYSTEP = 500.0; // meters
        static constexpr double UPPER_SURFACE_BOUND = 1.0 + (5000 / anita::POLAR_EARTH_RADIUS); // Roughly higher than Mt Vinson in terms of normalized radius
        static constexpr double LOWER_SURFACE_BOUND = 1.0 - (3000 / anita::POLAR_EARTH_RADIUS); // Roughly lower than Bentley Subglacial Trench in terms of normalized radius

        std::vector<float>* data;

    public:
        DataRaster(std::string filePath);
        ~DataRaster();
    };

    // Template code commented out due to memory allocation bug
    // Explicit implementation above

    // template <typename T> class DataRaster{
    //     T* values;
    //     unsigned int columns, rows;

    //     std::streamsize getSizeOfFile(const std::string &filePath){
    //         std::ifstream dataFile;
    //         dataFile.open(filePath, std::ifstream::in | std::ifstream::binary);            
    //         if(!dataFile.is_open()){
    //             std::cout << "Error opening data file: " << filePath << std::endl;
    //             return 0;
    //         }
    //         dataFile.ignore(std::numeric_limits<std::streamsize>::max());
    //         std::streamsize length = dataFile.gcount();
    //         return length;
    //     }

    //     public:
    //         DataRaster(const std::string &filePath){
    //             // std::streamsize byteCount = getSizeOfFile(filePath);
    //             // values = new T[byteCount/sizeof(T)];
    //             values = new T[6667*6667];
    //             // std::cout << byteCount/sizeof(T) << std::endl;
    //             std::ifstream dataFile(filePath, std::ifstream::in | std::ios::binary);
    //             dataFile.read(reinterpret_cast<char*>(values), 6667*6667*sizeof(T));
    //             std::cout << values[100] << std::endl;
    //         }
    //         ~DataRaster(){
    //             delete [] values;
    //         }
    //     };
}
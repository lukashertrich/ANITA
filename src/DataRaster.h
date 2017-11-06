#pragma once
#define _USE_MATH_DEFINES
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include "PREM.h"

namespace anita{

    template <typename T> class DataRaster{
        T* values = nullptr;
        unsigned int columns, rows;
        const int dataInterval = 1000;

        std::streamsize getSizeOfFile(const std::string &filePath){
            std::ifstream dataFile;
            dataFile.open(filePath, std::ifstream::in | std::ifstream::binary);            
            if(!dataFile.is_open()){
                std::cout << "Error opening data file: " << filePath << std::endl;
                return 0;
            }
            dataFile.ignore(std::numeric_limits<std::streamsize>::max());
            std::streamsize length = dataFile.gcount();
            return length;
        }

        public:
            void importData(const std::string& filePath){
                std::streamsize byteCount = getSizeOfFile(filePath);
                values = new T[byteCount/sizeof(T)];
                std::ifstream dataFile(filePath, std::ifstream::in | std::ios::binary);
                dataFile.read(reinterpret_cast<char*>(values), byteCount);
            }
            
            ~DataRaster(){
                if(values){
                    delete [] values;
                }
            }

            float getvalue(double x, double y){
                return 0.;
            }
        };
}
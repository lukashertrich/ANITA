#pragma once
#define _USE_MATH_DEFINES
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include "PREM.h"
#include "BEDMAP.h"
#include "Vector2.h"

namespace anita{

    template <typename T>
    class DataRaster{
        T* values = nullptr;

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
                std::ifstream dataFile(filePath, std::ifstream::in | std::ifstream::binary);
                dataFile.read(reinterpret_cast<char*>(values), byteCount);
            }
            
            DataRaster(){}

            ~DataRaster(){
                if(values){
                    delete[] values;
                }
            }

            std::vector<T> getCellValues(const anita::Vector2<double>& dataCoords) const {	
                int xFloor = (int16_t)std::min(std::max(floor(dataCoords.x/anita::DATA_INTERVAL), 0.), anita::DATA_COLUMNS - 1.);
                int yFloor = (int16_t)std::min(std::max(floor(dataCoords.y/anita::DATA_INTERVAL), 0.), anita::DATA_ROWS - 1.);
                std::vector<T> cellValues{
                    values[yFloor*(anita::DATA_COLUMNS - 1) + xFloor],
                    values[yFloor*(anita::DATA_COLUMNS - 1) + std::min(xFloor + 1, anita::DATA_COLUMNS - 1)],
                    values[std::min(yFloor + 1, anita::DATA_ROWS - 1)*(anita::DATA_COLUMNS - 1) + xFloor],
                    values[std::min(yFloor + 1, anita::DATA_ROWS - 1)*(anita::DATA_COLUMNS - 1) + std::min(xFloor + 1, anita::DATA_COLUMNS - 1)],
                };
                return cellValues;
            }

        };
}
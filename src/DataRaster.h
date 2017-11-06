#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <limits>

namespace anita{
    template <typename T> class DataRaster{
        T* values;
        unsigned int columns, rows;

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
            DataRaster(const std::string &filePath){
                // std::streamsize byteCount = getSizeOfFile(filePath);
                // values = new T[byteCount/sizeof(T)];
                values = new T[6667*6667];
                // std::cout << byteCount/sizeof(T) << std::endl;
                std::ifstream dataFile(filePath, std::ifstream::in | std::ios::binary);
                dataFile.read(reinterpret_cast<char*>(values), 6667*6667*sizeof(T));
                std::cout << values[100] << std::endl;
            }
            ~DataRaster(){
                delete [] values;
            }
        };
}
#pragma once
#include <string>
#include <fstream>
#include <limits>

namespace anita{
    template <typename T> class DataRaster{
        T values[];
        unsigned int columns, rows;

        std::streamsize getSizeOfFile(const std::string &filePath){
            std::ifstream dataFile(filePath, std::ifstream::in | std::ifstream::binary);
            if(!dataFile.is_open()){
                return 0;
            }
            dataFile.ignore(std::numeric_limits<std::streamsize>::max());
            std::streamsize length = dataFile.gcount();
            return length;
        }

        public:
            DataRaster(const std::string &filePath){
                std::streamsize byteCount = getSizeOfFile(filePath);
                values = new T[byteCount];
                std::ifstream dataFile(filePath, std::ios::binary);
                dataFile.read(reinterpret_cast<char*>(values), byteCount);
            }
            ~DataRaster(){
                delete[] values;
            }
        };
}
#pragma once
#include <string>
#include <fstream>
#include <limits>

namespace anita{
    template <typename T> class DataRaster{
        T values[];
        unsigned int columns, rows;

        std::streamsize getSizeOfFile(const std::string &filePath);

        public:
            DataRaster(const std::string &filePath);
            ~DataRaster();
        };
}
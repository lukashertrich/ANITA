#include "DataRaster.h"

std::streamsize anita::DataRaster::getSizeOfFile(const std::string &filePath){
    std::ifstream dataFile(filePath, std::ifstream::in | std::ifstream::binary);
    if(!dataFile.is_open()){
        return 0;
    }
    dataFile.ignore(std::numeric_limits<std::streamsize>::max());
    std::streamsize length = dataFile.gcount();
    return length;
}

anita::DataRaster::DataRaster(const std::string &filePath){    
    std::streamsize dataCount = anita::DataRaster::getSizeOfFile(filePath);
    values = new T[dataCount];
    std::ifstream dataFile(filePath, std::ios::binary);
    dataFile.read(reinterpret_cast<char*>(values), dataCount*sizeof(T));
}

anita::DataRaster::~DataRaster(){
    delete values;
}
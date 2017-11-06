// #include "DataRaster.h"
// #include <iostream>
// #include <fstream>
// #include <string>

// anita::DataRaster::DataRaster(std::string filePath){
// 	data = new std::vector<float>(anita::DataRaster::DATA_ROWS*anita::DataRaster::DATA_COLUMNS);
// 	std::ifstream dataFile(filePath, std::ios::binary | std::ios::in);
// 	if(!dataFile){
// 		std::cout << "Error reading: " << filePath << std::endl;
// 		return;
// 	}
// 	dataFile.read(reinterpret_cast<char*>(data->data()), data->size()*sizeof(float));
// 	std::cout << (*data)[0] << std::endl;
// 	dataFile.close();
// }

// anita::DataRaster::~DataRaster(){
// 	delete data;
// }
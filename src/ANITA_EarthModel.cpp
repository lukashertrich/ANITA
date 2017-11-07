/*
 * 		ANITA Earth Model
 */

#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <thread>
#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>

// Several constants, classes, and non-member functions
// are separated into their own headers to reduce clutter
// Everything lives under the 'anita' namespace

#include "Vector2.h"
#include "Vector3.h"
#include "DataRaster.h"
#include "Diagnostics.h"
#include "QuadraticSolver.h"
#include "ESS.h"
#include "Filepaths.h"
#include "ReverseFloat.h"
#include "BEDMAP.h"
#include "Raycasting.h"


// DataRasters
anita::DataRaster<float> geoidDataRaster;
anita::DataRaster<float> bedDataRaster;
anita::DataRaster<float> surfaceDataRaster;
anita::DataRaster<float> iceThicknessDataRaster;   

void importData(anita::DataRaster<float>& dataRaster, std::string filePath){
	dataRaster.importData(filePath);
}

/*
 * 		MAIN PROGRAM
 */

int main(int argc, char **argv)
{
	// loadData();
	// std::cout << surfaceData[0] << std::endl;
	
	// setDataRaster(filePathGeoid, geoidDataRaster);
	std::thread t1(importData, std::ref(geoidDataRaster), std::ref(anita::filePathGeoid));
	std::thread t2(importData, std::ref(bedDataRaster), std::ref(anita::filePathBed));
	std::thread t3(importData, std::ref(surfaceDataRaster), std::ref(anita::filePathSurface));
	std::thread t4(importData, std::ref(iceThicknessDataRaster), std::ref(anita::filePathIceThickness));
	t1.join();
	t2.join();
	t3.join();
	t4.join();
	std::cout << "Done..." << std::endl;
	std::cout << "Press any key to close." << std::endl;
	std::cin.get(); // Wait for user input to terminate
	return 0;
}

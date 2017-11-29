/*
 * 		ANITA Earth Model
 */

// Build from terminal in a project directory containing src folder, not inside src itself.
// g++ -m64 -Wall -std=c++11 -pthread -march=native src/Raycasting.cpp src/Diagnostics.cpp src/QuadraticSolver.cpp src/Neutrino.cpp src/ANITA_EarthModel.cpp -o bin/ANITA_EarthModel

#include <iostream>
#include <string>
#include <thread>
#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>

/*
 * Several constants, classes, and non-member functions
 * are separated into their own headers to reduce clutter
 * Everything lives under the 'anita' namespace
 * 
 * Rays are cast backwards along path of neutrino from potential point of interaction
 */

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

using namespace anita;

// DataRasters
DataRaster<float> geoidDataRaster;
DataRaster<float> bedDataRaster;
DataRaster<float> surfaceDataRaster;
DataRaster<float> iceThicknessDataRaster;   

void importData(DataRaster<float>& dataRaster, const std::string& filePath){
	dataRaster.importData(filePath);
}

// Loading data simultaneously with threads may be faster
void initialize(){
	std::thread t1(importData, std::ref(geoidDataRaster), std::cref(filePathGeoid));
	std::thread t2(importData, std::ref(bedDataRaster), std::cref(filePathBed));
	std::thread t3(importData, std::ref(surfaceDataRaster), std::cref(filePathSurface));
	std::thread t4(importData, std::ref(iceThicknessDataRaster), std::cref(filePathIceThickness));
	t1.join();
	t2.join();
	t3.join();
	t4.join();
}

/*
* 		MAIN PROGRAM
*/

int main(int argc, char **argv)
{	
	std::cout << "Loading BEDMAP 2 data into memory..." << std::endl;
	initialize();
	//std::cout << getDataValue(Vector3<double>(5000., -30000., -POLAR_EARTH_RADIUS), surfaceDataRaster)<< std::endl;
	//testInteractionLength();
	// testFluxMapOutput();
	testAngularTrace();
	std::cout << "Done..." << std::endl;
	std::cout << "Press any key to close." << std::endl;
	std::cin.get();
	return 0;
}


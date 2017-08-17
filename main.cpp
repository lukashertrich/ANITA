#include <iostream>

/*
 * 		GLOBAL CONSTANTS
 */

const double DISTANCE_TO_HORIZON = 600.0e3; // meters
const double MEAN_EARTH_RADIUS = 6371.0e3; // meters
const double POLAR_EARTH_RADIUS = 6356.775e3; // meters
const double MEAN_EARTH_DENSITY = 5.515e3; // kg per cubic meter

const double INDEX_OF_REFRACTION_ICE = 1.78;

const double SPEED_OF_LIGHT = 2.9979e8; // meters per second

// Coefficients for a polynomial fit to European Spallation Source over 1e17 to 1e21 eV.
const double COEFFICIENTS[6] = { 1.661142492611783e+04,
	-4.616049469872646e+03,
	5.104845489878782e+02,
	-2.812808939782252e+01,
	7.727397573928863e-01,
	-8.473959043935221e-03 };

/*
 * 		DIAGNOSTIC FUNCTIONS
 */

void printConstants(){
	std::cout << "Enumeration of physical constants used in code:" << std::endl;
	std::cout << "Distance to horizon:			" << DISTANCE_TO_HORIZON << " meters." << std::endl;
	std::cout << "Mean Earth radius:			" << MEAN_EARTH_RADIUS << " meters." << std::endl;
	std::cout << "Polar Earth radius:			" << POLAR_EARTH_RADIUS << " meters."<< std::endl;
	std::cout << "Mean Earth density:			" << MEAN_EARTH_DENSITY << " kilograms per cubic meter." << std::endl;
	std::cout << "Index of refraction of ice:		" << INDEX_OF_REFRACTION_ICE << std::endl;
	std::cout << "Vacuum speed of light:			" << SPEED_OF_LIGHT << " meters per second." << std::endl;
}

/*
 * 		MAIN PROGRAM
 */

int main(int argc, char **argv)
{
	printConstants();
	std::cout << "Done..." << std::endl;
	return 0;
}

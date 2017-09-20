#include <iostream>
#include <string>
#include <thread>
#include <random>
#include <vector>

/*
 * 		GLOBAL CONSTANTS
 */

const double DISTANCE_TO_HORIZON = 6.000e5; // meters
const double MEAN_EARTH_RADIUS = 6.3710e6; // meters
const double POLAR_EARTH_RADIUS = 6.356775e6; // meters
const double MEAN_EARTH_DENSITY = 5.515e3; // kg per cubic meter

const double INDEX_OF_REFRACTION_ICE = 1.78;

const double SPEED_OF_LIGHT = 2.9979e8; // meters per second

// Coefficients for a polynomial fit to European Spallation Source from 1e17 to 1e21 eV.
const double ESS_COEFFICIENTS[6] = {
	1.661142492611783e4,
	-4.616049469872646e3,
	5.104845489878782e2,
	-2.812808939782252e1,
	7.727397573928863e-1,
	-8.473959043935221e-3 };

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
	std::cout << std::endl;
}

void printUsage(){
	std::cout << "Usage: anita_earthmodel2 [Nevt][ spectrum='ESS' or 'Enu'][Enu(eV)][maxdepth(m)][cross section factor]"<< std::endl;
	std::cout << "If single energy 'Enu' then argv[3] is the neutrino energy, else ignored."<< std::endl;
	std::cout << std::endl;
}

/*
 * 		FUNCTIONS
 */

double getFirnDensity(double depth){
	return 0;
}

std::vector<double> solveQuadratic(double a, double b, double c){
	std::vector<double> solutions;
	double discriminant = b * b - 4 * a * c;
	if(discriminant < 0){
		return solutions;
	}
	solutions.push_back((-b+sqrt(discriminant))/ 2*a);
	solutions.push_back((-b+sqrt(discriminant))/ 2*a);
}

// Use piecewise analytical description of radially dependent density provided by Preliminary Reference Earth Model (PREM)
double getDensity(double radius, double angle){
	double normalizedRadius = radius / MEAN_EARTH_RADIUS;
	std::vector<double> intersections;
	bool testForIntersection = true;
	
	return 0;
}

/*
 * 		MAIN PROGRAM
 */

int main(int argc, char **argv)
{
	printConstants();
	if(argc < 3){
		printUsage();
	}
	else{
		const int Nevt = atoi(argv[1]);;
		bool ESS = false;
		
		std::string argv2 = argv[2];	// Create a string from command line argument to easily compare against text 
		if(argv2 == "ESS"){
			ESS = true;
			std::cout << "ESS specified." << std::endl;
			std::cout << std::endl;
		}
		double Enu = atof(argv[3]);
		double maxdepth = atof(argv[4]);
		double crossSectionFactor = atof(argv[5]);
		
	}
	std::cout << "Done..." << std::endl;
	return 0;
}

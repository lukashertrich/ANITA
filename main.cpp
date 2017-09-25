#include <iostream>
#include <string>
#include <thread>
#include <random>
#include <vector>

/*
 * 		GLOBAL CONSTANTS
 */

// Equatorial and polar Earth radii are defined by the WGS84 ellipsoid
const double DISTANCE_TO_HORIZON = 6.000e5; // meters
const double MEAN_EARTH_RADIUS = 6.3710e6; // meters
const double EQUATORIAL_EARTH_RADIUS = 6378137; // meters
const double EQUATORIAL_EARTH_RADIUS_SQR = EQUATORIAL_EARTH_RADIUS * EQUATORIAL_EARTH_RADIUS; // meters squared
const double INVERSE_FLATTENING = 298.257223563; // Used to determine polar radius
const double POLAR_EARTH_RADIUS = EQUATORIAL_EARTH_RADIUS * (1.0 - (1.0 / INVERSE_FLATTENING)); // meters
const double POLAR_EARTH_RADIUS_SQR = POLAR_EARTH_RADIUS * POLAR_EARTH_RADIUS; // meters squared
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

// Normalized radii of density profile shells listed from the surface inwards to the core. Obtained from PREM
const std::vector<double> DENSITY_PROFILE_RADII{
	1.0,																				// Ocean (Exclude ocean from antarctic side of ray traversal by using crust density up to sea level)
	6368000.0 / MEAN_EARTH_RADIUS,		// Crust
	6356000.0 / MEAN_EARTH_RADIUS,		// Crust
	6346600.0 / MEAN_EARTH_RADIUS,		// LVZ/LID
	6151000.0 / MEAN_EARTH_RADIUS,		// Transition Zone
	5971000.0 / MEAN_EARTH_RADIUS,		// Transition Zone
	5771000.0 / MEAN_EARTH_RADIUS,		// Transition Zone
	5701000.0 / MEAN_EARTH_RADIUS,		// Lower Mantle
	3480000.0 / MEAN_EARTH_RADIUS,		// Outer Core
	1221500.0 / MEAN_EARTH_RADIUS		// Inner Core
} ;

// x^0, x^1, x^2, x^3, etc.
const std::vector<std::vector<double>> DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS{	// Given in kg / m^3
	std::vector<double>{1020.0},															// Ocean
	std::vector<double>{2600.0},															// Crust
	std::vector<double>{2900.0},															// Crust
	std::vector<double>{2691.0, 692.4},												// LVZ/LID
	std::vector<double>{7108.9, -3804.5},											// Transition Zone
	std::vector<double>{11249.4, -8029.8},										// Transition Zone
	std::vector<double>{5319.7, -1483.6},											// Transition Zone
	std::vector<double>{7956.5, -6476.1, 5528.3, -3080.7}, 	// Lower Mantle
	std::vector<double>{12581.5, -1263.8, -3642.6, -5528.1},	// Outer Core
	std::vector<double>{13088.5, -8838.1}										// Inner Core
	};

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
	solutions.push_back(((-b+sqrt(discriminant))/ (2*a)));
	solutions.push_back(((-b+sqrt(discriminant))/ (2*a)));
}

std::vector<std::vector<double>> getIntersections(double a, double b, double c){
	std::vector<double> intersectionsIngoing;
	std::vector<double> intersectionsOutgoing;
	
	for(unsigned int i = 0; i < DENSITY_PROFILE_RADII.size(); i++){
		std::vector<double> solutions = solveQuadratic(a, b, c - (DENSITY_PROFILE_RADII[i] * DENSITY_PROFILE_RADII[i]));
		if(solutions.size() == 0){
			break;
		}
		intersectionsIngoing.push_back(solutions[0]);
		intersectionsOutgoing.push_back(solutions[1]);
	}
	return std::vector<std::vector<double>>{intersectionsIngoing, intersectionsOutgoing};
}

// Provide initial ray position and direction in order to calculate integral of traversed density across PREM profiles
double getDensityTraversed(std::vector<double> position, std::vector<double> direction){
	double x = position[0];
	double y = position[1];
	double z = position[2];
	double dirX = direction[0];
	double dirY = direction[1];
	double dirZ = direction[2];
	
	// Coefficients of quadratic equation to solve against squared PREM profile radii
	double a = (dirX * dirX + dirY * dirY) / EQUATORIAL_EARTH_RADIUS_SQR + (dirZ * dirZ) / POLAR_EARTH_RADIUS_SQR;
	double b = 2.0 * ((x*dirX + y*dirY) / EQUATORIAL_EARTH_RADIUS_SQR + (z*dirZ) / POLAR_EARTH_RADIUS_SQR);
	double c =  (x*x + y*y) / EQUATORIAL_EARTH_RADIUS_SQR + (z*z) / POLAR_EARTH_RADIUS_SQR; // subtract respective profile radius squared to complete coefficient
	
	// LEFT OFF HERE: Get interesctions <=================
	
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

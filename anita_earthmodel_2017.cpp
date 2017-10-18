/*
 * 		ANITA Earth Model
 */

#include <iostream>
#include <string>
#include <thread>
#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>

/*
 * 		GLOBAL CONSTANTS
 */

// Data filepaths
const std::string filePathBedElev = "bedmap2_bed.txt";
const std::string filePathIceElev = "bedmap2_surface.txt";
const std::string filePathGeoidToWGS84 = "gl04c_geoid_to_wgs84.txt"; // Note that BEDMAP 2 provides 'geiod' spelling. File must be renamed.

// Equatorial and polar Earth radii are defined by the WGS84 ellipsoid, and are assumed to define sea level
// Antarctic data will be mapped on top of crust extruded through sea level
const double DISTANCE_TO_HORIZON = 6.000e5; // meters
const double MEAN_EARTH_RADIUS = 6.3710e6; // meters
const double EQUATORIAL_EARTH_RADIUS = 6378137; // meters
const double EQUATORIAL_EARTH_RADIUS_SQR = EQUATORIAL_EARTH_RADIUS * EQUATORIAL_EARTH_RADIUS; // meters squared
const double INVERSE_FLATTENING = 298.257223563; // Used to determine polar radius
// WGS84 polar radius is defined by inverse flattening term
const double POLAR_EARTH_RADIUS = EQUATORIAL_EARTH_RADIUS * (1.0 - (1.0 / INVERSE_FLATTENING)); // meters
const double POLAR_EARTH_RADIUS_SQR = POLAR_EARTH_RADIUS * POLAR_EARTH_RADIUS; // meters squared
const double MEAN_EARTH_DENSITY = 5.515e3; // kg per cubic meter

const double SPEED_OF_LIGHT = 2.9979e8; // meters per second

// Coefficients for a polynomial fit to European Spallation Source from 1e17 to 1e21 eV.
const double ESS_COEFFICIENTS[6] = {
	1.661142492611783e4,
	-4.616049469872646e3,
	5.104845489878782e2,
	-2.812808939782252e1,
	7.727397573928863e-1,
	-8.473959043935221e-3 };

// Normalized radii of density profile shells listed from the surface inwards to the core.
// Obtained from spherical PREM to be mapped to normalized WGS84 ellipsoid equation
const std::vector<double> DENSITY_PROFILE_RADII{
	1.0,																				// Ocean (Exclude ocean from antarctic side of ray traversal by using constant crust density extruded  up to sea level)
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

// x^0, x^1, x^2, x^3, etc. from PREM
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
 *		GLOBAL VARIABLES
 */

// Serial packing of rows
std::vector<double> geoidToWGS84;
std::vector<double> bedRadNormSqr; 
std::vector<double> iceRadNormSqr;

/*
 * 		FUNCTIONS
 */

void precomputeRadNormSqr(){
	std::ifstream geoidToWGS84File(filePathGeoidToWGS84);
	std::ifstream bedElevFile(filePathBedElev);
	std::ifstream iceElevFile(filePathIceElev);
	// std::vector<double> bedRadNormSqr{
	// 	std::istream_iterator<double>(bedElevFile),
	// 	std::istream_iterator<double>(),
	// 	std::back_inserter(bedRadNormSqr)};
}

void normalizeVector(std::vector<double> v){
	double magsqr = 0;
	for(int i = 0; i < v.size(); i++){
		magsqr += v[i] * v[i];
	}
	if(magsqr == 0){
		return;
	}
	magsqr = sqrt(magsqr);
	for(int i = 0; i < v.size(); i++){
		v[i] /= magsqr;
	}
}

double getFirnDensity(double depth){
	//TODO: need reference materials to obtain good function. Does exponential packing suffice?
	return 0;
}

// Supply coefficients  in ax^2+bx+c form to obtain real solutions
std::vector<double> solveQuadratic(double a, double b, double c){
	std::vector<double> solutions;
	double discriminant = b * b - 4 * a * c;
//	std::cout << "Discriminant: " << discriminant << std::endl;
	if(discriminant < 0){ // No real solutions found
		return solutions; // Empty vector, size() == 0
	}
	// Solve real roots in a numerically stable way
	if(b < 0){
		solutions.push_back((2.0*c) / (-b + discriminant));
		solutions.push_back((-b+discriminant) / (2.0*a));
	}
	else{
		solutions.push_back((-b-discriminant) / (2.0*a));
		solutions.push_back((2.0*c) / (-b - discriminant));
	}
	
//	std::cout << "Solutions: " <<   solutions[0] << ", "<< solutions[1] << std::endl;
	return solutions;
}

std::vector<std::vector<double>> getIntersections(double a, double b, double c){
	std::vector<double> intersectionsIngoing;
	std::vector<double> intersectionsOutgoing;
	// Test for PREM density shell intersection from surface inwards to core
	// Loop will break out as soon as no intersection is found, meaning current and deeper shells aren't traversed.
	for(int i = 0; i < DENSITY_PROFILE_RADII.size(); i++){
		std::vector<double> solutions = solveQuadratic(a, b, c - (DENSITY_PROFILE_RADII[i] * DENSITY_PROFILE_RADII[i]));
		if(solutions.size() == 0){ // No real solutions found, indicating shell isn't traversed
			break;
		}
		else{
			// Restrict solutions to zero and greater to enforce traversal from point of interest forward along ray
			// Otherwise density traversed along a direction out from a point that lies BENEATH the Earth's
			// outer surface will include density behind the point (ray in other direction) as well.
			intersectionsIngoing.push_back(std::max(solutions[0],0.));
			intersectionsOutgoing.push_back(std::max(solutions[1],0.));
//			std::cout << "Intersections at: " << solutions[0] << ", " << solutions[1] << std::endl;
		}
	}	
	return std::vector<std::vector<double>>{intersectionsIngoing, intersectionsOutgoing};
}

// t is the traversal distance given by intersection; a, b, c are aggregate coefficients of the ray traversal
// The compiler will likely optimize recurring terms, so, for clarity, code is a direct implementation of indefinite integral
// provided by output from Wolfram Alpha
double computeLinearDensityTraversal(double a, double b, double c, double t){
	return ((2.0 * a * t + b) * sqrt(t * (a * t + b) + c)) / (4.0 * a) - (((b*b - 4.0*a*c) * log(2.0*sqrt(a)*sqrt(t*(a*t+b)+c) + 2.0*a*t+b)) / (8.0 * sqrt(a*a*a)));
}

double computeQuadraticDensityTraversal(double a, double b, double c, double t){
	return (a / 3.0) * t * t * t + 0.5 * b * t * t + c * t;
}

double computeCubicDensityTraversal(double a, double b, double c, double t){	
	return (1.0 / (128.0 * sqrt(a*a*a*a*a)))
	* (2.0 * sqrt(a) *(2.0*a*t+b)*sqrt(t*(a*t+b)+c)
	* (8.0*a*b*t + 4.0*a*(2.0*a*t*t+5.0*c)-3.0*b*b)
	+ 3.0 * (b*b-4.0*a*c)*(b*b-4.0*a*c) * log(2.0*sqrt(a)*sqrt(t*(a*t+b)+c)+2.0*a*t+b));
}

// Provide initial ray position in cartesian meters relative to center of Earth and a unit direction
// in order to calculate integral of traversed density across PREM profiles applied to WGS84 ellipsoid
double getDensityTraversed(std::vector<double> position, std::vector<double> direction){
	double x = position[0];
	double y = position[1];
	double z = position[2];
	double dirX = direction[0];
	double dirY = direction[1];
	double dirZ = direction[2];
	
	// Coefficients of quadratic equation to solve against squared normalized PREM profile radii, (ax^2 + bx + c)
	double a = (dirX * dirX + dirY * dirY) / EQUATORIAL_EARTH_RADIUS_SQR + (dirZ * dirZ) / POLAR_EARTH_RADIUS_SQR;
	double b = 2.0 * ((x*dirX + y*dirY) / EQUATORIAL_EARTH_RADIUS_SQR + (z*dirZ) / POLAR_EARTH_RADIUS_SQR);
	double c =  (x*x + y*y) / EQUATORIAL_EARTH_RADIUS_SQR + (z*z) / POLAR_EARTH_RADIUS_SQR; // subtract respective profile radius squared to complete coefficient
	
//	std::cout << "Ray coefficients: " << a << " x^2, " << b << " x, " << c << std::endl;
	
	auto intersections = getIntersections(a,b,c);
	
	double densityTraversed = 0.0;
	
	if(intersections[0].size() > 0){
	for(int i = 0; i < intersections[0].size() - 1; i++){ // Last shell must be integrated between last outbound and inbound
//		std::cout << intersections[0].size() << std::endl;
		// Integrate both parts of each shell traversal
		switch (DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i].size()){
			case 4: // Cubic term and all less
			densityTraversed += DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][3] * (
			computeCubicDensityTraversal(a,b,c,intersections[0][i+1]) - computeCubicDensityTraversal(a,b,c,intersections[0][i]) +
			computeCubicDensityTraversal(a,b,c,intersections[1][i]) - computeCubicDensityTraversal(a,b,c,intersections[1][i+1])
			);
			case 3: // Quadratic term and all less
			densityTraversed += DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][2] * (
			computeQuadraticDensityTraversal(a,b,c,intersections[0][i+1]) - computeQuadraticDensityTraversal(a,b,c,intersections[0][i]) +
			computeQuadraticDensityTraversal(a,b,c,intersections[1][i]) - computeQuadraticDensityTraversal(a,b,c,intersections[1][i+1])
			);
			case 2: // Linear term and all less
			densityTraversed += DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][1] * (
			computeLinearDensityTraversal(a,b,c,intersections[0][i+1]) - computeLinearDensityTraversal(a,b,c,intersections[0][i]) +
			computeLinearDensityTraversal(a,b,c,intersections[1][i]) - computeLinearDensityTraversal(a,b,c,intersections[1][i+1])
			);
			case 1: // Constant term
			densityTraversed += DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][0] * (intersections[0][i+1] - intersections[0][i] +
			intersections[1][i] - intersections[1][i+1]);
			break;
			default: // Shouldn't occur
			break;
		}
//		std::cout << "Density traversed: " << densityTraversed << std::endl;
	}	
		int i = intersections[0].size() - 1;
		switch (DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i].size()){
			case 4: // Cubic term and all less
			densityTraversed += DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][3] * (
			computeCubicDensityTraversal(a,b,c,intersections[1][i]) - computeCubicDensityTraversal(a,b,c,intersections[0][i])
			);
			case 3: // Quadratic term and all less
			densityTraversed += DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][2] * (
			computeQuadraticDensityTraversal(a,b,c,intersections[1][i]) - computeQuadraticDensityTraversal(a,b,c,intersections[0][i])
			);
			case 2: // Linear term and all less
			densityTraversed += DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][1] * (
			computeLinearDensityTraversal(a,b,c,intersections[1][i]) - computeLinearDensityTraversal(a,b,c,intersections[0][i])
			);
			case 1: // Constant term
			densityTraversed += DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][0] * (intersections[1][i] - intersections[0][i]);
			break;
			default: // Shouldn't occur
			break;
		}
//		std::cout << "Density traversed: " << densityTraversed << std::endl;
	}
	
	return densityTraversed;
}

/*
 * 		DIAGNOSTIC FUNCTIONS
 */

void printConstants(){
	std::cout << "Enumeration of physical constants used in code:" << std::endl;
	std::cout << "Distance to horizon:			" << DISTANCE_TO_HORIZON << " meters." << std::endl;
	std::cout << "Mean Earth radius:			" << MEAN_EARTH_RADIUS << " meters." << std::endl;
	std::cout << "Polar Earth radius:			" << POLAR_EARTH_RADIUS << " meters."<< std::endl;
	std::cout << "Mean Earth density:			" << MEAN_EARTH_DENSITY << " kilograms per cubic meter." << std::endl;
	std::cout << "Vacuum speed of light:			" << SPEED_OF_LIGHT << " meters per second." << std::endl;
	std::cout << std::endl;
}

void printUsage(){
	std::cout << "Usage: anita_earthmodel2 [Nevt][ spectrum='ESS' or 'Enu'][Enu(eV)][maxdepth(m)][cross section factor]"<< std::endl;
	std::cout << "If single energy 'Enu' then argv[3] is the neutrino energy, else ignored."<< std::endl;
	std::cout << std::endl;
}

void testDensityTraversal(){
	std::cout << "Generating test data of density traversal for swept hemisphere through Earth from point at south pole and saving to 'diagnostic.dat'." << std::endl;
	std::ofstream outFile;
	outFile.precision(20);
	std::string outfilePath = "diagnostic.dat";
	outFile.open(outfilePath);
	auto position = std::vector<double>{0,0,-POLAR_EARTH_RADIUS - 1000}; // geometric south pole 
	auto direction = std::vector<double>{0,0,1.00}; // +Z direction
//	double densityTraversed = getDensityTraversed(position, direction);
//	std::cout << "Density traversed: " << densityTraversed << " kg m / m^3" << std::endl;
	unsigned int n = 100;
	double phi;
	double theta;
	for(unsigned int i = 0; i <= n; i++){
		phi = (2.0 * M_PI / n) * i;
		for (unsigned int j = 0; j <=n; j++){			
			theta = j *((M_PI) * (1.0 / n));
			direction[0] = sin(theta)*cos(phi);
			direction[1] = sin(theta)*sin(phi);
			direction[2] = cos(theta);
			outFile << phi << "	" << theta << "	" << getDensityTraversed(position, direction) << std::endl;			
		}
		outFile << std::endl;
	}
	outFile.close();
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
	testDensityTraversal();
	std::cout << "Done..." << std::endl;
	return 0;
}

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
#include "Vector.h"
#include "DataRaster.h"
#include "Diagnostics.h"
#include "QuadraticSolver.h"


/*
 * 		GLOBAL CONSTANTS
 */

// Data filepaths
const std::string filePathBed = "bedmap2_bed.flt";
const std::string filePathSurface = "bedmap2_surface.flt";
const std::string filePathIceThickness = "bedmap2_thickness.flt";
const std::string filePathGeoid = "gl04c_geiod_to_wgs84.flt"; // Note that BEDMAP 2 provides incorrect 'geiod' spelling.

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
	std::vector<double>{1020.0},								// Ocean
	std::vector<double>{2600.0},								// Crust
	std::vector<double>{2900.0},								// Crust
	std::vector<double>{2691.0, 692.4},							// LVZ/LID
	std::vector<double>{7108.9, -3804.5},						// Transition Zone
	std::vector<double>{11249.4, -8029.8},						// Transition Zone
	std::vector<double>{5319.7, -1483.6},						// Transition Zone
	std::vector<double>{7956.5, -6476.1, 5528.3, -3080.7}, 		// Lower Mantle
	std::vector<double>{12581.5, -1263.8, -3642.6, -5528.1},	// Outer Core
	std::vector<double>{13088.5, -8838.1}						// Inner Core
	};
	
const int DATA_ROWS = 6667;
const int DATA_COLUMNS = DATA_ROWS;
const int DATA_INTERVAL = 1000; // meters
const int PROJECTION_PLANE_LAT = 71; // degrees

const double Z_PLANE = sin(PROJECTION_PLANE_LAT*(M_PI/180.)) * POLAR_EARTH_RADIUS;

const double EPSILON = 0.001; // 1mm for numerical gradient
const double RAYSTEP = 500.0; // meters
const double UPPER_SURFACE_BOUND = 1.0 + (5000 / POLAR_EARTH_RADIUS); // Roughly higher than Mt Vinson in normalized radius
const double LOWER_SURFACE_BOUND = 1.0 - (3000 / POLAR_EARTH_RADIUS); // Roughly lower than Bentley Subglacial Trench in normalized radius

/*
 *		GLOBAL VARIABLES
 */

// Serial packing of rows, initialize with zeros to be overwritten by direct file read
std::vector<float> geoidData(DATA_ROWS * DATA_COLUMNS);
std::vector<float> bedData(DATA_ROWS * DATA_COLUMNS);
std::vector<float> surfaceData(DATA_ROWS * DATA_COLUMNS);
std::vector<float> iceThicknessData(DATA_ROWS * DATA_COLUMNS);

///
anita::DataRaster<float> geoidDataRaster;
anita::DataRaster<float> bedDataRaster;
anita::DataRaster<float> surfaceDataRaster;
anita::DataRaster<float> iceThicknessDataRaster;

/*
 * 		================= FUNCTIONS =========================
 */

/*
 *		Misc.
 */

void print(std::string message){
	std::cout << message << std::endl;
}

// For some machines endianness is reversed
float reverseFloat( const float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) &inFloat;
   char *returnFloat = ( char* ) &retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

/*
 *		DATA POLLING
 */

std::vector<double> getPositionOfTau(const double x0, const double y0, const double z0,
	const double xDir, const double yDir, const double zDir, const double tau){
	std::vector<double> position{
		x0 + tau * xDir,
		y0 + tau * yDir,
		z0 + tau * zDir
	};
	return position;
}

std::vector<double> getProjectionCoordinates(const double x, const double y, const double z){
	double alpha = (-POLAR_EARTH_RADIUS - Z_PLANE)/(-POLAR_EARTH_RADIUS + z);
	std::vector<double> projCoords{
		alpha * x,
		alpha * y
	};
	return projCoords;
}

std::vector<double> getDataCoordinates(const double xProj, const double yProj){
	std::vector<double> dataCoordinates{
		xProj+(DATA_COLUMNS*DATA_INTERVAL/2),
		yProj+(DATA_ROWS*DATA_INTERVAL/2)
	};
	return dataCoordinates;
}

std::vector<double> getNormalizedUVCoordinates(const double xData, const double yData){
	double xFrac, yFrac;
	modf(xData/DATA_INTERVAL, &xFrac);
	modf(yData/DATA_INTERVAL, &yFrac);
	std::vector<double> normalizedUVs{
		xFrac, yFrac
	};
	return normalizedUVs;
}

std::vector<double> getCellValues(const double xDataCoord, const double yDataCoord, const std::vector<double> &dataVector){	
	int xFloor = (uint16_t)std::min(std::max(floor(xDataCoord/DATA_INTERVAL), 0.), DATA_COLUMNS - 1.);
	int yFloor = (uint16_t)std::min(std::max(floor(yDataCoord/DATA_INTERVAL), 0.), DATA_ROWS - 1.);
	std::vector<double> cellValues{
		dataVector[yFloor*(DATA_COLUMNS - 1) + xFloor],
		dataVector[yFloor*(DATA_COLUMNS - 1) + std::min(xFloor + 1, DATA_COLUMNS - 1)],
		dataVector[std::min(yFloor + 1, DATA_ROWS - 1)*(DATA_COLUMNS - 1) + xFloor],
		dataVector[std::min(yFloor + 1, DATA_ROWS - 1)*(DATA_COLUMNS - 1) + std::min(xFloor + 1, DATA_COLUMNS - 1)],
	};
	return cellValues;
}

double getDataValue(const double x, const double y, const double z, std::vector<double> &dataVector){	
	std::vector<double> projCoords = getProjectionCoordinates(x, y, z);
	std::vector<double> dataCoordinates = getDataCoordinates(projCoords[0], projCoords[1]);
	std::vector<double> cellValues = getCellValues(dataCoordinates[0], dataCoordinates[1], dataVector);
	std::vector<double> normalizedUVs = getNormalizedUVCoordinates(dataCoordinates[0], dataCoordinates[1]);
	// Return bilinear interpolation of cell value based on UVs
	return cellValues[0]
	+ normalizedUVs[0]*(cellValues[1]-cellValues[0])
	+ normalizedUVs[1]*(cellValues[2]-cellValues[0])
	+ normalizedUVs[0]*normalizedUVs[1]*(cellValues[3]-cellValues[2]-cellValues[1]-cellValues[0]);
}

void loadData(){
	print("Loading Antarctic bedrock and ice surface elevations into memory...");
	std::ifstream dataFile(filePathGeoid, std::ios::binary);
	dataFile.read(reinterpret_cast<char*>(geoidData.data()), geoidData.size()*sizeof(float));	
	dataFile.close();
	dataFile.open(filePathBed);
	dataFile.read(reinterpret_cast<char*>(bedData.data()), bedData.size()*sizeof(float));
	dataFile.close();
	dataFile.open(filePathSurface);
	dataFile.read(reinterpret_cast<char*>(surfaceData.data()), surfaceData.size()*sizeof(float));
	dataFile.close();
	dataFile.open(filePathIceThickness);
	dataFile.read(reinterpret_cast<char*>(iceThicknessData.data()), iceThicknessData.size()*sizeof(float));
	dataFile.close();

	// Add geoid data to bed and ice to convert to WGS84 reference
	for(unsigned long long i = 0; i < geoidData.size(); i++){
		#ifdef _MSBF // Most significant bit first
		geoidData[i] = reverseFloat(geoidData[i]);
		bedData[i] = reverseFloat(bedData[i]);
		surfaceData[i] = reverseFloat(surfaceData[i]);
		iceThicknessData[i] = reverseFloat(iceThicknessData[i]);
		#endif
		if(geoidData[i] != -9999.){
			bedData[i] += geoidData[i];
			surfaceData[i] += geoidData[i];
		}
		std::cout << geoidData[i] << std::endl;		
	}

	print("Loading complete.");
}



/*
 *		PREM DENSITY TRAVERSAL
 */

void normalizeVector(std::vector<double> &v){
	double magsqr = 0;
	for(unsigned int i = 0; i < v.size(); i++){
		magsqr += v[i] * v[i];
	}
	if(magsqr == 0){
		return;
	}
	magsqr = sqrt(magsqr);
	for(unsigned int i = 0; i < v.size(); i++){
		v[i] /= magsqr;
	}
}

double getRadiusSqr(const double x, const double y, const double z){
	return x*x + y*y + z*z;
}

double getRadius(const double x, const double y, const double z){
	return sqrt(getRadiusSqr(x,y,z));
}

double getEllipsoidalRadius(const double x, const double y, const double z){	
	return sqrt((EQUATORIAL_EARTH_RADIUS_SQR*(x*x + y*y))/(x*x + y*y + z*z) + (POLAR_EARTH_RADIUS_SQR*z*z)/(x*x + y*y + z*z));
}

double getDistanceToSurface(const double x, const double y, const double z, std::vector<double> &dataVector){
	return getDataValue(x, y, z, dataVector) + getEllipsoidalRadius(x, y, z) - getRadius(x, y, z);
}

double getFirnDensity(double depth){
	//TODO: need reference materials to obtain good function. Does exponential packing suffice?
	return 0;
}

// Supply coefficients  in ax^2+bx+c form to obtain real solutions
std::vector<double> solveQuadratic(double a, double b, double c){
	std::vector<double> solutions;
	double discriminant = b * b - 4 * a * c;
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
	return solutions;
}

std::vector<std::vector<double>> getIntersections(double a, double b, double c){
	std::vector<double> intersectionsIngoing;
	std::vector<double> intersectionsOutgoing;
	// Test for PREM density shell intersection from surface inwards to core
	// Loop will break out as soon as no intersection is found, meaning current and deeper shells aren't traversed.
	for(unsigned int i = 0; i < DENSITY_PROFILE_RADII.size(); i++){
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
		}
	}	
	return std::vector<std::vector<double>>{intersectionsIngoing, intersectionsOutgoing};
}

// t is the traversal distance given by intersection; a, b, c are aggregate coefficients of the ray traversal
// The compiler will likely optimize recurring terms, so, for clarity, code is a direct implementation of integral
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

std::vector<double> getQuadraticCoefficientsOfNormalizedEllipsoidalRay(const double x, const double y, const double z,
const double dirX, const double dirY, const double dirZ){
	std::vector<double> coefficients;
	coefficients.push_back((dirX * dirX + dirY * dirY) / EQUATORIAL_EARTH_RADIUS_SQR + (dirZ * dirZ) / POLAR_EARTH_RADIUS_SQR);
	coefficients.push_back(2.0 * ((x*dirX + y*dirY) / EQUATORIAL_EARTH_RADIUS_SQR + (z*dirZ) / POLAR_EARTH_RADIUS_SQR));
	coefficients.push_back((x*x + y*y) / EQUATORIAL_EARTH_RADIUS_SQR + (z*z) / POLAR_EARTH_RADIUS_SQR);
	return coefficients;
}

// Provide initial ray position in cartesian meters relative to center of Earth and a unit direction
// in order to calculate integral of traversed density across PREM profiles applied to WGS84 ellipsoid
double getDensityTraversed(std::vector<double> &position, std::vector<double> &direction){
	// Relabel vector parts for ease of readability
	// Compiler will likely optimize out dummy variables
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
		
	auto intersections = getIntersections(a,b,c);
	
	double densityTraversed = 0.0;
	
	if(intersections[0].size() > 0){
	for(unsigned int i = 0; i < intersections[0].size() - 1; i++){ // Last shell must be integrated between last outbound and inbound
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
	}
	
	return densityTraversed;
}

/*
*		QUADRATIC RAYSTEPPING
*/

double getGradientAlongRayAtPoint(const double x, const double y, const double z,
	const double xDir, const double yDir, const double zDir, std::vector<double> &dataVector){
		return (getDistanceToSurface(x+EPSILON*xDir, y+EPSILON*yDir, z+EPSILON*zDir, dataVector)
		- getDistanceToSurface(x-EPSILON*xDir, y-EPSILON*yDir, z-EPSILON*zDir, dataVector))
		/ (2 * EPSILON);
}

std::vector<double> getCoefficientsOfQuadraticRaystep(const double x, const double y, const double z,
	const double xDir, const double yDir, const double zDir, const double previousGradientAlongRay,
	const double previousDistanceToSurface, const double tau0, std::vector<double> &dataVector){
		std::vector<double> coefficients(3);
		double tau1 = tau0 + RAYSTEP;
		double nextX = x + RAYSTEP*xDir;
		double nextY = y + RAYSTEP*yDir;
		double nextZ = z + RAYSTEP*zDir;
		double nextGradientAlongRay = getGradientAlongRayAtPoint(nextX, nextY, nextZ, xDir, yDir, zDir, dataVector);
		coefficients[0] = (nextGradientAlongRay - previousGradientAlongRay) / (2 * RAYSTEP);
		double nextDistanceToSurface = getDistanceToSurface(nextX, nextY, nextZ, dataVector);
		coefficients[1] = (nextDistanceToSurface - previousDistanceToSurface - coefficients[0] * (tau1*tau1 - tau0*tau0)) / RAYSTEP;
		coefficients[2] = previousDistanceToSurface - coefficients[1] * tau0 - coefficients[0] * tau0 * tau0;
		return coefficients;
}

// // Returns a vector of x,y,z,tau vectors
// std::vector<std::vector<double>> intersectSurface(const double x, const double y, const double z,
// 	const double xDir, const double yDir, const double zDir, const double tau, std::vector<double> &dataVector){	
// 	// Compute taus of intersection with bounding ellipsoids
// 	// Terrain radii are bounded by lower and upper radius to determine valid surface polling areas
// 	std::vector<double> quadCoeff = getQuadraticCoefficientsOfNormalizedEllipsoidalRay(x,y,z,xDir,yDir,zDir);
// 	std::vector<std::vector<double>> boundIntersections;
// 	boundIntersections.push_back(solveQuadratic(quadCoeff[0],quadCoeff[1],quadCoeff[2]-(LOWER_SURFACE_BOUND*LOWER_SURFACE_BOUND)));
// 	boundIntersections.push_back(solveQuadratic(quadCoeff[0],quadCoeff[1],quadCoeff[2]-(UPPER_SURFACE_BOUND*UPPER_SURFACE_BOUND)));
	
// }

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

// void testData(){
// 	int n = 6667;
// 	int d = 66;
// 	std::ofstream outFile;
// 	outFile.open("map.dat");
// 	for(int y = 0; y < n; y += d){
// 		for(int x = 0; x < n; x += d){
// 			outFile << x * 1000 << "	" << y * 1000 << "	" << ((iceThicknessData[x + n*y] == -9999.0f) ? 0 : iceThicknessData[x + n*y]) << std::endl;
// 		}
// 		outFile << std::endl;
// 	}
// 	outFile.close();
// }

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
	std::thread t1(importData, std::ref(geoidDataRaster), std::ref(filePathGeoid));
	std::thread t2(importData, std::ref(bedDataRaster), std::ref(filePathBed));
	std::thread t3(importData, std::ref(surfaceDataRaster), std::ref(filePathSurface));
	std::thread t4(importData, std::ref(iceThicknessDataRaster), std::ref(filePathIceThickness));
	t1.join();
	t2.join();
	t3.join();
	t4.join();
	std::cout << geoidDataRaster.values[0] << std::endl;
	print("Done...");
	print("Press any key to close.");
	std::cin.get(); // Wait for user input to terminate
	return 0;
}

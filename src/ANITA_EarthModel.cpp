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
		return (getDistanceToSurface(x+anita::EPSILON*xDir, y+anita::EPSILON*yDir, z+anita::EPSILON*zDir, dataVector)
		- getDistanceToSurface(x-anita::EPSILON*xDir, y-anita::EPSILON*yDir, z-anita::EPSILON*zDir, dataVector))
		/ (2 * anita::EPSILON);
}

std::vector<double> getCoefficientsOfQuadraticRaystep(const double x, const double y, const double z,
	const double xDir, const double yDir, const double zDir, const double previousGradientAlongRay,
	const double previousDistanceToSurface, const double tau0, std::vector<double> &dataVector){
		std::vector<double> coefficients(3);
		double tau1 = tau0 + anita::RAYSTEP;
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
	std::cout << "Distance to horizon:			" << anita::DISTANCE_TO_HORIZON << " meters." << std::endl;
	std::cout << "Mean Earth radius:			" << anita::MEAN_EARTH_RADIUS << " meters." << std::endl;
	std::cout << "Polar Earth radius:			" << anita::POLAR_EARTH_RADIUS << " meters."<< std::endl;
	std::cout << "Mean Earth density:			" << anita::MEAN_EARTH_DENSITY << " kilograms per cubic meter." << std::endl;
	std::cout << "Vacuum speed of light:			" << anita::SPEED_OF_LIGHT << " meters per second." << std::endl;
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
	auto position = std::vector<double>{0,0,-anita::POLAR_EARTH_RADIUS - 1000}; // geometric south pole 
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
	std::thread t1(importData, std::ref(anita::geoidDataRaster), std::ref(anita::filePathGeoid));
	std::thread t2(importData, std::ref(anita::bedDataRaster), std::ref(anita::filePathBed));
	std::thread t3(importData, std::ref(anita::surfaceDataRaster), std::ref(anita::filePathSurface));
	std::thread t4(importData, std::ref(anita::iceThicknessDataRaster), std::ref(anita::filePathIceThickness));
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

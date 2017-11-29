#define _USING_MATH_DEFINES
#include "Raycasting.h"
#include "BEDMAP.h"
#include "QuadraticSolver.h"
#include <cmath>
namespace anita{

	Vector3<double> getPositionOfTau(const Vector3<double>& pos, const Vector3<double>& dir, const double tau){
		return pos + dir.norm() * tau;
	}

	Vector2<double> getProjectionCoordinates(const Vector3<double>& pos){
		const auto alpha = (-POLAR_EARTH_RADIUS - Z_PLANE)/(-POLAR_EARTH_RADIUS + pos.z);
		return Vector2<double>(pos.x * alpha, pos.y * alpha);
	}

	Vector2<double> getDataCoordinates(const Vector2<double>& projCoords){
		return Vector2<double>(
			projCoords.x + (DATA_COLUMNS*DATA_INTERVAL/2),
			projCoords.y + (DATA_ROWS*DATA_INTERVAL/2)
		);
	}

	Vector2<double> getNormalizedUVCoordinates(const Vector2<double>& dataCoords){
		double xFrac, yFrac, xInt, yInt;
		xFrac = modf(dataCoords.x / DATA_INTERVAL, &xInt);
		yFrac = modf(dataCoords.y / DATA_INTERVAL, &yInt);
		return Vector2<double>(xFrac, yFrac);
	}

	double getDataValue(const Vector3<double>& pos, const DataRaster<float>& dataRaster){	
		const auto projCoords = getProjectionCoordinates(pos);
		const auto dataCoords = getDataCoordinates(projCoords);
		const auto cellValues = dataRaster.getCellValues(dataCoords);
		const auto normalizedUVs = getNormalizedUVCoordinates(dataCoords);
		// Return bilinear interpolation of cell value based on UVs
		return cellValues[0]*((1.0-normalizedUVs.x)*(1.0-normalizedUVs.y))
			+ cellValues[1]*((normalizedUVs.x)*(1.0-normalizedUVs.y))
			+ cellValues[2]*((1.0-normalizedUVs.x)*(normalizedUVs.y))
			+ cellValues[3]*((normalizedUVs.x)*(normalizedUVs.y));
	}

	double getEllipsoidalRadius(const Vector3<double>& pos){	
		return sqrt(
			((EQUATORIAL_EARTH_RADIUS_SQR*(pos.x*pos.x + pos.y*pos.y))
			+ (POLAR_EARTH_RADIUS_SQR*pos.z*pos.z))
			/ pos.magSqr()
		);
	}

	double getDistanceToSurface(const Vector3<double>& pos, const DataRaster<float>& dataRaster){
		return getEllipsoidalRadius(pos) + getDataValue(pos, dataRaster) - pos.mag();
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
				intersectionsIngoing.push_back(std::max(solutions[0], 0.0));
				intersectionsOutgoing.push_back(std::max(solutions[1], 0.0));
			}
		}	
		return std::vector<std::vector<double>>{intersectionsIngoing, intersectionsOutgoing};
	}

	// t is the traversal distance given by intersection; a, b, c are aggregate coefficients of the ray traversal
	// The compiler will likely optimize recurring terms, so, for clarity, code is a direct implementation of integral
	// provided by output from Wolfram Alpha
	double computeLinearDensityTraversal(const double a, const double b, const double c, const double t){
		return ((2.0 * a * t + b) * sqrt(t * (a * t + b) + c)) / (4.0 * a) - (((b*b - 4.0*a*c) * log(2.0*sqrt(a)*sqrt(t*(a*t+b)+c) + 2.0*a*t+b)) / (8.0 * sqrt(a*a*a)));
	}

	double computeQuadraticDensityTraversal(const double a, const double b, const double c, const double t){
		return (a / 3.0) * t * t * t + 0.5 * b * t * t + c * t;
	}

	double computeCubicDensityTraversal(const double a, const double b, const double c, const double t){	
		return (1.0 / (128.0 * sqrt(a*a*a*a*a)))
		* (2.0 * sqrt(a) *(2.0*a*t+b)*sqrt(t*(a*t+b)+c)
		* (8.0*a*b*t + 4.0*a*(2.0*a*t*t+5.0*c)-3.0*b*b)
		+ 3.0 * (b*b-4.0*a*c)*(b*b-4.0*a*c) * log(2.0*sqrt(a)*sqrt(t*(a*t+b)+c)+2.0*a*t+b));
	}

	std::vector<double> getQuadraticCoefficientsOfNormalizedEllipsoidalRay(const Vector3<double>& pos, const Vector3<double>& direction){
		auto dir = direction.norm();
		std::vector<double> coefficients;
		coefficients.push_back((dir.x * dir.x + dir.y * dir.y) / EQUATORIAL_EARTH_RADIUS_SQR + (dir.z * dir.z) / POLAR_EARTH_RADIUS_SQR);
		coefficients.push_back(2.0 * ((pos.x*dir.x + pos.y*dir.y) / EQUATORIAL_EARTH_RADIUS_SQR + (pos.z*dir.z) / POLAR_EARTH_RADIUS_SQR));
		coefficients.push_back((pos.x*pos.x + pos.y*pos.y) / EQUATORIAL_EARTH_RADIUS_SQR + (pos.z*pos.z) / POLAR_EARTH_RADIUS_SQR);
		return coefficients;
	}

	// Provide initial ray position in cartesian meters relative to center of Earth and a unit direction
	// in order to calculate integral of traversed density across PREM profiles applied to WGS84 ellipsoid
	double getDensityTraversed(const Vector3<double>& pos, const Vector3<double>& direction){
		const auto dir = direction.norm();
		
		// Coefficients of quadratic equation to solve against squared normalized PREM profile radii, (ax^2 + bx + c)
		const double a = (dir.x * dir.x + dir.y * dir.y) / EQUATORIAL_EARTH_RADIUS_SQR + (dir.z * dir.z) / POLAR_EARTH_RADIUS_SQR;
		const double b = 2.0 * ((pos.x*dir.x + pos.y*dir.y) / EQUATORIAL_EARTH_RADIUS_SQR + (pos.z*dir.z) / POLAR_EARTH_RADIUS_SQR);
		const double c =  (pos.x*pos.x + pos.y*pos.y) / EQUATORIAL_EARTH_RADIUS_SQR + (pos.z*pos.z) / POLAR_EARTH_RADIUS_SQR; // subtract respective profile radius squared to complete coefficient
			
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

	double getGradientAlongRayAtPoint(const Vector3<double>& pos, const Vector3<double>& direction, const DataRaster<float>& dataRaster){
		auto dir = direction.norm();
		return (getDistanceToSurface(pos + dir * EPSILON, dataRaster)
			- getDistanceToSurface(pos - dir * EPSILON, dataRaster))
			/ (2 * EPSILON);
	}

	std::vector<double> getCoefficientsOfQuadraticRaystep(const Vector3<double>& pos, const Vector3<double>& direction, const double previousGradientAlongRay,
		const double previousDistanceToSurface, const double tau0, const DataRaster<float>& dataRaster){
			const auto dir = direction.norm();
			std::vector<double> coefficients(3);
			const double tau1 = tau0 + RAYSTEP;
			const auto nextPos = pos + dir * RAYSTEP;
			const double nextGradientAlongRay = getGradientAlongRayAtPoint(nextPos, dir, dataRaster);
			coefficients[0] = (nextGradientAlongRay - previousGradientAlongRay) / (2 * RAYSTEP);
			double nextDistanceToSurface = getDistanceToSurface(nextPos, dataRaster);
			coefficients[1] = (nextDistanceToSurface - previousDistanceToSurface - coefficients[0] * (tau1*tau1 - tau0*tau0)) / RAYSTEP;
			coefficients[2] = previousDistanceToSurface - coefficients[1] * tau0 - coefficients[0] * tau0 * tau0;
			return coefficients;
	}
		
	std::vector<double> getTauOfSurfaceIntersections(const Vector3<double>& pos, const Vector3<double> direction, const DataRaster<float>& dataRaster){
		const auto dir = direction.norm();
		const auto quadraticCoefficients = getQuadraticCoefficientsOfNormalizedEllipsoidalRay(pos, dir);
		const auto lower = solveQuadratic(quadraticCoefficients[0], quadraticCoefficients[1], quadraticCoefficients[2] - LOWER_SURFACE_BOUND * LOWER_SURFACE_BOUND);
		const auto upper = solveQuadratic(quadraticCoefficients[0], quadraticCoefficients[1], quadraticCoefficients[2] - UPPER_SURFACE_BOUND * UPPER_SURFACE_BOUND);
		
		// Surface intersections
		std::vector<double> intersections;

		// Does ray never potentially intersect Earth? (no intersections with upper bounds)
		// Or, does ray point away from Earth? (Both tau of intersections are negative, behind starting point)
		// (Short-circuit evaluation prevents accessing indices if size() == 0)
		if((upper.size() == 0) || ((upper[0] < 0) && (upper[1] < 0))){
			return intersections; // These are surface intersections, NOT upper bounds intersections
		}

		double tau;
		// double distanceToSurface;
		
		Vector3<double> position = pos;

		// Does ray intersect lower bounds, or are intersections behind starting point?
		// TODO: lower bounds considerations to optimize evaluation, ignoring for now with 'true'
		if(true || (lower.size() == 0) || ((lower[0] < 0) && (lower[1] < 0))){
			// ignore lower bounds intersections and iterate to intersect surface
			double gradientAlongRay;
			tau = std::max(0.0, upper[0]);
			while (tau < upper[1]){
				gradientAlongRay = getGradientAlongRayAtPoint(position, dir, dataRaster);
			}
		}
		else{
			// TODO
		}

		
		return intersections;
	}

	// double getProbabilityOfInteraction(const Neutrino& neutrino){
		
	// }
}
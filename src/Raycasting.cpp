#include "Raycasting.h"
#include "BEDMAP.h"
#include "QuadraticSolver.h"

anita::Vector3<double> anita::getPositionOfTau(const anita::Vector3<double>& pos, const anita::Vector3<double>& dir, const double tau){
    return pos + dir.norm() * tau;
}

anita::Vector2<double> anita::getProjectionCoordinates(const anita::Vector3<double>& pos){
    const auto alpha = (-anita::POLAR_EARTH_RADIUS - anita::Z_PLANE)/(-anita::POLAR_EARTH_RADIUS + pos.z);
    return anita::Vector2<double>(pos.x * alpha, pos.y * alpha);
}

anita::Vector2<double> anita::getDataCoordinates(const anita::Vector2<double>& projCoords){
    return anita::Vector2<double>(
        projCoords.x + (anita::DATA_COLUMNS*anita::DATA_INTERVAL/2),
        projCoords.y + (anita::DATA_ROWS*anita::DATA_INTERVAL/2)
    );
}

anita::Vector2<double> anita::getNormalizedUVCoordinates(const anita::Vector2<double>& dataCoords){
    double xFrac, yFrac;
    modf(dataCoords.x / anita::DATA_INTERVAL, &xFrac);
    modf(dataCoords.y / anita::DATA_INTERVAL, &yFrac);
    return anita::Vector2<double>(xFrac, yFrac);
}

double anita::getDataValue(const anita::Vector3<double>& pos, const anita::DataRaster<float>& dataRaster){	
    const auto projCoords = anita::getProjectionCoordinates(pos);
    const auto dataCoords = anita::getDataCoordinates(projCoords);
    const auto cellValues = dataRaster.getCellValues(dataCoords);
    const auto normalizedUVs = anita::getNormalizedUVCoordinates(dataCoords);
    // Return bilinear interpolation of cell value based on UVs
    return cellValues[0]
    + normalizedUVs.x*(cellValues[1]-cellValues[0])
    + normalizedUVs.y*(cellValues[2]-cellValues[0])
    + normalizedUVs.x*normalizedUVs.y*(cellValues[3]-cellValues[2]-cellValues[1]-cellValues[0]);
}

double anita::getEllipsoidalRadius(const anita::Vector3<double>& pos){	
	return sqrt(
		((anita::EQUATORIAL_EARTH_RADIUS_SQR*(pos.x*pos.x + pos.y*pos.y))
		+ (anita::POLAR_EARTH_RADIUS_SQR*pos.z*pos.z))
		/ pos.magSqr()
	);
}

double anita::getDistanceToSurface(const anita::Vector3<double>& pos, anita::DataRaster<float> dataRaster){
	return anita::getEllipsoidalRadius(pos) + anita::getDataValue(pos, dataRaster) - pos.mag();
}

std::vector<std::vector<double>> anita::getIntersections(double a, double b, double c){
	std::vector<double> intersectionsIngoing;
	std::vector<double> intersectionsOutgoing;
	// Test for PREM density shell intersection from surface inwards to core
	// Loop will break out as soon as no intersection is found, meaning current and deeper shells aren't traversed.
	for(unsigned int i = 0; i < anita::DENSITY_PROFILE_RADII.size(); i++){
		std::vector<double> solutions = anita::solveQuadratic(a, b, c - (anita::DENSITY_PROFILE_RADII[i] * anita::DENSITY_PROFILE_RADII[i]));
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
double anita::computeLinearDensityTraversal(const double a, const double b, const double c, const double t){
	return ((2.0 * a * t + b) * sqrt(t * (a * t + b) + c)) / (4.0 * a) - (((b*b - 4.0*a*c) * log(2.0*sqrt(a)*sqrt(t*(a*t+b)+c) + 2.0*a*t+b)) / (8.0 * sqrt(a*a*a)));
}

double anita::computeQuadraticDensityTraversal(const double a, const double b, const double c, const double t){
	return (a / 3.0) * t * t * t + 0.5 * b * t * t + c * t;
}

double anita::computeCubicDensityTraversal(const double a, const double b, const double c, const double t){	
	return (1.0 / (128.0 * sqrt(a*a*a*a*a)))
	* (2.0 * sqrt(a) *(2.0*a*t+b)*sqrt(t*(a*t+b)+c)
	* (8.0*a*b*t + 4.0*a*(2.0*a*t*t+5.0*c)-3.0*b*b)
	+ 3.0 * (b*b-4.0*a*c)*(b*b-4.0*a*c) * log(2.0*sqrt(a)*sqrt(t*(a*t+b)+c)+2.0*a*t+b));
}

std::vector<double> anita::getQuadraticCoefficientsOfNormalizedEllipsoidalRay(const anita::Vector3<double>& pos, const anita::Vector3<double>& direction){
	auto dir = direction.norm();
	std::vector<double> coefficients;
	coefficients.push_back((dir.x * dir.x + dir.y * dir.y) / anita::EQUATORIAL_EARTH_RADIUS_SQR + (dir.z * dir.z) / anita::POLAR_EARTH_RADIUS_SQR);
	coefficients.push_back(2.0 * ((pos.x*dir.x + pos.y*dir.y) / anita::EQUATORIAL_EARTH_RADIUS_SQR + (pos.z*dir.z) / anita::POLAR_EARTH_RADIUS_SQR));
	coefficients.push_back((pos.x*pos.x + pos.y*pos.y) / anita::EQUATORIAL_EARTH_RADIUS_SQR + (pos.z*pos.z) / anita::POLAR_EARTH_RADIUS_SQR);
	return coefficients;
}

// Provide initial ray position in cartesian meters relative to center of Earth and a unit direction
// in order to calculate integral of traversed density across PREM profiles applied to WGS84 ellipsoid
double anita::getDensityTraversed(const anita::Vector3<double>& pos, const anita::Vector3<double>& direction){
	auto dir = direction.norm();
	
	// Coefficients of quadratic equation to solve against squared normalized PREM profile radii, (ax^2 + bx + c)
	double a = (dir.x * dir.x + dir.y * dir.y) / anita::EQUATORIAL_EARTH_RADIUS_SQR + (dir.z * dir.z) / anita::POLAR_EARTH_RADIUS_SQR;
	double b = 2.0 * ((pos.x*dir.x + pos.y*dir.y) / anita::EQUATORIAL_EARTH_RADIUS_SQR + (pos.z*dir.z) / anita::POLAR_EARTH_RADIUS_SQR);
	double c =  (pos.x*pos.x + pos.y*pos.y) / anita::EQUATORIAL_EARTH_RADIUS_SQR + (pos.z*pos.z) / anita::POLAR_EARTH_RADIUS_SQR; // subtract respective profile radius squared to complete coefficient
		
	auto intersections = anita::getIntersections(a,b,c);
	
	double densityTraversed = 0.0;
	
	if(intersections[0].size() > 0){
	for(unsigned int i = 0; i < intersections[0].size() - 1; i++){ // Last shell must be integrated between last outbound and inbound
		// Integrate both parts of each shell traversal
		switch (anita::DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i].size()){
			case 4: // Cubic term and all less
			densityTraversed += anita::DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][3] * (
				anita::computeCubicDensityTraversal(a,b,c,intersections[0][i+1]) - anita::computeCubicDensityTraversal(a,b,c,intersections[0][i]) +
				anita::computeCubicDensityTraversal(a,b,c,intersections[1][i]) - anita::computeCubicDensityTraversal(a,b,c,intersections[1][i+1])
			);
			case 3: // Quadratic term and all less
			densityTraversed += anita::DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][2] * (
				anita::computeQuadraticDensityTraversal(a,b,c,intersections[0][i+1]) - anita::computeQuadraticDensityTraversal(a,b,c,intersections[0][i]) +
				anita::computeQuadraticDensityTraversal(a,b,c,intersections[1][i]) - anita::computeQuadraticDensityTraversal(a,b,c,intersections[1][i+1])
			);
			case 2: // Linear term and all less
			densityTraversed += anita::DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][1] * (
				anita::computeLinearDensityTraversal(a,b,c,intersections[0][i+1]) - anita::computeLinearDensityTraversal(a,b,c,intersections[0][i]) +
				anita::computeLinearDensityTraversal(a,b,c,intersections[1][i]) - anita::computeLinearDensityTraversal(a,b,c,intersections[1][i+1])
			);
			case 1: // Constant term
			densityTraversed += anita::DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][0] * (intersections[0][i+1] - intersections[0][i] +
			intersections[1][i] - intersections[1][i+1]);
			break;
			default: // Shouldn't occur
			break;
		}
	}	
		int i = intersections[0].size() - 1;
		switch (anita::DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i].size()){
			case 4: // Cubic term and all less
			densityTraversed += anita::DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][3] * (
				anita::computeCubicDensityTraversal(a,b,c,intersections[1][i]) - anita::computeCubicDensityTraversal(a,b,c,intersections[0][i])
			);
			case 3: // Quadratic term and all less
			densityTraversed += anita::DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][2] * (
				anita::computeQuadraticDensityTraversal(a,b,c,intersections[1][i]) - anita::computeQuadraticDensityTraversal(a,b,c,intersections[0][i])
			);
			case 2: // Linear term and all less
			densityTraversed += anita::DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][1] * (
				anita::computeLinearDensityTraversal(a,b,c,intersections[1][i]) - anita::computeLinearDensityTraversal(a,b,c,intersections[0][i])
			);
			case 1: // Constant term
			densityTraversed += anita::DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS[i][0] * (intersections[1][i] - intersections[0][i]);
			break;
			default: // Shouldn't occur
			break;
		}
	}
	
	return densityTraversed;
}

double anita::getGradientAlongRayAtPoint(const anita::Vector3<double>& pos, const anita::Vector3<double>& direction, const anita::DataRaster<float>& dataRaster){
	auto dir = direction.norm();
	return (anita::getDistanceToSurface(pos + dir * anita::EPSILON, dataRaster)
		- anita::getDistanceToSurface(pos - dir * anita::EPSILON, dataRaster))
		/ (2 * anita::EPSILON);
}

std::vector<double> anita::getCoefficientsOfQuadraticRaystep(const anita::Vector3<double>& pos, const anita::Vector3<double>& direction, const double previousGradientAlongRay,
	const double previousDistanceToSurface, const double tau0, const anita::DataRaster<float>& dataRaster){
		const auto dir = direction.norm();
		std::vector<double> coefficients(3);
		const double tau1 = tau0 + anita::RAYSTEP;
		const auto nextPos = pos + dir * anita::RAYSTEP;
		const double nextGradientAlongRay = anita::getGradientAlongRayAtPoint(nextPos, dir, dataRaster);
		coefficients[0] = (nextGradientAlongRay - previousGradientAlongRay) / (2 * anita::RAYSTEP);
		double nextDistanceToSurface = anita::getDistanceToSurface(nextPos, dataRaster);
		coefficients[1] = (nextDistanceToSurface - previousDistanceToSurface - coefficients[0] * (tau1*tau1 - tau0*tau0)) / anita::RAYSTEP;
		coefficients[2] = previousDistanceToSurface - coefficients[1] * tau0 - coefficients[0] * tau0 * tau0;
		return coefficients;
}
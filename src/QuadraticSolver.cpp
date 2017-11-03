#include "QuadraticSolver.h"

std::vector<double> anita::solveQuadratic(const double squareTerm, const double linearTerm, const double constantTerm){
	std::vector<double> solutions;
	// Compiler will optimize out dummy vars
    double a = squareTerm;
    double b = linearTerm;
    double c = constantTerm;
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
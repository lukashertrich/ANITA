#pragma once
#include <vector>
#include "Vector2.h"
#include "Vector3.h"
#include "DataRaster.h"
namespace anita{
    anita::Vector3<double> getPositionOfTau(const anita::Vector3<double>& pos, const anita::Vector3<double>& dir, const double tau);    
    anita::Vector2<double> getProjectionCoordinates(const anita::Vector3<double>& pos);    
    anita::Vector2<double> getDataCoordinates(const anita::Vector2<double>& projCoords);    
    anita::Vector2<double> getNormalizedUVCoordinates(const anita::Vector2<double>& dataCoords);
    double getDataValue(const anita::Vector3<double>& pos, const anita::DataRaster<float>& dataRaster);
    double getEllipsoidalRadius(const anita::Vector3<double>& pos);    
    double getDistanceToSurface(const anita::Vector3<double>& pos, anita::DataRaster<float> dataRaster);
    std::vector<std::vector<double>> getIntersections(double a, double b, double c);
    double computeLinearDensityTraversal(const double a, const double b, const double c, const double t);
    double computeQuadraticDensityTraversal(const double a, const double b, const double c, const double t);
    double computeCubicDensityTraversal(const double a, const double b, const double c, const double t);
    std::vector<double> getQuadraticCoefficientsOfNormalizedEllipsoidalRay(const anita::Vector3<double>& pos, const anita::Vector3<double>& direction);
    double getDensityTraversed(const anita::Vector3<double>& pos, const anita::Vector3<double>& direction);
    double getGradientAlongRayAtPoint(const anita::Vector3<double>& pos, const anita::Vector3<double>& direction, const anita::DataRaster<float>& dataRaster);
    std::vector<double> getCoefficientsOfQuadraticRaystep(const anita::Vector3<double>& pos, const anita::Vector3<double>& direction, const double previousGradientAlongRay,const double previousDistanceToSurface, const double tau0, const anita::DataRaster<float>& dataRaster);
}
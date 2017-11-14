#pragma once
#include <vector>
#include "Vector2.h"
#include "Vector3.h"
#include "DataRaster.h"
namespace anita{
    Vector3<double> getPositionOfTau(const Vector3<double>& pos, const Vector3<double>& dir, const double tau);    
    Vector2<double> getProjectionCoordinates(const Vector3<double>& pos);    
    Vector2<double> getDataCoordinates(const Vector2<double>& projCoords);    
    Vector2<double> getNormalizedUVCoordinates(const Vector2<double>& dataCoords);
    double getDataValue(const Vector3<double>& pos, const DataRaster<float>& dataRaster);
    double getEllipsoidalRadius(const Vector3<double>& pos);    
    double getDistanceToSurface(const Vector3<double>& pos, const DataRaster<float>& dataRaster);
    std::vector<std::vector<double>> getIntersections(double a, double b, double c);
    double computeLinearDensityTraversal(const double a, const double b, const double c, const double t);
    double computeQuadraticDensityTraversal(const double a, const double b, const double c, const double t);
    double computeCubicDensityTraversal(const double a, const double b, const double c, const double t);
    std::vector<double> getQuadraticCoefficientsOfNormalizedEllipsoidalRay(const Vector3<double>& pos, const Vector3<double>& direction);
    double getDensityTraversed(const Vector3<double>& pos, const Vector3<double>& direction);
    double getGradientAlongRayAtPoint(const Vector3<double>& pos, const Vector3<double>& direction, const DataRaster<float>& dataRaster);
    std::vector<double> getCoefficientsOfQuadraticRaystep(const Vector3<double>& pos, const Vector3<double>& direction, const double previousGradientAlongRay,const double previousDistanceToSurface, const double tau0, const DataRaster<float>& dataRaster);
    std::vector<double> getTauOfSurfaceIntersections(const Vector3<double>& pos, const Vector3<double> direction, const DataRaster<float>& dataRaster);
}
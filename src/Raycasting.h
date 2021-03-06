#pragma once
#include <vector>
#include "Vector2.h"
#include "Vector3d.h"
#include "DataRaster.h"
#include "Neutrino.h"
namespace anita{
    Vector3d getPositionOfTau(const Vector3d& pos, const Vector3d& dir, const double tau);    
    Vector2<double> getProjectionCoordinates(const Vector3d& pos);    
    Vector2<double> getDataCoordinates(const Vector2<double>& projCoords);    
    Vector2<double> getNormalizedUVCoordinates(const Vector2<double>& dataCoords);
    double getDataValue(const Vector3d& pos, const DataRaster<float>& dataRaster);
    double getEllipsoidalRadius(const Vector3d& pos);    
    double getDistanceToSurface(const Vector3d& pos, const DataRaster<float>& dataRaster);
    std::vector<std::vector<double>> getIntersections(double a, double b, double c);
    double computeLinearDensityTraversal(const double a, const double b, const double c, const double t);
    double computeQuadraticDensityTraversal(const double a, const double b, const double c, const double t);
    double computeCubicDensityTraversal(const double a, const double b, const double c, const double t);
    std::vector<double> getQuadraticCoefficientsOfNormalizedEllipsoidalRay(const Vector3d& pos, const Vector3d& direction);
    double getInteractionLength(const Vector3d& pos, const Vector3d& direction);
    double getGradientAlongRayAtPoint(const Vector3d& pos, const Vector3d& direction, const DataRaster<float>& dataRaster);
    std::vector<double> getCoefficientsOfQuadraticRaystep(const Vector3d& pos, const Vector3d& direction, const double previousGradientAlongRay,const double previousDistanceToSurface, const double tau0, const DataRaster<float>& dataRaster);
    std::vector<double> getTauOfSurfaceIntersections(const Vector3d& pos, const Vector3d direction, const DataRaster<float>& dataRaster);
    double getProbabilityOfInteraction(const Neutrino& neutrino);
    std::vector<double> getTransmittedFraction(const double energy, const double interactionLength);
    
}
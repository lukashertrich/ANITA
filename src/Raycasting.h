#pragma once
#include "Vector2.h"
#include "Vector3.h"
#include "DataRaster.h"
namespace anita{
    anita::Vector3<double> getPositionOfTau(const anita::Vector3<double>& pos, const anita::Vector3<double>& dir, const double tau);    
    anita::Vector2<double> getProjectionCoordinates(const anita::Vector3<double>& pos);    
    anita::Vector2<double> getDataCoordinates(const anita::Vector2<double>& projCoords);    
    anita::Vector2<double> getNormalizedUVCoordinates(const anita::Vector2<double>& dataCoords);
    double getDataValue(const anita::Vector3<double>& pos, anita::DataRaster<float>& dataRaster);
    double getEllipsoidalRadius(const anita::Vector3<double>& pos);    
    double getDistanceToSurface(const anita::Vector3<double>& pos, anita::DataRaster<float> dataRaster);
}
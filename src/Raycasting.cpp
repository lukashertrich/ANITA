#include "Raycasting.h"
#include "BEDMAP.h"

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

double anita::getDataValue(const anita::Vector3<double>& pos, anita::DataRaster<float>& dataRaster){	
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

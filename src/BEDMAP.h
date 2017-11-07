#pragma once
#include "DataRaster.h"
namespace anita{
    const int DATA_ROWS = 6667;
    const int DATA_COLUMNS = DATA_ROWS;
    const int DATA_INTERVAL = 1000; // meters
    const int PROJECTION_PLANE_LAT = 71; // degrees
    
    const double Z_PLANE = sin(PROJECTION_PLANE_LAT*(M_PI/180.)) * POLAR_EARTH_RADIUS;
    
    const double EPSILON = 0.001; // 1mm for numerical gradient
    const double RAYSTEP = 500.0; // meters
    const double UPPER_SURFACE_BOUND = 1.0 + (5000 / POLAR_EARTH_RADIUS); // Roughly higher than Mt Vinson in normalized radius
    const double LOWER_SURFACE_BOUND = 1.0 - (3000 / POLAR_EARTH_RADIUS); // Roughly lower than Bentley Subglacial Trench in normalized radius

    // DataRasters
    anita::DataRaster<float> geoidDataRaster;
    anita::DataRaster<float> bedDataRaster;
    anita::DataRaster<float> surfaceDataRaster;
    anita::DataRaster<float> iceThicknessDataRaster;    
}
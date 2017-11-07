#pragma once
#include "DataRaster.h"
#include <cmath>
namespace anita{
    constexpr int DATA_ROWS = 6667;
    constexpr int DATA_COLUMNS = DATA_ROWS;
    constexpr int DATA_INTERVAL = 1000; // meters
    constexpr int PROJECTION_PLANE_LAT = 71; // degrees
    
    constexpr double Z_PLANE = sin(PROJECTION_PLANE_LAT*(acos(-1.)/180.)) * POLAR_EARTH_RADIUS;
    
    constexpr double EPSILON = 0.001; // 1mm for numerical gradient
    constexpr double RAYSTEP = 500.0; // meters
    constexpr double UPPER_SURFACE_BOUND = 1.0 + (5000. / POLAR_EARTH_RADIUS); // Roughly higher than Mt Vinson in normalized radius
    constexpr double LOWER_SURFACE_BOUND = 1.0 - (3000. / POLAR_EARTH_RADIUS); // Roughly lower than Bentley Subglacial Trench in normalized radius
    
}
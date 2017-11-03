#pragma once
namespace anita{
    const double MEAN_EARTH_RADIUS = 6.3710e6; // meters
    const double EQUATORIAL_EARTH_RADIUS = 6378137; // meters
    const double EQUATORIAL_EARTH_RADIUS_SQR = EQUATORIAL_EARTH_RADIUS * EQUATORIAL_EARTH_RADIUS; // meters squared
    const double INVERSE_FLATTENING = 298.257223563; // Used to determine polar radius
    // WGS84 polar radius is defined by inverse flattening term
    const double POLAR_EARTH_RADIUS = EQUATORIAL_EARTH_RADIUS * (1.0 - (1.0 / INVERSE_FLATTENING)); // meters
    const double POLAR_EARTH_RADIUS_SQR = POLAR_EARTH_RADIUS * POLAR_EARTH_RADIUS; // meters squared
    
    // Normalized radii of density profile shells listed from the surface inwards to the core.
    // Obtained from spherical PREM to be mapped to normalized WGS84 ellipsoid equation
    const double DENSITY_PROFILE_RADII[]{
        1.0,								// Ocean (Exclude ocean from antarctic side of ray traversal by using constant crust density up to bedrock)
        6368000.0 / MEAN_EARTH_RADIUS,		// Crust
        6356000.0 / MEAN_EARTH_RADIUS,		// Crust
        6346600.0 / MEAN_EARTH_RADIUS,		// LVZ/LID
        6151000.0 / MEAN_EARTH_RADIUS,		// Transition Zone
        5971000.0 / MEAN_EARTH_RADIUS,		// Transition Zone
        5771000.0 / MEAN_EARTH_RADIUS,		// Transition Zone
        5701000.0 / MEAN_EARTH_RADIUS,		// Lower Mantle
        3480000.0 / MEAN_EARTH_RADIUS,		// Outer Core
        1221500.0 / MEAN_EARTH_RADIUS		// Inner Core
    } ;

    // x^0, x^1, x^2, x^3, etc. from PREM
    const std::vector<std::vector<double>> DENSITY_PROFILE_POLYNOMIAL_COEFFICIENTS{	// Given in kg / m^3
        std::vector<double>{1020.0},								// Ocean
        std::vector<double>{2600.0},								// Crust
        std::vector<double>{2900.0},								// Crust
        std::vector<double>{2691.0, 692.4},							// LVZ/LID
        std::vector<double>{7108.9, -3804.5},						// Transition Zone
        std::vector<double>{11249.4, -8029.8},						// Transition Zone
        std::vector<double>{5319.7, -1483.6},						// Transition Zone
        std::vector<double>{7956.5, -6476.1, 5528.3, -3080.7}, 		// Lower Mantle
        std::vector<double>{12581.5, -1263.8, -3642.6, -5528.1},	// Outer Core
        std::vector<double>{13088.5, -8838.1}						// Inner Core
	};
}
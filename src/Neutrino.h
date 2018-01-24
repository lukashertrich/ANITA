#pragma once
#include <cmath>
#include <vector>
#include "Vector3d.h"

namespace anita{
    class Neutrino{
        public:
            const Vector3d destination;
            const Vector3d reverseFlightDir;
            const double energy; // eV
            Neutrino(Vector3d& destination, Vector3d& reverseFlightDir, double& energy);
    };

    // Cross section constants from fermilab paper, converted to SI units (square meters & eV)
    
    constexpr double NEUTRINO_XSECTION_POW = 0.363; // Used by cross section function

    constexpr double NEUTRINO_XSECTION_CC = 5.53e-40; // meters squared
    constexpr double NEUTRINO_XSECTION_NC = 2.31e-40;
    constexpr double NEUTRINO_XSECTION_TOTAL = 7.84e-40;

    constexpr double ANTINEUTRINO_XSECTION_CC = 5.52e-40;
    constexpr double ANTINEUTRINO_XSECTION_NC = 2.29e-40;
    constexpr double ANTINEUTRINO_XSECTION_TOTAL = 7.80e-40;

    constexpr double NUCLEON_MASS = 1.67375e-27; // kg

    std::vector<double> getCrossSections(const double energy);
}
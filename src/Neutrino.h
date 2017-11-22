#pragma once
#include "Vector3.h"
namespace anita{
    class Neutrino{
        public:
            const Vector3<double> destination;
            const Vector3<double> reverseFlightDir;
            const double energy;
            Neutrino(Vector3<double>& destination, Vector3<double>& reverseFlightDir, double& energy);
    };
}
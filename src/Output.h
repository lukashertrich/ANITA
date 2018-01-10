#pragma once
#include "Vector3.h"
namespace anita{
    void outputFluxMap(const Vector3<double> position, const double energy, const unsigned long long resolution);
    void outputAngularTrace(const Vector3<double>& position, const double energy, const unsigned long long resolution);
    void outputEnergySpectrumAngularTrace(const Vector3<double>& position, const double minEnergy, const double maxEnergy, const unsigned long long resolution);
}
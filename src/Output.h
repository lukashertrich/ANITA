#pragma once
#include "Vector3d.h"
namespace anita{
    void outputFluxMap(const Vector3d& position, const double energy, const unsigned long long resolution);
    void outputAngularTrace(const Vector3d& position, const double energy, const unsigned long long resolution);
    void outputEnergySpectrumAngularTrace(const Vector3d& position, const double minEnergy, const double maxEnergy, const unsigned long long resolution);
}
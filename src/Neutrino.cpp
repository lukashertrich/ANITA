#include "Neutrino.h"

namespace anita{
    Neutrino::Neutrino(Vector3<double>& destination, Vector3<double>& reverseFlightDir, double& energy): destination(destination), reverseFlightDir(reverseFlightDir), energy(energy){}
    std::vector<double> getCrossSections(const double energy){
        std::vector<double> crossSections;
        crossSections.reserve(6);

        double factor = pow(energy / 1.0e9, NEUTRINO_XSECTION_POW);

        crossSections.push_back(NEUTRINO_XSECTION_CC * factor);
        crossSections.push_back(NEUTRINO_XSECTION_NC * factor);
        crossSections.push_back(NEUTRINO_XSECTION_TOTAL * factor);

        crossSections.push_back(ANTINEUTRINO_XSECTION_CC * factor);
        crossSections.push_back(ANTINEUTRINO_XSECTION_NC * factor);
        crossSections.push_back(ANTINEUTRINO_XSECTION_TOTAL * factor);

        return crossSections;
    }
}
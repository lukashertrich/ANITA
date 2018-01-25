#include "Neutrino.h"

namespace anita{
    Neutrino::Neutrino(Vector3d& destination, Vector3d& reverseFlightDir, double& energy): destination(destination), reverseFlightDir(reverseFlightDir), energy(energy){}
    std::vector<double> getCrossSections(const double energy){
        std::vector<double> crossSections;
        crossSections.reserve(6);

        double factor = pow(energy / 1.0e9, NEUTRINO_XSECTION_POW);

        // ***** DEPRECATED *****
        //crossSections.push_back(NEUTRINO_XSECTION_CC * factor);
        //crossSections.push_back(NEUTRINO_XSECTION_NC * factor);
        constexpr double logOfCentimeterSqrXSection = 36.0 * log(10.0);

        crossSections.push_back(1.0e-40 * exp(logOfCentimeterSqrXSection - 98.8 * pow(log(energy / 1.0e9), -0.0964)));
        crossSections.push_back(crossSections.back() * 2.39); // No direct model was given so a simple multiplier is used
        crossSections.push_back(NEUTRINO_XSECTION_TOTAL * factor); // No longer valid

        crossSections.push_back(ANTINEUTRINO_XSECTION_CC * factor); // No longer valid
        crossSections.push_back(ANTINEUTRINO_XSECTION_NC * factor); // No longer valid
        crossSections.push_back(ANTINEUTRINO_XSECTION_TOTAL * factor); // No longer valid

        return crossSections;
    }
}
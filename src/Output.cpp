#include "Output.h"
#include "Raycasting.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Vector3d.h"
#include "Constants.h"

namespace anita{
    void outputFluxMap(const Vector3d& position, const double energy, const unsigned long long resolution){
		std::ofstream outFile("fluxMap.dat");
		double theta, phi, interactionLength;
		for(unsigned long long xy = 0; xy < resolution; xy++){
			phi = 2 * M_PI * (1.0 * xy / resolution);
			for(unsigned long long z = 0; z < resolution; z++){
				theta = acos(2.0 * (1.0 * z / resolution) - 1.0);
				auto direction = Vector3d{sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
				interactionLength = getInteractionLength(position, direction);
				auto transmittedFractions = getTransmittedFraction(energy, interactionLength);
				outFile << std::setprecision(20) << phi << "	" << theta << "	" << std::scientific 
				<< transmittedFractions[2] << std::endl;
			}
			outFile << std::endl;
		}
		outFile.close();
	}

	void outputAngularTrace(const Vector3d& position, const double energy, const unsigned long long resolution){
		std::ofstream outFile("angularTrace.dat");
		outFile << std::fixed << std::setprecision(20) << std::scientific;
		double theta, interactionLength;
		for (unsigned long long i = 0; i < resolution; i++){
			theta = M_PI_2 * (1.0 * i / resolution);
			auto direction = Vector3d{cos(theta), 0.0, sin(theta)};
			interactionLength = getInteractionLength(position, direction);
			auto transmittedFractions = getTransmittedFraction(energy,interactionLength);
			outFile << 90.0 - (theta * 180.0 / M_PI) << "	";
			outFile << transmittedFractions[0] << "	";
			outFile << transmittedFractions[1] << "	";
			outFile << transmittedFractions[2] << "	";
			outFile << transmittedFractions[3] << "	";
			outFile << transmittedFractions[4] << "	";
			outFile << transmittedFractions[5] << std::endl;
		}
        outFile.close();
	}

    void outputEnergySpectrumAngularTrace(const Vector3d& position, const double minEnergy, const double maxEnergy, const unsigned long long resolution){
        // Enforce energy bounds
        auto minE = std::min(minEnergy, maxEnergy);
        auto maxE = std::max(minEnergy, maxEnergy);
        auto range = maxE - minE;
        auto interval = range / (double)resolution;
        std::ofstream outFile("energySpectrumAngularTrace.dat");
		outFile << std::fixed << std::setprecision(20) << std::scientific;
		double theta, interactionLength, energy;
        for(unsigned long long j = 0; j <= resolution; j++){
            energy = minE + j * interval;
            for (unsigned long long i = 0; i < resolution; i++){
                theta = M_PI_2 * (1.0 * i / resolution);
                auto direction = Vector3d{cos(theta), 0.0, sin(theta)};
                interactionLength = getInteractionLength(position, direction);
                auto transmittedFractions = getTransmittedFraction(energy,interactionLength);
                outFile << energy << "  ";
                outFile << theta << "	";
                outFile << transmittedFractions[0] << "	";
                outFile << transmittedFractions[1] << "	";
                outFile << transmittedFractions[2] << "	";
                outFile << transmittedFractions[3] << "	";
                outFile << transmittedFractions[4] << "	";
                outFile << transmittedFractions[5] << std::endl;
            }
        // Extra break for gnuplot isoline
        outFile << std::endl;
        }
        outFile.close();
    }
}
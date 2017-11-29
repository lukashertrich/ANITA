
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "Diagnostics.h"
#include "PREM.h"
#include "DataRaster.h"
#include "Vector3.h"
#include "Raycasting.h"

namespace anita{

    void printConstants(){
        std::cout << "Enumeration of physical constants used in code:" << std::endl;
        std::cout << "Mean Earth radius:			" << MEAN_EARTH_RADIUS << " meters." << std::endl;
        std::cout << "Polar Earth radius:			" << POLAR_EARTH_RADIUS << " meters."<< std::endl;
        std::cout << std::endl;
    }

    void printUsage(){
        std::cout << "Usage: anita_earthmodel2 [Nevt][ spectrum='ESS' or 'Enu'][Enu(eV)][maxdepth(m)][cross section factor]"<< std::endl;
        std::cout << "If single energy 'Enu' then argv[3] is the neutrino energy, else ignored."<< std::endl;
        std::cout << std::endl;
    }

    void testDataRaster(){
        // DataRaster testRaster(std::string("gl04c_geiod_to_wgs84.flt"));
    }

    void testInteractionLength(){
        constexpr double energy = 1.0e19;
        auto position = Vector3<double>(0.0, 1.0, -POLAR_EARTH_RADIUS - 1.0); // South pole with small y offset to avoid origin on raycast
        auto direction = Vector3<double>(1.0, 0.0, 0.01);
        double interactionLength = getInteractionLength(position, direction);
        auto transmittedFraction = getTransmittedFraction(energy, interactionLength);
        auto crossSections = getCrossSections(energy);
        
        std::cout << "Fraction of initial flux transmitted through Earth model at: " << energy << " eV," << std::endl;
        std::cout << "with interaction length: " << interactionLength << " kg / m^2" << std::endl;
        std::cout << std::fixed << std::setprecision(4) << std::scientific;
        std::cout << "Charged current transmittance:                 " << transmittedFraction[0] << std::endl;
        std::cout << "Neutral current transmittance:                 " << transmittedFraction[1] << std::endl;
        std::cout << "Total transmittance:                           " << transmittedFraction[2] << std::endl;
        std::cout << "Charged current transmittance [Antineutrino]:  " << transmittedFraction[3] << std::endl;
        std::cout << "Neutral current transmittance [Antineutrino]:  " << transmittedFraction[4] << std::endl;
        std::cout << "Total transmittance [Antineutrino]:            " << transmittedFraction[5] << std::endl;
        std::cout << std::endl;
        std::cout << "Cross sections [m^2]:" << std::endl;
        std::cout << "Charged current:                 " << crossSections[0] << std::endl;
        std::cout << "Neutral current:                 " << crossSections[1] << std::endl;
        std::cout << "Total:                           " << crossSections[2] << std::endl;
        std::cout << "Charged current [Antineutrino]:  " << crossSections[3] << std::endl;
        std::cout << "Neutral current [Antineutrino]:  " << crossSections[4] << std::endl;
        std::cout << "Total [Antineutrino]:            " << crossSections[5] << std::endl;
    }



    // void testDensityTraversal(){
    //     std::cout << "Generating test data of density traversal for swept hemisphere through Earth from point at south pole and saving to 'diagnostic.dat'." << std::endl;
    //     std::ofstream outFile;
    //     outFile.precision(20);
    //     std::string outfilePath = "diagnostic.dat";
    //     outFile.open(outfilePath);
    //     auto position = std::vector<double>{0,0,-POLAR_EARTH_RADIUS - 1000}; // geometric south pole 
    //     auto direction = std::vector<double>{0,0,1.00}; // +Z direction
    //     unsigned int n = 100;
    //     double phi;
    //     double theta;
    //     for(unsigned int i = 0; i <= n; i++){
    //         phi = (2.0 * M_PI / n) * i;
    //         for (unsigned int j = 0; j <=n; j++){			
    //             theta = j *((M_PI) * (1.0 / n));
    //             direction[0] = sin(theta)*cos(phi);
    //             direction[1] = sin(theta)*sin(phi);
    //             direction[2] = cos(theta);
    //             outFile << phi << "	" << theta << "	" << getDensityTraversed(position, direction) << std::endl;			
    //         }
    //         outFile << std::endl;
    //     }
    //     outFile.close();
    // }

}

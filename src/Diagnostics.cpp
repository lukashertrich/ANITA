
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "Diagnostics.h"
#include "PREM.h"
#include "DataRaster.h"

void anita::printConstants(){
    std::cout << "Enumeration of physical constants used in code:" << std::endl;
    std::cout << "Mean Earth radius:			" << anita::MEAN_EARTH_RADIUS << " meters." << std::endl;
    std::cout << "Polar Earth radius:			" << anita::POLAR_EARTH_RADIUS << " meters."<< std::endl;
    std::cout << std::endl;
}

void anita::printUsage(){
    std::cout << "Usage: anita_earthmodel2 [Nevt][ spectrum='ESS' or 'Enu'][Enu(eV)][maxdepth(m)][cross section factor]"<< std::endl;
    std::cout << "If single energy 'Enu' then argv[3] is the neutrino energy, else ignored."<< std::endl;
    std::cout << std::endl;
}

void anita::testDataRaster(){
    anita::DataRaster<float> testRaster("gl04c_geiod_to_wgs84.flt");
}

// void anita::testDensityTraversal(){
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
//             outFile << phi << "	" << theta << "	" << anita::getDensityTraversed(position, direction) << std::endl;			
//         }
//         outFile << std::endl;
//     }
//     outFile.close();
// }



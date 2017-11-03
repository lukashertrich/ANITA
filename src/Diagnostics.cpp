#include "Diagnostics.h"
void printConstants(){
    std::cout << "Enumeration of physical constants used in code:" << std::endl;
    std::cout << "Distance to horizon:			" << anita::DISTANCE_TO_HORIZON << " meters." << std::endl;
    std::cout << "Mean Earth radius:			" << anita::MEAN_EARTH_RADIUS << " meters." << std::endl;
    std::cout << "Polar Earth radius:			" << anita::POLAR_EARTH_RADIUS << " meters."<< std::endl;
    std::cout << "Mean Earth density:			" << anita::MEAN_EARTH_DENSITY << " kilograms per cubic meter." << std::endl;
    // std::cout << "Vacuum speed of light:			" << SPEED_OF_LIGHT << " meters per second." << std::endl;
    std::cout << std::endl;
}

void printUsage(){
    std::cout << "Usage: anita_earthmodel2 [Nevt][ spectrum='ESS' or 'Enu'][Enu(eV)][maxdepth(m)][cross section factor]"<< std::endl;
    std::cout << "If single energy 'Enu' then argv[3] is the neutrino energy, else ignored."<< std::endl;
    std::cout << std::endl;
}

void testDensityTraversal(){
    std::cout << "Generating test data of density traversal for swept hemisphere through Earth from point at south pole and saving to 'diagnostic.dat'." << std::endl;
    std::ofstream outFile;
    outFile.precision(20);
    std::string outfilePath = "diagnostic.dat";
    outFile.open(outfilePath);
    auto position = std::vector<double>{0,0,-POLAR_EARTH_RADIUS - 1000}; // geometric south pole 
    auto direction = std::vector<double>{0,0,1.00}; // +Z direction
    unsigned int n = 100;
    double phi;
    double theta;
    for(unsigned int i = 0; i <= n; i++){
        phi = (2.0 * M_PI / n) * i;
        for (unsigned int j = 0; j <=n; j++){			
            theta = j *((M_PI) * (1.0 / n));
            direction[0] = sin(theta)*cos(phi);
            direction[1] = sin(theta)*sin(phi);
            direction[2] = cos(theta);
            outFile << phi << "	" << theta << "	" << getDensityTraversed(position, direction) << std::endl;			
        }
        outFile << std::endl;
    }
    outFile.close();
}

void testData(){
    int n = 6667;
    int d = 66;
    std::ofstream outFile;
    outFile.open("map.dat");
    for(int y = 0; y < n; y += d){
        for(int x = 0; x < n; x += d){
            outFile << x * 1000 << "	" << y * 1000 << "	" << ((iceThicknessData[x + n*y] == -9999.0f) ? 0 : iceThicknessData[x + n*y]) << std::endl;
        }
        outFile << std::endl;
    }
    outFile.close();
}

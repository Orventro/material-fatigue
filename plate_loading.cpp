#include "mfem.hpp"
#include <string>
#include <iostream>
#include <filesystem>
#include "Environment.hpp"

int main(int argc, char* argv[]) {
    // Defining computation parameters
    #ifdef MFEM_USE_CUDA
        mfem::Device device("cuda");
    #else
        mfem::Device device("cpu");
    #endif

    if (argc < 4) {
        std::cout << "Provide run num, mesh and simulation directories.\n";
        return 0;
    }
    
    mfem::StopWatch timer;
    timer.Start();
    
    std::filesystem::create_directory(argv[3]);
    int runNum = std::stoi(argv[1]);
    std::string psiFile = "";
    if (std::stoi(argv[1]) > 0) 
        psiFile = std::string(argv[3]) + "/psi_" + std::to_string(runNum-1);
    Environment env(argv[2], psiFile);

    const int numOfSteps = 50000;
    const int saveStep = 1000;
    const double dt = 3e-9 / 16.0;

    for (uint i = 0; i < numOfSteps; i++) {
        env.step(dt);

        if ( !((i+1) % saveStep) ) 
            env.saveDisplacement(std::string(argv[3]) + "/plateMesh_" + 
                                 std::to_string(runNum) + "_" + 
                                 std::to_string(i+1) + ".vtk");
    }
    
    env.savePsi(std::string(argv[3]) + "/psi_" + std::to_string(runNum));

    timer.Stop();
    std::cout << "\nProgram took " << timer.RealTime() << " seconds." << std::endl;
}

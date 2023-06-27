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
    // if (std::stoi(argv[1]) > 0) 
    //     psiFile = std::string(argv[3]) + "/psi_" + std::to_string(runNum-1);

    const int numOfSteps = 1000;
    const int saveStep = 100;
    const double dt = 3e-9 / 16.0;
    const int fatigueCycles = 300, fatigueSteps = 800;
    bool saveSteps = false;
    
    std::ofstream runInfo(std::string(argv[3]) + "/runInfo");
    runInfo << "mesh=" << argv[2] << std::endl;
    runInfo << "numOfSteps=" << numOfSteps << std::endl;
    runInfo << "saveStep=" << saveStep << std::endl;
    runInfo << "dt=" << dt << std::endl;
    runInfo << "fatigueCycles=" << fatigueCycles << std::endl;
    runInfo << "fatigueSteps=" << fatigueSteps << std::endl;

    Environment env(argv[2], runInfo);//, psiFile);
    env.savePsi(std::string(argv[3]) + "/psi_0");

    runInfo.close();

    std::string snapshotPath = std::string(argv[3]) + "/envSnapshot";
    
    if (runNum == 0) {
        for (int i = 0; i < numOfSteps; i++) {
            env.step(dt);

            if ( !((i+1) % saveStep) ) {
                std::string filename = std::string(argv[3]) + "/plateMesh_" + 
                                    std::to_string(i+1) + ".vtk";
                std::cout << "\rSaving data to " << filename << std::endl;
                std::ofstream vtkFile(filename);
                env.saveDisplacement(vtkFile);
                env.saveSigma(vtkFile);
                env.calcPsi();
                env.saveVars(vtkFile);
                vtkFile.close();
            }
        }
        env.saveEnv(snapshotPath+"_0");
        // env.savePsi(std::string(argv[3]) + "/psi_" + std::to_string(runNum));
    } 
    // else {
    //     env.loadPsi(std::string(argv[3]) + "/psi_" + std::to_string(runNum-1), 0);
    //     for (int i = 0; i < numOfSteps; i++) {
    //         env.step(dt);

    //         if ( (i+1) == numOfSteps ) {
    //             std::string filename = std::string(argv[3]) + "/mesh_iter_" + 
    //                                 std::to_string(runNum) + ".vtk";
    //             std::cout << "\rSaving data to " << filename << std::endl;
    //             std::ofstream vtkFile(filename);
    //             env.saveDisplacement(vtkFile);
    //             env.saveSigma(vtkFile);
    //             env.calcPsi();
    //             env.saveVars(vtkFile);
    //             vtkFile.close();
    //         }
    //     }
    //     env.savePsi(std::string(argv[3]) + "/psi_" + std::to_string(runNum));
    // }
    // return 0;
    
    for (int i = runNum+1; i <= fatigueCycles; i++) {
        env.loadEnv(snapshotPath+"_0");//+std::to_string(i-1));
        env.loadPsi(std::string(argv[3]) + "/psi_" + std::to_string(i-1), numOfSteps*dt);

        // if(i % 5 == 0) saveSteps = true;
        // else saveSteps = false;
        for (int j = 0; j < fatigueSteps; j++) {
            env.step(dt, true);
            if (saveSteps && j%saveStep==(saveStep-1)) {
                std::string filename = std::string(argv[3]) + "/mesh_nn_" + std::to_string(i) + "_" + std::to_string(j+1) + ".vtk";
                std::cout << "\rSaving data to " << filename << std::endl;
                std::ofstream vtkFile(filename);
                env.saveDisplacement(vtkFile);
                env.saveSigma(vtkFile);
                env.calcPsi();
                env.saveVars(vtkFile);
                vtkFile.close();
            }
        }
        std::string filename = std::string(argv[3]) + "/mesh_fc_" + std::to_string(i) + ".vtk";
        std::cout << "\rSaving data to " << filename << std::endl;
        std::ofstream vtkFile(filename);
        env.saveDisplacement(vtkFile);
        env.saveSigma(vtkFile);
        env.calcPsi();
        env.saveVars(vtkFile);
        vtkFile.close();
        env.savePsi(std::string(argv[3]) + "/psi_" + std::to_string(i));
    }
    

    timer.Stop();
    std::cout << "\nProgram took " << timer.RealTime() << " seconds." << std::endl;
}

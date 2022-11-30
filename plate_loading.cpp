#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>
#include <string>
#include "ElasticWaveOperator.hpp"
#include "SigmaCoefficient.hpp"
#include "FatigueCoefficient.hpp"

using namespace mfem;

void displacementOnBoundary(const Vector& x, double time, Vector& d) {
    const double freq = 3e6;
    const double amplitude = 4e-5; // 4 mm
    //d(0) = 0.0;
    // d(1) = amplitude * tanh(time * freq * 2. * M_PI);
    d(1) = amplitude * sin(time * freq * 2. * M_PI);
}

int main(int argc, char* argv[]) {
    // Defining computation parameters
    #ifdef MFEM_USE_CUDA
        Device device("cuda");
    #else
        Device device("cpu");
    #endif

    if (argc < 3) {
        std::cout << "No mesh file and save directory provided.\n";
        return 0;
    }
    
    StopWatch timer;
    timer.Start();
    const int ORDER = 4;
    double calculationTime = 0.0;
    Mesh mesh(argv[1]);
    std::cout << "Mesh loaded.\n";

    /*
     * Boundary attributes:
     *    top:      1
     *    left:     2
     *    bottom:   3
     *    right:    4
     *    circle:   5
     */
    Array<int> boundaryAttribute(mesh.bdr_attributes.Size());
    boundaryAttribute = 0;
    boundaryAttribute[0] = 1;

    double dt = 3e-9 / 16.0; // 9e-9 for 2nd order // order = 2, Newmark scheme => k = 0.577, k_order = 2, dt < k * h_min / Cp / k_order = 1.9e-8
    const uint NumberOfTimeSteps = 1000000; // 2e-6 [s] takes the P-wave to move across the plate

    // Defining material parameters
    const double E0 = 116e9; // Young's modulus [Pa]
    const double nu0 = 0.32; // Poisson's ratio
    const double rho0 = 4370; // density [kg/m^3]
    const double lambda0 = E0 * nu0 / ((1 + nu0) * (1 - 2 * nu0));
    const double mu0 = E0 / (2 * (1 + nu0));

    // Define FE space
    H1_FECollection fe_collection(ORDER, mesh.Dimension(), BasisType::GaussLobatto);
    FiniteElementSpace fespace1(&mesh, &fe_collection, 1);
    FiniteElementSpace fespace2(&mesh, &fe_collection, 2);
    FiniteElementSpace fespace3(&mesh, &fe_collection, 3);
    GridFunction displacementGF(&fespace2);
    GridFunction sigmaGF(&fespace3);
    sigmaGF = 0;
    GridFunction lambdaGF(&fespace1);
    GridFunction muGF(&fespace1);

    ConstantCoefficient rhoCoef(4500.0);

    SigmaCoefficient sigmaCoeff(displacementGF);
    FatigueCoefficient fatigueCoeff(sigmaGF, fespace1);
    LambdaCoefficient lambdaCoeff(fatigueCoeff);
    MuCoefficient muCoeff(fatigueCoeff);
    sigmaCoeff.setMaterial(&lambdaCoeff, &muCoeff);

    ElasticWaveOperator2D oper(fespace2, rhoCoef, lambdaCoeff, muCoeff);
    oper.SetTime(calculationTime);
    std::cout << "Operator created.\n";
    mfem::CentralDifferenceSolver ode_solver;
    ode_solver.Init(oper);

    Vector u(fespace2.GetVSize()), du(fespace2.GetVSize());
    u  = 0.0;
    du = 0.0;

    // std::cout << fespace3.GetVSize() << std::endl;
    // return 0;

    VectorFunctionCoefficient displacementCoeff(2, displacementOnBoundary);
    displacementCoeff.SetTime(calculationTime);

    std::cout << "Preparation took " << timer.RealTime() << " seconds.\nStarting forward simulation.\n";

    std::filesystem::create_directory(argv[2]);
    const uint saveStep = 100;
    const uint fatigueStep = 10;
    const int REF = 4;
    for (uint i = 0; i < NumberOfTimeSteps; i++) {
        // Setting correct values on boundary
        ode_solver.Step(u, du, calculationTime, dt);

        // du *= 1-2e-4; // only for conevrging simualtions

        displacementCoeff.SetTime(calculationTime);
        displacementGF.SetFromTrueDofs(u);
        displacementGF.ProjectBdrCoefficient(displacementCoeff, boundaryAttribute);
        displacementGF.GetTrueDofs(u);

        if ( !((i+1) % fatigueStep) ) {
            sigmaGF.ProjectCoefficient(sigmaCoeff);

            fatigueCoeff.Update();
            
            // fatigueGF.ProjectCoefficient(fatigueCoeff);
            // std::cout << fatigueCoeff.GetCalls() << ' ' << fatigueGF.Size() << ' ' << fatigueCoeff.GetMaxNo() << std::endl;
            // break;
        }

        if ( !((i+1) % saveStep) ) {
            std::string filename = std::string(argv[2]) + "/plateMesh_" + std::to_string(i+1) + ".vtk";
            std::cout << "\rSaving data to " << filename << std::endl;

            std::ofstream vtkFile(filename);
            mesh.PrintVTK(vtkFile, REF);
            sigmaGF.SaveVTK(vtkFile, "sigma", REF);
            displacementGF.SaveVTK(vtkFile, "displacement", REF);
            lambdaGF.ProjectCoefficient(lambdaCoeff);
            lambdaGF.SaveVTK(vtkFile, "lambda", REF);
            muGF.ProjectCoefficient(muCoeff);
            muGF.SaveVTK(vtkFile, "mu", REF);
        }
        
    }

    timer.Stop();
    std::cout << "\nProgram took " << timer.RealTime() << " seconds." << std::endl;
}

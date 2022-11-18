#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>
#include <string>
#include "ElasticWaveOperator.hpp"
#include "SigmaCoefficient.hpp"

using namespace mfem;

void displacementOnBoundary(const Vector& x, double time, Vector& d) {
    const double freq = 3e6;
    const double amplitude = 4e-5;
    //d(0) = 0.0;
    d(1) = amplitude * tanh(time * freq * 2. * M_PI);
    // d(1) = amplitude * sin(time * freq * 2. * M_PI);
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
    std::cout << "Loading mesh...\n";
    Mesh mesh(argv[1]);
    std::cout << "Mesh loaded.\n";

    // mesh.attributes.Print();
    /*
     * Boundary attributes:
     *    top:      1
     *    left:     2
     *    bottom:   3
     *    right:    4
     *    circle:   5
     */
    double dt = 3e-9 / 16.0; // 9e-9 for 2nd order // order = 2, Newmark scheme => k = 0.577, k_order = 2, dt < k * h_min / Cp / k_order = 1.9e-8
    const uint NumberOfTimeSteps = 1000000; // 2e-6 [s] takes the P-wave to move across the plate

    // Define FE space
    H1_FECollection fe_collection(ORDER, mesh.Dimension(), BasisType::GaussLobatto);
    FiniteElementSpace fespace(&mesh, &fe_collection, mesh.Dimension());
    GridFunction displacement(&fespace);

    // Defining material parameters
    ConstantCoefficient rhoCoef(4500.0);
    Vector lambda(2);
    lambda(0) = 7.7e10;
    lambda(1) = 7.7e8;
    PWConstCoefficient lambdaCoef(lambda);
    Vector mu(2);
    mu(0) = 4.4e10;
    mu(1) = 4.4e8;
    PWConstCoefficient muCoef(mu);

    std::cout << "Preparation took " << timer.RealTime() << " seconds.\nStarting forward simulation.\n";

    // Defining operator for spatial discretization and ODE solver for time discretization
    ElasticWaveOperator2D oper(fespace, rhoCoef, lambdaCoef, muCoef);
    printf("Operator done!\n");
    oper.SetTime(calculationTime);
    mfem::CentralDifferenceSolver ode_solver;
    ode_solver.Init(oper);

    // Initializing (to zeros) displacement u and velocity du
    Vector u(fespace.GetVSize()), du(fespace.GetVSize());
    u  = 0.0;
    du = 0.0;
    // u.UseDevice(true);
    // du.UseDevice(true);

    FiniteElementSpace scalarFESpace(&mesh, &fe_collection, 1);
    VectorFunctionCoefficient displacementCoeff(2, displacementOnBoundary);
    displacementCoeff.SetTime(calculationTime);

    
    
    Array<int> boundaryAttribute(mesh.bdr_attributes.Size());
    boundaryAttribute = 0;
    boundaryAttribute[0] = 1; // top boundary

    SigmaCoefficient sigmaCoef(displacement, lambdaCoef, muCoef);
    GridFunction sigmaGF(&fespace);


    std::filesystem::create_directory(argv[2]);
    const uint saveStep = 100;
    const uint fatigueStep = 10;
    const int REF = 4;
    for (uint i = 0; i < NumberOfTimeSteps; i++) {
        // Setting correct values on boundary
        ode_solver.Step(u, du, calculationTime, dt);

        du *= 1-2e-4; // only for conevrging simualtions

        displacementCoeff.SetTime(calculationTime);
        displacement.SetFromTrueDofs(u);
        displacement.ProjectBdrCoefficient(displacementCoeff, boundaryAttribute);
        displacement.GetTrueDofs(u);

        if ( !((i+1) % saveStep) ) {
            std::string filename = std::string(argv[2]) + "/plateMesh_" + std::to_string(i+1) + ".vtk";
            std::cout << "\rSaving data to " << filename << std::endl;

            std::ofstream vtkFile(filename);
            mesh.PrintVTK(vtkFile, REF);
            std::cout << "Mesh saved.\n";
            sigmaCoef.setComponent(0);
            sigmaGF.ProjectCoefficient(sigmaCoef);
            std::cout << "projected sigma" << std::endl;
            sigmaGF.SaveVTK(vtkFile, "sigma", REF);
            std::cout << "saved sigma" << std::endl;
            displacement.SaveVTK(vtkFile, "displacement", REF);
            
        }

        std::cout << "\rFinished time step " << i << std::flush;
    }

    timer.Stop();
    std::cout << "\nProgram took " << timer.RealTime() << " seconds." << std::endl;
}

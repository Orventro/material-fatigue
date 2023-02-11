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

#define p2(a) ((a)*(a))

/*
* Boundary attributes:
*    top:      1
*    left:     2
*    bottom:   3
*    right:    4
*    circle:   5
*/
void displacementOnBoundary(const Vector& x, double time, Vector& d) {
    const double freq = 3e6;
    const double amplitude = 4e-5; // 4 mm
    //d(0) = 0.0;
    // d(1) = amplitude * tanh(time * freq * 2. * M_PI);
    d(1) = amplitude * sin(time * freq * 2. * M_PI);
    // d(2) = -amplitude * sin(time * freq * 2. * M_PI);
}

int main(int argc, char* argv[]) {
    // Defining computation parameters
    #ifdef MFEM_USE_CUDA
        Device device("cuda");
    #else
        Device device("cpu");
    #endif

    if (argc < 4) {
        std::cout << "Provide run num, mesh and simulation directories.\n";
        return 0;
    }
    
    StopWatch timer;
    timer.Start();
    const int ORDER = 4;
    double calculationTime = 0.0;
    int numRun = std::stoi(argv[1]);
    Mesh mesh(argv[2]);

    std::cout << "Mesh loaded.\n";

    Array<int> boundaryAttribute(mesh.bdr_attributes.Size());
    boundaryAttribute = 0;
    boundaryAttribute[0] = 1;
    boundaryAttribute[2] = 1;

    double dt = 3e-9 / 16.0; // 9e-9 for 2nd order // order = 2, Newmark scheme => k = 0.577, k_order = 2, dt < k * h_min / Cp / k_order = 1.9e-8
    const uint NumberOfTimeSteps = 50'000; // 2e-6 [s] takes the P-wave to move across the plate

    // Defining material parameters
    const double E0 = 116e9; // Young's modulus [Pa]
    const double nu0 = 0.32; // Poisson's ratio
    const double beta = 0.31, sigma_u = 340e6, sigma_B = 1160e6, psi_star = 0.98, gamma = 0.5;

    // Define FE space
    H1_FECollection fe_collection(ORDER, mesh.Dimension(), BasisType::GaussLobatto);
    FiniteElementSpace fespace1(&mesh, &fe_collection, 1);
    FiniteElementSpace fespace2(&mesh, &fe_collection, 2);
    FiniteElementSpace fespace3(&mesh, &fe_collection, 3);
    GridFunction displacementGF(&fespace2);
    GridFunction sigmaGF(&fespace1);
    sigmaGF = 0;
    GridFunction lambdaGF(&fespace1);
    GridFunction muGF(&fespace1);

    Vector mu, lambda, psi;
    mu.SetSize(mesh.GetNE());
    lambda.SetSize(mesh.GetNE());
    psi.SetSize(mesh.GetNE());
    if(numRun > 0) {
        std::string psiFilename = "psi_" + std::to_string(numRun-1);
        std::ifstream psiFile(std::string(argv[3]) + "/" + psiFilename);
        psi.Load(psiFile);
        for(int i = 0; i < psi.Size(); i++) {
            double E = E0*(1 - 0.5*psi(i))*((psi_star > psi(i)) + 1e-3);
            lambda(i) = E * nu0 / ((1 + nu0) * (1 - 2 * nu0));
            mu(i) = E / (2 * (1 + nu0));
        }
        printf("Successfully loaded '%s'\n", psiFilename.c_str());
    } else {
        for(int i = 0; i < mu.Size(); i++) {
            lambda(i) = E0 * nu0 / ((1 + nu0) * (1 - 2 * nu0));
            mu(i) = E0 / (2 * (1 + nu0));
        }
    }
    
    PWConstCoefficient muCoeff(mu);
    PWConstCoefficient lambdaCoeff(lambda);
    lambdaGF.ProjectCoefficient(lambdaCoeff);
    muGF.ProjectCoefficient(muCoeff);
    

    ConstantCoefficient rhoCoef(4370.0);

    // FatigueCoefficient fatigueCoeff(sigmaGF, fespace1);
    // LambdaCoefficient lambdaCoeff(fatigueCoeff);
    // MuCoefficient muCoeff(fatigueCoeff);
    // ConstantCoefficient lambdaCoeff(lambda0);
    // ConstantCoefficient muCoeff(mu0);
    SigmaCoefficient sigmaCoeff(displacementGF, lambdaCoeff, muCoeff);
   

    ElasticWaveOperator2D *oper = new ElasticWaveOperator2D(fespace2, rhoCoef, lambdaCoeff, muCoeff);
    oper->SetTime(calculationTime);
    std::cout << "Operator created.\n";
    mfem::CentralDifferenceSolver ode_solver;
    ode_solver.Init(*oper);

    Vector u(fespace2.GetVSize()), du(fespace2.GetVSize());
    u  = 0.0;
    du = 0.0;

    // std::cout << fespace3.GetVSize() << std::endl;
    // return 0;

    VectorFunctionCoefficient displacementCoeff(2, displacementOnBoundary);
    displacementCoeff.SetTime(calculationTime);

    std::cout << "Preparation took " << timer.RealTime() << " seconds.\nStarting forward simulation.\n";

    std::filesystem::create_directory(argv[3]);
    const uint saveStep = 1000;
    const int REF = 4;

    Vector sigmaMin(mesh.GetNE()), sigmaMax(mesh.GetNE());
    sigmaMin = 0.;
    sigmaMax = 0.;

    for (uint i = 0; i < NumberOfTimeSteps; i++) {
        // Setting correct values on boundary
        ode_solver.Step(u, du, calculationTime, dt);

        // du *= 1-2e-4; // only for conevrging simulations

        displacementCoeff.SetTime(calculationTime);
        displacementGF.SetFromTrueDofs(u);
        displacementGF.ProjectBdrCoefficient(displacementCoeff, boundaryAttribute);
        displacementGF.GetTrueDofs(u);

        for(int i = 0; i < mesh.GetNE(); i++) {
            Vector center(2);
            mesh.GetElementCenter(i, center);
            IntegrationPoint ip;
            ip.Set2(center.GetData());
            double s = sigmaCoeff.Eval(*mesh.GetElementTransformation(i), ip);
            sigmaMin(i) = std::min(s, sigmaMin(i));
            sigmaMax(i) = std::max(s, sigmaMax(i));
        }

        if ( !((i+1) % saveStep) ) {
            std::string filename = std::string(argv[3]) + "/plateMesh_" + 
                                   std::to_string(numRun) + "_" + 
                                   std::to_string(i+1) + ".vtk";
            std::cout << "\rSaving data to " << filename << std::endl;

            std::ofstream vtkFile(filename);
            mesh.PrintVTK(vtkFile, REF);
            sigmaGF.ProjectCoefficient(sigmaCoeff);
            sigmaGF.SaveVTK(vtkFile, "sigma", REF);
            displacementGF.SaveVTK(vtkFile, "displacement", REF);
            // lambdaGF.ProjectCoefficient(lambdaCoeff);
            lambdaGF.SaveVTK(vtkFile, "lambda", REF);
            // muGF.ProjectCoefficient(muCoeff);
            muGF.SaveVTK(vtkFile, "mu", REF);
            // fatigueCoeff.getPsi().SaveVTK(vtkFile, "psi", REF);
        }
        
        
    }

    // calculating damage
    double deltaN = 1e18;
    // Vector psi(mesh.GetNE());
    psi = 0.;

    for (int i = 0; i < mesh.GetNE(); i++) {
        double sx = sigmaMax(i);
        if (sx <= 0)
            continue;
        double sn = sigmaMin(i);
        double deltaSigma = sx - sn;
        double sigma_eq = std::sqrt((sx*deltaSigma)/2);
        if (sigma_eq <= sigma_u)
            continue;

        Vector elmentCenter(2);
        mesh.GetElementCenter(i, elmentCenter);
        if (elmentCenter.Norml2() > 5e-3)
            continue;

        double sigmaFrac = (sigma_eq - sigma_u)/(sigma_B - sigma_u);
        double B = 1e-3 * (std::pow(sigmaFrac, 1/beta)) / (2*(1-gamma));
        double x = 0.5/B * (0.5/(1-gamma) - std::pow(psi(i), 1-gamma)/(1-gamma) + std::pow(psi(i), 2-2*gamma)/(2-2*gamma));
        deltaN = std::min(x, deltaN);
    }

    for (int i = 0; i < mesh.GetNE(); i++) {
        double sx = sigmaMax(i);
        if (sx <= 0)
            continue;
        double sn = sigmaMin(i);
        double deltaSigma = sx - sn;
        double sigma_eq = std::sqrt((sx*deltaSigma)/2);
        if (sigma_eq <= sigma_u)
            continue;
        
        Vector elmentCenter(2);
        mesh.GetElementCenter(i, elmentCenter);
        if (elmentCenter.Norml2() > 5e-3)
            continue;

        double sigmaFrac = (sigma_eq - sigma_u)/(sigma_B - sigma_u);
        double B = 1e-3 * (std::pow(sigmaFrac, 1/beta)) / (2*(1-gamma));
        psi(i) = std::pow(1 - std::sqrt(p2(1-std::pow(psi(i), 1-gamma)) - 2*(1-gamma)*deltaN*B), 1/(1-gamma));
    }

    std::ofstream psiFile(std::string(argv[3]) + "/psi_" + std::to_string(numRun));
    psiFile << psi.Size() << std::endl;
    psi.Print(psiFile, 1);
    psi *= 100;
    psiFile.close();

    timer.Stop();
    std::cout << "\nProgram took " << timer.RealTime() << " seconds." << std::endl;
}

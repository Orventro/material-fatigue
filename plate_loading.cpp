#include "mfem.hpp"
#include <iostream>
#include <cmath>
#include <filesystem>
#include <string>
#include "ElasticWaveOperator.hpp"

using namespace mfem;


void displacementOnBoundary(const Vector& x, double time, Vector& d) {
    const double freq = 3e6;
    const double amplitude = 0.04;
    //d(0) = 0.0;
    // d(1) = amplitude * tanh(time * freq * 2. * M_PI);
    d(1) = amplitude * sin(time * freq * 2. * M_PI);
}

Vector calcSigma(Vector &dudx, Vector &dudy, double mu, double lambda){
    Vector l(dudx);
    l += dudy;
    l *= lambda;
    
    Vector sx(dudx);
    sx *= 2*mu;
    sx += l;
    sx *= sx;
    
    Vector sy(dudy);
    sy *= 2*mu;
    sy += l;
    sy *= sy;
    
    sx += sy;
    return sx;
}

int main(int argc, char* argv[]) {
    // Defining computation parameters
    #ifdef MFEM_CUDA
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
    GridFunction velocity(&fespace);

    // Defining material parameters
    ConstantCoefficient rhoCoef(4500.0);

    // ConstantCoefficient lambdaCoef(7.7e10);
    // ConstantCoefficient muCoef(4.4e10);

    Vector lambda(mesh.attributes.Max());
    lambda = 7.7e10;
    lambda(1) = lambda(0)*0.01;
    PWConstCoefficient lambdaCoef(lambda);
    Vector mu(mesh.attributes.Max());
    mu = 4.4e10;
    mu(1) = mu(0)*0.01;
    PWConstCoefficient muCoef(mu);
    // PWConstCoefficient lambdaCoef(2);
    // lambdaCoef(0) = 7.7e10;
    // lambdaCoef(1) = 7.7e10;
    // PWConstCoefficient muCoef(2);
    // muCoef(0) = 4.4e10;
    // muCoef(1) = 4.4e10;

    // Defining operator for spatial discretization and ODE solver for time discretization
    ElasticWaveOperator2D oper(fespace, rhoCoef, lambdaCoef, muCoef);
    printf("Operator done!\n");
    oper.SetTime(calculationTime);
    auto ode_solver = new CentralDifferenceSolver();
    ode_solver->Init(oper);

    // Initializing (to zeros) displacement u and velocity du
    Vector u(fespace.GetVSize()), du(fespace.GetVSize());
    u  = 0.0;
    du = 0.0;

    std::cout << "Preparation took " << timer.RealTime() << " seconds.\nStarting forward simulation.\n";

    VectorFunctionCoefficient displacementCoeff(2,displacementOnBoundary);
    displacementCoeff.SetTime(calculationTime);
    Array<int> attr(mesh.bdr_attributes.Size());
    attr = 0;
    attr[0] = 1; // top boundary

    // Displacement derivatives
    FiniteElementSpace scalarFESpace(&mesh, &fe_collection, 1);
    GridFunction duxdx(&scalarFESpace), duydy(&scalarFESpace), sigma(&scalarFESpace);
    GridFunction duydx(&scalarFESpace), duxdy(&scalarFESpace), tau(&scalarFESpace);
    Vector  duxdx_v(scalarFESpace.GetVSize()),
            duydy_v(scalarFESpace.GetVSize()),
            duxdy_v(scalarFESpace.GetVSize()),
            duydx_v(scalarFESpace.GetVSize()),
            sigma_v(scalarFESpace.GetVSize()),
            tau_v  (scalarFESpace.GetVSize());
    // Stresses
    // GridFunction sigma_x, sigma_y, sigma_abs;
    
    // time loop
    const uint saveStep = 1000;
    for (int i = 0; i < NumberOfTimeSteps; i++) {
        // Setting correct values on boundary
        ode_solver->Step(u, du, calculationTime, dt);

        // du *= 1-1e-3; // ONLY for conevrging simualtions

        displacementCoeff.SetTime(calculationTime);
        displacement.SetFromTrueDofs(u);
        displacement.ProjectBdrCoefficient(displacementCoeff, attr);
        displacement.GetTrueDofs(u);
        

        // Saving results
        std::ostringstream filename;
        std::filesystem::create_directory(argv[2]);

        if ( !((i+1) % saveStep) ) {
            filename.str("");
            std::string filename = std::string(argv[2]) + "/plateMesh_" + std::to_string(i+1) + ".vtk";
            std::cout << "\rSaving data to " << filename << std::endl;
            std::ofstream vtkFile(filename);
            mesh.PrintVTK( vtkFile, 4);
            const int ref = 4;
            displacement.SetFromTrueDofs(u);
            displacement.SaveVTK( vtkFile, "displacement", ref);
            velocity.SetFromTrueDofs(du);
            velocity.SaveVTK( vtkFile, "velocity", ref);
            

            // Calculate sigmas
            displacement.GetDerivative(1, 0, duxdx);
            duxdx.SaveVTK(vtkFile, "duxdx", ref);
            displacement.GetDerivative(2, 1, duydy);
            duydy.SaveVTK(vtkFile, "duydy", ref);

            displacement.GetDerivative(1, 1, duxdy);
            duxdy.SaveVTK(vtkFile, "duxdy", ref);
            displacement.GetDerivative(2, 0, duydx);
            duydx.SaveVTK(vtkFile, "duydx", ref);

            duxdx.GetTrueDofs(duxdx_v);
            duydy.GetTrueDofs(duydy_v);
            duxdy.GetTrueDofs(duxdy_v);
            duydx.GetTrueDofs(duydx_v);
            sigma_v = calcSigma(duxdx_v, duydy_v, muCoef.constant, lambdaCoef.constant);
            // tau_v = duxdy_v;
            // tau_v += duydx_v;
            // tau_v *= muCoef.constant;
            // sigma.SetFromTrueDofs(sigma_v);
            // sigma.SaveVTK(vtkFile, "sigma", 4);
            // tau.SetFromTrueDofs(tau_v);
            // tau.SaveVTK(vtkFile, "tau", 4);
            vtkFile.close();

            char vishost[] = "localhost";
            int  visport   = 8080;
            socketstream sol_sock(vishost, visport);
            sol_sock.precision(8);
            sol_sock << "solution\n" << mesh << displacement << std::flush;
        }

        std::cout << "\rFinished time step " << i << std::flush;
    }

    timer.Stop();
    std::cout << "\nProgram took " << timer.RealTime() << " seconds." << std::endl;
}

#pragma once

#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>
#include <string>
#include "ElasticWaveOperator.hpp"
#include "SigmaCoefficient.hpp"

void displacementOnBoundary(const mfem::Vector& x, double time, mfem::Vector& d) {
    const double freq = 3e6;
    const double amplitude = 4e-5; // 4 mm
    d(0) = 0.0;
    d(1) = amplitude * sin(time * freq * 2. * M_PI);
}

class Environment{
    mfem::Mesh mesh;

    mfem::H1_FECollection feCollection;
    mfem::FiniteElementSpace fespace1;
    mfem::FiniteElementSpace fespace2;
    mfem::GridFunction displacementGF;
    mfem::CentralDifferenceSolver odeSolver;
    mfem::VectorFunctionCoefficient displacementCoeff;
    mfem::CentralDifferenceSolver solver;
    mfem::Vector mu, lambda, psi, u, du, sigmaMin, sigmaMax;
    mfem::Array<int> boundaryAttribute;
    mfem::ConstantCoefficient rhoCoeff;
    mfem::PWConstCoefficient muCoeff;
    mfem::PWConstCoefficient lambdaCoeff;

    ElasticWaveOperator2D* oper;
    SigmaCoefficient* sigmaCoeff;

    double t = 0;
    const int REF = 4;
    
public:

    const double E0 = 116e9,
                 nu0 = 0.32,
                 beta = 0.31,
                 sigma_u = 340e6,
                 sigma_b = 1160e6,
                 psi_star = 0.98,
                 gamma = 0.5;

    Environment(std::string meshPath, std::string psiFilename="") :
        mesh(meshPath.c_str()),
        feCollection(4, mesh.Dimension(), mfem::BasisType::GaussLobatto),
        fespace1(&mesh, &feCollection, 1),
        fespace2(&mesh, &feCollection, 2),
        displacementGF(&fespace2),
        displacementCoeff(2, displacementOnBoundary),
        u(fespace2.GetVSize()),
        du(fespace2.GetVSize()),
        sigmaMin(mesh.GetNE()),
        sigmaMax(mesh.GetNE()),
        boundaryAttribute(mesh.bdr_attributes.Size()),
        rhoCoeff(4370.)
    {
        
        // init boundary displacement
        boundaryAttribute = 0;
        boundaryAttribute[0] = 1;
        boundaryAttribute[2] = 1;

        displacementCoeff.SetTime(0);

        // load material
        sigmaMin = 0.;
        sigmaMax = 0.;
        mu.SetSize(mesh.GetNE());
        lambda.SetSize(mesh.GetNE());
        if (psiFilename.size() > 1) {
            std::ifstream psiFile(psiFilename);
            psi.Load(psiFile);
            std::cout << "psi loaded\n";
            MFEM_ASSERT(psi.Size() == mesh.GetNE(), "Loaded psi from file '" + 
                            psiFile + "', it's size is not equal to mesh size");
            for(int i = 0; i < psi.Size(); i++) {
                double E = E0*(1 - 0.5*psi(i))*((psi_star > psi(i)) + 1e-3);
                lambda(i) = E * nu0 / ((1 + nu0) * (1 - 2 * nu0));
                mu(i) = E / (2 * (1 + nu0));
            }
        } else {
            psi.SetSize(mesh.GetNE());

            mu = E0 / (2 * (1 + nu0));
            lambda = E0 * nu0 / ((1 + nu0) * (1 - 2 * nu0));
            psi = 0.;
        }
        
        muCoeff = mfem::PWConstCoefficient(mu);
        lambdaCoeff = mfem::PWConstCoefficient(lambda);

        

        // init solver
        oper = new ElasticWaveOperator2D(fespace2, rhoCoeff, lambdaCoeff, muCoeff);
        oper->SetTime(0);
        solver.Init(*oper);
        sigmaCoeff = new SigmaCoefficient(displacementGF, lambdaCoeff, muCoeff);

        std::cout << "Init done\n";
    }

    void step(double dt) {
        solver.Step(u, du, t, dt);

        displacementCoeff.SetTime(t);
        displacementGF.SetFromTrueDofs(u);
        displacementGF.ProjectBdrCoefficient(displacementCoeff, boundaryAttribute);
        displacementGF.GetTrueDofs(u);

        for(int i = 0; i < mesh.GetNE(); i++) {
            mfem::Vector center(2);
            mesh.GetElementCenter(i, center);
            mfem::IntegrationPoint ip;
            ip.Set2(center.GetData());
            double s = sigmaCoeff->Eval(*mesh.GetElementTransformation(i), ip);
            sigmaMin(i) = std::min(s, sigmaMin(i));
            sigmaMax(i) = std::max(s, sigmaMax(i));
        }
    }

    void saveDisplacement(std::string filename) {
        std::cout << "\rSaving data to " << filename << std::endl;

        std::ofstream vtkFile(filename);
        mesh.PrintVTK(vtkFile, REF);
        displacementGF.SaveVTK(vtkFile, "displacement", REF);
    }

    void savePsi(std::string filename) {
        double deltaN = 1e18;
        psi = 0.;
        mfem::Vector elmentCenter(2);

        for (int i = 0; i < mesh.GetNE(); i++) {
            double sx = sigmaMax(i);
            if (sx <= 0)
                continue;
            double sn = sigmaMin(i);
            double deltaSigma = sx - sn;
            double sigma_eq = std::sqrt((sx*deltaSigma)/2);
            if (sigma_eq <= sigma_u)
                continue;

            mesh.GetElementCenter(i, elmentCenter);
            if (elmentCenter.Norml2() > 5e-3)
                continue;

            double sigmaFrac = (sigma_eq - sigma_u)/(sigma_b - sigma_u);
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

            mesh.GetElementCenter(i, elmentCenter);
            if (elmentCenter.Norml2() > 5e-3)
                continue;

            double sigmaFrac = (sigma_eq - sigma_u)/(sigma_b - sigma_u);
            double B = 1e-3 * (std::pow(sigmaFrac, 1/beta)) / (2*(1-gamma));
            psi(i) = std::pow(1 - std::sqrt(p2(1-std::pow(psi(i), 1-gamma)) - 2*(1-gamma)*deltaN*B), 1/(1-gamma));
        }

        std::ofstream psiFile(filename);
        psiFile << psi.Size() << std::endl;
        psi.Print(psiFile, 1);
        psiFile.close();
    }

};



#pragma once

#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>
#include <string>
#include <stdexcept>
#include "ElasticWaveOperator.hpp"
#include "SigmaCoefficient.hpp"

double amplitude = 0.8e-6;
const double freq = 2e6;

double gentle_sin(double t){
    return (exp(3*t)-1) / exp(3*t) * sin(t);
}

void displacementOnTop(const mfem::Vector& x, double time, mfem::Vector& d) {
    d(0) = 0;
    d(1) = amplitude * gentle_sin(time * freq * 2. * M_PI);
}

void displacementOnBottom(const mfem::Vector& x, double time, mfem::Vector& d) {
    d(0) = 0;
    d(1) = -amplitude * gentle_sin(time * freq * 2. * M_PI);
}

class Environment{
    mfem::Mesh mesh;

    mfem::H1_FECollection feCollection;
    mfem::FiniteElementSpace fespace1;
    mfem::FiniteElementSpace fespace2;
    mfem::GridFunction displacementGF;
    mfem::GridFunction sigmaGF;
    mfem::VectorFunctionCoefficient displTopCoeff, displBottomCoeff;
    mfem::CentralDifferenceSolver solver;
    mfem::Vector mu, lambda, psi0, psi, u, du, sigmaMin, sigmaMax, sigma_eq, B, dn;
    mfem::Array<int> bdrAtrTop, bdrAtrBottom;
    mfem::ConstantCoefficient rhoCoeff;
    mfem::PWConstCoefficient muCoeff;
    mfem::PWConstCoefficient lambdaCoeff;

    ElasticWaveOperator2D* oper=0;
    SigmaCoefficient* sigmaCoeff=0;

    double t = 0;
    
public:
    const double E0 = 116e9,
                 nu0 = 0.32,
                 beta = 0.31,
                 sigma_u = 3.7e8,
                 sigma_b = 1.16e9,
                 psi_star = 0.98,
                 gamma = 0.5,
                 rho_0 = 4500,
                 dPsi0 = 0.1,
                 kappa = 1.0,
                 resid = 1e-3;

    Environment(std::string meshPath, std::string psiFilename="", double t=0) :
        mesh(meshPath.c_str()),
        feCollection(2, mesh.Dimension(), mfem::BasisType::GaussLobatto),
        fespace1(&mesh, &feCollection, 1),
        fespace2(&mesh, &feCollection, 2),
        displacementGF(&fespace2),
        sigmaGF(&fespace1),
        displTopCoeff(2, displacementOnTop),
        displBottomCoeff(2, displacementOnBottom),
        u(fespace2.GetVSize()),
        du(fespace2.GetVSize()),
        sigmaMin(mesh.GetNE()),
        sigmaMax(mesh.GetNE()),
        sigma_eq(mesh.GetNE()),
        dn(mesh.GetNE()),
        B(mesh.GetNE()),
        bdrAtrTop(mesh.bdr_attributes.Size()),
        bdrAtrBottom(mesh.bdr_attributes.Size())
    {
        bdrAtrTop = 0;
        bdrAtrTop[0] = 1;
        displTopCoeff.SetTime(0);
        bdrAtrBottom = 0;
        bdrAtrBottom[2] = 1;
        displBottomCoeff.SetTime(0);

        sigmaMin = 0.;
        sigmaMax = 0.;
        sigma_eq = 0.;
        B = 0.;
        // u = 0.;
        // du = 0.;
        mu.SetSize(mesh.GetNE());
        lambda.SetSize(mesh.GetNE());
        if (psiFilename.size() > 1) {
            loadPsi(psiFilename);
        } else {
            psi0.SetSize(mesh.GetNE());

            mu = E0 / (2 * (1 + nu0));
            lambda = E0 * nu0 / ((1 + nu0) * (1 - 2 * nu0));
            psi0 = 0.;
        }
        psi = mfem::Vector(psi0);
        
        initMaterial(t);

        std::cout << "Init done\n";
    }

    void initMaterial(double t){
        muCoeff = mfem::PWConstCoefficient(mu);
        lambdaCoeff = mfem::PWConstCoefficient(lambda);
        rhoCoeff = mfem::ConstantCoefficient(rho_0);

        if (oper) {
            delete oper;
            oper = 0;
        }

        oper = new ElasticWaveOperator2D(fespace2, rhoCoeff, lambdaCoeff, muCoeff);
        oper->SetTime(t);
        solver.Init(*oper);
        if (sigmaCoeff) {
            delete sigmaCoeff;
            sigmaCoeff = 0;
        }
        sigmaCoeff = new SigmaCoefficient(displacementGF, lambdaCoeff, muCoeff);
    }

    void saveEnv(std::string envPath) {
        std::ofstream envFile(envPath);
        if (envFile.fail())
            throw std::invalid_argument("Could not write to " + envPath);
        envFile << u.Size() << '\n';
        u.Print(envFile);
        du.Print(envFile);
        envFile.close();
    }

    void loadEnv(std::string envPath) {
        std::ifstream envFile(envPath);
        if (envFile.fail())
            throw std::invalid_argument("Could not read from " + envPath);
        u.Load(envFile);
        du.Load(envFile, u.Size());
        sigmaMin = 0.0;
        sigmaMax = 0.0;
        B = 0.0;
        dn = 0.0;
    }

    void step(double dt) {
        solver.Step(u, du, t, dt);

        // u *= 1-2e-5; // energy dissipation
        
        displTopCoeff.SetTime(t);
        displBottomCoeff.SetTime(t);
        displacementGF.SetFromTrueDofs(u);
        // if (t >= dt*4000) amplitude *= decay;
        // if (t <= dt*9000) {
            displacementGF.ProjectBdrCoefficient(displTopCoeff, bdrAtrTop);
            displacementGF.ProjectBdrCoefficient(displBottomCoeff, bdrAtrBottom);
        // }
        displacementGF.GetTrueDofs(u);
        sigmaGF.ProjectCoefficient(*sigmaCoeff);
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


    void saveDisplacement(std::ostream& out) {
        mesh.PrintVTK(out, 1);
        displacementGF.SaveVTK(out, "displacement", 1);
    }

    void saveSigma(std::ostream& out) {
        sigmaGF.ProjectCoefficient(*sigmaCoeff);
        sigmaGF.SaveVTK(out, "sigma", 1);
    }

    void saveVtk(std::ostream& out, mfem::Vector& vec, std::string name) {
        out << "SCALARS " << name << " double 1\n";
        out << "LOOKUP_TABLE default\n";
        vec.Print(out, 1);
    }

    void saveVars(std::ostream& out) {
        out << "CELL_DATA " << mesh.GetNE() << '\n';
        saveVtk(out, sigmaMin, "sigma_min");
        saveVtk(out, sigmaMax, "sigma_max");
        saveVtk(out, sigma_eq, "sigma_eq");
        saveVtk(out, B, "B");
        saveVtk(out, dn, "deltaN");
        saveVtk(out, psi, "psi");
    }

    void calcPsi() {
        mfem::Vector elmentCenter(2);
        double deltaN = 1e18;
        B = 0.;
        psi = mfem::Vector(psi0); // copy
        dn = 0.;

        // for(int j = 0; j < 4; j++) { // speed up sim
            for (int i = 0; i < mesh.GetNE(); i++) {
                double sx = sigmaMax(i);
                if (sx <= 0)
                    continue;
                double sn = sigmaMin(i);
                double deltaSigma = sx - sn;
                sigma_eq(i) = std::sqrt((sx*deltaSigma)/2);
                if (sigma_eq(i) <= sigma_u)
                    continue;

                double sigmaFrac = (sigma_eq(i) - sigma_u)/(sigma_b - sigma_u);
                B(i) = 1e-3 * std::pow(sigmaFrac, 1/beta) / (2*(1-gamma));
                
                if (psi(i) < 1)
                    dn(i) = 0.5/((1-gamma)*B(i)) * (p2(1-std::pow(psi(i), 1-gamma)) - p2(1-std::pow(std::min(1.0, psi(i) + dPsi0), 1-gamma)));
                    // dn(i) = 0.5/((1-gamma)*B(i)) * (p2(1-std::pow(psi(i), 1-gamma)) - p2(1-std::pow(std::min(psi_star, 1 - (1- psi(i))*0.5), 1-gamma)));
                else 
                    dn(i) = 1e18;
                
                // if(psi(i) < 1)
                //     dn(i) = 0.5/((1-gamma)*B(i)) * p2(1-std::pow(psi(i), 1-gamma));
                // else dn(i) = 1e18;
                deltaN = std::min(deltaN, dn(i));
            }
            for (int i = 0; i < mesh.GetNE(); i++) {
                if (B(i) == 0)
                    continue;
                psi(i) = std::pow(1 - std::sqrt(p2(1-std::pow(psi(i), 1-gamma)) - 2*(1-gamma)*deltaN*B(i)), 1/(1-gamma));
                psi(i) = std::min(1.0, psi(i));
                dn(i) = std::min(dn(i), deltaN*10);
            }
        // }
    }

    void savePsi(std::string filename) {
        double deltaN = 1e18;
        std::ofstream psiFile(filename);
        psiFile << psi.Size() << std::endl;
        psi.Print(psiFile, 1);
        psiFile.close();
    }

    void loadPsi(std::string filename, double time=0) {
        std::ifstream psiFile(filename);
        if (psiFile.fail())
            throw std::invalid_argument("Psi file does not exist");
        psi0.Load(psiFile);
        std::cout << "Psi loaded from " << filename << std::endl;
        if (psi0.Size() != mesh.GetNE())
            throw std::invalid_argument("Loaded psi from file '" + 
                        filename + "', it's size (" + 
                        std::to_string(psi0.Size()) + 
                        ") is not equal to mesh size (" + 
                        std::to_string(mesh.GetNE()) + ")");
        for(int i = 0; i < psi0.Size(); i++) {
            double E = E0*(1 - kappa*psi0(i))*((psi_star > psi0(i)) + 1e-3);
            lambda(i) = E * nu0 / ((1 + nu0) * (1 - 2 * nu0));
            mu(i) = E / (2 * (1 + nu0));
        }

        initMaterial(t);
        t = time;
    }

    float getTime(){
        return oper->GetTime();
    }

    void setTime(double t){
        oper->SetTime(t);
    }

};

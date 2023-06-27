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

double amplitude = 0.75e-6;
const double freq = 5e6;

double gentle_sin(double t){
    return (exp(3*t)-1) / exp(3*t) * sin(t);
}

double single_wave(double t){
    return t < 2*M_PI ? (1-cos(t))*0.5 : 0;
}

void displacementOnTop(const mfem::Vector& x, double time, mfem::Vector& d) {
    d(0) = 0;
    d(1) = amplitude * single_wave(time * freq * 2. * M_PI);
}

void displacementOnBottom(const mfem::Vector& x, double time, mfem::Vector& d) {
    d(0) = 0;
    d(1) = -amplitude * single_wave(time * freq * 2. * M_PI);
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
                 beta_l = 0.31,
                 beta_v = 0.27,
                 sigma_u = 3.37e8,
                 sigma_uw = 2.5e8,      
                 sigma_b = 1.16e9,
                 psi_star = 0.98,
                 gamma = 0.5,
                 alpha = 1-gamma,
                 rho_0 = 4500,
                 dPsi0 = 0.1,
                 kappa = 1.0,
                 resid = 1e-3;

    Environment(std::string meshPath, std::ostream &runInfo, std::string psiFilename="", double t=0) :
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
        runInfo << "E0=" << E0 << std::endl;
        runInfo << "nu0=" << nu0 << std::endl; 
        runInfo << "beta_l=" << beta_l << std::endl; 
        runInfo << "beta_v=" << beta_v << std::endl; 
        runInfo << "sigma_u=" << sigma_u << std::endl; 
        runInfo << "sigma_uw=" << sigma_uw << std::endl; 
        runInfo << "sigma_b=" << sigma_b << std::endl; 
        runInfo << "psi_star=" << psi_star << std::endl; 
        runInfo << "gamma=" << gamma << std::endl; 
        runInfo << "rho_0=" << rho_0 << std::endl; 
        runInfo << "dPsi0=" << dPsi0 << std::endl; 
        runInfo << "kappa=" << kappa << std::endl; 
        runInfo << "resid=" << resid << std::endl; 
        runInfo << "amplitude=" << amplitude << std::endl; 
        runInfo << "freq=" << freq << std::endl; 

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

    void step(double dt, bool projectDispl=true) {
        solver.Step(u, du, t, dt);

        // u *= 1-2e-5; // energy dissipation
        
        displTopCoeff.SetTime(t);
        displBottomCoeff.SetTime(t);
        displacementGF.SetFromTrueDofs(u);
        if (projectDispl) {
            displacementGF.ProjectBdrCoefficient(displTopCoeff, bdrAtrTop);
            displacementGF.ProjectBdrCoefficient(displBottomCoeff, bdrAtrBottom);
        }
        displacementGF.GetTrueDofs(u);
        sigmaGF.ProjectCoefficient(*sigmaCoeff);
        for(int i = 0; i < mesh.GetNE(); i++) {
            mfem::Array<int> vert(4);
            mesh.GetElementVertices(i, vert);

            for (int j = 0; j < 4; j++) {
                for (int k = j; k < 4; k++) {
                    mfem::IntegrationPoint ip1, ip2;
                    ip1.Set2(mesh.GetVertex(vert[j]));
                    ip2.Set2(mesh.GetVertex(vert[k]));
                    ip1.x = (ip1.x + ip2.x) / 2;
                    ip1.y = (ip1.y + ip2.y) / 2;
                    mesh.GetElementTransformation(i)->SetIntPoint(&ip1);
                    double s = sigmaCoeff->Eval(*mesh.GetElementTransformation(i), ip1);
                    sigmaMin(i) = std::min(s, sigmaMin(i));
                    sigmaMax(i) = std::max(s, sigmaMax(i));
                }
            }

            // mfem::Vector center(2);
            // mesh.GetElementCenter(i, center);
            // mfem::IntegrationPoint ip;
            // ip.Set2(center.GetData());
            // double s = sigmaCoeff->Eval(*mesh.GetElementTransformation(i), ip);
            // sigmaMin(i) = std::min(s, sigmaMin(i));
            // sigmaMax(i) = std::max(s, sigmaMax(i));
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

        for (int i = 0; i < mesh.GetNE(); i++) {
            double sx = sigmaMax(i);
            if (sx <= 0)
                continue;
            double sn = sigmaMin(i);
            double deltaSigma = sx - sn;
            sigma_eq(i) = std::sqrt((sx*deltaSigma)/2);

            double sigmaFrac_l = std::max(sigma_eq(i) - sigma_u, 0.0)/(sigma_b - sigma_u);
            double sigmaFrac_v = std::max(sigma_eq(i) - sigma_uw, 0.0)/(sigma_u - sigma_uw);    
            B(i) = std::max(1e-3 * std::pow(sigmaFrac_l, 1/beta_l) / (2*alpha),
                            1e-8 * std::pow(sigmaFrac_v, 1/beta_v) / (2*alpha));
            
            if (B(i) == 0 || psi(i) == 1)
                continue;

            dn(i) = 0.5/(alpha*B(i)) * (p2(1-std::pow(psi(i), alpha)) - p2(1-std::pow(std::min(1.0, psi(i) + dPsi0), alpha)));
            // dn(i) = 0.5/(alpha*B(i)) * (p2(1-std::pow(psi(i), alpha)) - p2(1-std::pow(std::min(psi_star, 1 - (1- psi(i))*0.5), alpha)));
            // dn(i) = 0.25/(alpha*B(i)) * p2(1-std::pow(psi(i), alpha));
            
            // else dn(i) = 1e18;
            deltaN = std::min(deltaN, dn(i));
        }

        for (int i = 0; i < mesh.GetNE(); i++) {
            if (B(i) == 0)
                continue;
            psi(i) = std::pow(1 - std::sqrt(p2(1-std::pow(psi(i), alpha)) - 2*alpha*deltaN*B(i)), 1/alpha);
            psi(i) = std::min(1.0, psi(i));
            dn(i) = std::min(dn(i), deltaN*10);
        }
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

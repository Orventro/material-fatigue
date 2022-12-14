#pragma once

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include "mfem.hpp"

#define p2(a) ((a)*(a))

class FatigueCoefficient
{
private:
    mfem::GridFunction &sigma, psi;
    const double beta = 0.31, sigma_u = 340e6, sigma_B = 1160e6, psi_star = 0.98, gamma = 0.5, k = 0.5, E0 = 116e9, nu0 = 0.32, rho0 = 4370;

public:
    FatigueCoefficient(mfem::GridFunction &_sigma,
                       mfem::FiniteElementSpace &space)
            : sigma(_sigma), psi(&space) { 
        psi = 0.;
        MFEM_ASSERT(psi.Size()*3 == sigma.Size(), "invalid input"); 
        std::cout << "FatigueCoefficient: " << psi.Size() << std::endl;
    }

    virtual void Update() {
        double deltaN = 1e18; // just a big number

        for(int i = 0; i < psi.Size(); i++) {
            double deltaSigma = sigma(i*3+2) - sigma(i*3+1);
            double sigma_eq = std::sqrt((sigma(i*3+2)*deltaSigma)/2);
            if (sigma_eq > sigma_u) {
                double sigmaFrac = (sigma_eq - sigma_u)/(sigma_B - sigma_u);
                double B = 1e-3 * (std::pow(sigmaFrac, 1/beta)) / (2*(1-gamma));
                
                double newN = 0.5/B * (0.5/(1-gamma) - std::pow(psi(i), 1-gamma)/(1-gamma) + std::pow(psi(i), 2-2*gamma)/(2-2*gamma));
                // printf(" gamma %f\n pow1 %f\n pow2 %f\n psi %f\n B %f\n powB %f\n sigmaFrac %f\n sEq %f\n sU %f\n sB %f\n N1 %f\n", 
                //     gamma, std::pow(psi(i), 1-gamma), std::pow(psi(i), 2-2*gamma), psi(i), B, std::pow(sigmaFrac, 1/beta), sigmaFrac,
                //     sigma_eq, sigma_u, sigma_B, 1/B);
                deltaN = std::min(deltaN, newN);
            }
        }

        // std::cout << "deltaN = " << deltaN << std::endl;

        for(int i = 0; i < psi.Size(); i++) {
            double deltaSigma = sigma(i*3+2) - sigma(i*3+1);
            double sigma_eq = std::sqrt((sigma(i*3+2)*deltaSigma)/2);
            if (sigma_eq > sigma_u) {
                double sigmaFrac = (sigma_eq - sigma_u)/(sigma_B - sigma_u);
                double B = 1e-3 * (std::pow(sigmaFrac, 1/beta)) / (2*(1-gamma));

                psi(i) = std::pow(1 - std::sqrt(p2(1-std::pow(psi(i), 1-gamma)) - 2*(1-gamma)*deltaN*B), 1/(1-gamma));
                // printf("psi %f\n sqrt2 %f\n", psi(i), p2(1-std::pow(psi(i), 1-gamma)) - 2*(1-gamma)*deltaN*B);
            }
        }
    }

    mfem::GridFunction &getPsi() { return psi;}

    friend class LambdaCoefficient;
    friend class MuCoefficient;
};


class LambdaCoefficient : public mfem::Coefficient
{
private:
    FatigueCoefficient &fatigue;
public:

    LambdaCoefficient(FatigueCoefficient &c) : fatigue(c) { }

    virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
    {   
        double p = fatigue.psi.GetValue(T, ip);
        double E = fatigue.E0 * (1 - fatigue.k*p) * (((fatigue.psi_star > p) ? 1:0) + 1e-3);
        return E*fatigue.nu0/((1+fatigue.nu0)*(1-2*fatigue.nu0));
    }
};


class MuCoefficient : public mfem::Coefficient
{
private:
    FatigueCoefficient &fatigue;
public:

    MuCoefficient(FatigueCoefficient &c) : fatigue(c) { }

    virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
    {   
        double p = fatigue.psi.GetValue(T, ip);
        double E = fatigue.E0 * (1 - fatigue.k*p) * (((fatigue.psi_star > p) ? 1:0) + 1e-3);
        return E / (2 * (1 + fatigue.nu0));
    }
};

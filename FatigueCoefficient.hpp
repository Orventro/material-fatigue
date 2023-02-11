#pragma once

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include "mfem.hpp"

#define p2(a) ((a)*(a))

class FatigueCoefficient : public mfem::Coefficient
{
private:
    mfem::GridFunction &sigma, psi, fatigue;
    const double beta = 0.31, sigma_u = 340e6, sigma_B = 1160e6, psi_star = 0.98, gamma = 0.5, k = 0.5, E0 = 116e9, nu0 = 0.32, rho0 = 4370;
    double deltaN = 1e18;
    bool calcMode = true;

public:
    FatigueCoefficient(mfem::GridFunction &_sigma,
                       mfem::FiniteElementSpace &space)
            : sigma(_sigma), psi(&space), fatigue(&space)
    { 
        psi = 0.;
        MFEM_ASSERT(psi.Size()*3 == sigma.Size(), "invalid input"); 
    }

    // virtual void Update() 
    // {

    //     for(int i = 0; i < psi.Size(); i++) {
    //         double deltaSigma = sigma(i*3+2) - sigma(i*3+1);
    //         double sigma_eq = std::sqrt((sigma(i*3+2)*deltaSigma)/2);
    //         double sigmaFrac = (sigma_eq - sigma_u)/(sigma_B - sigma_u);
    //         double B = 1e-3 * (std::pow(sigmaFrac, 1/beta)) / (2*(1-gamma));
            
    //         deltaN = std::min(deltaN, 0.5/B * (0.5/(1-gamma) - std::pow(psi(i), 1-gamma)/(1-gamma) + std::pow(psi(i), 2-2*gamma)/(2-2*gamma)));
    //     }


    //     for(int i = 0; i < psi.Size(); i++) {
    //         double deltaSigma = sigma(i*3+2) - sigma(i*3+1);
    //         double sigma_eq = std::sqrt((sigma(i*3+2)*deltaSigma)/2);
    //         double sigmaFrac = (sigma_eq - sigma_u)/(sigma_B - sigma_u);
    //         double B = 1e-3 * (std::pow(sigmaFrac, 1/beta)) / (2*(1-gamma));

    //         psi(i) = std::pow(1 - std::sqrt(p2(1-std::pow(psi(i), 1-gamma)) - 2*(1-gamma)*deltaN*B), 1/(1-gamma));
    //     }
    // }

    virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
    {
        if(calcMode) {
            double s2 = sigma.GetValue(T, ip, 2);
            double s3 = sigma.GetValue(T, ip, 3);

            double deltaSigma = s3 - s2;
            double sigma_eq = std::sqrt((s3*deltaSigma)/2);
            if (sigma_eq <= sigma_u)
                return deltaN;
            double sigmaFrac = (sigma_eq - sigma_u)/(sigma_B - sigma_u);
            double B = 1e-3 * (std::pow(sigmaFrac, 1/beta)) / (2*(1-gamma));

            double x = 0.5/B * (0.5/(1-gamma) - std::pow(psi.GetValue(T, ip), 1-gamma)/(1-gamma) + std::pow(psi.GetValue(T, ip), 2-2*gamma)/(2-2*gamma));

            // if (x != x) {
            //     double ps = psi.GetValue(T, ip);
            //     std::cout << "x " << x << ' ' << B << ' ' << ps << ' ' << sigmaFrac << ' ' << (0.5/B * 0.5/(1-gamma)) << std::endl;
            // }

            deltaN = std::min(x, deltaN);//std::min(deltaN, 1/2/B * (0.5/(1-gamma) - std::pow(psi.GetValue(T, ip), 1-gamma)/(1-gamma) + std::pow(psi.GetValue(T, ip), 2-2*gamma)/(2-2*gamma)));
            return deltaN;
        } else {
            double s2 = sigma.GetValue(T, ip, 2);
            double s3 = sigma.GetValue(T, ip, 3);
            
            double deltaSigma = s3 - s2;
            double sigma_eq = std::sqrt((s3*deltaSigma)/2);
            if (sigma_eq <= sigma_u)
                return psi.GetValue(T, ip);
            double sigmaFrac = (sigma_eq - sigma_u)/(sigma_B - sigma_u);
            double B = 1e-3 * (std::pow(sigmaFrac, 1/beta)) / (2*(1-gamma));

            if(p2(1-std::pow(psi.GetValue(T, ip), 1-gamma)) < 2*(1-gamma)*deltaN*B) 
                return psi.GetValue(T, ip);

            double ps = std::pow(1 - std::sqrt(p2(1-std::pow(psi.GetValue(T, ip), 1-gamma)) - 2*(1-gamma)*deltaN*B), 1/(1-gamma));

            // std::cout << "ps " << psi.GetValue(T, ip) << ' ' << std::pow(psi.GetValue(T, ip), 1-gamma) << ' ' << 2*(1-gamma)*deltaN*B << ' ' << (p2(1-std::pow(psi.GetValue(T, ip), 1-gamma)) - 2*(1-gamma)*deltaN*B) << std::endl;

            return ps;
        }
    }

    void calcDelta() {
        deltaN = 1e18;
        calcMode = true;
        fatigue.ProjectCoefficient(*this);
        // std::cout << "deltan " << deltaN << std::endl;
    }

    void calcPsi() {
        calcMode = false;
        psi.ProjectCoefficient(*this);
    }

    mfem::GridFunction& getPsi() {
        return psi;
    }

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
        double E = fatigue.E0*(1 - fatigue.k*p)*(((fatigue.psi_star > p) ? 1:0) + 1e-3);
        // return E*fatigue.nu0/((1+fatigue.nu0)*(1-2*fatigue.nu0));
        double l = E*fatigue.nu0/((1+fatigue.nu0)*(1-2*fatigue.nu0));
        // if (l != l) {
        //     std::cout << p << std::endl;
        // }
        return l;
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
        double E = fatigue.E0*(1 - fatigue.k*p)*(((fatigue.psi_star > p) ? 1:0) + 1e-3);
        return E / (2 * (1 + fatigue.nu0));
    }
};

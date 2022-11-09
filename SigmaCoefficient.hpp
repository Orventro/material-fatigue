#pragma once

#include <fstream>
#include <stdexcept>
#include "mfem.hpp"

class SigmaCoefficient : public mfem::Coefficient
{
private:
    mfem::GridFunction &u;
    mfem::Coefficient &lambda, &mu;
    mfem::DenseMatrix eps, sigma;

public:
    SigmaCoefficient(mfem::GridFunction &_u, mfem::Coefficient &_lambda, mfem::Coefficient &_mu)
            : u(_u), lambda(_lambda), mu(_mu) { component = -1; }

    virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
    {
        u.GetVectorGradient(T, eps);  // eps = grad(u)
        eps.Symmetrize();             // eps = (1/2)*(grad(u) + grad(u)^t)
        double l = lambda.Eval(T, ip);
        double m = mu.Eval(T, ip);
        sigma.Diag(l*eps.Trace(), eps.Size()); // sigma = lambda*trace(eps)*I
        sigma.Add(2*m, eps);          // sigma += 2*mu*eps

        switch (component)
        {
            case 0:
                return sigma(0, 0); // Sxx
            case 1:
                return sigma(0, 1); // Sxy
            case 2:
                return sigma(1, 1); // Syy
            default:
                throw std::runtime_error("Only 0(Sxx), 1(Sxy) and 2(Syy) are supported.");
        }
    }
    
    virtual void Read(std::istream &in) { }
    virtual ~SigmaCoefficient() { }
    void setComponent (char _component) { component = _component; }

private:
    char component;
};

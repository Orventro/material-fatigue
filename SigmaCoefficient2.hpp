#pragma once

#include <fstream>
#include <stdexcept>
#include "mfem.hpp"

class SigmaCoefficient : public mfem::VectorCoefficient
{
private:
    mfem::GridFunction &u;
    mfem::Coefficient &lambda, &mu;
    mfem::DenseMatrix eps, sigma;

public:
    SigmaCoefficient(mfem::GridFunction &_u, mfem::Coefficient &_lambda, mfem::Coefficient &_mu)
            : VectorCoefficient(3), u(_u), lambda(_lambda), mu(_mu) {  }

    virtual void Eval(mfem::Vector &V, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
    {
        u.GetVectorGradient(T, eps);  // eps = grad(u)
        eps.Symmetrize();             // eps = (1/2)*(grad(u) + grad(u)^t)
        double l = lambda.Eval(T, ip);
        double m = mu.Eval(T, ip);
        sigma.Diag(l*eps.Trace(), eps.Size()); // sigma = lambda*trace(eps)*I
        sigma.Add(2*m, eps);          // sigma += 2*mu*eps

        V(0) = sigma(0,0);
        V(1) = sigma(0,1);
        V(2) = sigma(1,1);
    }
    
    virtual void Read(std::istream &in) { }
    virtual ~SigmaCoefficient() { }

private:
};

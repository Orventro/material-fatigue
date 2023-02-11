#pragma once

#include <fstream>
#include <stdexcept>
#include "mfem.hpp"

#define p2(a) ((a)*(a))

class SigmaCoefficient : public mfem::Coefficient
{
private:
    mfem::GridFunction &u;
    mfem::Coefficient &lambda, &mu;
    mfem::DenseMatrix eps, sigma;
    bool firstEval = 1;

public:
    SigmaCoefficient(mfem::GridFunction &_u, mfem::Coefficient &_lambda, mfem::Coefficient &_mu)
            :u(_u), lambda(_lambda), mu(_mu) {  }

    virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
    {
        u.GetVectorGradient(T, eps);  // eps = grad(u)
        eps.Symmetrize();             // eps = (1/2)*(grad(u) + grad(u)^t)
        double l = lambda.Eval(T, ip);
        double m = mu.Eval(T, ip);
        sigma.Diag(l*eps.Trace(), eps.Size()); // sigma = lambda*trace(eps)*I
        sigma.Add(2*m, eps);          // sigma += 2*mu*eps

        return (sigma(0,0) + sigma(1,1))/2 + 0.5*std::sqrt(p2(sigma(0,0) - sigma(1,1)) + 4*p2(sigma(0,1)));
    }
    
    void hadEval() {firstEval = 0;}

    virtual void Read(std::istream &in) { }
    virtual ~SigmaCoefficient() { }
};

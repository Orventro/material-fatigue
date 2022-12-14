#pragma once

#include <fstream>
#include <stdexcept>
#include "mfem.hpp"

#define p2(a) ((a)*(a))

class SigmaCoefficient : public mfem::VectorCoefficient
{
private:
    mfem::GridFunction &u;
    mfem::Coefficient *lambda, *mu;
    mfem::DenseMatrix eps, sigma;
    bool firstEval = 1;

public:
    SigmaCoefficient(mfem::GridFunction &_u, mfem::Coefficient *_lambda=0, mfem::Coefficient *_mu=0)
            :u(_u), lambda(_lambda), mu(_mu), mfem::VectorCoefficient(3) {  }

    void setMaterial(mfem::Coefficient *_lambda, mfem::Coefficient *_mu) {
        lambda = _lambda;
        mu = _mu;
    }

    virtual void Eval(mfem::Vector &V, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
    {
        MFEM_ASSERT(lambda && mu, "material not set");

        u.GetVectorGradient(T, eps);  // eps = grad(u)
        eps.Symmetrize();             // eps = (1/2)*(grad(u) + grad(u)^t)
        double l = lambda->Eval(T, ip);
        double m = mu->Eval(T, ip);
        sigma.Diag(l*eps.Trace(), eps.Size()); // sigma = lambda*trace(eps)*I
        sigma.Add(2*m, eps);          // sigma += 2*mu*eps

        if(firstEval) {
            V(1) = 1e18;
            V(2) = -1e18;
        }

        V(0) = (sigma(0,0) + sigma(1,1))/2 + 1/2*std::sqrt(p2(sigma(0,0) - sigma(1,1)) + 4*p2(sigma(0,1)));
        V(1) = std::min(V(1), V(0));
        V(2) = std::max(V(2), V(0));
    }
    
    void hadEval() {firstEval = 0;}

    virtual void Read(std::istream &in) { }
    virtual ~SigmaCoefficient() { }
};

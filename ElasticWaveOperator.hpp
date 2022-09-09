#pragma once

#include "mfem.hpp"

class ElasticWaveOperator2D : public mfem::SecondOrderTimeDependentOperator {
protected:
    mfem::FiniteElementSpace fespace;
    mfem::Array<int> listOfEssentialDOFs;
    mfem::BilinearForm M, K;
    mfem::LinearForm   * f;
    mfem::SparseMatrix * KMat;
    mfem::SparseMatrix * MMatInv;

public:
    ElasticWaveOperator2D(mfem::FiniteElementSpace &fspace, mfem::Coefficient &rhoCoef, mfem::Coefficient &lambdaCoef,
                          mfem::Coefficient &muCoef);
    
    void Mult(const mfem::Vector& u, const mfem::Vector &du, mfem::Vector &d2u) const override;
    void ImplicitSolve(const double beta, const double, const mfem::Vector& u,const mfem::Vector &du, mfem::Vector &d2u) override;
    // void 
    void SetTime(const double t) override;
    ~ElasticWaveOperator2D() override;
};

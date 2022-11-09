#pragma once

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include "mfem.hpp"

class FatigueCalculator {
protected:

    mfem::GridFunction displacement,
                       duxdxGf,
                       duxdyGf,
                       duydxGf,
                       duydyGf,
                       sigmaGf;

    mfem::Vector u,
                 duxdx,
                 duydx,
                 duxdy,
                 duydy,
                 sigmaxx,
                 sigmaxy,
                 sigmayy,
                 mu,
                 lambda,
                 elementCenters,
                 maxSigma;

    SigmaCoefficient sigmaCoeff;

    // mfem::Array<double> maxSigma;
    mfem::Array<int> fatigue;
    mfem::Array<bool> stressed;

    bool duxdxActual,
         duxdyActual,
         duydxActual,
         duydyActual,
         sigmaActual,
         fatigueActual;

public:

    FatigueCalculator(mfem::FiniteElementSpace &fspace, 
                      mfem::FiniteElementSpace &scalarSpace, 
                      mfem::Mesh &mesh, mfem::Vector &mu, mfem::Vector &lambda);
    void update(mfem::Vector &newU);
    void calcDuxdx();
    void calcDuxdy();
    void calcDuydx();
    void calcDuydy();
    void calcSigma();
    void updateFatigue();
    void printSigma(std::ostream &out);
    void printDuxdx(std::ostream &out);
    void printDuxdy(std::ostream &out);
    void printDuydx(std::ostream &out);
    void printDuydy(std::ostream &out);
    void printFatigue(std::ostream &out);
    void printMaxSigma(std::ostream &out);
};


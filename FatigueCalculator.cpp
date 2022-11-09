#include "FatigueCalculator.hpp"

FatigueCalculator::FatigueCalculator(mfem::FiniteElementSpace &fspace,
                                     mfem::FiniteElementSpace &scalarSpace,
                                     mfem::Mesh &mesh,
                                     mfem::GridFunction &u,
                                     mfem::Coefficient &mu,
                                     mfem::Coefficient &lambda) : 
        displacement(&fspace),
        duxdxGf(&scalarSpace), 
        duxdyGf(&scalarSpace), 
        duydxGf(&scalarSpace), 
        duydyGf(&scalarSpace),
        sigmaGf(&scalarSpace),
        u(fspace.GetNDofs()*2),
        duxdx(mesh.GetNE()), 
        duydx(mesh.GetNE()), 
        duxdy(mesh.GetNE()), 
        duydy(mesh.GetNE()), 
        sigmaxx(mesh.GetNE()),
        fatigue(mesh.GetNE()),
        stressed(mesh.GetNE()),
        maxSigma(mesh.GetNE()),
        elementCenters(mesh.GetNE()*mesh.Dimension()),
        sigmaCoeff(mu, lambda) {
    
    duxdxActual = false;
    duxdyActual = false;
    duydxActual = false;
    duydyActual = false;
    sigmaActual = false;
    fatigueActual = false;

    fatigue = 0;
    stressed = 0;
    maxSigma = 0.0;

    // initialize coefficients
    auto tab = fspace.GetElementToDofTable();
    mu = mfem::Vector(fspace.GetNDofs());
    lambda = mfem::Vector(fspace.GetNDofs());
    mfem::Array<int> row;
    for (int i = 0; i < tab.Size(); i++) {
        int attr = mesh.GetAttribute(i) - 1;
        tab.GetRow(i, row);
        for (int j = 0; j < row.Size(); j++) {
            lambda(row[j]) = lambdaV(attr);
            mu(row[j]) = muV(attr);
        }
    }
    
    mfem::Vector center(mesh.Dimension());
    for (int i = 0; i < mesh.GetNE(); i++) {
        mesh.GetElementCenter(i, center);
        for (int j = 0; j < mesh.Dimension(); j++) {
            elementCenters(i*mesh.Dimension() + j) = center(j);
        }
    }

}

void FatigueCalculator::update(mfem::Vector &newU) {
    u = newU;
    displacement.SetFromTrueDofs(u);
    // displacement.GetTrueDofs(u);
    duxdxActual = false;
    duxdyActual = false;
    duydxActual = false;
    duydyActual = false;
    sigmaActual = false;
    fatigueActual = false;
}

void FatigueCalculator::calcDuxdx() {
    if (!duxdxActual) {
        displacement.GetDerivative(1, 0, duxdxGf);
        mfem::IntegrationPoint point;
        for(int i = 0; i < duxdx.Size(); i++) {
            point.Set2(elementCenters(i*2), elementCenters(i*2 + 1));
            duxdx(i) = duxdxGf.GetValue(i, point);
        }
        duxdxActual = true;
    }
}

void FatigueCalculator::calcDuxdy() {
    if (!duxdyActual) {
        displacement.GetDerivative(1, 1, duxdyGf);
        duxdyGf.GetTrueDofs(duxdy);
        duxdyActual = true;
    }
}

void FatigueCalculator::calcDuydx() {
    if (!duydxActual) {
        displacement.GetDerivative(2, 0, duydxGf);
        duydxGf.GetTrueDofs(duydx);
        duydxActual = true;
    }
}

void FatigueCalculator::calcDuydy() {
    if (!duydyActual) {
        displacement.GetDerivative(2, 1, duydyGf);
        mfem::IntegrationPoint point;
        for(int i = 0; i < duxdx.Size(); i++) {
            point.Set2(elementCenters(i*2), elementCenters(i*2 + 1));
            duydy(i) = duydyGf.GetValue(i, point);
        }
        duydyActual = true;
    }
}

void FatigueCalculator::calcSigma() {
    if (!sigmaActual) {
        calcDuxdx();
        calcDuydy();
        int mi=0;
        double ms=0;
        for (int i = 0; i < duxdx.Size(); i++){
            double l = (duxdx(i) + duydy(i)) * lambda(i);
            double sx = duxdx(i) * mu(i) * 2 + l;
            double sy = duydy(i) * mu(i) * 2 + l;
            sigma(i) = sqrt(sx*sx + sy*sy);
        }
        sigmaGf.SetFromTrueDofs(sigma);
        sigmaActual = true;
    }
}

void FatigueCalculator::updateFatigue() {
    if(!fatigueActual) {
        calcSigma();
        for (int i = 0; i < fatigue.Size(); i++) {
            if (stressed[i]) {
                if (sigma(i) < 0.8e13) {
                    stressed[i] = 0;
                }
            } else {
                if (sigma(i) > 1.2e13) {
                    stressed[i] = 1;
                    fatigue[i] += 1;
                }
            }
            maxSigma(i) = std::max(maxSigma(i), sigma(i));
            fatigueActual = true;
        }
    }
}

void FatigueCalculator::printDuxdx(std::ostream &out) {
    calcDuxdx();
    duxdxGf.SaveVTK(out, "duxdx", 4);
}

void FatigueCalculator::printDuxdy(std::ostream &out) {
    calcDuxdy();
    duxdyGf.SaveVTK(out, "duxdy", 4);
}

void FatigueCalculator::printDuydx(std::ostream &out) {
    calcDuydx();
    duydxGf.SaveVTK(out, "duydx", 4);
}

void FatigueCalculator::printDuydy(std::ostream &out) {
    calcDuydy();
    duydyGf.SaveVTK(out, "duydy", 4);
}

void FatigueCalculator::printSigma(std::ostream &out) {
    calcSigma();
    sigmaGf.SaveVTK(out, "sigma", 4);
}

void FatigueCalculator::printFatigue(std::ostream &out) {
    updateFatigue();
    fatigue.Print(out, 1);
}

void FatigueCalculator::printMaxSigma(std::ostream &out) {
    updateFatigue();
    maxSigma.Print(out, 1);
}
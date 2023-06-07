#include "mfem.hpp"

bool checkDiag(mfem::SparseMatrix m){
    for(int i = 0; i < m.Height(); i++) 
        for(int j = m.GetI()[i]; j < m.GetI()[i+1]; j++) 
            if ((m.GetJ()[j] != i) & (abs(m.GetData()[j]) > 1e-15)) 
                return 0;
    return 1;
}

class ElasticWaveOperator2D : public mfem::SecondOrderTimeDependentOperator {
protected:
    mfem::FiniteElementSpace fespace;
    mfem::Array<int> listOfEssentialDOFs;
    mfem::BilinearForm M, K;
    mfem::LinearForm   * f;
    mfem::SparseMatrix * KMat;
    mfem::SparseMatrix * MMatInv;

public:
    ElasticWaveOperator2D(
                        mfem::FiniteElementSpace  &fspace, 
                        mfem::ConstantCoefficient         &rhoCoef, 
                        mfem::Coefficient         &lambdaCoef,
                        mfem::Coefficient         &muCoef):
        SecondOrderTimeDependentOperator(fspace.GetTrueVSize(), 0.0),
        fespace(fspace),
        M(&fespace),
        K(&fespace),
        f(nullptr),
        KMat(nullptr),
        MMatInv(nullptr) {

    const int order = fespace.GetOrder(0);
    mfem::IntegrationRule i_rule = mfem::IntegrationRules(0, mfem::Quadrature1D::GaussLobatto).Get(mfem::Geometry::SQUARE,2*order-1);

    // Setting boundary markers for the essential (Dirichlet) boundary conditions
    mfem::Array<int> essentialBoundaryMarkers(fespace.GetMesh()->bdr_attributes.Max());
    mfem::Array<int> subListOfEssDOFs;
    // essentialBoundaryMarkers = 0;
    // essentialBoundaryMarkers[2] = 1; // bottom
    // fespace.GetEssentialTrueDofs(essentialBoundaryMarkers, subListOfEssDOFs, 1);
    // listOfEssentialDOFs.Append(subListOfEssDOFs);
    essentialBoundaryMarkers = 0;
    essentialBoundaryMarkers[3] = 1; // right
    essentialBoundaryMarkers[1] = 1; // left
    fespace.GetEssentialTrueDofs(essentialBoundaryMarkers, subListOfEssDOFs, 0); // x
    listOfEssentialDOFs.Append(subListOfEssDOFs);
    // essentialBoundaryMarkers = 0;
    // essentialBoundaryMarkers[4] = 1; // hole x
    // fespace.GetEssentialTrueDofs(essentialBoundaryMarkers, subListOfEssDOFs, 0);
    // listOfEssentialDOFs.Append(subListOfEssDOFs);
    // essentialBoundaryMarkers = 0;
    // essentialBoundaryMarkers[4] = 1; // hole y
    // fespace.GetEssentialTrueDofs(essentialBoundaryMarkers, subListOfEssDOFs, 1);
    // listOfEssentialDOFs.Append(subListOfEssDOFs);
    

    // Compute inverse mass matrix MMatInv
    M.AddDomainIntegrator(new mfem::VectorMassIntegrator(rhoCoef, &i_rule));
    M.Assemble();
    M.Finalize();

    mfem::SparseMatrix Mmat;
    M.FormSystemMatrix(listOfEssentialDOFs, Mmat);
    if (!checkDiag(Mmat)) {
        std::cout << "Mass matrix is not diagonal!\n\
                      Make sure your mesh contains only squares!\n";
        exit(1);
    }
    mfem::Vector diagonal(fespace.GetTrueVSize());
    Mmat.GetDiag(diagonal);
    Mmat.~SparseMatrix();
    for (int i = 0; i < diagonal.Size(); i++) {
        if (diagonal(i) < 1e-12) {
            std::cout << "Zero encountered" << std::endl;
            exit(1);
        } else
            diagonal(i) = 1.0 / diagonal(i); // FIXME: possible division by zero
    }

    MMatInv = new mfem::SparseMatrix(diagonal);
    MMatInv->Finalize();
    diagonal.Destroy(); // freeing memory

    // Compute stiffness matrix Kmat
    KMat = new mfem::SparseMatrix();
    auto elasticityIntegrator = new mfem::ElasticityIntegrator(lambdaCoef, muCoef);
    elasticityIntegrator->SetIntegrationRule(i_rule);
    K.AddDomainIntegrator(elasticityIntegrator);
    K.Assemble();
    K.Finalize(0); // TODO: what does skip_zeros mean?
    K.FormSystemMatrix(listOfEssentialDOFs, *KMat);

    // Form right-hand side f
    f = new mfem::LinearForm(&fespace);
    f->Assemble();
}
    
    void Mult(const mfem::Vector& u, const mfem::Vector &du, mfem::Vector &d2u) const {
        KMat->AddMult(u, *f, -1.0);
        MMatInv->Mult(*f, d2u);
    }

    void ImplicitSolve(const double beta, const double, const mfem::Vector& u,const mfem::Vector &du, mfem::Vector &d2u) {
        MFEM_ASSERT((abs(beta) < 1e-12), "ImplicitSove implemented only for the explicit case with beta=0.")
        Mult(u, du, d2u);
    }

    void SetTime(const double t) {
        SecondOrderTimeDependentOperator::SetTime(t);
        *f = 0.0;
    }

    ~ElasticWaveOperator2D() {
        delete MMatInv;
        delete KMat;
    }
};

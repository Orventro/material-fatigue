#include <iostream>

#include"ElasticWaveOperator.hpp"

/**
 * Class implementing right-hand side of the equation u_tt = M^{-1} [ f - K u ].
 * Using this class is required for standard MFEM time-stepping algorithms.
 */
ElasticWaveOperator2D::ElasticWaveOperator2D(
                        mfem::FiniteElementSpace  &fspace, 
                        mfem::Coefficient         &rhoCoef, 
                        mfem::Coefficient         &lambdaCoef,
                        mfem::Coefficient         &muCoef):
        SecondOrderTimeDependentOperator(fspace.GetTrueVSize(), 0.0),
        fespace(fspace),
        M(&fespace),
        K(&fespace),
        //f_coef(rhs),
        f(nullptr),
        KMat(nullptr),
        MMatInv(nullptr) {

    const int order = fespace.GetOrder(0);
    mfem::IntegrationRule i_rule = mfem::IntegrationRules(0, mfem::Quadrature1D::GaussLobatto).Get(mfem::Geometry::SQUARE,2*order-1);

    // Setting boundary markers for the essential (Dirichlet) boundary conditions
    mfem::Array<int> essentialBoundaryMarkers(fespace.GetMesh()->bdr_attributes.Max());
    mfem::Array<int> subListOfEssDOFs;
    essentialBoundaryMarkers = 0;
    essentialBoundaryMarkers[0] = 1; // top
    fespace.GetEssentialTrueDofs(essentialBoundaryMarkers, subListOfEssDOFs, 1);
    essentialBoundaryMarkers = 0;
    essentialBoundaryMarkers[2] = 1; // bottom
    fespace.GetEssentialTrueDofs(essentialBoundaryMarkers, subListOfEssDOFs, 1);
    listOfEssentialDOFs.Append(subListOfEssDOFs);
    essentialBoundaryMarkers = 0;
    essentialBoundaryMarkers[3] = 1; // right
    essentialBoundaryMarkers[1] = 1; // left
    fespace.GetEssentialTrueDofs(essentialBoundaryMarkers, subListOfEssDOFs, 0);
    listOfEssentialDOFs.Append(subListOfEssDOFs);

    printf("Bounds done!\n");

    // Compute inverse mass matrix MMatInv
    M.AddDomainIntegrator(new mfem::VectorMassIntegrator(rhoCoef, &i_rule));
    M.Assemble();
    M.Finalize();
    
    mfem::SparseMatrix Mmat;
    M.FormSystemMatrix(listOfEssentialDOFs, Mmat);
    mfem::Vector diagonal(fespace.GetTrueVSize());
    Mmat.GetDiag(diagonal);
    Mmat.~SparseMatrix();
    // delete Mmat; // freeing memory, for some reason Destroy() method is protected
    for (int i = 0; i < diagonal.Size(); i++) {
        if (diagonal(i) < 1e-12)
            std::cout << "Zero encountered" << std::endl;
        else
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
    //f_coef.SetTime(0.0);
    //IntegrationRule integrationRuleEdge = IntegrationRules(0, Quadrature1D::GaussLobatto).Get(Geometry::SEGMENT,2*order-1);
    //f->AddBoundaryIntegrator(new VectorBoundaryFluxLFIntegrator(f_coef, 1e3/*, &integrationRuleEdge*/));
    f->Assemble();
}

    /**
     * Solve d2u = M^{-1} [ f - K u] with respect to d2u.
     * @param u   - displacement, known at current time
     * @param du  - velocity du/dt, known at current time; not used here
     * @param d2u - acceleration d^2u/dt^2, found at current time
     */
void ElasticWaveOperator2D::Mult(const mfem::Vector& u, const mfem::Vector &du, mfem::Vector &d2u) const {
    KMat->AddMult(u, *f, -1.0);
    MMatInv->Mult(*f, d2u);
};

    /**
     * Solve d2u = M^{-1} [ f - K (u + beta * d2u) ] with respect to d2u.
     * Assuming beta == 0!
     * @param beta  - parameter of Newmark scheme beta * dt^2, 0 for explicit scheme
     * @param gamma - parameter of Newmark scheme gamma * dt, not used here
     * @param u     - displacement at given time
     * @param du    - velocity du/dt at given time; not used here
     * @param d2u   - acceleration d^2u/dt^2, found at given time
     */
void ElasticWaveOperator2D::ImplicitSolve(const double beta, const double, const mfem::Vector& u,const mfem::Vector &du, mfem::Vector &d2u) {
    MFEM_ASSERT((abs(beta) < 1e-15), "ImplicitSove implemented only for the explicit case with beta=0.")
    Mult(u, du, d2u);
}

// void Elas

    /**
     * Set current time and update time-dependent parts of the operator (here only right-hand side).
     * @param t - current calculation time
     */
void ElasticWaveOperator2D::SetTime(const double t) {
    SecondOrderTimeDependentOperator::SetTime(t);
    //f_coef.SetTime(t);
    //f->Assemble();
    *f = 0.0;
}

ElasticWaveOperator2D::~ElasticWaveOperator2D() {
    // std::cout << "Deleting EWO" << std::endl;
    delete MMatInv;
    delete KMat;
    // delete f;
    // std::cout << "Deleted EWO" << std::endl;
};
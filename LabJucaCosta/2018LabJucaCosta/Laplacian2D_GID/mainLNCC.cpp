//
//  PZ
//

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"
#include "TPZMatLaplacian.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"


#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "pzgengrid.h"
#include "pzfunction.h"

#include "pzlog.h"

#include <iostream>
#include <math.h>
#include "pzskylstrmatrix.h"
#include "TPZSpStructMatrix.h"
#include <time.h>
#include <stdio.h>

#include "TPZSSpStructMatrix.h"

// To Identifier of the material into the domain
int ElementIDMat = 5;
int DimProblem = 2;

// Exact solution of the differential equation  (Laplacian)
static void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);

// Boundary conditions
static void NeumannEsquerda(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
static void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
static void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result);

// Permeability function to create material
void PermeabilityFunc(const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &perm);

// To create geometrical and computational meshes.
void CreatingGeometricMesh(TPZGeoMesh* gMesh,int elementIdMat);
void CreatingComputationalMesh(TPZCompMesh* cMesh, int elementIdMat, bool DGFEM);

// To apply uniform refinement on previous geometrical mesh
void UniformRefinement(TPZGeoMesh* gmesh,int ndiv);

// Utilitaries
void SetMaterialIdForElements(TPZGeoMesh *gmesh,int matid);
void SetMatBCIdForElements(TPZGeoMesh *gmesh,int matid, TPZVec<REAL> &P0,TPZVec<REAL> &P1);


// PROGRAM
int main(int argc, char *argv[]) {

    /// PRIORITARY -> RELATIVE TO GEOMETRICAL MESH
    // Type of elements
    int typeel = 3; // 3 - triangles, 4 - quadrilaterals
    // To create geometric mesh from GID file or not
    bool fromgid = true;
    TPZGeoMesh * gmesh = new TPZGeoMesh();
    if(!fromgid)
        CreatingGeometricMesh(gmesh,ElementIDMat);
    else {
        // IMPORTANDO MALHA GEOMETRICA DESDE ARQUIVO GID
        TPZReadGIDGrid myreader;
        std::string GeoGridFile;
        GeoGridFile="Curso1_Q.dump";
        gmesh = myreader.GeometricGIDMesh(GeoGridFile);
        // Inserting material id
        SetMaterialIdForElements(gmesh,ElementIDMat);
        TPZManVector<REAL> P0(3,0.);
        TPZManVector<REAL> P1(3,0.);
        P0[0] = -0.5;
        SetMatBCIdForElements(gmesh, -3, P0,P1);
        P0 = P1;
        P1[0] = 0.5;
        SetMatBCIdForElements(gmesh, -2, P0,P1);
        P0 = P1;
        P1[1] = 1.;
        SetMatBCIdForElements(gmesh, -4, P0,P1);
        P0 = P1;
        P1[0] = -0.5;
        SetMatBCIdForElements(gmesh, -6, P0,P1);
        P0 = P1;
        P1[1] = 0.;
        SetMatBCIdForElements(gmesh, -5, P0,P1);
    }
    
    // Refining mesh (uniform)
    int nref = 3;
    UniformRefinement(gmesh,nref);
    
    /// PRIORITARY -> RELATIVE TO COMPUTATIONAL MESH
    ///Indicacao de malha DGFem. Se false, vamos criar malha H1
    bool DGFEM = true;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    CreatingComputationalMesh(cmesh,ElementIDMat,DGFEM);
    
    // PRIORITARY -> RELATIVE TO SOLVING THE PROBLEM
    ///Analysis : construção do problema algébrico e inversão do sistema
    TPZAnalysis an(cmesh);

    ///Criando structural matrix - skyline não simétrica por causa do DGFEM usado
    TPZSkylineNSymStructMatrix matriz(cmesh);
    an.SetStructuralMatrix(matriz);
    
    ///Decomposição LU
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    
    ///Assemble da matriz de rigidez e vetor de carga
    an.Assemble();
    
    ///Resolução do sistema
    an.Solve();
    
    // UNNECESSARY -> RELATIVE TO CALCULATING ERROR AFTER SOLVE PROCESS
    ///Calculando erro da aproximacao
    an.SetExact(SolExataSteklov);///definindo solucao exata do problema
    TPZVec<REAL> erro;
    std::ofstream anPostProcessFile("postprocess.txt");
    an.PostProcess(erro,anPostProcessFile);///calculando erro
    std::cout << "\nErro de aproximação:\n";
    std::cout << "Norma H1 = " << erro[0] << "\nNorma L2 = " << erro[1]
    << "\nSeminorma H1 = " << erro[2] << "\n\n";
    
    // PRIORITARY -> RELATIVE TO PRINT SOLUTION TO VISUALIZATION BY PARAVIEW
    ///Exportando para Paraview
    TPZVec<std::string> scalarVars(1), vectorVars(0);
    std::stringstream sout;
    sout << "Laplacian" << DimProblem << "D_MESHINIT_E" << typeel << "H" << std::setprecision(2) << nref << "_cont" << DGFEM << ".vtk";

    scalarVars[0] = "Solution";
 
    an.DefineGraphMesh(DimProblem,scalarVars,vectorVars,sout.str());
    int resolution = 0;
    an.PostProcess(resolution);

    // Cleaning created meshes
    if(!DGFEM)
        delete cmesh;
    delete gmesh;

    return 0;
}

// To create simple geometrical mesh using 6 nodes
void CreatingGeometricMesh(TPZGeoMesh*gMesh,int elementIdMat) {
    
    ///Criando nós
    const int nnodes = 6;
    double coord[nnodes][2] = { {-0.5,0},{0,0},{0.,0.5},{-0.5,0.5},{0.5,0},{0.5,0.5} };
    for (int i = 0; i < nnodes; i++) {
        int nodind = gMesh->NodeVec().AllocateNewElement();
        TPZManVector<REAL, 3> nodeCoord(3);
        nodeCoord[0] = coord[i][0];
        nodeCoord[1] = coord[i][1];
        nodeCoord[2] = 0.;
        gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
    }
    
    ///Criando elementos
    const int nel = 2;
    int els[nel][4] = { {0,1,2,3},{1,4,5,2} };
    for (int iel = 0; iel < nel; iel++) {
        TPZManVector<int64_t, 4> nodind(4);
        int64_t index;
        nodind[0] = els[iel][0];
        nodind[1] = els[iel][1];
        nodind[2] = els[iel][2];
        nodind[3] = els[iel][3];
        gMesh->CreateGeoElement(EQuadrilateral, nodind, elementIdMat, index);
    }
    
    ///Criando elementos de contorno
    const int nelbc = 6;
    int bcels[nelbc][3] = { {0,1,-3},{1,4,-2},{4,5,-4},{5,2,-6},{2,3,-6},{3,0,-5} };
    for (int iel = 0; iel < nelbc; iel++) {
        TPZManVector<int64_t, 4> nodind(2);
        int64_t index;
        nodind[0] = bcels[iel][0];
        nodind[1] = bcels[iel][1];
        int matid = bcels[iel][2];
        gMesh->CreateGeoElement(EOned, nodind, matid, index);
    }
    
    ///Construindo conectividade da malha
    gMesh->BuildConnectivity();
    
    
}

///Refinamento uniforme da malha
void UniformRefinement(TPZGeoMesh* gmesh,int ndiv) {
    ///Inicializando padrões de refinamento uniforme
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    for (int i = 0; i < ndiv; i++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl * gel = gmesh->ElementVec()[iel];
            if (!gel) continue;
            if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
            TPZVec<TPZGeoEl*> filhos;
            gel->Divide(filhos);
        }///iel
    }///i
}
// Inserting material id to elements from GID mesh
void SetMaterialIdForElements(TPZGeoMesh *gmesh,int matid) {
    int nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        if (!gel) continue;
        if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
        if(gel->Dimension()!=DimProblem) continue;
        gel->SetMaterialId(matid);
    }///iel
}
// Inserting material id to elements from GID mesh
void SetMatBCIdForElements(TPZGeoMesh *gmesh,int matid, TPZVec<REAL> &P0,TPZVec<REAL> &P1) {
    int count = 0, nel = gmesh->NElements();
    TPZManVector<REAL> Center(3,0.), Point(3,0.);
    REAL t = -1., r=-1., p=-1.;
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        if (!gel) continue;
        if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
        if(gel->Dimension()!=DimProblem-1) continue;
        gel->CenterPoint(gel->NSides()-1,Center);
        gel->X(Center,Point);
        if(IsZero(P1[0]-P0[0])) {
            if(!IsZero(Point[0]-P0[0])) continue;
            t = 0.;
        } else
            t=(Point[0]-P0[0])/(P1[0]-P0[0]);
        if(IsZero(P1[1]-P0[1])) {
            if(!IsZero(Point[1]-P0[1])) continue;
            r = 0.;
        }
        else
            r=(Point[1]-P0[1])/(P1[1]-P0[1]);
        if(IsZero(P1[2]-P0[2])) {
            if(!IsZero(Point[2]-P0[2])) continue;
            p = 0.;
        }
        else
            p=(Point[2]-P0[2])/(P1[2]-P0[2]);
        
        if(r<0. || t<0. || p<0.) continue;
        if(r>1. || t> 1. || p>1.) continue;
        if((IsZero(t) && IsZero(r)) || (IsZero(t) && IsZero(p)) || (IsZero(r) && IsZero(p))) {
            gel->SetMaterialId(matid); count++;
            continue;
        }
        if((IsZero(t) && IsZero(r-p)) || (IsZero(r) && IsZero(t-p)) ||(IsZero(p) && IsZero(r-t))) {
            gel->SetMaterialId(matid);
            count++;
        }
    }///iel
    std::cout << std::endl << "BC Condition "  << matid << " com " << count << "elements." << std::endl;
}

// To computational mesh
void CreatingComputationalMesh(TPZCompMesh*cMesh,int elementIdMat,bool DGFEM) {
    /// criar materiais
    const int dim = 2;
    TPZMatLaplacian * mat = new TPZMatLaplacian(elementIdMat, dim);
    const STATE K = 1., F = 0.;
    mat->SetParameters(K, F);
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(PermeabilityFunc, 5);
    dummy->SetPolynomialOrder(0);
    TPZAutoPointer<TPZFunction<STATE> > func(dummy);
    mat->SetPermeabilityFunction(func);
    
    if (DGFEM) {
        ///Formulacao nao-simetria de Baumann, Oden e Babuska sem penalizacao
        mat->SetNoPenalty();
        mat->SetNonSymmetric();
    }
    cMesh->InsertMaterialObject(mat);
    
    const int pOrder = 4;
    
    ///Condições de contorno
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    TPZBndCond * BCondDirichletNulo = mat->CreateBC(mat, -3, 0, val1, val2);//0 = Dirichlet
    cMesh->InsertMaterialObject(BCondDirichletNulo);
    
    TPZBndCond * BCondNeumannZero = mat->CreateBC(mat, -2, 1, val1, val2);//1 = Neumann
    cMesh->InsertMaterialObject(BCondNeumannZero);
    
    TPZMaterial * BCondNeumannEsq = mat->CreateBC(mat, -5, 1, val1, val2);//1 = Neumann
    BCondNeumannEsq->SetForcingFunction(NeumannEsquerda, pOrder);
    cMesh->InsertMaterialObject(BCondNeumannEsq);
    
    TPZMaterial * BCondNeumannDir = mat->CreateBC(mat, -4, 1, val1, val2);//1 = Neumann
    BCondNeumannDir->SetForcingFunction(NeumannDireita, pOrder);
    cMesh->InsertMaterialObject(BCondNeumannDir);
    
    TPZMaterial * BCondNeumannAcima = mat->CreateBC(mat, -6, 1, val1, val2);//1 = Neumann
    BCondNeumannAcima->SetForcingFunction(NeumannAcima, pOrder);
    cMesh->InsertMaterialObject(BCondNeumannAcima);
    
    
    cMesh->SetDefaultOrder(pOrder);
    cMesh->SetDimModel(dim);//dim = 2
    
    if (DGFEM) {
        ///Criando malha de Galerkin descontínuo
        cMesh->SetAllCreateFunctionsDiscontinuous();
    }
    else {
        ///cria malha H1
        cMesh->SetAllCreateFunctionsContinuous();
    }

    ///Criando elementos computacionais
    cMesh->AutoBuild();

    if (DGFEM) {
        ///Cria elementos de interface
        TPZCreateApproximationSpace::CreateInterfaces(*cMesh);
    }
}

static void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du) {
    
    const REAL n = 0;
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    const REAL sol = pow((REAL)2,0.25 + n/2.)*pow(r,0.5 + n)*cos((0.5 + n)*t);
    u[0] = sol;
    
    du(0,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
    du(1,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));
    
}

static void NeumannEsquerda(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {   ///Jorge 2017, TPZFMatrix<STATE> &du){
    REAL normal[2] = {-1,0};
    
    TPZManVector<STATE> u(1);
    TPZFMatrix<STATE> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

static void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {   ///Jorge 2017, TPZFMatrix<STATE> &du){
    REAL normal[2] = {+1,0};
    
    TPZManVector<STATE> u(1);
    TPZFMatrix<STATE> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

static void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {   ///Jorge 2017, TPZFMatrix<STATE> &du){
    REAL normal[2] = {0,+1};
    
    TPZManVector<STATE> u(1);
    TPZFMatrix<STATE> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void PermeabilityFunc(const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &perm) {
    int i, j;
    int n = perm.Rows();
    int m = perm.Cols();
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            if (i==j || j==i-n)
                perm(i, j) = 1.;
            else
                perm(i, j) = 0.;
        }
    }
}


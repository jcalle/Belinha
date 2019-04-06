#include "myheader.h"
#include "pzerror.h"

#include "TPZCompElDisc.h"
#include "pztransientanalysis.h"
#include "pzpoisson3d.h"

#include "linlaw.h"
#include "linlawbi.h"
#include "conelaw.h"
#include "burger1d.h"
#include "burgerbi.h"
#include "euler.h"
#include "euler4c.h"
#include "eulerbi.h"

#include "pzfstrmatrix.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgnode.h"
#include "pzsolve.h"
#include "pzstepsolver.h"
#include "pzelmat.h"
#include "pzgraphmesh.h"

#include "filedat1d.h"
#include "filedat2d.h"
#include "conslaw.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "ttimeanalysis.h"
#include "tarungek.h"
#include "tgdimplicit.h"

#include "hadaptive.h"

#include <fstream>
#include <time.h>

#define MAXERROR 2.

using namespace std;

int m_rodada = -1;

void ExactSolutionBurger(const TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix<STATE> &deriv);
void ExactSolutionLinear(const TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix<STATE> &deriv);
void ExactSolutionAngle(const TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix<STATE> &deriv);
void ExactSolutionCone(const TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix<STATE> &deriv);
int AdjustNamePostProcess(std::string name, int from, int continuo, int order, int rank);

int main(int argc, char *argv[])
{
	int type, rank = 0, nuzero = 1;
	//  if(argc!=1)
	//    rank = atoi(argv[1]);

	  //Input file
	char name[128], filein[32], fileout[32];
	std::string lawname;
	strncpy(filein, "data00.dat", 16);
	ifstream input(filein);
	int i, j, n, order, posdim;
	posdim = 3;

	GetDataCommented(input, name, 64);  // Problem title
	//Constructing name of output file
	for (i = 0; i < 6; i++) {
		name[i] = (char)tolower(name[i]);
		if (name[i] == ' ' || name[i] == '.') break;
	}
	name[i] = '\0';
	lawname = name;
	strncpy(fileout, name, 8);
	strncat(fileout, "out.dat", 8);
	ofstream dataout(fileout, ios::app);
	int timedif = 0;
	time_t t1 = time(0);             // Get current time
	dataout << endl << "   *******       *******   " << endl << "HYPERBOLIC PROBLEM"
		<< endl << "Current time - STARTING " << ctime(&t1) << endl;
	char comando[32];
	strncpy(comando, "hostname >> ", 14);
	strncat(comando, fileout, 16);
	system(comando);
	dataout << endl;

	//Geometric mesh name
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	GetDataCommented(input, name, 128);
	gmesh->SetName(name);

	//Determining spatial dimension of the problem
	GetDataCommented(input, posdim);   // Problem geometrical dimension

	int matindex = -1;

	if (posdim == 2) {
		filein[6] = '\0';
		strncat(filein, "bi.dat", 7);
		input.close();
		input.open(filein);
		GetDataCommented(input, name, 64);  // Problem title
		GetDataCommented(input, name, 128);
	}

	// **To computational mesh
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	GetDataCommented(input, name, 128);
	cmesh->SetName(name);

	//Information to linear conservation law
	GetDataCommented(input, order);
	TConservationLaw *law;
	if (posdim == 1) {
		if (!strncmp(lawname.c_str(), "linear", 7)) law = new TLinearLaw(matindex, order);
		else if (!strncmp(lawname.c_str(), "burger", 7)) law = new TBurgerLaw1D(matindex);
		else if (!strncmp(lawname.c_str(), "euler", 6)) {
			if (order == 4) law = new TEulerLaw4C(matindex);
			else law = new TEulerLaw1D(matindex);
		}
		else { cout << "Undefined conservation law." << endl; return -10; }
	}
	else {
		if (!strncmp(lawname.c_str(), "linear", 7)) law = new TLinearLaw2D(matindex, order);
		else if (!strncmp(lawname.c_str(), "circle", 7)) law = new TLinearLaw2DCircle(matindex);
		else if (!strncmp(lawname.c_str(), "cone", 5)) law = new TConeLaw(matindex);
		else if (!strncmp(lawname.c_str(), "burger", 7)) law = new TBurgerLaw2D(matindex);
		else if (!strncmp(lawname.c_str(), "euler", 6)) law = new TEulerLaw2D(matindex);
		else { cout << "Undefined conservation law." << endl; return -23; }
	}
	order = law->NStateVariables();

	//Flux type to use as numerical flux for high and lower resolution
	GetDataCommented(input, type);
	int typelower;
	GetDataCommented(input, typelower);
	law->SetNumericalFluxType(type, typelower); // Numerical flux type
	//Interpolation order of the finite elements
	int ordercmesh;
	GetDataCommented(input, ordercmesh);
	dataout << "TPZCompEl::gOrder = " << ordercmesh << endl;
	cmesh->SetDefaultOrder(ordercmesh);
	//Construction of the computational mesh
	GetCommentary(input);
	n = law->Dimension() + 1;
	TPZVec<int> types(n);
	// Type of computational element to point, onedimen e/ou twodimen 
	for (i = 0; i < n; i++) {
		input >> types[i];
		SetCreateFunctionToElement(i, types[i]);   //SetCreateFunctionToElement(i,type);
	}

	if (posdim == 1) {
		TDatafile1d file1d(input, filein);
		file1d.InitializeMesh(cmesh, law);
	}
	else {
		TDatafile2d file2d(input, filein);
		file2d.InitializeMesh(cmesh, law);
	}

	//Recupering the Dirichlet value to UZero function
	int ncond = NumberOfDirichletConditions(*cmesh);
	ncond *= order;
	for (j = 0; j < ncond; j++)
		for (i = 0; i < order; i++)
			u0[j*order + i] = DetectDirichletValue(cmesh, i);
	//Making refinement in first mesh
	int level, continuo;
	while (strncmp(name, "Number of Cycle", 15)) input.getline(name, 128);
	input >> level;
	GetCommentary(input);
	// Type of computational element to point, onedimen e/ou twodimen 
	for (i = 0; i < n; i++)
		input >> types[i];
	if (level)
		Refinements(*cmesh, level, types, 0);

	GetDataCommented(input, continuo);
	if (continuo) {
		cmesh->SetAllCreateFunctionsContinuous();
		int64_t indexdisc, newindex;
		int64_t ii, nelems = cmesh->NElements();
		for(indexdisc=0;indexdisc<nelems;indexdisc++)
			cmesh->Discontinuous2Continuous(indexdisc, newindex);
	}
	AdjustNamePostProcess(lawname, 3, continuo, cmesh->GetDefaultOrder(), rank);

	// Analysis information
	TTimeAnalysis *an = 0;
	int m, aux, nerros, severalplots, explic;
	GetDataCommented(input, explic);
	if (explic > 0)
		an = new TTimeRungeKutta(input, dataout, cmesh, level);
	else
		an = new TGDImplicit(input, dataout, cmesh, level);

	//Applying initial function in T = 0 
	PutUZero(nuzero, an, input, posdim);

	//Inserting the stiffness matrix, it is band.
  //	TPZBandStructMatrix stiff(cmesh);
	TPZFStructMatrix stiff(cmesh);
	//	TPZBlockDiagonalStructMatrix stiff(cmesh);
	an->SetStructuralMatrix(stiff);
	TPZStepSolver<STATE> sol;
	sol.SetDirect(ELU);
	an->SetSolver(sol);

	//Stablizing initial function and defining variables to plot mesh
	TPZVec<std::string> scalnames(order);
	TPZVec<std::string> vecnames(0);
	law->VariablesName(scalnames);   // Names of the variables

	//Solve and print result.
	an->GetSchemeType(lawname);
	GetDataCommented(input, nerros);     // numero de rodadas para diferentes CFLs
	int nerroinit = 0;
	TPZVec<REAL> errorvec(3 * abs(nerros), 0.);   // abs -> valor absoluto de inteiros
	GetDataCommented(input, aux);        // armazena os primeros (severalplots) arquivos
	if (aux < 0) severalplots = nerros;   // de pos-processamento.
	else severalplots = aux;                    // quando aux<0 armazena todos!
	TPZVec<int> numberplot(severalplots, -1);
	int countplot = 0;
	int resolution = 1;
	if (aux < 0) for (m = 0; m < abs(nerros); m++) numberplot[m] = m;
	else for (m = 0; m < severalplots; m++) input >> numberplot[m];
	GetDataCommented(input, resolution);
	if (nerros < 0) {
		GetOldErros(errorvec, nerroinit, fileout);
		for (i = 0; i < nerroinit; i++) GetCommentary(input, 16);
		nerros *= -1;
	}
	for (m = nerroinit; m < nerros; m++) {
		m_rodada = m;
		// Limpando os dados no objeto analise que independem da malha e do material
		an->CleanToStartRun();
		// Completing the information for analysis object 
		an->ReadData(input);
		/*
			virtual void DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile);
    /** @brief Define GrapMesh as VTK with tensorial names depending on extension of the file 
		virtual void DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensnames, const std::string &plotfile);
		*/
		if (severalplots && countplot < severalplots && m == numberplot[countplot]) {
			countplot++;
			lawname = plotfile;
			an->DefineGraphMesh(posdim, scalnames, vecnames, plotfile);
			an->GraphMesh(posdim)->SetResolution(resolution);
		}
		else {
			plotfile[0] = '\0';
			an->DefineGraphMesh(posdim, scalnames, vecnames, plotfile);
		}
		t1 = time(0);
		timedif = t1;
		dataout << "\nRun " << m << "   Start " << ctime(&t1) << "   End ";

		//Solve the system 
		an->Run(input, dataout);

		//Compute execution time 
		t1 = time(0);
		timedif = t1 - timedif;
		dataout << ctime(&(t1 = time(0))) << "   Dif = " << timedif << endl;

		// To evaluate error 
		if (!strncmp(lawname.c_str(), "cir", 3))
			cmesh->EvaluateError(ExactSolutionAngle, true, errorvec);
		else if (!strncmp(lawname.c_str(), "bur", 3))
			cmesh->EvaluateError(ExactSolutionBurger,true, errorvec);
		else if (!strncmp(lawname.c_str(), "lin", 3))
			cmesh->EvaluateError(ExactSolutionLinear,true, errorvec);
		else if (!strncmp(lawname.c_str(), "con", 3))
			cmesh->EvaluateError(ExactSolutionCone,true, errorvec);
		for (i = 0; i < 3; i++) if (!(errorvec[3 * m + i] < BIGNUMBER)) errorvec[3 * m + i] = MAXERROR;
		dataout << "\nEstimated errors :\nTrue_error = " << errorvec[3 * m];
		dataout << "\nL2_error = " << errorvec[3 * m + 1] << "\nEstimate Sum = " << errorvec[3 * m + 2] << endl;
	}
	GetDataCommented(input, m);           // number of errors agrouped.
	PrintErrorsDX(errorvec, m, lawname );   // Print evaluated errors and time of the last evaluation.
	dataout << "FINISH -   Current time " << ctime(&(t1 = time(0))) << endl << endl;
	dataout << "   *******       *******   " << endl << endl;
	cout << "FINISH  -  File " << lawname << endl;

	// Clear analysis, and geometrical and computational mesh 
	input.close();
	dataout.close();
	delete an;
	delete cmesh;
	delete gmesh;

	return 0;
}

// To create geometrical and computational meshes.
void CreatingGeometricMesh(TPZGeoMesh* gMesh, int elementIdMat);
void CreatingComputationalMesh(TPZCompMesh* cMesh, int elementIdMat, bool DGFEM);

// To apply uniform refinement on previous geometrical mesh
void UniformRefinement(TPZGeoMesh* gmesh, int ndiv);

// Utilitaries
void SetMaterialIdForElements(TPZGeoMesh *gmesh, int matid);
void SetMatBCIdForElements(TPZGeoMesh *gmesh, int matid, TPZVec<REAL> &P0, TPZVec<REAL> &P1);

#include "belinha.h"

#include "TPZReadGIDGrid.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZMatLaplacian.h"

// To Identifier of the material into the domain
int ElementIDMat = 5;
int DimProblem = 2;

int mainpz(int argc, char *argv[])
{

	/// PRIORITARY -> RELATIVE TO GEOMETRICAL MESH
	// Type of elements
	int typeel = 3; // 3 - triangles, 4 - quadrilaterals
	// To create geometric mesh from GID file or not
	bool fromgid = true;
	TPZGeoMesh * gmesh = new TPZGeoMesh();
	if (!fromgid)
		CreatingGeometricMesh(gmesh, ElementIDMat);
	else {
		// IMPORTANDO MALHA GEOMETRICA DESDE ARQUIVO GID
		TPZReadGIDGrid myreader;

		std::string GeoGridFile = BELSOURCEDIR;
		if (typeel == 3)
			GeoGridFile += "/Projects_LabJC/Projects_LabJC/Laplacian2D_GID/Lapl2DT.dump";
		else
			GeoGridFile += "/Projects_LabJC/Projects_LabJC//Laplacian2D_GID/Lapl2DQ.dump";
		gmesh = myreader.GeometricGIDMesh(GeoGridFile);
		// Inserting material id
		SetMaterialIdForElements(gmesh, ElementIDMat);
		TPZManVector<REAL> P0(3, 0.);
		TPZManVector<REAL> P1(3, 0.);
		P0[0] = -0.5;
		SetMatBCIdForElements(gmesh, -3, P0, P1);
		P0 = P1;
		P1[0] = 0.5;
		SetMatBCIdForElements(gmesh, -2, P0, P1);
		P0 = P1;
		P1[1] = 1.;
		SetMatBCIdForElements(gmesh, -4, P0, P1);
		P0 = P1;
		P1[0] = -0.5;
		SetMatBCIdForElements(gmesh, -6, P0, P1);
		P0 = P1;
		P1[1] = 0.;
		SetMatBCIdForElements(gmesh, -5, P0, P1);
	}

	// Refining mesh (uniform)
	int nref = 0;
	UniformRefinement(gmesh, nref);

	/// PRIORITARY -> RELATIVE TO COMPUTATIONAL MESH
	///Indicacao de malha DGFem. Se false, vamos criar malha H1
	bool DGFEM = true;
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	CreatingComputationalMesh(cmesh, ElementIDMat, DGFEM);

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
	an.SetExact(ExactSolutionLinear);///definindo solucao exata do problema
	TPZVec<REAL> erro;
	std::ofstream anPostProcessFile("postprocess.txt");
	an.PostProcess(erro, anPostProcessFile);///calculando erro
	std::cout << "\nErro de aproximação:\n";
	std::cout << "Norma H1 = " << erro[0] << "\nNorma L2 = " << erro[1]
		<< "\nSeminorma H1 = " << erro[2] << "\n\n";

	// PRIORITARY -> RELATIVE TO PRINT SOLUTION TO VISUALIZATION BY PARAVIEW
	///Exportando para Paraview
	TPZVec<std::string> scalarVars(1), vectorVars(0);
	std::stringstream sout;
	sout << "Laplacian" << DimProblem << "D_MESHINIT_E" << typeel << "H" << std::setprecision(2) << nref << "_cont" << DGFEM << ".vtk";

	scalarVars[0] = "Solution";

	an.DefineGraphMesh(DimProblem, scalarVars, vectorVars, sout.str());
	int resolution = 0;
	an.PostProcess(resolution);

	// Cleaning created meshes
	if (!DGFEM)
		delete cmesh;
	delete gmesh;

	return 0;
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

// To create simple geometrical mesh using 6 nodes
void CreatingGeometricMesh(TPZGeoMesh*gMesh, int elementIdMat) {

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
void UniformRefinement(TPZGeoMesh* gmesh, int ndiv) {
	///Inicializando padrões de refinamento uniforme

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
void SetMaterialIdForElements(TPZGeoMesh *gmesh, int matid) {
	int nel = gmesh->NElements();
	for (int iel = 0; iel < nel; iel++) {
		TPZGeoEl * gel = gmesh->ElementVec()[iel];
		if (!gel) continue;
		if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
		if (gel->Dimension() != DimProblem) continue;
		gel->SetMaterialId(matid);
	}///iel
}
// Inserting material id to elements from GID mesh
void SetMatBCIdForElements(TPZGeoMesh *gmesh, int matid, TPZVec<REAL> &P0, TPZVec<REAL> &P1) {
	int count = 0, nel = gmesh->NElements();
	TPZManVector<REAL> Center(3, 0.), Point(3, 0.);
	REAL t = -1., r = -1., p = -1.;
	for (int iel = 0; iel < nel; iel++) {
		TPZGeoEl * gel = gmesh->ElementVec()[iel];
		if (!gel) continue;
		if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
		if (gel->Dimension() != DimProblem - 1) continue;
		gel->CenterPoint(gel->NSides() - 1, Center);
		gel->X(Center, Point);
		if (IsZero(P1[0] - P0[0])) {
			if (!IsZero(Point[0] - P0[0])) continue;
			t = 0.;
		}
		else
			t = (Point[0] - P0[0]) / (P1[0] - P0[0]);
		if (IsZero(P1[1] - P0[1])) {
			if (!IsZero(Point[1] - P0[1])) continue;
			r = 0.;
		}
		else
			r = (Point[1] - P0[1]) / (P1[1] - P0[1]);
		if (IsZero(P1[2] - P0[2])) {
			if (!IsZero(Point[2] - P0[2])) continue;
			p = 0.;
		}
		else
			p = (Point[2] - P0[2]) / (P1[2] - P0[2]);

		if (r < 0. || t < 0. || p < 0.) continue;
		if (r > 1. || t > 1. || p > 1.) continue;
		if ((IsZero(t) && IsZero(r)) || (IsZero(t) && IsZero(p)) || (IsZero(r) && IsZero(p))) {
			gel->SetMaterialId(matid); count++;
			continue;
		}
		if ((IsZero(t) && IsZero(r - p)) || (IsZero(r) && IsZero(t - p)) || (IsZero(p) && IsZero(r - t))) {
			gel->SetMaterialId(matid);
			count++;
		}
	}///iel
	std::cout << std::endl << "BC Condition " << matid << " com " << count << "elements." << std::endl;
}
void PermeabilityFunc(const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &perm) {
	int i, j;
	int n = perm.Rows();
	int m = perm.Cols();
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			if (i == j || j == i - n)
				perm(i, j) = 1.;
			else
				perm(i, j) = 0.;
		}
	}
}


// To computational mesh
void CreatingComputationalMesh(TPZCompMesh*cMesh, int elementIdMat, bool DGFEM) {
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
//	BCondNeumannEsq->SetForcingFunction(NeumannEsquerda, pOrder);
	cMesh->InsertMaterialObject(BCondNeumannEsq);

	TPZMaterial * BCondNeumannDir = mat->CreateBC(mat, -4, 1, val1, val2);//1 = Neumann
//	BCondNeumannDir->SetForcingFunction(NeumannDireita, pOrder);
	cMesh->InsertMaterialObject(BCondNeumannDir);

	TPZMaterial * BCondNeumannAcima = mat->CreateBC(mat, -6, 1, val1, val2);//1 = Neumann
//	BCondNeumannAcima->SetForcingFunction(NeumannAcima, pOrder);
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

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//void ComoChamarUZero() {
//	TPZFMatrix<REAL> InitialSol;
//	TPZTransientAnalysis<TPZMatPoisson3d> tranAnalysis;
//	TPZCompMesh *cmesh;
//	InitialSolutionLinearConvection(InitialSol, cmesh);

//	tranAnalysis.SetInitialSolution(InitialSol);

//}
void InitialSolutionLinearConvection(TPZFMatrix<REAL> &InitialSol, TPZCompMesh * cmesh) {
	InitialSol.Redim(cmesh->NEquations(), 1);
	InitialSol.Zero();
	for (int iel = 0; iel < cmesh->NElements(); iel++) {
		TPZCompEl * cel = cmesh->ElementVec()[iel];
		if (!cel) continue;
		TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cel);
		if (!disc) continue;
		if (disc->NConnects() == 0) continue;
		int bl = disc->Connect(0).SequenceNumber();
		int blpos = cmesh->Block().Position(bl);
		int blocksize = cmesh->Block().Size(bl);

		TPZGeoEl * gel = cel->Reference();
		TPZVec<REAL> xi(3), xVec(3);
		gel->CenterPoint(gel->NSides() - 1, xi);
		gel->X(xi, xVec);
		double x = xVec[0];
		double y = xVec[1];
		double u = 0.125;

		double xCircle = 0.25;
		double yCircle = 0.5;
		double R = 0.1;
		if ((x - xCircle)*(x - xCircle) + (y - yCircle)*(y - yCircle) <= R * R) u = 1.;

		InitialSol(blpos + blocksize - 20 + 0, 0) = u;
		InitialSol(blpos + blocksize - 20 + 1, 0) = 0.;
		InitialSol(blpos + blocksize - 20 + 2, 0) = 0.;
		InitialSol(blpos + blocksize - 20 + 3, 0) = 0.;
		InitialSol(blpos + blocksize - 20 + 4, 0) = 0.;

	}//for iel

	TPZVec<REAL> celerity(3, 0.);
	celerity[0] = 1.;
#ifdef LinearConvection
	TPZEulerEquation::SetLinearConvection(cmesh, celerity);
#endif

}//method

void ExactSolutionBurger(const TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix<STATE> &deriv) {
	if (loc[0] < .5) val[0] = 1.;
	else val[0] = -1.;
}
void ExactSolutionLinear(const TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix<STATE> &deriv) {
	//  if(loc[0] < -0.5) 
	//	  val[0] = 1.;
	//  else val[0] = 0.;

	/*  if(loc[0] < 0.25) val[0] = 1.;
	  else if(loc[0] < 0.75) val[0] = 0.;
	  else val[0] = -1.;
	  if(EndTime>.4) PZError << "Burger - SolucaoExata. Onda rarefacao e onda shock formaram um novo shock.\n";
	  REAL ul = 0., uc = 1., ur = 0.;
	  REAL x0 = .3, x1 = .5;
	  REAL x0fim = uc*EndTime+x0, x1fim = .5*(uc+ur)*EndTime+x1;
	  if(loc[0]<x0) val[0] = 0.;
	  else if(loc[0]<x0fim) val[0] = (1./EndTime)*(loc[0]-x0);
	  else if(loc[0]<x1fim) val[0] = 1.;
	  else val[0] = 0.;
	*/
	val[0] = sin(2 * PI_Value*loc[0]);

}
void ExactSolutionAngle(const TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix<STATE> &deriv) {
	if (loc[0] > 0. && loc[1] < 0.)
		val[0] = sin(2 * PI_Value*loc[0])*sin(-2.*PI_Value*loc[1]);
	else
		val[0] = 0.;
}
void ExactSolutionCone(const TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix<STATE> &deriv) {
	TPZVec<REAL> local = loc;
	UZeroCone(local, val, 0);
}

int AdjustNamePostProcess(std::string name, int from, int continuo, int order, int rank) {
	for (int i = 0; i < from; i++)
		if (name[i] == '\0') PZError << "AdjustNamePostProcess has string previously finishing" << endl;
	name[from] = '\0';
	char aux[8];
	int index = 0;
	aux[index++] = Itoa(rank);
	aux[index++] = Itoa(order);
	if (order > 9) aux[index++] = Itoa(order % 10);
	if (continuo) aux[index++] = 'c';
	else aux[index++] = 'd';
	aux[index++] = '\0';
	name += aux;
	return index;
}

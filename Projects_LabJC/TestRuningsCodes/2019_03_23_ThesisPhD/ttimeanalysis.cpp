/*******       File : tarungek.c

This file contains the method definitions for class TTimeAnalysis.

*******              *******/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "ttimeanalysis.h"
#include "pzsolve.h"
#include "pzstepsolver.h"
#include "pzbndmat.h"
#include "pzmvmesh.h"
#include "pzdxmesh.h"
#include "pzv3dmesh.h"
#include "conslaw.h"
#include "pzgnode.h"
#include "pzcompel.h"
#include "pzvec.h"
#include "pzgeoel.h"
#include "pzintel.h"
#include "interface.h"
#include "pzquad.h"

#include "pzbndcond.h"
#include "pzstrmatrix.h"

#include "pzelmat.h"
#include "myheader.h"
#include "hadaptive.h"

int ReAssembling = 0;

/*******       TTimeAnalysis, class derived from TPZAnalysis       *******/
TTimeAnalysis::TTimeAnalysis(std::istream &input, std::ostream &out,TPZCompMesh *mesh,int level)
     : TPZAnalysis(mesh,true,out) {
  fUk = 0;
  fName[0] = '\0';
  strncpy(fName,"Time_Analysis",15);
/*  GetDataCommented(input,fOrder);
  if(fOrder < 1  || fOrder > 4) {
   std::cout << "TTimeAnalysis constructor. Bad parameter fOrder. fOrder = 3.\n";
    fOrder = 1;
  }
  int i, nelem = mesh->NMaterials();
  for(i=0;i<4;i++)
    fCoef[i][0] = fCoef[i][1] = 0.;
  fCoef[0][1] = 1.;
  switch(fOrder) {
  case 1:
    break;
  case 2: {
    fCoef[1][0] = 1./2.;
    fCoef[1][1] = 1./2.;
    break;
  }
  case 3: {
    fCoef[1][0] = fCoef[1][1] = .25;
    fCoef[2][0] = 1./6.;
    fCoef[2][1] = 2./3.;
    break;
  }
  case 4: {
    fCoef[1][0] = 0.;
    fCoef[1][1] = 0.5;
    fCoef[2][0] = 0.;
    fCoef[2][1] = 1.;
    fCoef[3][0] = 1./3.;
    fCoef[3][1] = 1./6.;
    break;
  }
  }
  for(i=0;i<nelem;i++) {
    if(mesh->MaterialVec()[i]->Id() < 0) continue;
    fLaw = (TConservationLaw *)mesh->MaterialVec()[i];
    break;
  }
  GetDataCommented(input,fStartTime);
  GetDataCommented(input,fNEndTimes);
	fEndTime = new REAL[fNEndTimes];
	for(i=0;i<fNEndTimes;i++) input >> fEndTime[i];
  GetDataCommented(input,MINIMETIMESTEP);
  GetDataCommented(input,MAXROWS);
  /**Adaptive or non-adaptive scheme
  GetDataCommented(input,i);        // i=0 nonadaptive, i!=0 adaptive
	if(i) {
		/** Criando objeto para adaptatividade h-p ou do fluxo numerico
		fAdaptive = new TAdaptive(input,fLaw->Dimension(),1);
		fAdaptive->MinLevel = level;
	}
	else {
		GetCommentary(input,14);
		fAdaptive = 0;
	}
*/
  fCFL = 1.;
//  fCFLDiffussion = 0.;
  AppendMethodName(Name());
  fUZero = 0;
  fLevel = level+1;
//  fUseLimiter = 0;
/*  PartialStiff = 0;
  PartialRhs = 0;
  PartialSol = 0;*/
#ifdef PARALLELVERSION
  fRank = 0;
  fSize = 1;
#endif
}

TTimeAnalysis::TTimeAnalysis(std::istream &input, std::ostream &out,TPZCompMesh *cmesh,REAL TEnd) :
    TPZAnalysis(cmesh,true,out) {

  int i, nelem = cmesh->NMaterials();
  std::map<int, TPZMaterial * >::const_iterator mit;
  for (mit = fCompMesh->MaterialVec().begin(); mit != fCompMesh->MaterialVec().end(); mit++) {
	  TPZMaterial *mat = mit->second;
    if(mat->Id() < 0) continue;
    fLaw = (TConservationLaw *)mat;
    break;
  }
  fUk = 0;

/*  GetDataCommented(input,fName,2);
	GetDataCommented(input,fOrder);
  if(fOrder < 1  || fOrder > 3) {
   std::cout << "TTimeAnalysis constructor. Bad parameter fOrder. fOrder = 3.\n";
    fOrder = 1;
  }
  for(i=0;i<4;i++)
    fCoef[i][0] = fCoef[i][1] = 0.;
  fCoef[0][1] = 1.;
  switch(fOrder) {
  case 1:
    break;
  case 2: {
    fCoef[1][0] = fCoef[1][1] = 1./2.;
    break;
  }
  case 3: {
    fCoef[1][0] = fCoef[1][1] = .25;
    fCoef[2][0] = 1./6.;
    fCoef[2][1] = 2./3.;
    break;
  }
  case 4: {
    fCoef[1][0] = 0.;
    fCoef[1][1] = 0.5;
    fCoef[2][0] = 0.;
    fCoef[2][1] = 1.;
    fCoef[3][0] = 1./3.;
    fCoef[3][1] = 1./6.;
    break;
  }
  }
  GetDataCommented(input,fStartTime);
  GetDataCommented(input,fNEndTimes);
	fEndTime = new REAL[fNEndTimes];
	for(i=0;i<fNEndTimes;i++) input >> fEndTime[i];
  GetDataCommented(input,MINIMETIMESTEP);
  GetDataCommented(input,MAXROWS);
  /**Adaptive or non-adaptive scheme
  GetDataCommented(input,i);        // i=0 nonadaptive, i!=0 adaptive
	if(i) {
		/** Criando objeto para adaptatividade h-p ou do fluxo numerico
		fAdaptive = new TAdaptive(input,fLaw->Dimension(),1);
		fAdaptive->MinLevel = 1;
	}
	else fAdaptive = 0;
  if(fEndTime[0]<TEnd) {
    fStartTime = fEndTime[0];
    fEndTime[0] = TEnd;
  }
  GetDataCommented(input,fCFL);
  GetDataCommented(input,fCFLDiffussion);
  GetDataCommented(input,fRhsAdaptNorm);
  AppendMethodName(Name());
  fUZero = 0;
/*  PartialStiff = 0;
  PartialRhs = 0;
  PartialSol = 0;*/
#ifdef PARALLELVERSION
  fRank = 0;
  fSize = 1;
#endif
}

TTimeAnalysis::TTimeAnalysis(std::ostream &out) : TPZAnalysis(0,1,out) {
  fLaw = 0;
//  fOrder = 1;
	fRhsAdaptNorm = 0.75;
/*  for(i=0;i<4;i++)
    fCoef[i][0] = fCoef[i][1] = 0.;
  fCoef[0][1] = 1.;
  switch(fOrder) {
  case 1:
    break;
  case 2: {
    fCoef[1][0] = fCoef[1][1] = 1./2.;
    break;
  }
  case 3: {
    fCoef[1][0] = fCoef[1][1] = .25;
    fCoef[2][0] = 1./6.;
    fCoef[2][1] = 2./3.;
    break;
  }
  case 4: {
    fCoef[1][0] = 0.;
    fCoef[1][1] = 0.5;
    fCoef[2][0] = 0.;
    fCoef[2][1] = 1.;
    fCoef[3][0] = 1./3.;
    fCoef[3][1] = 1./6.;
    break;
  }
  }*/
  fStartTime = 0.;
	fNEndTimes = 1;
	fEndTime = new REAL(1.);
  MAXROWS = 700;
  fUk = 0;
  fUZero = 0;
	fAdaptive = 0;
//	fCFLDiffussion = 0.;
  MINIMETIMESTEP =  1.e-7;
/*  PartialStiff = 0;
  PartialRhs = 0;
  PartialSol = 0;*/
#ifdef PARALLELVERSION
  fRank = 0;
  fSize = 1;
#endif
}

TTimeAnalysis::~TTimeAnalysis() {
  if(fUk) delete fUk;
//  if(PartialStiff) delete PartialStiff;
//  if(PartialRhs) delete PartialRhs;
//  if(PartialSol) delete PartialSol;
	if(fAdaptive) delete fAdaptive;
	if(fEndTime) delete fEndTime;
}

/**Imprime estado atual da malha geometrica, computacional, solucao, etc*/
void TTimeAnalysis::PrintCurrentState() {
  /**Imprime objeto analysis corrente*/
	std::ofstream anout("analysis.dat");
  anout << "DATA TIME_ANALYSIS" <<std::endl;
  anout << "Name analysis" <<std::endl;
  anout << fName<< std::endl;
//  anout << "Order Runge-Kutta" <<std::endl;
//  anout << fOrder <<std::endl;
  anout << "Initial Time" <<std::endl;
  anout << fStartTime <<std::endl;
  anout << "End Time" <<std::endl;
  anout << fEndTime <<std::endl;
  anout << "CFL_Condition" <<std::endl;
  anout << fCFL <<std::endl;
//  anout << "Coeficiente de difussividade" <<std::endl;
//  anout << fCFLDiffussion <<std::endl;
  anout << "Refinement Level" <<std::endl;
  anout << fLevel <<std::endl;
	if(fAdaptive) {
	  anout << "Minime Refinement Level" <<std::endl;
  	anout << fAdaptive->MinLevel <<std::endl;
	}
//  anout << "Use Limiter" <<std::endl;
//  anout << fUseLimiter <<std::endl;
  anout << "Current Time (law)" <<std::endl;
  anout << CurrentTime();
  anout.close();
//  anout.open("gmesh.dat");
//  fCompMesh->Reference()->PrintData(anout);
//  anout.close();
//  anout.open("cmesh.dat");
//  fCompMesh->PrintData(anout);
//  anout.close();
}

int TTimeAnalysis::OrderLaw() { return fLaw->NStateVariables(); }

double TTimeAnalysis::CurrentTime() {
  if(fLaw) return fLaw->Time();
  return 0.;
}

/**To adequate the stiffness matrix and rhs depending on the fSolution vector.
   It is necessary when the mesh was refined and fSolution was interpolated*/
void TTimeAnalysis::AdequateMatrix() {
  int numeq = fCompMesh->NEquations();
  if(!fCompMesh || !fStructMatrix || !fSolver) {
    PZError << "TTimeAnalysis::AdequateMatrix. Hasn't stiffness or rhs matrix, or mesh.\n";
    return;
  }
	if(fRhs.Rows() != numeq) {
    fRhs.Redim(numeq,1);
    fSolution.Redim(numeq,1);
    fUk->Redim(numeq,1);
    fCompMesh->InitializeBlock();

    if(fSolution.Rows() != numeq)
      PZError << "TTimeAnalysis::AdequateMatrix has incompatible solution.\n";
  }
}

/**Assemble only the vector on the right hand*/
void TTimeAnalysis::AssembleRhs() {

  if(!fCompMesh || !fStructMatrix.operator->() || !fSolver) return;

	fRhs.Redim(fCompMesh->NEquations(),1);
  fRhs.Zero();
	fSolver->SetMatrix(0);
//	fSolver->SetMatrix((fStructMatrix.operator->()).CreateAssemble(fRhs)); //aqui TPZFMatrix não é nula
  fRhsNorm = Norm(fRhs);
}

void TTimeAnalysis::Solve() {
	int numeq = fCompMesh->NEquations();
	if(fRhs.Rows() != numeq) return;
  fSolution.Zero();
	TPZFMatrix<STATE> residual(fRhs);
	TPZFMatrix<STATE> delu(numeq,1,0.);
	  /** Prints only by tests */
//	ofstream outmat("matrizes.txt",ios::app);
//	  fSolver->Matrix()->Print("fMat_3",outmat);
//	  fUk->Print("fUk_3",outmat);
//	  fRhs.Print("fRhs_3",outmat);
//	  fSolution.Print("fSol_1",outmat);
//	  outmat.flush();
	if(fSolution.Rows() != numeq) {
	  fSolution.Redim(numeq,1);
	} else {
	  fSolver->Matrix()->Residual(fSolution,fRhs,residual);
	}
//	  fSolver->Matrix()->Print("fMat_4",outmat);
//	  fUk->Print("fUk_4",outmat);
//	  fRhs.Print("fRhs_4",outmat);
//	  fSolution.Print("fSol_2",outmat);
//	  outmat.close();
//	REAL normres  = Norm(residual);
	fSolver->Solve(residual, delu);
	fSolution += delu;
  ApplyBCDirichlet();
  LoadSolution();
}

int TTimeAnalysis::MatrixTruncated(int n) {
  int i, band = ((TPZFBMatrix<STATE> *)(fSolver->Matrix()).operator->())->GetBand();
  int row = n*MAXROWS-1;
  int nvar = OrderLaw(), numeq = fCompMesh->NEquations();
  if((row+nvar*band+1) >= numeq) return numeq;
  int j, k;
  char mask;
	TPZMatrix<STATE> *matrix = fSolver->Matrix().operator->();
  for(i=0;i<(nvar+1)*band;i++) {
    mask = 'n';
    /**Achando linha cujo ultimo valor nao zero esteja na diagonal */
    for(j=0;j<band;j++) {
      if(matrix->GetVal(row,row+j+1)!=0.) { mask = 'y'; break; }
    }
    /**Verificando se estamos na linha final do bloco achado para truncamento */
    k = 1;
    while(mask=='n' && k<band) {
      for(j=0;j<band-k;j++) {
        if(matrix->GetVal(row-k,row+1+j)!=0.) { mask = 'y'; break; }
      }
      k++;
    }
    row++;
    /**Verificando que a proxima linha pertence a um outro bloco */
    if(mask=='n') return row;
  }
 std::cout << "TTimeAnalysis::MatrixTruncated. Insecure partition matrix." <<std::endl;
  return numeq;
}

void TTimeAnalysis::SolveToImplicitChunks() {
/*	TPZMatrixSolver *solver = new TPZStepSolver(*((TPZStepSolver *)fSolver));
  int neq = fCompMesh->NEquations();
  int ngroups = (neq/MAXROWS)+1;
	TPZMatrix *matrix = fSolver->Matrix();
  int band = ((TPZFBMatrix *)matrix)->GetBand();
  int i, j, k, nn = 0;   //posicao da primeira linha em fStiffness para ser trabalhada
  int nrowi, nrows = 0;
	
//	solver->SetMatrix(PartialStiff);

  /**Redimensiona as matrices e vetor auxiliares para particionar as matrices grandes
  for(i=0;i<ngroups;i++) {
    /**Averigua linha perto da linha n*mil onde eh possivel truncar a matrix
    nrows = MatrixTruncated(i+1);
    nrowi = nrows - nn;

    PartialRhs->Redim(nrowi,1);
    PartialSol->Redim(nrowi,1);
    PartialStiff->Redim(nrowi,nrowi);
    /**Passando os valores das matrices grandes as sub-matrices
    for(j=0;j<nrowi;j++) {
      (*PartialStiff)(j,j) = matrix->GetVal(nn+j,nn+j);
      (*PartialRhs)(j,0) = fRhs.GetVal(nn+j,0);
    }
    for(j=0;j<band;j++) {
      for(k=j+1;k<nrowi;k++) {
        (*PartialStiff)(k-j-1,k) = matrix->GetVal(nn+k-j-1,nn+k);
        (*PartialStiff)(k,k-j-1) = matrix->GetVal(nn+k,nn+k-j-1);
      }
    }
    /**Resolvendo sub-matriz sequencialmente por emquanto
    PartialStiff->SetIsDecomposed(0);
//		fSolver->SetMatrix(0);
	  fSolver->SetMatrix(PartialStiff);
    fSolver->Solve(*PartialRhs,*PartialSol);
    for(j=0;j<nrowi;j++)
      fSolution(nn+j,0) = (*PartialSol)(j,0);
    nn = nrows;
  }   // Ateh aqui pode ser nomeada de funcao AdequateMatrixToParallel
	matrix->SetIsDecomposed(1);
	fSolver->SetMatrix(matrix);
	/** Delete temporary solver 
//	if(solver) delete solver;
*/
}

void TTimeAnalysis::SolveToImplicitParallel() {
#ifdef PARALLELVERSION
if(!fRank) {
  int neq = fCompMesh->NEquations();
	TPZMatrix *matrix = fStructMatrix->Matrix();
  int ngroups = (neq/MAXROWS)+1;
  int band = ((TPZFBMatrix *)matrix)->GetBand();
  int i, j, k, nn = 0;   //posicao da primeira linha em fStiffness para ser trabalhada
  int nrowi, nrows=0;
//  fSolver->SetMatrix(PartialStiff);

  /**Redimensiona as matrices e vetor auxiliares para particionar as matrices grandes*/
  for(i=0;i<ngroups;i++) {
    /**Averigua linha perto da linha n*mil onde eh possivel truncar a matrix*/
    nrows = MatrixTruncated(i+1);
    nrowi = nrows - nn;
    
    PartialRhs->Redim(nrowi,1);
    PartialSol->Redim(nrowi,1);
    PartialStiff->Redim(nrowi,nrowi);
    /**Passando os valores das matrices grandes as sub-matrices*/
    for(j=0;j<nrowi;j++) {
      (*PartialStiff)(j,j) = matrix->GetVal(nn+j,nn+j);
      (*PartialRhs)(j,0) = fRhs.GetVal(nn+j,0);
    }
    for(j=0;j<band;j++) {
      for(k=j+1;k<nrowi;k++) {
        (*PartialStiff)(k-j-1,k) = matrix->GetVal(nn+k-j-1,nn+k);
        (*PartialStiff)(k,k-j-1) = matrix->GetVal(nn+k,nn+k-j-1);
      }
    }
    /**Resolvendo sub-matriz sequencialmente por emquanto*/
    PartialStiff->SetIsDecomposed(0);
    fSolver->Solve(*PartialRhs,*PartialSol);
    for(j=0;j<nrowi;j++)
      fSolution(nn+j,0) = (*PartialSol)(j,0);
    nn = nrows;
  }   // Ateh aqui pode ser nomeada de funcao AdequateMatrixToParallel
}
#else
  return;
#endif
}

void TTimeAnalysis::SolveToImplicit() {
  fSolution.Zero();

// if(fExplicit < 0) {
//#ifdef PARALLELVERSION
//    SolveToImplicitParallel();
//#else
//    SolveToImplicitChunks();
//#endif
//  }
//  else {
    fSolver->Matrix()->SetIsDecomposed(0);
    fSolver->Solve(fRhs,fSolution);
//  }

  /**Aplicando condicoes Dirichlet se tiver e transferir a solucao para a malha*/
  ApplyBCDirichlet();
  LoadSolution();
}

void TTimeAnalysis::GetSchemeType(std::string &filename) {
  /**We can to choose implicit or explicit scheme to diffusion*/
  char aux[16];
  int index=0;
  aux[index++] = 'R'; aux[index++] = 'K';
  aux[index++] = '1';
 std::cout << "Time Analysis Adaptive - Parallel and Implicit scheme.\n";
  aux[index++] = 'p'; aux[index++] = 'i';
	aux[index++] = '\0';
  filename = aux;
}

void TTimeAnalysis::ReDraw(int stepgraph,int ending) {
	for(int k=0;k<3;k++) {
		if(fGraphMesh[k]) {
			if(fAdaptive)
				ReDefineGraphMesh(k+1);

			if(ending) {
				if(fAdaptive) fGraphMesh[k]->DrawMesh(stepgraph+1);
				((TPZDXGraphMesh *)fGraphMesh[k])->SetNumCases(stepgraph+1);
			}
			else {
				if(fAdaptive) fGraphMesh[k]->DrawMesh(stepgraph+300);
			}
			fGraphMesh[k]->DrawSolution(stepgraph,CurrentTime());  //,fScalarNames[k],fVectorNames[k]);
			if(ending) ((TPZDXGraphMesh *)fGraphMesh[k])->Close();
		}
	}
	std::cout << "eq: " << fCompMesh->NEquations() << " el: ";
	CountElements(*fCompMesh);
}

/**Put cfldiff as diffusion coefficient into the computational elements of the vector 
void TTimeAnalysis::SetFactor(TPZAdmChunkVector<TPZCompEl *> &elvec,REAL cfldiff) {
	int i, n = elvec.NElements();
	TPZCompEl *cel;
	for(i=0;i<n;i++) {
		cel = elvec[i];
		if(!cel || !cel->IsInterpolated()) continue;
		((TPZInterpolatedElement *)cel)->SetDiffusionCoefficient(cfldiff);
	}
}*/

void TTimeAnalysis::AdjustBeforeAssemble(REAL time) {
  fLaw->SetTime(time,fTimeStep);
  fLaw->SetMaxEigen(0.);
  ReAssembling = 0;
}
void TTimeAnalysis::ReAssemble(REAL time) {
  fTimeStep *= 0.5;
  /*Incrementa-se o termo diffussivo para aproximar com volumes finitos*/
  AdjustBeforeAssemble(time);
  Assemble();
}

void TTimeAnalysis::Run(std::istream &input, std::ostream &out) {
  double ctime = 0.;
  int numitera=0, refine = 1;
//  int niterations;
  std::string lawname;
  lawname = plotfile;
  if(!fCompMesh || !fLaw) {
    PZError << "TTimeAnalysis::Run, hasn't computational mesh or equation associated.\n";
    return;
  }

//  int l, k, count;
  int stepgraph = 0;
  int onlyrhs = 0;
  /** Put initial time in fTimes and into conservation law*/
  ctime = fStartTime;

  /** Initializating the solution at initial time 
  for(k=0;k<3;k++) 
		if(fGraphMesh[k]) {
			fGraphMesh[k]->SetNames(fScalarNames[k],fVectorNames[k]);
//			fGraphMesh[k]->DrawMesh(stepgraph+300);
		}

  ComputeAreas(fCompMesh,&minarea);
  fLaw->SetMinDeltaX(minarea);

  AdequateMatrix();
  ApplyUZero();   // Discretization of UZero

  /** Adapting the numerical flux: high aprox. into discontinuous region
     and lower aprox. into the regular regions
	if(fAdaptive) fAdaptive->SteadyState(fSteadyState);
  Adapting(stepgraph,0);

  /**Here: We can to change the mesh and acerting stiff matrix, solver and fUk
  ApplyBCNonFlux();
  AdequateMatrix();
  /**To temporary storage of L(u0)+L(u1)+L(u2)+...
  TPZFMatrix RhsSum;
	int end = 0;
do {
  /** Numero de plots para o arquivo de post-processo 
	count = 0;
  TPZVec<REAL> TimePost(fNPlots);
  REAL timepost = (fEndTime[end]-ctime)/fNPlots;
  TimePost[count] = fStartTime + ctime;
  for(l=1;l<fNPlots-1;l++) TimePost[l] = TimePost[l-1]+timepost;
  TimePost[l] = fEndTime[end];

for(k=0;k<3;k++)
	if(fGraphMesh[k])
		fGraphMesh[k]->DrawMesh(stepgraph+300);
stepgraph = 0;

  /** Iterating while current time is lower than end time or it achieved maxime iterations
  while(ctime < fEndTime[end] && numitera < MAXITERATIONS) {
    /**Incrementa numero de iteracoes
    numitera++;

    /**Compute time step from CFL condition (and adjust next time at end time)
       and compute maxime area of the 
    ComputeNextTime(ctime,fEndTime[end]);

    /**Adapting the numerical flux: high aprox. into discontinuous region
       and lower aprox. into the regular regions
    AdjustBeforeAssemble(ctime);
    /** Drawing to Post-processing 
    if(ctime>=TimePost[count]) {
      count++;
      ReDraw(stepgraph++);
    }

    /**Fixing the solution at time tk and adjusting solver
    fSolutionNext = fCompMesh->Solution();
    if(fExplicit>0) {
      /** Pode-se trabalhar com sub-malhas, assim sobre ellas pode dar-se
          um valor diferente ao tipo de fluxo a utilizar.
      /**Iterating over the order of the Runge-Kutta solver
      for(l=0;l<fOrder;l++) {
				if(fExplicit==1) CleanDiffusion();
        if(ReAssembling && TimePost[count]!=fEndTime[end]) {
          ReAssemble(ctime);
          onlyrhs = 0;
        }
        else {
          if(!onlyrhs) {
            Assemble();
            onlyrhs = 1;
          }
          else AssembleRhs();
        }

        /**Computing the sumatory of operators c(l,0)*[H(u0)+H(u1)+...] 
        if(!l) { RhsSum = fRhs; RhsSum.Zero(); }
        (*fUk) = (fCoef[l][0]*fTimeStep)*RhsSum;
        RhsSum += fRhs;
        fRhs *= (fCoef[l][1]*fTimeStep);
        fRhs += (*fUk);
        if(!onlyrhs) {
          ApplyBC();
          onlyrhs = 1;
        }
        else ApplyRhsBC();
				fPreviousSol = fSolution;
        /**Solving K(u(l+1)-u(k)) = H and finding u(l+1)
        Solve();
				fPreviousSol -= fSolution;
				fRhsNorm = Norm(fPreviousSol);
        /**It is necessary to adjust the time data for next intermediate computation
      }
    }
    else {  //Case implicito
			// construindo matrix com igual dimensao que fRhs
      TPZFMatrix RhsOnlyIntegral(fRhs);
			// Calcula o termo produto interno das funcoes teste com as variaveis de estado
      AssembleOnlyIntegral(RhsOnlyIntegral);
      /**Iterating over the order of the Runge-Kutta solver
			// seta o coeficiente para o metodo Runge-Kutta (Euler explicito)
      fLaw->SetCoef(1.);
			//Parametro Reassembling vira 1 se acontece um gradiente muito grande, caso a malha seja modificada
      if(ReAssembling && TimePost[count]!=fEndTime[end]) ReAssemble(ctime);
      else Assemble();
			
      /**Computing the sumatory of operators c(l,0)*[H(u0)+H(u1)+...] 
      fRhs *= (fCoef[0][1]*fTimeStep);
      fRhs += RhsOnlyIntegral;
      /**Solving K(u(l+1)-u(k)) = H and finding u(l+1)

      ApplyBC();
			fPreviousSol = fSolution;
			fSolution.Zero();
			Solve(); // SolveToImplicit();
			if(fSteadyState) {
				fPreviousSol -= fSolution;
  	    fRhsNorm = Norm(fPreviousSol);
			}
    }
		/**Adapting the numerical flux: high aprox. into discontinuous region
			and lower aprox. into the regular regions and hp-adaptivity
		refine = Adapting(stepgraph,0);
		if(refine) {
			numitera++;
			onlyrhs = 0;
		}
		/** Imprime resultados nao para todo passo de tempo. Pode depender de parametro 
		if(!(numitera%30)) {
			out << "t: " << CurrentTime() << " eq: " << fCompMesh->NEquations();
			if(fSteadyState) out << " fRhsNorm = " << fRhsNorm;
			out <<std::endl;
			CountElements(*fCompMesh);
		}
		if(fUseLimiter) {
			for(int nn=0;nn<fLaw->NStateVariables();nn++)
				CockburnLimiter(fCompMesh->ElementVec(),nn);
		}

		/** Acerting current time  and the CFL value
		ctime += fTimeStep;
	}

	/** Asking to continue process, can to put new end time 
	out << "\nLast End Time = " << fEndTime[end];
	end++;

  /**Imprime estado atual da malha geometrica, computacional, solucao, etc
  AdjustBeforeAssemble(ctime);
	CountElements(*fCompMesh);
  /**To post-processing at last time -> (ending=1) 
  ReDraw(stepgraph,1);
  for(k=0;k<3;k++)
		if(fGraphMesh[k] && end < fNEndTimes) {
			ofstream *pout = (ofstream *)(fGraphMesh[k]->Out());
			pout->close();
			delete pout;
			PlotFileName(lawname,end,plotfile);
			pout = new ofstream(plotfile);
			fGraphMesh[k]->SetOutFile(*pout);
    }

}while(!IsZero(ctime-fEndTime[fNEndTimes-1]));
*/
}

int TTimeAnalysis::Adapting(int &stepgraph,int onlyflux) {
	/** Clean object to adaptive work if is the case */
	if(!fAdaptive) return 0;
  REAL minarea;
  fAdaptive->CleanVectors();
	
  /** Search wavelet coefficients to refine or coarse the elements */
	if((fSteadyState && fRhsAdaptNorm < fRhsNorm) || !fAdaptive->IsNecessaryToRefine(*fCompMesh,fLevel,onlyflux))
	  return 0;
		
  /** Making adaptive refinement */
  fAdaptive->WaveletHPAdaptive(*fCompMesh);

  ComputeAreas(fCompMesh,&minarea);
  fLaw->SetMinDeltaX(minarea);

  /**Here: We can to change the mesh and acerting stiff matrix, solver and fUk*/
  AdequateMatrix();
	return 1;
}

void TTimeAnalysis::ComputeNextTime(double ctime,REAL EndTime) {
  /*Computing time step from CFL Condition*/
  fTimeStep = StabilityCondition();
  /*Adjusting dt whether it is very small*/
  if(fTimeStep<MINIMETIMESTEP) fTimeStep = MINIMETIMESTEP;

  /**Acerting time step near of end time and store current time*/
  if(ctime+fTimeStep > EndTime)
    fTimeStep = EndTime - ctime;
}

void TTimeAnalysis::Print(char *name, std::ostream &out,int all) {
  out <<std::endl << name <<std::endl <<std::endl << "TTimeAnalysis Analysis" <<std::endl;
//  out << "\tRungeKutta method order    = " << fOrder <<std::endl;
  out << "\tSpatial Analysis used name = " << Name() <<std::endl;
  out << "\tInitial Time               = " << fStartTime <<std::endl;
  out << "\tFinish Time                = " << fEndTime <<std::endl;
  if(all) PrintSolution(out);
}

void TTimeAnalysis::PrintSolution(std::ostream &out) {
  if(!fCompMesh) {
   std::cout << "TTimeAnalysis::PrintSolution hasn't mesh associated.\n";
    return;
  }
  if(LastTimeReached()==fEndTime[fNEndTimes-1]) {
    out <<std::endl << "Time reached = " << fEndTime <<std::endl << "Connect Solution :" <<std::endl;
    fCompMesh->ConnectSolution(out);
    return;
  }
  out <<std::endl << "End time is not reached" <<std::endl;
  out << "current time = " << LastTimeReached() <<std::endl;
  fCompMesh->ConnectSolution(out);
}

double TTimeAnalysis::StabilityCondition() {
  REAL area = fLaw->MinDeltaX();
//  if(fLaw->Dimension()==2) area = sqrt(area);
  REAL maxeig = fLaw->MaxEigen();
	if(maxeig < 0.) PZError << "TTimeAnalysis::StabilityCondition. Bad maxeig.\n";
  if(IsZero(maxeig*.1)) return MINIMETIMESTEP;
  if(maxeig > BIGMAXEIG) return fEndTime[fNEndTimes-1];
  return (fCFL*area)/maxeig;
}

void TTimeAnalysis::ApplyUZero() {
  if(!fUZero) {
    PZError << "TTimeAnalysis::ApplyUZero. Function UZero() is undefined.\n";
    return;
  }
  int i, dim, index, nelem = fCompMesh->NMaterials();
  std::map<int, TPZMaterial * >::const_iterator mit;
  for (mit = fCompMesh->MaterialVec().begin(); mit != fCompMesh->MaterialVec().end(); mit++) {
	  TPZMaterial *mat = mit->second;
	    if(mat->Id()>-1) {
			dim = mat->Dimension();
			break;
		}
  }
	TPZFMatrix<REAL> InitialSol;

  InitialSol.Redim(Mesh()->NEquations(), 1);
  InitialSol.Zero();
  for (int iel = 0; iel < Mesh()->NElements(); iel++) {
	  TPZCompEl * cel = Mesh()->ElementVec()[iel];
	  if (!cel) continue;
	  TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cel);
	  if (!disc) continue;
	  if (disc->NConnects() == 0) continue;
	  int bl = disc->Connect(0).SequenceNumber();
	  int blpos = Mesh()->Block().Position(bl);
	  int blocksize = Mesh()->Block().Size(bl);

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

  SetInitialSolution(InitialSol);

  ApplyBCDirichlet();
}

  /*
  int elord;
  TPZBlock<STATE> &block = fCompMesh->Block();
  block.SetMatrix(&fCompMesh->Solution());
  block.Matrix()->Zero();
  nelem = fCompMesh->NElements();
  TPZInterpolatedElement *cel;
  TPZGeoEl *gel;
  int jn, ncon, dfseq, dfvar = fLaw->NStateVariables();
  TPZVec<REAL> u(dfvar,0.);
  TPZVec<int> MaskConnects(fCompMesh->NConnects(),0);
  REAL maxeig;
  TPZVec<REAL> normal(2,0.);
  normal[0] = 1.;
  
  for(i=0;i<nelem;i++) {
    TPZCompEl *el = fCompMesh->ElementVec()[i];
    if(!el || el->IsInterface()) continue;
    cel = (TPZInterpolatedElement *)el;
    gel = cel->Reference();
    dfvar = cel->Material()->NStateVariables();
    u.Resize(dfvar,0.);
    elord = cel->GetPreferredOrder();
    if(!elord) {
		continue;   // ???
//      if(!cel->CanBeDiscontinuous() || cel->NShapeF()!=1) {
//        PZError << "TTimeAnalysis::ApplyUZero.Bad continuity or shape function.\n";
//        continue;
//      }
//      cel->SetDiffusionCoefficient(0.);
      
      REAL area = gel->Volume();
      int geldim = gel->Dimension();
      if(area==0.) {
        if(!cel->Dimension())
          area = 1.;
        else
				  area = gel->Volume();
      }
      ncon = cel->NConnects();
      dfseq = cel->Connect(ncon-1).SequenceNumber();
      if(dfvar!=block.Size(dfseq)) PZError << "TTimeAnalysis::ApplyUZero. Bad size.\n";
      Integral(cel,fUZero,u);
      for(jn=0;jn<dfvar;jn++) {
        u[jn] /= area;
        block(dfseq,0,jn,0) = u[jn];
      }

      //** To compute the maxime eigenvalue of jacobian to CFL 
      // **Center of geometrical element
      TPZVec<REAL> center(geldim);
      if(!geldim) center.Resize(1,0.);
      TPZVec<REAL> x(3,0.);
	  gel->CenterPoint(gel->NSides() - 1, center);
      gel->X(center,x);
      fLaw->SetPoint(x);
      maxeig = fLaw->MaxEigJacob(u,normal);
      if(fLaw->MaxEigen() < maxeig) fLaw->SetMaxEigen(maxeig);
    }
    else {
	    int in, idisc = 0;
      ncon = cel->NCornerConnects();
      for(in=0;in<ncon;in++) {
//		  fSolution = fCompMesh->Solution();   // Apenas para identificar os connect no vetor solucao
        TPZVec<REAL> x(3,0.);
        for(jn=0;jn<3;jn++)
          x[jn]=gel->NodePtr(in)->Coord(jn);

        index = cel->ConnectIndex(in);
        fUZero(x,u,cel);

        //** To compute the maxime eigenvalue of jacobian to CFL 
        fLaw->SetPoint(x);
        maxeig = fLaw->MaxEigJacob(u,normal);
        if(fLaw->MaxEigen() < maxeig) fLaw->SetMaxEigen(maxeig);

        //** Fill value into adequated block, depend of the connect continuity 
        if(MaskConnects[index]) continue;

		//   if(cel->IsConnectContinuous(in)) {
     //     MaskConnects[index] = 1;
        //  dfseq = cel->Connect(in).SequenceNumber();
      //    if(dfvar!=block.Size(dfseq))
    //        PZError << "TTimeAnalysis::ApplyUZero. Bad size." <<std::endl;
  //        for(jn=0;jn<dfvar;jn++)
//            block(dfseq,0,jn,0) = u[jn];
			
   //     }
 //       else {
     
	 dfseq = cel->Connect(cel->Reference()->NSides()-1).SequenceNumber();
          for(jn=0;jn<dfvar;jn++)
            block(dfseq,0,idisc*dfvar+jn,0) = u[jn];
          idisc++;
       // }
      }
    }
  }*/

void TTimeAnalysis::Integral(TPZInterpolatedElement *el,
     void (*fp)(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel),TPZVec<REAL> &integral) {
  int i,dim = el->Dimension();
  TPZFMatrix<STATE> axes(3,3,0.);
  TPZFMatrix<STATE> jacobian(dim,dim);
  TPZFMatrix<STATE> jacinv(dim,dim);
  REAL detjac;
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  double weight = 0.;

  TPZIntPoints *rule;
  rule = &(el->GetIntegrationRule());

  int order = el->Material()->NStateVariables();
  TPZVec<REAL> u(order,0.);
  for(i=0;i<order;i++) integral[i] = 0.;

  for(i=0;i<rule->NPoints();++i) {
    rule->Point(i,intpoint,weight);
    el->Reference()->Jacobian(intpoint,jacobian,axes,detjac,jacinv);
    el->Reference()->X(intpoint, x);
    weight *= detjac;

    fp(x,u,el);
    for(i=0;i<order;i++) integral[i] += (u[i]*weight);
  }
}

/**Return maxime area and minime area using the pointer paramenter*/
REAL TTimeAnalysis::ComputeAreas(TPZCompMesh *cmesh,REAL *areamin) {
  int i,nelems = cmesh->NElements();
  REAL area, areamax = 0.;
  if(areamin) (*areamin) = 1.e12;
  TPZCompEl *el;
  int dim = fLaw->Dimension();

  for(i=0;i<nelems;i++) {
    el = fCompMesh->ElementVec()[i];
    if(!el || el->IsInterface()) continue;
    TPZInterpolatedElement *intel = (TPZInterpolatedElement *)el;
    if(intel->Dimension() != dim) continue;
    area = intel->Reference()->Volume();
    areamax = (area < areamax) ? areamax : area;
    (*areamin) = (area < (*areamin)) ? area : (*areamin);
  }
  return areamax;
}

void TTimeAnalysis::ApplyBC() {
  int i,j,k,l,nelem,nrows;
int RhsRows = fRhs.Rows();
  nelem = fCompMesh->NElements();
  TPZAdmChunkVector<TPZCompEl *> &elvec = fCompMesh->ElementVec();
  TPZCompEl *cel;
  TPZInterpolatedElement *el;
  for(i=0;i<nelem;i++) {
    cel = elvec[i];
    if(!cel || cel->IsInterface()) continue;
    int bcindex = cel->Reference()->MaterialId();
    el = (TPZInterpolatedElement *)cel;
    if(bcindex < 0) {
    TPZBndCond *bc = (TPZBndCond *)el->Material();
    int nvar = bc->Val1().Rows();
    int type = bc->Type();
    int ncon, icon, maskdisc = 0;
    switch(type) {
    case 1:
    case 2:
    case 7:
    case 8:
      {
        int seqnum, firsteq;
        if(!el->PreferredSideOrder(el->Reference()->NSides()-1)) {
				// para ordem zero -> volumes finitos
          TPZConnect &connect = el->Connect(el->Reference()->NSides());
          seqnum = connect.SequenceNumber();
          firsteq = fCompMesh->Block().Position(seqnum);
					nrows = nvar;
          for(k=firsteq;k<firsteq+nvar;k++) {
            for(l=1;l<nrows;l++) {
		if((k+l)<RhsRows) {
              fSolver->Matrix()->PutVal(k,l,0.);
              fSolver->Matrix()->PutVal(l,k,0.);
		}
            }
            fSolver->Matrix()->PutVal(k,k,1.);
            fRhs.PutVal(k,0,0.);
          }
        }
        else {
				// NOTA : O seguinte nao deve(ria) ser feito para connects continuos.
				// para interpolacao maior ou igual a um
          ncon = el->NCornerConnects();
				nrows = nvar * ncon;
          for(j=0;j<ncon;j++) {
//            if(el->IsConnectContinuous(j))
  //            icon = j;
    //        else {
              if(j==1 && !maskdisc) maskdisc = 0;
              else maskdisc = nvar;
              icon = el->Reference()->NSides();
      //      }
						// zera matriz tangente e vetor de carga, com diagonal um
            TPZConnect &connect = el->Connect(icon);
            seqnum = connect.SequenceNumber();
            firsteq = fCompMesh->Block().Position(seqnum);
            if(j==1) firsteq += maskdisc;
            for(k=firsteq;k<firsteq+nvar;k++) {
              for(l=1;l<nrows;l++) {
                if((k+l)<RhsRows) {
		fSolver->Matrix()->PutVal(k,l,0.);
                fSolver->Matrix()->PutVal(l,k,0.);
		}
              }
              fSolver->Matrix()->PutVal(k,k,1.);
              fRhs.PutVal(k,0,0.);
            }
          }
        }
        break;
      }
    case 0:
    case 4:
      {
        int seqnum, firsteq;
        if(!el->PreferredSideOrder(el->Reference()->NSides()-1)) {
          TPZConnect &connect = el->Connect(el->Reference()->NSides());
          seqnum = connect.SequenceNumber();
          firsteq = fCompMesh->Block().Position(seqnum);
					nrows = nvar;
          for(k=firsteq;k<firsteq+nvar;k++) {
            for(l=1;l<nrows;l++) {
		if((k+l)<RhsRows) {
              fSolver->Matrix()->PutVal(k,l,0.);
              fSolver->Matrix()->PutVal(l,k,0.);
		}
            }
            fSolver->Matrix()->PutVal(k,k,1.);
						fRhs.PutVal(k,0,0.);
          }
        }
        else {
          ncon = el->NCornerConnects();
					nrows = nvar * ncon;
          for(j=0;j<ncon;j++) {
        //    if(el->IsConnectContinuous(j))
         //     icon = j;
          //  else {
              if(j==1 && !maskdisc) maskdisc = 0;
              else maskdisc = nvar;
              icon = el->Reference()->NSides();
            //}
            TPZConnect &connect = el->Connect(icon);
            seqnum = connect.SequenceNumber();
            firsteq = fCompMesh->Block().Position(seqnum);
            if(j==1) firsteq += maskdisc;
            for(k=firsteq;k<firsteq+nvar;k++) {
              for(l=1;l<nrows;l++) {
		if((k+l)<RhsRows) {
                fSolver->Matrix()->PutVal(k,l,0.);
                fSolver->Matrix()->PutVal(l,k,0.);
		}
              }
              fSolver->Matrix()->PutVal(k,k,1.);
              fRhs.PutVal(k,0,0.);
            }
          }
        }
        break;
      }
//      case 1:   //Neumann -> Fluxo dado
      //case 2:
      //case 4:
      case 3:
      {
        ncon = el->Reference()->NSides();
        int nn = el->PreferredSideOrder(ncon-1);
        if(nn<0) nn = 0;
        if(!nn) {
          TPZCompElSide thisside(el,ncon-1);
          TPZStack<TPZCompElSide> elvec;
          thisside.EqualLevelElementList(elvec,0,0);
          TPZConnect &con = elvec[0].Element()->Connect(elvec[0].Element()->Reference()->NSides());
          int seqnumneigh = con.SequenceNumber();
          int firsteqneigh = fCompMesh->Block().Position(seqnumneigh);
          int seqnum = el->Connect(el->Reference()->NSides()).SequenceNumber();
          int firsteq = fCompMesh->Block().Position(seqnum);
					nrows = nvar * ncon;
          for(k=firsteq;k<firsteq+nvar;k++) {
            for(l=1;l<nrows;l++) {
		if((k+l)<RhsRows) {
              fSolver->Matrix()->PutVal(k,l,0.);
              fSolver->Matrix()->PutVal(l,k,0.);
		}
            }
            fSolver->Matrix()->PutVal(k,k,1.);
            fRhs.PutVal(k,0,fRhs.GetVal(firsteqneigh+k-firsteq,0));
          }
        }
        else {
					nrows = nvar * ncon;
          for(j=0;j<ncon;j++) {
            /** Find neighboard elements in the same level of the current element 
            TPZCompElSide thisside(el,j);
            TPZStack<TPZCompElSide> elvec(3);
            thisside.EqualLevelElementList(elvec,0,0);
            ((TPZInterpolatedElement *)elvec[0].Element())->MakeConnectContinuous(elvec[0].Side());
            ((TPZInterpolatedElement *)el)->MakeConnectDiscontinuous(j);*/
          }
        }
        break;
      }
      case 5:
        break;
      default:
        PZError << "TTimeAnalysis::ApplyBC. Boundary condition " << bcindex <<
          ", type " << type <<" is not implemented.\n";
        break;
      }
    }
  }
  for(k=0;k<fRhs.Rows();k++)
    if((*(fSolver->Matrix()))(k,k)==0.)
      fSolver->Matrix()->PutVal(k,k,1.);
}

/*void TTimeAnalysis::ApplyBCNonFlux() {
  int i, nelem;
  nelem = fCompMesh->NElements();
  TPZAdmChunkVector<TPZCompEl *> &elvec = fCompMesh->ElementVec();
  TPZCompEl *cel;
  TPZInterpolatedElement *el;
  for(i=0;i<nelem;i++) {
    cel = elvec[i];
    if(!cel || !cel->IsInterpolated()) continue;
    int bcindex = cel->Reference()->MaterialId();
    el = (TPZInterpolatedElement *)cel;
    if(!(bcindex < 0)) continue;
    TPZBndCond *bc = (TPZBndCond *)el->Material();
    int type = bc->Type();
    if(type == 4) {
//      if(!el->Dimension())
//        ((TInterfaceElement *)(el->Interface(0)))->NonContribution();
//      else 
//        ((TInterfaceElement *)(el->Interface(el->Reference()->NSides()-1)))->NonContribution();
      fCompMesh->ElementVec()[el->Index()] = 0;
      delete el;
    }
/*    else if(type==8) {
      int mask, nvar = el->Material()->NStateVariables();
      int seqnum, firsteq, ncon = el->NSides();
      for(j=0;j<ncon;j++) {
        mask = 0;
        if(!el->IsConnectContinuous(j)) {
          el->MakeConnectContinuous(j);
          mask = 1;
        }
        TPZConnect &connect = el->Connect(j);
        seqnum = connect.SequenceNumber();
        firsteq = fCompMesh->Block().Position(seqnum);
        for(int k=firsteq+1;k<firsteq+nvar-1;k++)
          fSolution->PutVal(k,0,(bc->Val2()(k-firsteq,0)));
        if(mask) el->MakeConnectDiscontinuous(j);
      }
    }
  }
//  fCompMesh->CleanInterfaces();
//  fCompMesh->InitializeBlock();
}
*/
void TTimeAnalysis::ApplyBCDirichlet() {
  int i,j,k,nelem;
  nelem = fCompMesh->NElements();
  TPZAdmChunkVector<TPZCompEl *> &elvec = fCompMesh->ElementVec();
  TPZCompEl *cel;
  TPZInterpolatedElement *el;
  for(i=0;i<nelem;i++) {
    cel = elvec[i];
    if(!cel || cel->IsInterface()) continue;
    int bcindex = cel->Reference()->MaterialId();
    el = (TPZInterpolatedElement *)cel;
    if(!(bcindex < 0)) continue;
    TPZBndCond *bc = (TPZBndCond *)el->Material();
    int nvar = bc->Val1().Rows(), nsides = el->Reference()->NSides();
    int type = bc->Type();
    int ncon, icon, maskdisc = 0;
    switch(type) {                                         // Only Dirichlet
    case 5:
      break;
    case 0: {
      int seqnum, firsteq;
      if(!el->PreferredSideOrder(el->Reference()->NSides()-1)) {
        TPZConnect &connect = el->Connect(nsides);
        seqnum = connect.SequenceNumber();
        firsteq = fCompMesh->Block().Position(seqnum);
        for(k=firsteq;k<firsteq+nvar;k++) {
          fCompMesh->Solution().PutVal(k,0,(bc->Val2()(k-firsteq,0)));
		  fSolution.PutVal(k,0,(bc->Val2()(k-firsteq,0)));
		}
      }
      else {
        ncon = el->NCornerConnects();
        for(j=0;j<ncon;j++) {
//          if(el->IsConnectContinuous(j)) icon = j;
  //        else {
            if(j==1 && !maskdisc) maskdisc = 0;
            else maskdisc = nvar;
            icon = el->Reference()->NSides();
    //      }
          TPZConnect &connect = el->Connect(icon);
          seqnum = connect.SequenceNumber();
          firsteq = fCompMesh->Block().Position(seqnum);
          if(j==1) firsteq += maskdisc;
          for(k=firsteq;k<firsteq+nvar;k++) {
            fCompMesh->Solution().PutVal(k,0,(bc->Val2()(k-firsteq,0)));
			fSolution.PutVal(k,0,(bc->Val2()(k-firsteq,0)));
		  }
        }
      }
      /** Melhorar para caso de ordem maior que um e para caso 1-D na fronteira*/
      break;
    }
    case 1:   //Neumann -> Fluxo dado
    case 2:
    case 3:
    case 4:
    case 6:
    case 7:
    case 8:
      break;
/*      int side = nsides-1;
      TPZGeoElSide gelside(el->Reference(),side);
      gelside = gelside.Neighbour();
      side = gelside.Side();
      TPZCompEl *compel = gelside.Reference().Element();
      if(!compel) {
      nsides = gelside.Element()->NSideSubElements(side);
      for(k=0;k<nsides;k++) {
        TPZCompEl *compel = gelside.Reference().Element();
        while(!compel) {
          gelside = 
      ncon = el->NSides();
      int nn = el->SideOrder(ncon-1);
      if(nn<0) nn = 0;
      if(!nn) {
        TPZCompElSide thisside(el,ncon-1);
        TPZStack<TPZCompElSide> elvec(2);
        thisside.EqualLevelElementList(elvec,0,0);
        TPZConnect &con = elvec[0].Element()->Connect(elvec[0].Element()->NSides());
        int seqnumneigh = con.SequenceNumber();
        int firsteqneigh = fCompMesh->Block().Position(seqnumneigh);
        int seqnum = el->Connect(el->NSides()).SequenceNumber();
        int firsteq = fCompMesh->Block().Position(seqnum);
        for(k=firsteq;k<firsteq+nvar;k++) {
          for(l=0;l<nrows;l++) {
            fStiffness->PutVal(k,l,0.);
            fStiffness->PutVal(l,k,0.);
          }
          fStiffness->PutVal(k,k,1.);
          fRhs->PutVal(k,0,fRhs->GetVal(firsteqneigh+k-firsteq,0));
        }
      }
      else {
        TPZCompEl *cel = el->Reference()->Neighbour(ncon-1).Element()->Reference();
        for(j=0;j<ncon;j++) {
          /** Find neighboard elements in the same level of the current element 
          TPZCompElSide thisside(el,j);
          TPZStack<TPZCompElSide> elvec(3);
          thisside.EqualLevelElementList(elvec,0,0);
          int number = elvec.NElements();
          for(k=0;k<number-1;k++)
            if(cel == elvec[k].Element()) break;
          ((TPZInterpolatedElement *)elvec[k].Element())->MakeConnectContinuous(elvec[k].Side());
          ((TPZInterpolatedElement *)el)->MakeConnectDiscontinuous(j);
        }
      }*/
      break;
    default:
      PZError << "TTimeAnalysis::ApplyBCDirichlet. Boundary condition " << bcindex <<
      ", type " << type <<" is not implemented.\n";
      break;
    }
  }
}

void TTimeAnalysis::ApplyRhsBC() {
  int i,k,nelem;
  nelem = fCompMesh->NElements();
  TPZAdmChunkVector<TPZCompEl *> &elvec = fCompMesh->ElementVec();
  TPZCompEl *el;
  TPZInterpolatedElement *cel;
  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el || el->IsInterface()) continue;
    cel = (TPZInterpolatedElement *)el;
    int bcindex = cel->Reference()->MaterialId();
    if(bcindex < 0) {
      TPZBndCond *bc = (TPZBndCond *)cel->Material();
      int nvar = bc->Val1().Rows();
      int type = bc->Type();
      int ncon, maskdisc = 0;
      switch(type) {
      case 0:                               // Dirichlet
        break;
      case 1:   //Neumann -> Fluxo dado
      case 2:
      case 3:
      case 5:
			case 6:
			case 7:
			case 8:
      case 4:
        break;
      {
        ncon = cel->Reference()->NSides();
        int nn = cel->PreferredSideOrder(ncon-1);
        if(nn<0) nn = 0;
        if(!nn) {
          TPZCompElSide thisside(cel,ncon-1);
          TPZStack<TPZCompElSide> elvec;
          thisside.EqualLevelElementList(elvec,0,0);
          TPZConnect &con = elvec[0].Element()->Connect(elvec[0].Element()->Reference()->NSides());
          int seqnumneigh = con.SequenceNumber();
          int firsteqneigh = fCompMesh->Block().Position(seqnumneigh);
          int seqnum = cel->Connect(ncon).SequenceNumber();
          int firsteq = fCompMesh->Block().Position(seqnum);
          for(k=firsteq;k<firsteq+nvar;k++)
            fRhs.PutVal(k,0,fRhs.GetVal(firsteqneigh+k-firsteq,0));
        }
        else {    // Cuidado com a matrix stiffness para o caso ApplyBC()
//          for(j=0;j<ncon;j++) {
            /** Find neighboard elements in the same level of the current element */
//            TPZCompElSide thisside(cel,j);
//            TPZStack<TPZCompElSide> elvec(3);
//            thisside.EqualLevelElementList(elvec,0,0);
//            ((TPZInterpolatedElement *)elvec[0].Element())->MakeConnectContinuous(elvec[0].Side());
//            cel->MakeConnectDiscontinuous(j);
//          }
        }
        break;
      }
      default:
        PZError << "TTimeAnalysis::ApplyRhsBC. Boundary condition " << bcindex <<
          ", type " << type <<" is not implemented.\n";
        break;
      }
    }
  }
}

void TTimeAnalysis::CleanToStartRun() {
  fLaw->Clean();
  fLaw->SetTime(fStartTime,fTimeStep);
  fTimeStep = fCFL = 0.;
  if(!fUk) fUk = new TPZFMatrix<STATE>(fSolution);
  else (*fUk) = fSolution;
}

void TTimeAnalysis::ReDefineGraphMesh(int dim) {
  int dim1 = dim-1;
  TPZGraphMesh *graphmesh = fGraphMesh[dim1];
  if(!dim || dim>3 || !graphmesh) {
    PZError << "TTimeAnalysis::ReDefineGraphMesh. Bad paramenter.\n";
    return;
  }

  TPZDrawStyle style = graphmesh->Style();
  TPZGraphMesh *newgraph;
/*  switch(style) {
  case 0:
    newgraph = new TPZDXGraphMesh(fCompMesh,dim,(TPZDXGraphMesh *)graphmesh);
    break;
  case 1:
    newgraph = new TPZMVGraphMesh(fCompMesh,dim,(TPZMVGraphMesh *)graphmesh);
    break;
  case 2:
    newgraph = new TPZV3DGraphMesh(fCompMesh,dim,(TPZV3DGraphMesh *)graphmesh,fCompMesh);
    break;
  case EVTKStyle:
	  newgraph = new TPZVTKGraphMesh(fCompMesh, dim, (TPZV3DGraphMesh *)graphmesh);
	  break;
  default:
    PZError << "TTimeAnalysis::ReDefineGraphMesh. GraphMesh style is undefined.\n";
    delete graphmesh;
    fGraphMesh[dim1] = 0;
    return;
  }
  */
  newgraph->SetResolution(graphmesh->Res());

  delete fGraphMesh[dim1];
  fGraphMesh[dim1] = newgraph;
}

void TTimeAnalysis::ReadData(std::ifstream &input) {
	GetCommentary(input,1);
	int limiter;
	double diff;
  GetDataCommented(input,limiter);
	GetDataCommented(input,fSteadyState);
  GetDataCommented(input,fCFL);
  GetDataCommented(input,diff);
  GetDataCommented(input,fRhsAdaptNorm);
  GetDataCommented(input,fNPlots);
}

/*
void TTimeAnalysis::DefineGraphMesh(int dim,TPZVec<char *> &scalnames,
			       TPZVec<char *> &vecnames,TPZGraphMesh *graph) {
  int dim1 = dim-1;
  TPZDrawStyle style = graph->Style();

  if(fGraphMesh[dim1]) delete fGraphMesh[dim1];
  fScalarNames[dim1] = scalnames;
  fVectorNames[dim1] = vecnames;

  switch(style) {
  case 0:
    fGraphMesh[dim1] = new TPZDXGraphMesh(fCompMesh,dim,(TPZDXGraphMesh *)graph,fLaw);
    break;
  case 1:
    fGraphMesh[dim1] = new TPZMVGraphMesh(fCompMesh,dim,(TPZMVGraphMesh *)graph,fLaw);
    break;
  case 2:
    fGraphMesh[dim1] = new TPZV3DGraphMesh(fCompMesh,dim,(TPZV3DGraphMesh *)graph,fLaw);
    break;
  default:
   std::cout << "grafgrid was not created\n";
    fGraphMesh[dim1] = 0;
    break;
  }
}
*/
void TTimeAnalysis::SetInitialSolution(TPZFMatrix<STATE> & InitialSol) {
	const int nrows = Mesh()->Solution().Rows();
	const int ncols = Mesh()->Solution().Cols();
	if ((InitialSol.Rows() != nrows) || (InitialSol.Cols() != ncols)) {
		PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
	}
	else {
		fSolution = InitialSol;
	}
}

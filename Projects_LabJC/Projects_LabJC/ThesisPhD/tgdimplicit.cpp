/*******       File : tarungek.c

This file contains the method definitions for class TGDImplicit.

*******              *******/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "pzgraphmesh.h"

#include "tgdimplicit.h"
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

#include "pzstrmatrix.h"

#include "pzelmat.h"
#include "myheader.h"
#include "hadaptive.h"

/*******       TGDImplicit, class derived from TPZAnalysis       *******/
TGDImplicit::TGDImplicit(std::istream &input, std::ostream &out,TPZCompMesh *mesh,int level)
     : TTimeAnalysis(input,out,mesh,level) {
	int order;
  GetDataCommented(input,order);
  int i, nelem;
  std::map<int, TPZMaterial * >::const_iterator mit;
  for (mit = mesh->MaterialVec().begin(); mit != mesh->MaterialVec().end(); mit++) {
	  TPZMaterial *mat = mit->second;
	  if (mat->Id() < 0) continue;
    fLaw = (TConservationLaw *)mat;
  }
  GetDataCommented(input,fStartTime);
  GetDataCommented(input,fNEndTimes);
	fEndTime = new REAL[fNEndTimes];
	for(i=0;i<fNEndTimes;i++) input >> fEndTime[i];
  GetDataCommented(input,MINIMETIMESTEP);
  GetDataCommented(input,MAXROWS);
  /**Adaptive or non-adaptive scheme*/
  GetDataCommented(input,i);        // i=0 nonadaptive, i!=0 adaptive
	if(i) {
		/** Criando objeto para adaptatividade h-p ou do fluxo numerico*/
		fAdaptive = new TAdaptive(input,fLaw->Dimension(),1);
		fAdaptive->MinLevel = level;
	}
	else {
		GetCommentary(input,14);
		fAdaptive = 0;
	}

  fCFL = 1.;
  fCFLDiffussion = 0.;
  AppendMethodName(Name());
  fUZero = 0;
  fLevel = level+1;
#ifdef PARALLELVERSION
  fRank = 0;
  fSize = 1;
#endif
}

void TGDImplicit::AdjustBeforeAssemble(REAL time) {
	TTimeAnalysis::AdjustBeforeAssemble(time);
	((TConservationLaw *)fLaw)->SetFactor(fCFLDiffussion);
}

void TGDImplicit::CleanDiffusion() {
//	SetFactor(fCompMesh->ElementVec(),0.);
	fCFLDiffussion = 0.;
}

int TGDImplicit::Adapting(int &stepgraph,int onlyflux) {
	/** All finite elements receive a especified difusion coefficient */
//  SetFactor(fCompMesh->ElementVec(),fCFLDiffussion);
	((TConservationLaw *)fLaw)->SetFactor(fCFLDiffussion);

	return TTimeAnalysis::Adapting(stepgraph,onlyflux);
}

void TGDImplicit::CleanToStartRun() {
	TTimeAnalysis::CleanToStartRun();
	fCFLDiffussion = 0.;
}

void TGDImplicit::ReadData(std::ifstream &input) {
	GetCommentary(input,1);
	int limiter;
  GetDataCommented(input,limiter);
	GetDataCommented(input,fSteadyState);
  GetDataCommented(input,fCFL);
  GetDataCommented(input,fCFLDiffussion);
  GetDataCommented(input,fRhsAdaptNorm);
  GetDataCommented(input,fNPlots);
}

void TGDImplicit::GetSchemeType(char *filename) {
  /**We can to choose implicit or explicit scheme to diffusion*/
  char aux[16];
  int index=0;
  aux[index++] = 'G'; aux[index++] = 'D';
  std::cout << "Galerkin Discontinuous - Implicit scheme.\n";
  aux[index++] = 'i'; aux[index++] = 'm';
  aux[index++] = '\0';
  strncat(filename,aux,index);
}

void TGDImplicit::Run(std::istream &input, std::ostream &out) {
  double ctime = 0.;
  int numitera=0, niterations, refine = 1;
  std::string lawname;
  lawname = plotfile;
  if(!fCompMesh || !fLaw) {
    PZError << "TTimeRungeKutta::Run, hasn't computational mesh or equation associated.\n";
    return;
  }

  int l, k, count, stepgraph = 0;
  int onlyrhs = 0;
  /** Put initial time in fTimes and into conservation law*/
  ctime = fStartTime;
  REAL minarea;

  /** Initializating the solution at initial time */
  for(k=0;k<3;k++) 
		if(fGraphMesh[k]) {
			fGraphMesh[k]->SetNames(fScalarNames[k],fVectorNames[k]);
//			fGraphMesh[k]->DrawMesh(stepgraph+300);
		}

  ComputeAreas(fCompMesh,&minarea);
  fLaw->SetMinDeltaX(minarea);

  AdequateMatrix();
  ApplyUZero();   // Discretization of UZero
  fSolution = fCompMesh->Solution();
//  ApplyBCDirichlet();

  /** Adapting the numerical flux: high aprox. into discontinuous region
     and lower aprox. into the regular regions*/
	if(fAdaptive) fAdaptive->SteadyState(fSteadyState);
  Adapting(stepgraph,0);

  /**Here: We can to change the mesh and acerting stiff matrix, solver and fUk*/
  AdequateMatrix();
  /**To temporary storage of L(u0)+L(u1)+L(u2)+...*/
  TPZFMatrix<STATE> RhsSum;
	int end = 0;
	do {
  	/** Numero de plots para o arquivo de post-processo */
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

	  /** Iterating while current time is lower than end time or it achieved maxime iterations*/
  	while(ctime < fEndTime[end] && numitera < MAXITERATIONS) {
    	/**Incrementa numero de iteracoes*/
	    numitera++;

  	  /**Compute time step from CFL condition (and adjust next time at end time)
    	   and compute maxime area of the */
	    ComputeNextTime(ctime,fEndTime[end]);

	    /**Adapting the numerical flux: high aprox. into discontinuous region
  	     and lower aprox. into the regular regions*/
    	AdjustBeforeAssemble(ctime);
	    /** Drawing to Post-processing */
  	  if(ctime>=TimePost[count]) {
	      count++;
  	    ReDraw(stepgraph++);
	    }

	    /**Fixing the solution at time tk and adjusting solver*/
  	  fSolutionNext = fCompMesh->Solution();

      // construindo matrix com igual dimensao que fRhs
      TPZFMatrix<STATE> RhsOnlyIntegral(fRhs);

	  /** Prints only by tests */
//	  ofstream outmat("matrizes.txt",ios::app);
//	  fSolver->Matrix()->Print("fMat_1",outmat);
//	  fUk->Print("fUk_1",outmat);
//	  fRhs.Print("fRhs_1",outmat);
//	  outmat.flush();
			// Calcula o termo produto interno das funcoes teste com as variaveis de estado
      AssembleOnlyIntegral(RhsOnlyIntegral);
//	  RhsOnlyIntegral.Print("RhsOnlyInt",outmat);
//	  outmat.flush();
      /**Iterating over the order of the Runge-Kutta solver*/
			// seta o coeficiente para o metodo Runge-Kutta (Euler explicito)
      fLaw->SetCoef(.5);
			//Parametro Reassembling vira 1 se acontece um gradiente muito grande, caso a malha seja modificada
      if(ReAssembling && TimePost[count]!=fEndTime[end]) ReAssemble(ctime);
      else Assemble();
			
	  /** Prints only by tests */
//	  fRhs.Print("fRhs after assemble",outmat);
//	  outmat.flush();

      /**Computing the sumatory of operators c(l,0)*[H(u0)+H(u1)+...] */
      fRhs *= fTimeStep;
      fRhs += RhsOnlyIntegral;

      /**Solving K(u(l+1)-u(k)) = H and finding u(l+1)*/
      ApplyBC();
//	  fSolver->Matrix()->Print("fStiff after ApplyBC",outmat);
//	  fRhs.Print("fRhs after ApplyBC",outmat);
//	  outmat.flush();
//	  outmat.close();
	  fPreviousSol = fSolution;
	  fSolution.Zero();
	  Solve();   //SolveToImplicit();  // Solve();
	  /** Prints only by tests */
//	  fSolver->Matrix()->Print("fStiff after solve",outmat);
//	  fRhs.Print("fRhs after solve",outmat);
//	  outmat.close();
	  
	  if(fSteadyState) {
		  fPreviousSol -= fSolution;
		  fRhsNorm = Norm(fPreviousSol);
	  }

	  /**Adapting the numerical flux: high aprox. into discontinuous region
		and lower aprox. into the regular regions and hp-adaptivity*/
	  refine = Adapting(stepgraph,0);
	  if(refine) {
		  numitera++;
		  onlyrhs = 0;
	  }
	  /** Imprime resultados nao para todo passo de tempo. Pode depender de parametro */
	  if(!(numitera%30)) {
		  out << "t: " << CurrentTime() << " eq: " << fCompMesh->NEquations();
		  if(fSteadyState) out << " fRhsNorm = " << fRhsNorm;
		  out << std::endl;
		  CountElements(*fCompMesh);
	  }
	  
	  /** Acerting current time  and the CFL value*/
	  ctime += fTimeStep;
	}

		/** Asking to continue process, can to put new end time */
		out << "\nLast End Time = " << fEndTime[end];
		end++;

  	/**Imprime estado atual da malha geometrica, computacional, solucao, etc*/
	  AdjustBeforeAssemble(ctime);
		CountElements(*fCompMesh);
	  /**To post-processing at last time -> (ending=1) */
  	ReDraw(stepgraph,1);
	  for(k=0;k<3;k++)
			if(fGraphMesh[k] && end < fNEndTimes) {
				PlotFileName(lawname,end,plotfile);
				fGraphMesh[k]->SetFileName(plotfile);
  	  }

	}while(!IsZero(ctime-fEndTime[fNEndTimes-1]));
}

void TGDImplicit::AssembleOnlyIntegral(TPZFMatrix<STATE> &rhs) {
  rhs.Zero();
  int iel;
  int nelem = fCompMesh->NElements();
  TPZElementMatrix ek,ef;
  TPZManVector<int64_t> destinationindex(100);
  TPZManVector<int64_t> sourceindex(100);
  TPZBlock<STATE> &Block = fCompMesh->Block();

  fCompMesh->Block().SetMatrix(&fSolutionNext);
  for(iel=0; iel < nelem; iel++) {
	  TPZInterpolatedElement *el = (TPZInterpolatedElement *)fCompMesh->ElementVec()[iel];
    if(!el || el->IsInterface()) continue;
    el->CalcIntegral(ef);
    if(!ef.fMat) continue;

    if(!el->HasDependency()) {
      destinationindex.Resize(ef.fMat.Rows());
      int destindex = 0;
      int numnod = ef.NConnects();
      for(int in=0; in<numnod; in++) {
         int npindex = ef.ConnectIndex(in);
         TPZConnect &np = fCompMesh->ConnectVec()[npindex];
         int blocknumber = np.SequenceNumber();
         int firsteq = Block.Position(blocknumber);
         int ndf = Block.Size(blocknumber);
         for(int idf=0; idf<ndf; idf++) {
           destinationindex[destindex++] = firsteq+idf;
         }
      }
      rhs.AddFel(ef.fMat,destinationindex);
    }
    else {
      // the element has dependent nodes
		el->InitializeElementMatrix(ek, ef);
      int destindex = 0;
      int fullmatindex = 0;
      destinationindex.Resize(ef.fConstrMat.Rows());
      sourceindex.Resize(ef.fConstrMat.Rows());
      int numnod = ef.fConstrConnect.NElements();
      for(int in=0; in<numnod; in++) {
         int npindex = ef.fConstrConnect[in];
         TPZConnect &np = fCompMesh->ConnectVec()[npindex];
         int blocknumber = np.SequenceNumber();
         int firsteq = Block.Position(blocknumber);
      	int ndf = Block.Size(blocknumber);
         if(np.HasDependency()) {
           fullmatindex += ndf;
           continue;
         }
         for(int idf=0; idf<ndf; idf++) {
           sourceindex[destindex] = fullmatindex++;
           destinationindex[destindex++] = firsteq+idf;
         }
      }
      sourceindex.Resize(destindex);
      destinationindex.Resize(destindex);
      rhs.AddFel(ef.fConstrMat,sourceindex,destinationindex);
    }
  }
  fCompMesh->Block().SetMatrix(&fCompMesh->Solution());
}

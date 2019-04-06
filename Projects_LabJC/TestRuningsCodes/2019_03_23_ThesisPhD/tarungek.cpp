/*******       File : tarungek.c

This file contains the method definitions for class TTimeRungeKutta.

*******              *******/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "tarungek.h"
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

/*******       TTimeRungeKutta, class derived from TPZAnalysis       *******/
TTimeRungeKutta::TTimeRungeKutta(std::istream &input, std::ostream &out,TPZCompMesh *mesh,int level)
     : TTimeAnalysis(input,out,mesh,level) {
  GetDataCommented(input,fOrder);
  if(fOrder < 1  || fOrder > 4) {
	  std::cout << "TTimeRungeKutta constructor. Bad parameter fOrder. fOrder = 3.\n";
    fOrder = 1;
  }
  int i, nelem;
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
  std::map<int, TPZMaterial * >::const_iterator mit;
  for (mit = mesh->MaterialVec().begin(); mit != mesh->MaterialVec().end(); mit++) {
	  TPZMaterial *mat = mit->second;
	  if (mat->Id() < 0) continue;
	  fLaw = (TConservationLaw *)mat;
	  break;
  }

//    if(mesh->MaterialVec()[i]->Id() < 0) continue;
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

  fUseLimiter = 0;
}

TTimeRungeKutta::TTimeRungeKutta(std::istream &input, std::ostream &out,TPZCompMesh *cmesh,REAL TEnd) :
    TTimeAnalysis(input,out,cmesh,TEnd) {

  int i;

  GetDataCommented(input,fName,2);
	GetDataCommented(input,fOrder);
  if(fOrder < 1  || fOrder > 3) {
	  std::cout << "TTimeRungeKutta constructor. Bad parameter fOrder. fOrder = 3.\n";
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
  /**Adaptive or non-adaptive scheme */
  GetDataCommented(input,i);        // i=0 nonadaptive, i!=0 adaptive
	if(i) {
		/** Criando objeto para adaptatividade h-p ou do fluxo numerico */
		fAdaptive = new TAdaptive(input,fLaw->Dimension(),1);
		fAdaptive->MinLevel = 1;
	}
	else fAdaptive = 0;
  if(fEndTime[0]<TEnd) {
    fStartTime = fEndTime[0];
    fEndTime[0] = TEnd;
  }
  GetDataCommented(input,fCFL);
	double diff;
  GetDataCommented(input,diff);
  GetDataCommented(input,fRhsAdaptNorm);
  AppendMethodName(Name());
  fUZero = 0;
}

TTimeRungeKutta::TTimeRungeKutta(std::ostream &out) : TTimeAnalysis(out) {
  fLaw = 0;
  fOrder = 1;
	fRhsAdaptNorm = 0.75;
  int i;
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
  fStartTime = 0.;
	fNEndTimes = 1;
	fEndTime = new REAL(1.);
  MAXROWS = 700;
  fUk = 0;
  fUZero = 0;
	fAdaptive = 0;
//	fCFLDiffussion = 0.;
  MINIMETIMESTEP =  1.e-7;
}

void TTimeRungeKutta::GetSchemeType(char *filename) {
  /**We can to choose implicit or explicit scheme to diffusion*/
  char aux[16];
  int index=0;
  aux[index++] = 'R'; aux[index++] = 'K';
  aux[index++] = Itoa(fOrder);
 std::cout << "RKGD - Celular Diffusion scheme.\n";
  aux[index++] = 'g'; aux[index++] = 'd';
  aux[index++] = '\0';
  strncat(filename,aux,index);
}

void TTimeRungeKutta::Run(std::istream &input, std::ostream &out) {
  double ctime = 0.;
  int numitera=0, refine = 1;
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

  	  /** Pode-se trabalhar com sub-malhas, assim sobre ellas pode dar-se
	        um valor diferente ao tipo de fluxo a utilizar.*/
    	/**Iterating over the order of the Runge-Kutta solver*/
  	  for(l=0;l<fOrder;l++) {
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

    	  /**Computing the sumatory of operators c(l,0)*[H(u0)+H(u1)+...] */
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
    	  /**Solving K(u(l+1)-u(k)) = H and finding u(l+1)*/
  	    Solve();
				fPreviousSol -= fSolution;
				fRhsNorm = Norm(fPreviousSol);
    	  /**It is necessary to adjust the time data for next intermediate computation*/
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
			if(fUseLimiter) {
				for(int nn=0;nn<fLaw->NStateVariables();nn++)
					CockburnLimiter(fCompMesh->ElementVec(),nn);
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

void TTimeRungeKutta::CockburnLimiter(TPZAdmChunkVector<TPZCompEl *> &elvec,int var) {
  if(fLaw->Dimension()==1) CockburnLimiter1d(elvec,var);
  else if(fLaw->Dimension()==2) CockburnLimiter2d(elvec,var);
  else PZError << "TTimeRungeKutta::CockburnLimiter. Bad dimension law.\n";
}

void TTimeRungeKutta::CockburnLimiter1d(TPZAdmChunkVector<TPZCompEl *> &elvec,int var) {
  int nelem = elvec.NElements();
  int i, j, elorder;
  TPZStack<TPZCompElSide> elvecneigh;
  TPZBlock<STATE> &block = fCompMesh->Block();
  block.SetMatrix(&(fCompMesh->Solution()));
  TPZCompEl *el;
  TPZConnect connect;
  TPZVec<REAL> qsi(2,0.);   // To triangular element: qsi(2,(1./3.));

  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el || el->IsInterface() || el->Dimension()!=1) continue;
    TPZInterpolatedElement *intel = (TPZInterpolatedElement *)el;
    elorder = intel->PreferredSideOrder(intel->Reference()->NSides()-1);
    if(elorder!=1) return;
    TPZInterpolatedElement *neighel;
    REAL mean, meanleft, meanright;
    int bound = 0;
    mean = LinearApproximation(qsi,intel,var);
    //Sera que tambem precisa: mean = intel->MeanSolution(var);  to meanleft e meanright?
    for(j=0;j<2;j++) {
      elvecneigh.Resize(0);
      TPZCompElSide thisside(intel,j);
      thisside.EqualLevelElementList(elvecneigh,0,0);
      /** Asking if exist neighbour and if neighbour has dimension one */
      if(!elvecneigh.NElements() || ((TPZInterpolatedElement *)elvecneigh[0].Element())->Dimension()!=1) {
        if(!j) {
          bound = 1;
          continue;
        }
        else {
          if(bound) PZError << "TTimeRungeKutta::CockburnLimiter1d. Element has no any neighboards.\n";
          meanright = meanleft;
          break;
        }
      }
      neighel = (TPZInterpolatedElement *)elvecneigh[0].Element();
      if(!j) {
        meanleft = mean - neighel->MeanSolution(var);
        continue;
      }
      meanright = neighel->MeanSolution(var) - mean;
      if(bound) meanleft = meanright;
    }
    int mask = 0;
//    if(intel->IsConnectContinuous(0)) connect = intel->Connect(0);
//    else {
      connect = intel->Connect(3);
      mask = 1;
 //   }
    int numseq = connect.SequenceNumber();
    REAL minmod = MinMod(mean - block(numseq,0,var,0), meanleft, meanright);
    block(numseq,0,var,0) = mean - minmod;
 //   if(intel->IsConnectContinuous(1)) connect = intel->Connect(1);
   // else {
      connect = intel->Connect(3);
      mask++;
   // }
    numseq = connect.SequenceNumber();
    int nvar = var;
    if(mask) nvar += intel->Material()->NStateVariables();
    block(numseq,0,nvar,0) = mean + minmod;
  }
}

void TTimeRungeKutta::CockburnLimiter2d(TPZAdmChunkVector<TPZCompEl *> &elvec,int var) {
}

REAL TTimeRungeKutta::MinMod(REAL first,REAL second,REAL third) {
  int mask = (first<0) ? -1 : 1;
  if(((second<0)?-1:1)!=mask || ((third<0)?-1:1)!=mask) return 0.;
  REAL vafirst = (first<0) ? -first : first;
  REAL vasecond = (second<0) ? -second : second;
  REAL vathird = (third<0) ? -third : third;
  if(vafirst>vasecond)
    vafirst = vasecond;
  return mask*((vafirst<vathird) ? vafirst : vathird);
}

void TTimeRungeKutta::ReadData(std::ifstream &input) {
	GetCommentary(input,1);
  GetDataCommented(input,fUseLimiter);
	GetDataCommented(input,fSteadyState);
  GetDataCommented(input,fCFL);
	double diff;
  GetDataCommented(input,diff);
  GetDataCommented(input,fRhsAdaptNorm);
  GetDataCommented(input,fNPlots);
}

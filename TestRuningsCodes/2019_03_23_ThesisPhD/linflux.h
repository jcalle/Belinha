
/*******   FILE :  LINFLUX.H   *******

Header file for class TLinearFlux, derived from class TNumericalFlux.

Classe usada no estudo de leis de conservacao lineares.
Dada uma lei de conservacacao linear Ut + A Ux + B Uy + C Uz = s
precissa-se apenas de fornecer os dados das matrizes A, B e C, em um
arquivo de dados. (p.e. "Amatrix.dat")
Os dados da matrix devem ser ingressados na seguinte ordem :
	Elementos da matrix A, B e C, por linha;
	Os autovalores de A, B e C, preferencialmente na ordem crescente;
	Os autovetores de A, B a C na matrix R, na mesma ordem dos autovetores;
	Elementos da matrix inversa de R, por linha;

*******/


#ifndef LINFLUXH
#define LINFLUXH

#include "numflux.h"

class TLinearLaw;

class TLinearFlux : public TNumericalFlux {

  REAL *fA;
  REAL *fValA;
  REAL *fEigVal;
  REAL *fEigVect;
  REAL *fInvEigVect;

  void (TLinearFlux::*fFlux)(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void FluxBackwardEuler(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void FluxLaxFriedrichs(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void FluxGodunov(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void FluxLaxWendroff(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void FluxRitchmyer(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void FluxBeamWarming(TPZVec<REAL> &U,TPZVec<REAL> &Up1p2,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void FluxGodunovLaxWendroff(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void FluxLaxFriedrichsLaxWendroff(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void FluxGodunovBeamWarming(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void FluxLaxFriedrichsBeamWarming(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void FluxSlopeLimiter(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);

 public:
  TLinearFlux(TLinearLaw *L,int dim);
  TLinearFlux(TLinearLaw *L,int t,int dim);
  TLinearFlux(TConservationLaw *L,REAL *A,REAL *EigVal,REAL *EigVect,REAL *InvEigVect,int t);
  TLinearFlux(TConservationLaw *L,REAL *A,int t);
  TLinearFlux(TConservationLaw *L,REAL *ValA);
  ~TLinearFlux() {
  }

  void NumericalFlux(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);

  void ShowFluxType();
  void SetFluxType(int t);
};


#endif



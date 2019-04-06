/*******       File : linlaw.h

Header file for class TLinearLaw1D derived from TConservationLaw to
linear conservation laws.

*******              *******/

#ifndef CONELAWH
#define CONELAWH

#include "conslaw.h"

template<class T>
class TPZVec;

class TConeLaw : public TConservationLaw {

 public:
	 virtual int NStateVariables() { return 1; }

  void Flux(TPZVec<REAL> &u,TPZVec<REAL> &flux);
  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &jacob);
  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &);
  void ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob);
  int ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta);
  void RoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &Roe);
  void EigRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &EigRoe);
  void ValRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &ValRoe);

  //Construtores e destrutor
  TConeLaw(std::istream &input);
  TConeLaw(int id);
  TConeLaw(int id,char *name,int type);
  ~TConeLaw() { }

  REAL MaxEigJacob(TPZVec<REAL> &u,TPZVec<REAL> &normal);
  REAL ValEigJacob(TPZVec<REAL> &u,int order,int dim=2);

  /** @brief Returns the integrable dimension of the material */
  virtual int Dimension() const { return 2; }

  TPZMaterial *NewMaterial();

  /**compute the value of the flux function to be used by ZZ error estimator*/
  void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux) { }

  int VariableIndex(char *name);
  int NSolutionVariables(int index);
  /**returns the solution associated with the var index based on the finite element approximation*/
  void Solution(TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
		TPZVec<REAL> &Solout);
  void VariablesName(TPZVec<char *> &names);
  void IncrementDiffusion(TPZVec<REAL> &,TPZFMatrix<STATE> &,TPZVec<REAL> &,TPZFMatrix<STATE> &,REAL ,
        REAL ,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,REAL );

  /**To evaluate true-error, L2-error and estimate error*/
  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<REAL> &flux,
        TPZVec<REAL> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val);

  int IdBC(double *x);
};

inline TPZMaterial *TConeLaw::NewMaterial() {
  TPZMaterial *newmat = new TConeLaw(Id());
  int fluxtype = fNumericalFlux.FluxType();
  ((TConservationLaw *)newmat)->SetNumericalFluxType(fluxtype,fluxtype);
  return newmat;
}

#endif

/*******       File : burger1d.h

Header file for class TBurgerLaw1D, to one-dimensional Burger conservation law.
It is used to test Riemann problems.

*******              *******/

#ifndef BURGERLAW1DH
#define BURGERLAW1DH

#include "conslaw.h"
template<class T>
class TPZVec;

class TBurgerLaw1D : public TConservationLaw {

  REAL fSonicPoint;

 public:
  TBurgerLaw1D(int id);
  TBurgerLaw1D(int id,char *name,int type=1);
  TBurgerLaw1D(std::istream &input);
  TBurgerLaw1D(TBurgerLaw1D &burger);
  ~TBurgerLaw1D() {
  }
  
  /** @brief Returns the integrable dimension of the material */
  virtual int Dimension() const { return 1; }

  virtual int NStateVariables() { return 1; }

  void Flux(TPZVec<REAL> &U,TPZVec<REAL> &flux);
  void JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob);
  void JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &Beta);
  void ValJacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &valjacob);
  void RoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &Roe);
  void EigRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &EigRoe);
  void ValRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &ValRoe);
  int ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta);
  REAL MaxEigJacob(TPZVec<REAL> &U,TPZVec<REAL> &normal);
  REAL ValEigJacob(TPZVec<REAL> &u,int order,int dim=1);

  REAL SonicPoint() { return fSonicPoint; }

  virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,
	  TPZVec<REAL> &flux,TPZVec<REAL> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val);

  /**compute the value of the flux function to be used by ZZ error estimator*/
  virtual void Flux(TPZVec<REAL> &x,TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
		    TPZVec<REAL> &flux);

  /**To create another material of the same type*/
  virtual TPZMaterial *NewMaterial();
  /**Read data of the material from a istream (file data)*/
  virtual void SetData(std::istream &data);

  /** To postprocessing */
  int VariableIndex(char *name);
  int NSolutionVariables(int index);
  /**returns the solution associated with the var index based on the finite element approximation*/
  void Solution(TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
		TPZVec<REAL> &Solout);
  void VariablesName(TPZVec<std::string> &names);

  /* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
  int IdBC(double *x);
};

inline TPZMaterial *TBurgerLaw1D::NewMaterial() {
  TPZMaterial *newmat = new TBurgerLaw1D(Id());
	int fluxtype = fNumericalFlux.FluxType();
  ((TConservationLaw *)newmat)->SetNumericalFluxType(fluxtype,fluxtype);
  return newmat;
}
inline void TBurgerLaw1D::SetData(std::istream &data) {
  TConservationLaw::SetData(data);
}

#endif

/*	friend void UZero1(TPZVec<REAL> &coord,TPZVec<REAL> &result);   // Riemann problem, shock
	friend void UZero2(TPZVec<REAL> &coord,TPZVec<REAL> &result);   // Riemann problem, rarefaction
	friend void UZero3(TPZVec<REAL> &coord,TPZVec<REAL> &result);   // shock - continuous initial function
	friend void UZero4(TPZVec<REAL> &coord,TPZVec<REAL> &result);   // equal sin function
	void SetUZero(int t=0);
*/

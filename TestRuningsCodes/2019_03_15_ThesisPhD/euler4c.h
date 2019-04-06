#ifndef EULERLAW4CHH
#define EULERLAW4CHH

#include "conslaw.h"


class TEulerLaw4C : public TConservationLaw {

  double fGamma;

 public :
  TEulerLaw4C(std::istream &input) : TConservationLaw(input) {
    fGamma=1.4;
    fNumericalFlux.SetOrder(NStateVariables());
    fNumericalFlux.SetDimension(Dimension());
  }
  TEulerLaw4C(int id) : TConservationLaw(id,"Euler4C-1D Conservation Law") {
    fGamma=1.4;
    fNumericalFlux.SetOrder(NStateVariables());
    fNumericalFlux.SetDimension(Dimension());
  }
  TEulerLaw4C(int id,char *name,int type) : TConservationLaw(id,name) {
    fGamma=1.4;
    fNumericalFlux.SetOrder(NStateVariables());
    SetNumericalFluxType(type,type);
    fNumericalFlux.SetDimension(Dimension());
  }
  TEulerLaw4C(TEulerLaw4C &law) : TConservationLaw(law) {
    fGamma = law.Gamma();
  }
  ~TEulerLaw4C() {
  }

  virtual int NStateVariables() { return 4; }

  void Flux(TPZVec<REAL> &u,TPZVec<REAL> &flux);
  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &jacob);
  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &);
  void ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob);
  int ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) { return 0; }
  REAL MaxEigJacob(TPZVec<REAL> &u,TPZVec<REAL> &normal);
  REAL ValEigJacob(TPZVec<REAL> &u,int order,int dim=1);

  void RoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &Roe);
  // EigRoe must to have Roe matrix, Roe eigenvalues, Roe matrix R and R inverse
  void EigRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &EigRoe);
  void ValRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &ValRoe);

  void Print(std::ostream &out= std::cout);
  /** @brief Returns the integrable dimension of the material */
  virtual int Dimension() const { return 1; }
  void FluxGodunov(double *Ui,double *Fluxi);
	
	REAL Pression(TPZVec<REAL> &U);
	
  TPZMaterial *NewMaterial();
  /**Read data of the material from a istream (file data)*/
  void SetData(std::istream &data);

  int VariableIndex(char *name);
  int NSolutionVariables(int index);
  /**returns the solution associated with the var index based on the finite element approximation*/
  void Solution(TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
		TPZVec<REAL> &Solout);
  void VariablesName(TPZVec<std::string> &names);

  REAL Gamma() { return fGamma; }

  /**compute the value of the flux function to be used by ZZ error estimator*/
  void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux) { }
};

inline TPZMaterial *TEulerLaw4C::NewMaterial() {
  TPZMaterial *newmat = new TEulerLaw4C(Id());
	int fluxtype = fNumericalFlux.FluxType();
  ((TConservationLaw *)newmat)->SetNumericalFluxType(fluxtype,fluxtype);
  return newmat;
}

#endif

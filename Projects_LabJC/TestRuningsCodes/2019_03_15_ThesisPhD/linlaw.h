/*******       File : linlaw.h

Header file for class TLinearLaw1D derived from TConservationLaw to
linear conservation laws.

*******              *******/

#ifndef LINLAWH
#define LINLAWH

#include "conslaw.h"

template<class T> 
class TPZFMatrix;
template<class T>
class TPZVec;

class TLinearLaw : public TConservationLaw {

 protected:
  int fOrder;                   // Order of the equations system

  REAL fA[MAXORDER*MAXORDER];   // A, B and/or C -> convection matrix
  REAL fValA[MAXORDER*MAXORDER];// |A| -> |A|=R|EigVal|(InvR)
  REAL fMaxEigVal;              // Valor do maximo autovalor da matrix A
  REAL fEigVal[MAXORDER];         // Vetor de autovalores
  REAL fEigVect[MAXORDER*MAXORDER];        // R -> matrix de autovetores
  REAL fInvEigVect[MAXORDER*MAXORDER];     // Inversa de R

  REAL fInverseQrt[MAXORDER*MAXORDER];
  
 public:
	 virtual int NStateVariables() { return fOrder; }

  REAL A(int i){ return fA[i]; }
  REAL ValA(int i) { return fValA[i]; }
  REAL EigVal(int i) { return fEigVal[i]; }
  REAL EigVect(int i) { return fEigVect[i]; }
  REAL InvEigVect(int i) { return fInvEigVect[i]; }
  double *A() { return fA; }
  double *ValA() { return fValA; }
  double *EigVal(){ return fEigVal; }
  double *EigVect() { return fEigVect; }
  double *InvEigVect() { return fInvEigVect; }

  void Flux(TPZVec<REAL> &u,TPZVec<REAL> &flux);
  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &jacob);
  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &);
  void ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob);
  int ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta);
  void RoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &Roe);
  void EigRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &EigRoe);
  void ValRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &ValRoe);
	/** To compute the inverse jacobian matrix */
	void InverseJacob(TPZFMatrix<STATE> &mat);

  //Construtores e destrutor
  TLinearLaw(std::istream &input);
  TLinearLaw(int id,int order);
  TLinearLaw(int id,int order,char *name,int type);
  ~TLinearLaw() { }

  virtual void SetMatrix(TPZFMatrix<STATE> &A,TPZFMatrix<STATE> &EigVect,TPZFMatrix<STATE> &InvEigVect,TPZVec<REAL> &EigVal);
  virtual void SetMatrix(std::istream &input);
  virtual void RequireMatrix();
  virtual void Print(std::ostream &out= std::cout);
  REAL MaxEigJacob(TPZVec<REAL> &u,TPZVec<REAL> &normal);
  REAL ValEigJacob(TPZVec<REAL> &u,int order,int dim=1);

  /** @brief Returns the integrable dimension of the material */
  virtual int Dimension() const { return 1; }

  TPZMaterial *NewMaterial();
  /**Read data of the material from a istream (file data)*/
  void SetData(std::istream &data);

  /**compute the value of the flux function to be used by ZZ error estimator*/
  void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<STATE> &axes, TPZVec<REAL> &flux) { }

  int VariableIndex(char *name);
  int NSolutionVariables(int index);
  /**returns the solution associated with the var index based on the finite element approximation*/
  void Solution(TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
		TPZVec<STATE> &Solout);
  void VariablesName(TPZVec<std::string> &names);
  /**Compute contribution to the stiffness matrix and right hand side at an integration point */
//  void IncrementDiffusion(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix &dsol,
//        REAL weight,REAL area,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);

  void InverseQrt(TPZVec<REAL> &sol,TPZFMatrix<STATE> &C);
  /**Compute contribution to the stiffness matrix and right hand side at an integration point
     void ContributeRhs(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ef);*/

  /**To evaluate true-error, L2-error and estimate error*/
  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<REAL> &flux,
        TPZVec<REAL> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val);

  void ApplyPeriodic(TPZCompMesh &cmesh);
  /* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
  int IdBC(double *x);
  /**Return number of variables to print into the postprocessing file*/
  int AllVariablesPost() { return 1; }
  /**Return number of variables to print into the postprocessing file*/
  int AllVariablesImplement() { return fOrder; }
};

inline TPZMaterial *TLinearLaw::NewMaterial() {
  TPZMaterial *newmat = new TLinearLaw(Id(),fOrder);
	int fluxtype = fNumericalFlux.FluxType();
  ((TConservationLaw *)newmat)->SetNumericalFluxType(fluxtype,fluxtype);
  return newmat;
}

#endif

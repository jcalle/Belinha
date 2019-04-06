#ifndef EULERBIHH
#define EULERBIHH

#include "conslaw.h"

class TEulerLaw2D : public TConservationLaw {

	double fGamma;
  // Valor final do dominio no eixo X
  double fXend;

   public :
	TEulerLaw2D(std::istream &input);
	TEulerLaw2D(int id,char *name,int type);
	TEulerLaw2D(int id);
	TEulerLaw2D(TEulerLaw2D &law);
	~TEulerLaw2D() {
	}

	virtual int NStateVariables() { return 4; }
  /**Return number of variables to print into the postprocessing file*/
  int AllVariablesPost() { return 9; }
  /**Return possible number of variables to print into the postprocessing file*/
  int AllVariablesImplemented() { return 9; }

  void Flux(TPZVec<REAL> &u,TPZVec<REAL> &flux);
  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &jacob);
  REAL MaxEigJacob(TPZVec<REAL> &u,TPZVec<REAL> &normal);
  REAL ValEigJacob(TPZVec<REAL> &, int, int);
	
  void InverseJacob(TPZFMatrix<STATE> &mat);
  int ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal);
  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal);

  void RoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &Roe);
  // EigRoe must to have Roe matrix, Roe eigenvalues, Roe matrix R and R inverse
  void EigRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &EigRoe);
  void ValRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &ValRoe);

  void Print(std::ostream &out= std::cout);
  /** @brief Returns the integrable dimension of the material */
  virtual int Dimension() const { return 2; }
  void FluxGodunov(double *Ui,double *Fluxi);
	
	REAL Pressure(TPZVec<REAL> &U);
	
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

	/** To compute the tau matrix to diffusive term */
	void Tau(TPZFMatrix<STATE> &jacinv,TPZFMatrix<STATE> &valjacob,TPZFMatrix<STATE> &tau);
	
  /** To determine the flux type for interface element when are subdivided the cells */
  int FluxType() { return fFluxTypeLower; }
  /* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
  int IdBC(double *x);

  /**compute the value of the flux function to be used by ZZ error estimator*/
  void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux) { }

  void ContributeOverInterface(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZVec<REAL> &up1,REAL weight,
     REAL area,int type,TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &phi0,TPZFMatrix<STATE> &phi1,TPZVec<REAL> &normal,
     TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

  void IncrementDiffusion(TPZVec<REAL> &x,TPZFMatrix<STATE> &,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,REAL weight,
        TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,REAL );

	/** Para calculo do Tau em matrices separadas */
	void MatrixDiff(TPZVec<REAL> &sol,TPZFMatrix<REAL> &axes, TPZFMatrix<STATE> &jacinv,TPZFMatrix<STATE>
	&ATauA,TPZFMatrix<STATE> &ATauB,TPZFMatrix<STATE> &BTauA,TPZFMatrix<STATE> &BTauB);
	void InvJacob2d(TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &jacinv);
	void JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &Ajacob,TPZFMatrix<STATE>
	&Bjacob);

};

inline TPZMaterial *TEulerLaw2D::NewMaterial() {
  TPZMaterial *newmat = new TEulerLaw2D(Id());
	int fluxtype = fNumericalFlux.FluxType();
  ((TConservationLaw *)newmat)->SetNumericalFluxType(fluxtype,fluxtype);
  return newmat;
}


#endif

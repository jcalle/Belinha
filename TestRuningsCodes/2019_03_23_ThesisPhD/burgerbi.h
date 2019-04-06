#ifndef BURGERLAWBIH
#define BURGERLAWBIH

#include "conslaw.h"

template<class T>
class TPZVec;
class TCompMesh;

class TBurgerLaw2D : public TConservationLaw {
  
 public:
  TBurgerLaw2D(int id);
  TBurgerLaw2D(std::istream &input);
  TBurgerLaw2D(int id,char *file,int type=1);
  TBurgerLaw2D(TBurgerLaw2D &burger);
  ~TBurgerLaw2D() { }

  /** @brief Returns the integrable dimension of the material */
  virtual int Dimension() const { return 2; }
  int VariableIndex(char *name);
  int NSolutionVariables(int var);
  virtual int NStateVariables() { return 1; }
  void VariablesName(TPZVec<std::string> &names);

  void Flux(TPZVec<REAL> &U,TPZVec<REAL> &funcao);
  void JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob);
  void JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &);
  void ValJacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &valjacob);
  int ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta);
  REAL MaxEigJacob(TPZVec<REAL> &U,TPZVec<REAL> &normal);
  REAL ValEigJacob(TPZVec<REAL> &u,int order,int dim);

  virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,
	  TPZVec<REAL> &flux,TPZVec<REAL> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val);
  
  /**compute the value of the flux function to be used by ZZ error estimator*/
  virtual void Flux(TPZVec<REAL> &x,TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
		    TPZVec<REAL> &flux);

  /**returns the solution associated with the var index based on the finite element approximation*/
  void Solution(TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
		TPZVec<REAL> &Solout);

  /**To create another material of the same type*/
  TPZMaterial *NewMaterial();
  /**Read data of the material from a istream (file data)*/
  void SetData(std::istream &data);

  /* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
  int IdBC(double *x) { return 0; }
};

inline int TBurgerLaw2D::VariableIndex(char *name) {
//  if(!strcmp(name,"velocity")) return 0;
  return 0;
}

inline int TBurgerLaw2D::NSolutionVariables(int index) {
  if(index==0) return 1;
  return TPZMaterial::NSolutionVariables(index);
}
inline TPZMaterial *TBurgerLaw2D::NewMaterial() {
  return new TBurgerLaw2D(*this);
}
inline void TBurgerLaw2D::SetData(std::istream &data) {
  TConservationLaw::SetData(data);
}


#endif

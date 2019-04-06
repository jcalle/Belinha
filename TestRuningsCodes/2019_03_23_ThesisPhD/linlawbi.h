/*******       File : LinLaw2.h       *******/

#ifndef LINLAWBIH
#define LINLAWBIH

#include "linlaw.h"

class TLinearLaw2D : public TLinearLaw {
 protected:
  REAL fMaxEigValB;   // Valor do maximo autovalor da matriz B

 public:
  REAL B(int i) { return fA[fOrder*fOrder+i]; }
  REAL ValB(int i) { return fValA[fOrder*fOrder+i]; }
  REAL EigValB(int i) { return fEigVal[fOrder+i]; }
  REAL EigVectB(int i) { return fEigVect[fOrder*fOrder+i]; }
  REAL InvEigVectB(int i) { return fInvEigVect[fOrder*fOrder+i]; }
  double *B() { return fA+fOrder*fOrder; }
  double *ValB() { return fValA+fOrder*fOrder; }
  double *EigValB(){ return fEigVal+fOrder; }
  double *EigVectB() { return fEigVect+fOrder*fOrder; }
  double *InvEigVectB() { return fInvEigVect+fOrder*fOrder; }

  void Flux(TPZVec<REAL> &u,TPZVec<REAL> &flux);
  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &jacob);
	void JacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal);
  void ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob);
  int ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta);

  TLinearLaw2D(std::istream &input);
  TLinearLaw2D(int id,int ord);
  TLinearLaw2D(int id,int order,char *name,int type);
  ~TLinearLaw2D() { }

  void SetMatrix(TPZFMatrix<STATE> &A,TPZFMatrix<STATE> &EigVect,TPZFMatrix<STATE> &InvEigVect,TPZVec<REAL> &EigVal);
  void SetMatrix(std::istream &input);
  void RequireMatrix();
  void Print(std::ostream &out= std::cout);
  REAL MaxEigJacob(TPZVec<REAL> &u,TPZVec<REAL> &normal);
	REAL ValEigJacob(TPZVec<REAL> &u,int order,int dim);
	/** @brief Returns the integrable dimension of the material */
	virtual int Dimension() const { return 2; }

  TPZMaterial *NewMaterial();
  /**Read data of the material from a istream (file data)*/
  void SetData(std::istream &data);
  
  /* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
  int IdBC(double *x);
};

inline TPZMaterial *TLinearLaw2D::NewMaterial() {
  TPZMaterial *newmat = new TLinearLaw2D(Id(),fOrder);
	int fluxtype = fNumericalFlux.FluxType();
  ((TConservationLaw *)newmat)->SetNumericalFluxType(fluxtype,fluxtype);
  return newmat;
}

class TLinearLaw2DCircle : public TLinearLaw2D {
  REAL   fAngle;   // Angulo de giro solicitado pelo usuario
 public:
  TLinearLaw2DCircle(std::istream &input);
  TLinearLaw2DCircle(int id);
  TLinearLaw2DCircle(int id,REAL angle,char *name,int type);
  ~TLinearLaw2DCircle() { }

  void SetData(std::istream &input);
  void SetData(REAL angle);
  
  TPZMaterial *NewMaterial();
};

inline TPZMaterial *TLinearLaw2DCircle::NewMaterial() {
  TPZMaterial *newmat = new TLinearLaw2DCircle(Id());
	int fluxtype = fNumericalFlux.FluxType();
  ((TConservationLaw *)newmat)->SetNumericalFluxType(fluxtype,fluxtype);
  return newmat;
}

#endif


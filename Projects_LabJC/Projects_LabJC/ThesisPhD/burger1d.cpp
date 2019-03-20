/*******       File : NLinLaw1.c       *******/

#include "burger1d.h"
#include "pzvec.h"
#include "pzfmatrix.h"

/*******       TBurgerLaw1D       *******/

TBurgerLaw1D::TBurgerLaw1D(int id) : TConservationLaw(id,"Burger Conservation Law 1D") {
  fSonicPoint = 0.;
  fNumericalFlux.SetOrder(NStateVariables());
  fNumericalFlux.SetDimension(Dimension());
}
TBurgerLaw1D::TBurgerLaw1D(std::istream &input) : TConservationLaw(input) {
  fNumericalFlux.SetOrder(NStateVariables());
  fNumericalFlux.SetDimension(Dimension());
  fSonicPoint = 0.;
}
TBurgerLaw1D::TBurgerLaw1D(int id,char *name,int type) : TConservationLaw(id,name) {
  fNumericalFlux.SetOrder(NStateVariables());
  SetNumericalFluxType(type,type);
  fNumericalFlux.SetDimension(Dimension());
  fSonicPoint = 0.;
}
TBurgerLaw1D::TBurgerLaw1D(TBurgerLaw1D &burger) : TConservationLaw(burger) {
  fSonicPoint = 0.;
}

void TBurgerLaw1D::Flux(TPZVec<REAL> &U,TPZVec<REAL> &flux) {
  flux[0] = .5*U[0]*U[0];
}
void TBurgerLaw1D::JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob) {
  jacob(0,0) = U[0];
}
void TBurgerLaw1D::JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
  jacob(0,0) = normal[0]*U[0];
}
void TBurgerLaw1D::ValJacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &valjacob) {
  valjacob(0,0) = fabs(U[0]);
}
int TBurgerLaw1D::ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta) {
	valjacob(0,0) = fabs(Beta[0]*u[0]);
	return 0;
}
void TBurgerLaw1D::RoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &/*Up1*/,TPZFMatrix<STATE> &Roe) {
  Roe(0,0)=U[0];
}
void TBurgerLaw1D::EigRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &/*Up1*/,TPZVec<REAL> &EigRoe) {
  EigRoe[0]=U[0];
}
void TBurgerLaw1D::ValRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &/*Up1*/,TPZFMatrix<STATE> &ValRoe) {
  ValRoe(0,0)=fabs(U[0]);
}
REAL TBurgerLaw1D::MaxEigJacob(TPZVec<REAL> &U,TPZVec<REAL> &/*normal*/) {
  return fabs(U[0]);
}
REAL TBurgerLaw1D::ValEigJacob(TPZVec<REAL> &u,int order,int dim) {
#ifndef NOTDEBUG
  if(dim!=1) PZError << "TBurgerLaw1D::ValEigJacob. Bad parameter dim.\n";
#endif
  return u[0];
}

/**compute the value of the flux function to be used by ZZ error estimator*/
void TBurgerLaw1D::Flux(TPZVec<REAL> &/*x*/,TPZVec<REAL> &/*Sol*/,TPZFMatrix<STATE> &/*DSol*/,
                        TPZFMatrix<REAL> &/*axes*/,TPZVec<REAL> &/*flux*/) {

}

int TBurgerLaw1D::VariableIndex(char *name) {
  return 0;
}

int TBurgerLaw1D::NSolutionVariables(int index) {
  if(index==0) return 1;
  return TPZMaterial::NSolutionVariables(index);
}
void TBurgerLaw1D::VariablesName(TPZVec<std::string> &names) {
  names[0] = "state_x";
}

/**returns the solution associated with the var index based on the finite element approximation*/
void TBurgerLaw1D::Solution(TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
                 int var,TPZVec<REAL> &Solout){
  if(var == 0) {
    Solout.Resize(1);
    Solout[0] = Sol[0];
    return;
  }
  TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

/**To evaluate true-error, L2-error and estimate error*/
void TBurgerLaw1D::Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<REAL> &flux,
     TPZVec<REAL> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {

  REAL udif;
  //int i;

  values.Fill(0.);
  if(x[0]<0. || x[0]>1.) return;

  udif = sol[0] - u_exact[0];
  values[0] = udif*udif; 
  values[1] = udif*udif;
  values[2] = fabs(udif);
}

/* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
int TBurgerLaw1D::IdBC(double *x) {
// return 0;    // Para caso periodico nao deve ser negativo
  if(x[0]<.1) return -1;
  else return -2;
}

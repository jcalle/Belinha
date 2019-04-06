/*******   File :   burgerbi.c   *******/

#include "burgerbi.h"
#include "pzvec.h"
#include "pzfmatrix.h"

/*******   class : TBurgerLaw2D   *******/

TBurgerLaw2D::TBurgerLaw2D(int id) : TConservationLaw(id,"Burger Conservation Law 2D") {
  fNumericalFlux.SetOrder(NStateVariables());
  fNumericalFlux.SetDimension(Dimension());
}
TBurgerLaw2D::TBurgerLaw2D(std::istream &input) : TConservationLaw(input) {
  fNumericalFlux.SetOrder(NStateVariables());
  fNumericalFlux.SetDimension(Dimension());
}
TBurgerLaw2D::TBurgerLaw2D(int id,char *name,int type) : TConservationLaw(id,name) {
  fNumericalFlux.SetOrder(NStateVariables());
  SetNumericalFluxType(type,type);
  fNumericalFlux.SetDimension(Dimension());
}

TBurgerLaw2D::TBurgerLaw2D(TBurgerLaw2D &burger) : TConservationLaw(burger) {
}

void TBurgerLaw2D::Flux(TPZVec<REAL> &U,TPZVec<REAL> &flux) {
	flux[0] = flux[1] = .5*U[0]*U[0];
}
void TBurgerLaw2D::JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob) {
	jacob(0,0) = jacob(1,0) = U[0];
}
void TBurgerLaw2D::JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
	jacob(0,0) = (normal[0]+normal[1])*U[0];
}
void TBurgerLaw2D::ValJacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &valjacob) {
	valjacob(0,0) = valjacob(1,0) = fabs(U[0]);
}
int TBurgerLaw2D::ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) {
	valjacob(0,0) = fabs((normal[0]+normal[1])*u[0]);
	return 0;
}	

REAL TBurgerLaw2D::MaxEigJacob(TPZVec<REAL> &U,TPZVec<REAL> &normal) {
//	return fabs(U[0]) + fabs(U[0]);
  return fabs(U[0]);
}

REAL TBurgerLaw2D::ValEigJacob(TPZVec<REAL> &u,int order,int dim) {
#ifndef NOTDEBUG
  if(order) PZError << "TBurgerLaw2D::ValEigJacob. Bad parameter dim.\n";
#endif
  return u[0];
}

/**To evaluate true-error, L2-error and estimate error*/
void TBurgerLaw2D::Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<REAL> &flux,
     TPZVec<REAL> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {

  REAL udif;
//  int i;
  udif = sol[0] - u_exact[0];

  values.Fill(0.);

  values[0] += udif*udif; 
  values[1] += udif*udif;
  values[2] += fabs(udif);
}
  
void TBurgerLaw2D::Flux(TPZVec<REAL> &x,TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
		    TPZVec<REAL> &flux) {
}

/**returns the solution associated with the var index based on the finite element approximation*/
void TBurgerLaw2D::Solution(TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
                 int var,TPZVec<STATE> &Solout){
  if(var == 0) {
    Solout.Resize(1);
    Solout[0] = Sol[0];
    return;
  }
  TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

void TBurgerLaw2D::VariablesName(TPZVec<std::string> &names) {
  names[0] = "state_x";
}


/**
void TBurgerLaw2D::SetUZero(int t) {
  if(t==0) {
      cout << endl << "Choose UZero definition" << endl;
      cout << "\tDefinition 1 : Riemann Problem  Uleft=1.,  Uright=0." << endl;
      cout << "\tDefinition 2 : Riemann Problem  Uleft=0.,  Uright=1." << endl;
      cout << "\tDefinition 3 : Oblique initial data (vide JiET)" << endl;
      cout << "\tChoose = ";
      cin >> t;
   }
   if(t==1) fUZero=UZero1;
   else if(t==2) fUZero= UZero2;
   else if(t==3) fUZero= UZero3;
   else fUZero=NULL;
}

TBC *TBurgerLaw2D::CreateBC(int id,int type,TFMatrix &val1,TFMatrix &val2) {
	cout << "TBurgerLaw2D::CreateBC is called." << endl;
   TBC *bc=new TBC(id,type,val1,val2);
   return bc;
}

void UZero1(TPZVec<REAL> &coord,TPZVec<REAL> &result) {
	double x=coord[0],y=coord[1];
   if(x<0. && y<0.) result[0]=.5;
   else if(x<0. && y>0.) result[0]=-.2;
   else if(x>0. && y<0.) result[0]=.8;
   else result[0]=-1.;
}
void UZero2(TPZVec<REAL> &coord,TPZVec<REAL> &result) {
   if(coord[0]<0.) result[0] = 0.;
	else result[0] = 1.;
}
void UZero3(TPZVec<REAL> &coord,TPZVec<REAL> &result) {
	if(coord[0]<0.) result[0] = 1.;
   else result[0] = 0.;
}

*/

#include "stdio.h"
#include "conelaw.h"
#include "linflux.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include <fstream>
#include "myheader.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"

TConeLaw::TConeLaw(int id) : TConservationLaw(id,"Cone Conservation Law 2D") {
  fNumericalFlux.SetOrder(1);
  fNumericalFlux.SetDimension(Dimension());
}
TConeLaw::TConeLaw(std::istream &input) : TConservationLaw(input) {
  fNumericalFlux.SetOrder(1);
  fNumericalFlux.SetDimension(Dimension());
}
TConeLaw::TConeLaw(int id,char *name,int type) : TConservationLaw(id,name) {
  fNumericalFlux.SetOrder(1);
  fNumericalFlux.SetDimension(Dimension());
  SetNumericalFluxType(type,type);
}

void TConeLaw::Flux(TPZVec<REAL> &Ui,TPZVec<REAL> &funcao) {
  funcao[0] = -1.*Point[1]*Ui[0];
  funcao[1] = Point[0]*Ui[0];
}

int TConeLaw::VariableIndex(char *name) {
  return 0;
}

void TConeLaw::VariablesName(TPZVec<char *> &names) {
  names[0] = "state_cone";
}

int TConeLaw::NSolutionVariables(int index) {
  if(index < 1) return 1;
  return TPZMaterial::NSolutionVariables(index);
}

/**returns the solution associated with the var index based on the finite element approximation*/
void TConeLaw::Solution(TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
                 int var,TPZVec<STATE> &Solout){
  if(!var) {
    Solout.Resize(1);
    Solout[0] = Sol[var];
    return;
  }
  TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

void TConeLaw::JacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &jacob) {
  jacob(0,0) = -1.*Point[1];
  jacob(1,0) = Point[0];
}
void TConeLaw::JacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
  jacob(0,0) = -1.*Point[1]*normal[0] + normal[1]*Point[0];
}

REAL TConeLaw::MaxEigJacob(TPZVec<REAL> &/*U*/,TPZVec<REAL> &normal) {
  return Max(fabs(Point[0]),fabs(Point[1]));
}

REAL TConeLaw::ValEigJacob(TPZVec<REAL> &u,int order,int dim) {
#ifndef NOTDEBUG
  if(order) PZError << "TConeLaw::ValEigJacob, undefined order.\n";
#endif
  if(dim==1) 
    return -1.*Point[1];
  return Point[0];
}

void TConeLaw::EigRoeMatrix(TPZVec<REAL> &/*Ul*/,TPZVec<REAL> &/*Ur*/,TPZVec<REAL> &EigRoe) {
  EigRoe[0] = -1.*Point[1];
  EigRoe[1] = Point[0];
}

void TConeLaw::RoeMatrix(TPZVec<REAL> &/*Ul*/,TPZVec<REAL> &/*Ur*/,TPZFMatrix<STATE> &Roe) {
  Roe(0,0) = -1.*Point[1];
  Roe(1,0) = Point[0];
}

void TConeLaw::ValRoeMatrix(TPZVec<REAL> &/*Ul*/,TPZVec<REAL> &/*Ur*/,TPZFMatrix<STATE> &ValRoe) {
  ValRoe(0,0) = fabs(Point[1]);
  ValRoe(1,0) = fabs(Point[0]);
}

void TConeLaw::ValJacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &valjacob) {
  valjacob(0,0) = fabs(Point[1]);
  valjacob(1,0) = fabs(Point[0]);
}
int TConeLaw::ValJacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) {
  valjacob(0,0) = fabs((-1.)*normal[0]*Point[1]+normal[1]*Point[0]);
	return 0;
}

/**To evaluate true-error, L2-error and estimate error*/
void TConeLaw::Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<REAL> &flux,
     TPZVec<REAL> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {

//  SetPoint(x);
  REAL udif = sol[0] - u_exact[0];

  values[0] = fabs(udif);
  values[1] = udif * udif;
  values[2] = fabs(udif);
}

int TConeLaw::IdBC(double *x) {
  return -1;
}
/*
void TConeLaw::IncrementDiffusion(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,
        REAL weight,REAL area,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,
        TPZFMatrix &ek,TPZFMatrix &ef,REAL cfldiff) {
  SetPoint(x);
//  area = sqrt(area);
  int phr = phi.Rows();
  int i,j,k,p,c, nvar = NStateVariables(), dim=Dimension();
	if(nvar!=1 || dim!=2) {
	  PZError << "TConeLaw::IncrementDiffusion. Bad number of variables parameter.\n";
		return;
	}
	
  TPZFMatrix dphixy(dim,phr,0.);
  for(i=0;i<phr;i++)
    for(j=0;j<dim;j++)
      for(k=0;k<dim;k++) dphixy(j,i) += dphi(k,i)*axes(k,j);

  InvJacob2d(axes,jacinv);   // Retorna em jacinv a inv do jacobiano da transf do elemento mestre e o elemento deformado

  int dsolr = dsol.Rows(), dsolc = dsol.Cols();
  REAL module, alpha, beta;  // com o determinante da matriz L0(phi(i)), para normalizar
  //Jacobiano da funcao fluxo
  TPZFMatrix jacob(dim,1);
  JacobFlux(sol,jacob);
  //Valores para L0(phi(i))
	REAL tau = 0., prod, soma;
	for(c=0;c<dim;c++) {
	  soma = 0;
	  for(k=0;k<dim;k++)
		  soma += jacob(k,0)*jacinv(c,k);
		tau += soma*soma;
	}
	tau = 1./sqrt(tau);
	tau *= (weight*fDeltaT*fCoef*cfldiff*area);
	
 if(fExplicit<1) {
  for(i=0;i<phr;i++) {
    prod = 0;
    for(c=0;c<dim;c++)
      prod += jacob(c,0)*dphixy(c,i);

    for(j=0;j<phr;j++) {
		  soma = 0.;
      for(c=0;c<dim;c++)
        soma += jacob(c,0)*dphixy(c,j);
			ek(i,j) += (prod*tau*soma);
    }
  }
 }
 else if(fExplicit==1) {
  TPZVec<REAL> divflux(nvar,0.);
	TPZVec<REAL> Beta(nvar);
  for(i=0;i<dim;i++) {
    p = i*nvar;
    for(j=0;j<nvar;j++) {
//      Beta[p+j] = ValEigJacob(sol,j,i+1);
      for(k=0;k<nvar;k++)
        divflux[j] += jacob(j+p,k)*dsol(i,k);
    }
  }
  for(i=0;i<phr;i++) {
    p = i*nvar;
    for(j=0;j<nvar;j++) {
      module = 0.;
      for(c=0;c<dim;c++) module += Beta[c*nvar+j]*Beta[c*nvar+j];
      module = sqrt(module);
      if(IsZero(module)) break;
      for(c=0;c<phr;c++) {
        ek(p+j,c*nvar+j) += weight*cfldiff*area*phi(c,0);
      }
      ef(p+j,0) -= weight*cfldiff*area*divflux[j];
    }
  }
 }
 else {
  TPZVec<REAL> divflux(nvar,0.);
	TPZVec<REAL> Beta(nvar);
  for(i=0;i<dim;i++) {
    p = i*nvar;
    for(j=0;j<nvar;j++) {
//      Beta[p+j] = ValEigJacob(sol,j,i+1);
      for(k=0;k<nvar;k++)
        divflux[j] += jacob(j+p,k)*dsol(i,k);
    }
  }
  for(i=0;i<phr;i++) {
    p = i*nvar;
    for(j=0;j<nvar;j++) {
      module = 0.;
      for(c=0;c<dim;c++) module += Beta[c*nvar+j]*Beta[c*nvar+j];
      module = sqrt(module);
      if(IsZero(module)) break;
      ef(p+j,0) -= weight*cfldiff*area*divflux[j];
    }
  }
 }
}
*/

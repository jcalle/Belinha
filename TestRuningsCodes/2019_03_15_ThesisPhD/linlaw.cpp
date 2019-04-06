#include <stdio.h>
#include "linlaw.h"
#include "linflux.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include <fstream>
#include "myheader.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"

TLinearLaw::TLinearLaw(int id,int order) : TConservationLaw(id,"Linear Conservation Law 1D") {
  while(order < 1) {
    PZError << "TLinearLaw constructor has bad order " << order << std::endl;
    PZError << "New order (Quit=0) ";
	std::cin >> order;
    if(!order) return;
  }
  fOrder = order;
  fNumericalFlux.SetOrder(order);
  fNumericalFlux.SetDimension(Dimension());
}
TLinearLaw::TLinearLaw(std::istream &input) : TConservationLaw(input) {
  char c;
  GetDataCommented(input,fOrder);
  while(fOrder < 1) {
    PZError << "TLinearLaw constructor has bad order " << fOrder << std::endl;
    PZError << "New order (Quit=0) ";
	std::cin >> fOrder;
    if(!fOrder) return;
  }
  input >> c;
  SetMatrix(input);
  fNumericalFlux.SetOrder(fOrder);
  fNumericalFlux.SetDimension(Dimension());
}
TLinearLaw::TLinearLaw(int id,int order,char *name,int type) : TConservationLaw(id,name) {
  while(order < 0) {
    PZError << "TLinearLaw constructor has bad order " << order << std::endl;
    PZError << "New order (Quit=0) ";
	std::cin >> order;
    if(!order) return;
  }
  if(order>0) {
    fOrder = order;
    RequireMatrix();
  }
  else fOrder = 1;
  fNumericalFlux.SetOrder(fOrder);
  fNumericalFlux.SetDimension(Dimension());
  SetNumericalFluxType(type,type);
}

int TLinearLaw::IdBC(double *x) { 
//	if(x[0] < .5) 
		return -4;
//	else 
//		return -2;
}

/*void TLinearLaw::IncrementDiffusion(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix &dsol,
        REAL weight,REAL area,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {
  int phr = phi.Rows();
  int i,j,k,p,c, nvar = NStateVariables(), dim=Dimension();
  if(dim==2) area = sqrt(area);
  TPZFMatrix dphixy(dim,phr,0.);
  for(i=0;i<phr;i++)
    for(j=0;j<dim;j++)
      for(k=0;k<dim;k++) dphixy(j,i) += dphi(k,i)*axes(k,j);

  int dsolr = dsol.Rows(), dsolc = dsol.Cols();
  REAL module;

  TPZVec<REAL> divflux(nvar,0.);
  TPZVec<REAL> Beta(dim*nvar,0.);   // Beta = Jacob(F(U)) em soma
  TPZVec<REAL> diffvec(nvar,0.);
  TPZFMatrix jacob(nvar*dim,nvar);
  JacobFlux(sol,jacob);
  for(i=0;i<dim;i++) {
    p = i*nvar;
    for(j=0;j<nvar;j++) {
      Beta[p+j] = ValEigJacob(sol,j,i+1);
      for(k=0;k<nvar;k++)
        divflux[j] += jacob(j+p,k)*dsol(i,k);
    }
  }
  for(i=0;i<phr;i++) {
    p = i*nvar;
    for(j=0;j<nvar;j++) diffvec[j] = 0.;
    for(c=0;c<dim;c++)
      for(j=0;j<nvar;j++)
        diffvec[j] += Beta[c*nvar+j]*dphixy(c,i);
    for(j=0;j<nvar;j++) {
      module = 0.;
      for(c=0;c<dim;c++) module += Beta[c*nvar+j]*Beta[c*nvar+j];
      module = sqrt(module);
      if(IsZero(module)) break;
//      for(c=0;c<phr;c++)
//        ek(p+j,c*nvar+j) += weight*fFactor*area*(diffvec[j]/module)*phi(c,0);
      ef(p+j,0) -= weight*fFactor*area*divflux[j]*(diffvec[j]/module);
    }
  }
}*/

void TLinearLaw::Flux(TPZVec<REAL> &Ui,TPZVec<REAL> &funcao) {
  int i,j;
  for(i=0;i<fOrder;i++) {
    funcao[i]=0.;
    for(j=0;j<fOrder;j++)
      funcao[i]+=(fA[i*fOrder+j]*Ui[j]);
  }
}

int TLinearLaw::VariableIndex(char *name) {
  int index = atoi(name+6) - 1;
  if(index<fOrder) return index;
  return TPZMaterial::VariableIndex(name);
}

void TLinearLaw::VariablesName(TPZVec<std::string> &names) {
  int i, n = AllVariablesPost();
  for(i=0;i<n;i++) {
    switch(i) {
    case 0:
      names[0] = "state_1";
      break;
    case 1:
      names[1] = "state_2";
      break;
    case 2:
      names[2] = "state_3";
      break;
    default:
      names[i] = "state_nonumber";
      break;
    }
  }
}

int TLinearLaw::NSolutionVariables(int index) {
  if(index < AllVariablesImplemented()) return 1;
  return TPZMaterial::NSolutionVariables(index);
}

/**returns the solution associated with the var index based on the finite element approximation*/
void TLinearLaw::Solution(TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
                 int var,TPZVec<STATE> &Solout){
  if(var<fOrder) {
    Solout.Resize(1);
    Solout[0] = Sol[var];
    return;
  }
  TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

void TLinearLaw::ApplyPeriodic(TPZCompMesh &cmesh) {
  int i, nelem = cmesh.NElements();
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  TPZCompEl *cel;
  TPZInterpolatedElement *el;
	int seqnum1, seqnum2;

  for(i=0;i<nelem;i++) {
    cel = elvec[i];
    if(!cel || cel->IsInterface()) continue;
    el = (TPZInterpolatedElement *)cel;
    if(el->Material()->Id()>-1) continue;
    if(IsZero(el->Reference()->NodePtr(0)->Coord(0)-1.)) {
 //     if(el->IsConnectContinuous(0))
   //     seqnum2 = el->Connect(0).SequenceNumber();
     // else 
		seqnum2 = el->Connect(1).SequenceNumber();
    }
    else {
     // if(el->IsConnectContinuous(0))
       // seqnum1 = el->Connect(0).SequenceNumber();
     // else 
		seqnum1 = el->Connect(1).SequenceNumber();
    }
  }
  int size = cmesh.Block().Size(seqnum2);
  for(i=0;i<size;i++)
    cmesh.Block()(seqnum1,0,i,0) = cmesh.Block()(seqnum2,0,i,0);
}

void TLinearLaw::Print(std::ostream &out) {
  int i,j;
  //	PrintData1D();
  out << "\nMatrix A :";
  for(i=0;i<fOrder;i++) {
    out << std::endl;
    for(j=0;j<fOrder;j++)
      out << fA[i] << "\t";
  }
  out << "\nAuto-Valores de A :\n";
  for(i=0;i<fOrder;i++)
    out << fEigVal[i] << "\t";
  out << "\nMatrix dos auto-vetores direitos de A :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fEigVect[i] << "\t";
  }
  out << "\nInversa da matriz dos auto-vetores de A:";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fInvEigVect[i] << "\t";
  }
  out << "\nMaximo autovalor = \n" << fMaxEigVal;
  out << "\nMatrix |A| :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fValA[i] << "\t";
  }
}

void TLinearLaw::JacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &jacob) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      jacob(i,j)=fA[i*fOrder+j];
}
void TLinearLaw::JacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      jacob(i,j) = normal[0]*fA[i*fOrder+j];
}

REAL TLinearLaw::MaxEigJacob(TPZVec<REAL> &/*U*/,TPZVec<REAL> &/*normal*/) {
  return fMaxEigVal;
}
REAL TLinearLaw::ValEigJacob(TPZVec<REAL> &u,int order,int dim) {
#ifndef NOTDEBUG
  if(dim!=1) PZError << "TLinearLaw::ValEigJacob Bad parameter dim.\n";
#endif
  return fEigVal[order];
}

void TLinearLaw::EigRoeMatrix(TPZVec<REAL> &/*Ul*/,TPZVec<REAL> &/*Ur*/,TPZVec<REAL> &EigRoe) {
  int i,j,k=fOrder*fOrder;
  for(i=0;i<k;i++) EigRoe[i] = fA[i];
  for(i=0;i<fOrder;i++) EigRoe[k+i]=fEigVal[i];
  for(i=0;i<k;i++) EigRoe[k+fOrder+i]=fEigVect[i];
  j=fOrder*(1+2*fOrder);
  for(i=0;i<k;i++) EigRoe[j+i]=fInvEigVect[i];
}

void TLinearLaw::RoeMatrix(TPZVec<REAL> &/*Ul*/,TPZVec<REAL> &/*Ur*/,TPZFMatrix<STATE> &Roe) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      Roe(i,j)=fA[i*fOrder+j];
}

void TLinearLaw::ValRoeMatrix(TPZVec<REAL> &/*Ul*/,TPZVec<REAL> &/*Ur*/,TPZFMatrix<STATE> &ValRoe) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      ValRoe(i,j)=fValA[i*fOrder+j];
}

void TLinearLaw::ValJacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &valjacob) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      valjacob(i,j)=fValA[i*fOrder+j];
}
int TLinearLaw::ValJacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) {
  int i,j;
	if(fOrder>4) return 1;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      valjacob(i,j)=fabs(normal[0])*fValA[i*fOrder+j];
	return 0;
}

void TLinearLaw::InverseJacob(TPZFMatrix<STATE> &jac) {
	REAL a = jac(0,0);
	REAL det, b, c, d, e, f, g, h, i;
	REAL EIHF, FGDI, HDGE;
	switch(fOrder) {
	case 1:
		jac(0,0) = 1./a;
		break;
	case 2:
		det = a*jac(1,1)-jac(0,1)*jac(1,0);
		det = 1./det;
	  if(IsZero(det)) {
  	  PZError << "TLinearLaw::InverseJacob. Determinante zero, matriz singular.\n";
	  	exit(1);
	  }
		jac(0,0) = jac(1,1)*det;
		jac(1,1) = a*det;
		det *= -1.;
		jac(0,1) *= det;
		jac(1,0) *= det;
		break;
	case 3:
		b=jac.GetVal(0,1); c=jac.GetVal(0,2);
		d=jac.GetVal(1,0); e=jac.GetVal(1,1); f=jac.GetVal(1,2);
		g=jac.GetVal(2,0); h=jac.GetVal(2,1); i=jac.GetVal(2,2);
	  EIHF = e*i-h*f; FGDI = f*g-d*i; HDGE = h*d-g*e;
  	det = a*EIHF + b*FGDI + c*HDGE;
	  if(IsZero(det)) {
  	  PZError << "TLinearLaw::InverseJacob. Determinante zero, matriz singular.\n";
	  	exit(1);
	  }
  
  	det = 1./det;
	  jac(0,0) = det*EIHF;
  	jac(1,0) = det*FGDI;
	  jac(2,0) = det*HDGE;
	  jac(0,1) = det*(c*h-b*i);
  	jac(1,1) = det*(a*i-g*c);
	  jac(2,1) = det*(g*b-a*h);
  	jac(0,2) = det*(b*f-e*c);
	  jac(1,2) = det*(d*c-a*f);
	  jac(2,2) = det*(a*e-d*b);
		break;
	case 4:
	{
		b=jac.GetVal(0,1); c=jac.GetVal(0,2); d=jac.GetVal(0,3);
		e=jac.GetVal(1,0); f=jac.GetVal(1,1); g=jac.GetVal(1,2); h=jac.GetVal(1,3);
		i=jac.GetVal(2,0);
		REAL j=jac.GetVal(2,1), k=jac.GetVal(2,2), l=jac.GetVal(2,3);
		REAL m=jac.GetVal(3,0), n=jac.GetVal(3,1), r=jac.GetVal(3,2), s=jac.GetVal(3,3);
	  REAL BHDF = b*h-d*f, BKCJ = b*k-c*j, BSDN = b*s-d*n, CFBG = c*f-b*g;
		REAL CLDK = c*l-d*k, CNBR = c*n-b*r, DGCH = d*g-c*h, DJBL = d*j-b*l;
		REAL DRCS = d*r-c*s, FRGN = f*r-g*n, GJFK = g*j-f*k, GSHR = g*s-h*r;
		REAL FSHN = f*s-h*n, LNJS = l*n-j*s, HJFL = h*j-f*l, HKGL = h*k-g*l;
		REAL KNJR = k*n-j*r, LRKS = l*r-k*s;
	  det = DJBL*(g*m-e*r)+BHDF*(k*m-i*r)+BSDN*(g*i-e*k);
  	det += (FSHN*(a*k-c*i)+LNJS*(a*g-c*e)+HJFL*(a*r-c*m));
	 	if(IsZero(det)) {
  		PZError << "TEulerLaw2D::InverseJacob. Determinante zero, matriz singular.\n";
	 		exit(1);
	  }
  
  	det = 1./det;
	  jac(0,0) = det*(k*FSHN+g*LNJS+r*HJFL);
  	jac(1,0) = det*(m*HKGL+i*GSHR+e*LRKS);
	  jac(2,0) = (-det)*(m*HJFL+i*FSHN+e*LNJS);
  	jac(3,0) = det*(m*GJFK+i*FRGN+e*KNJR);
	  jac(0,1) = det*(d*KNJR+b*LRKS-c*LNJS);
	  jac(1,1) = det*(m*CLDK+i*DRCS-a*LRKS);
  	jac(2,1) = det*(m*DJBL+i*BSDN+a*LNJS);
	  jac(3,1) = det*(m*BKCJ+i*CNBR-a*KNJR);
  	jac(0,2) = det*(d*FRGN-c*FSHN+b*GSHR);
	  jac(1,2) = det*(m*DGCH-e*DRCS-a*GSHR);
  	jac(2,2) = det*(m*BHDF-e*BSDN+a*FSHN);
	  jac(3,2) = det*(m*CFBG-e*CNBR-a*FRGN);
  	jac(0,3) = det*(d*GJFK-c*HJFL+b*HKGL);
	  jac(1,3) = (-det)*(i*DGCH+e*CLDK+a*HKGL);
  	jac(2,3) = det*(a*HJFL-i*BHDF-e*DJBL);
	  jac(3,3) = (-det)*(i*CFBG+e*BKCJ+a*GJFK);
		break;
	}
	default:
		jac.Identity();
	}
}

void TLinearLaw::SetData(std::istream &input) {
  TConservationLaw::SetData(input);
  SetMatrix(input);
}

void TLinearLaw::SetMatrix(std::istream &input) {
  int i,j,k,index;
  double temp;
  fMaxEigVal=0.;
  //Limpando para armazenar
  for(i=0;i<MAXORDER*MAXORDER;i++) fA[i]=fEigVect[i]=fValA[i]=fInvEigVect[i]=0.;
  for(i=0;i<MAXORDER;i++) fEigVal[i]=0.;
  //Ingressa dados da matriz A
  GetCommentary(input);   // Matrix A, Eigenvalues, Eigenvectors and InverseEigenvectors
  for(i=0;i<fOrder*fOrder;i++)
    input >> fA[i];
  //Ingresso dos autovalores, autovetores e a inversa dos autovetores
  for(i=0;i<fOrder;i++) {
    input >> temp;
    fMaxEigVal=Max(fMaxEigVal,fabs(temp));
    fEigVal[i]=temp;
  }
  for(i=0;i<fOrder*fOrder;i++)
    input >> fEigVect[i];
  for(i=0;i<fOrder*fOrder;i++)
    input >> fInvEigVect[i];

  REAL prod, eigval;
  for(i=0;i<fOrder;i++) {
    prod = 0.;
    index = i*fOrder;
    for(j=0;j<fOrder;j++) {
      for(k=0;k<fOrder;k++) {
        eigval = fEigVal[k]/fabs(fEigVal[k]);
        prod += (fEigVect[index+k]*eigval*fInvEigVect[k*fOrder+j]);
      }
      fInverseQrt[index+j] = prod;
    }
  }
  double aux;
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
      aux=0.;
      for(k=0;k<fOrder;k++) {
        temp = ValAbs(fEigVal[k])*fInvEigVect[k*fOrder+j];
        temp *= fEigVect[index+k];
        aux += temp;
      }
      fValA[index+j] = aux;
    }
  }

  std::cout << "\nVerifique :\n Matrix A \t EigVector\t InvEigVect\t EigValue";
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
		std::cout << "\n" << fA[index+j] << "\t" << fEigVect[index+j];
		std::cout << "\t" << fInvEigVect[index+j];
    }
	std::cout << "\t" << fEigVal[i];
  }
  std::cout << "\n";
}

void TLinearLaw::SetMatrix(TPZFMatrix<STATE> &A,TPZFMatrix<STATE> &EigVect,TPZFMatrix<STATE> &InvEigVect,
			   TPZVec<REAL> &EigVal) {
  int i,j,k,index;
  double temp,eigval;
  fMaxEigVal=0.;
#ifdef DEBUG
  if(fOrder!=A.Rows() || fOrder!=EigVect.Rows() || fOrder!=InvEigVect.Rows() ||
     fOrder!=EigVal.NElements()) {
    PZError << "TLinearLaw1D::SetMatrix. Some matrix is out of the dimension.\n";
    return;
  }
#endif
  for(i=0;i<fOrder;i++) {
    index = i*fOrder;
    for(j=0;j<fOrder;j++)
      fA[index+j] = A(i,j);
  }

  //Ingresso dos autovalores, autovetores e a inversa dos autovetores
  for(i=0;i<fOrder;i++) {
    temp = EigVal[i];
    fMaxEigVal=Max(fMaxEigVal,fabs(temp));
    fEigVal[i]=temp;
  }
  for(i=0;i<fOrder;i++) {
    index = i*fOrder;
    for(j=0;j<fOrder;j++)
      fEigVect[index+j] = EigVect(i,j);
  }
  for(i=0;i<fOrder;i++) {
    index = i*fOrder;
    for(j=0;j<fOrder;j++)
      fInvEigVect[index+j] = InvEigVect(i,j);
  }
  /**Criando a matriz C = A*|A|¯1 */
  for(i=0;i<fOrder;i++) {
    index = i*fOrder;
    temp = 0.;
    for(j=0;j<fOrder;j++) {
      for(k=0;k<fOrder;k++) {
        eigval = fEigVal[k]/fabs(fEigVal[k]);
        temp += (fEigVect[index+k]*eigval*fInvEigVect[k*fOrder+j]);
      }
      fInverseQrt[index+j] = temp;
    }
  }
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
      temp=0.;
      for(k=0;k<fOrder;k++) {
	     eigval = ValAbs(fEigVal[k])*fInvEigVect[k*fOrder+j];
	     eigval *= fEigVect[index+k];
	     temp += eigval;
      }
      fValA[index+j] = temp;
    }
  }
}

void TLinearLaw::RequireMatrix() {
  int i;
  std::cout << "\nIntroduz a matriz A pelo teclado? (y/n) ";
  char c[5];
  std::cin >> c[0];
  if(c[0]=='y') {
    fMaxEigVal=0.;
	std::cout << "\nINGRESE DADOS POR LINHA\n";
	std::cout << "1- Matrix A,\n";
    for(i=0;i<fOrder*fOrder;i++)
		std::cin >> fA[i];
	std::cout << "2- Autovalores de A,\n";
    for(i=0;i<fOrder;i++) {
		std::cin >> fEigVal[i];
      fMaxEigVal=Max(fMaxEigVal,fabs(fEigVal[i]));
    }
	std::cout << "3- Matriz R dos autovetores de A,\n";
    for(i=0;i<fOrder*fOrder;i++)
		std::cin >> fEigVect[i];
	std::cout << "4- Matriz inversa de R,\n";
    for(i=0;i<fOrder*fOrder;i++)
		std::cin >> fInvEigVect[i];
	std::cout << "Dados Completos.\n";
  }
  else {
	  std::cout << "\nNome do arquivo com dados relacionados a matriz A ";
    char Amatrix[20];
	std::cin >> Amatrix;
	std::ifstream in(Amatrix);
    SetMatrix(in);
    in.close();
  }
  int j,k,index;
  REAL temp,eigval;
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
      temp=0.;
      for(k=0;k<fOrder;k++) {
	     eigval = ValAbs(fEigVal[k])*fInvEigVect[k*fOrder+j];
	     eigval *= fEigVect[index+k];
	     temp += eigval;
      }
      fValA[index+j] = temp;
    }
  }
}

void TLinearLaw::InverseQrt(TPZVec<REAL> &sol,TPZFMatrix<STATE> &C) {
  int i, j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      C(i,j) = fInverseQrt[i*fOrder+j];
}

/**To evaluate true-error, L2-error and estimate error*/
void TLinearLaw::Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<REAL> &flux,
     TPZVec<REAL> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {

  TPZVec<REAL> udif(sol);
  int i, nvar = NStateVariables();
  for(i=0;i<nvar;i++) udif[i] -= u_exact[i];

  values.Fill(0.);

  for(i=0;i<nvar;i++) {
    values[0] += udif[i]*udif[i];
    values[1] += udif[i]*udif[i];
    values[2] += fabs(udif[i]);
  }
}

/*
void TLinearLaw::IncrementDiffusion(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix &dsol,
        REAL weight,REAL area,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {

  TPZVec<REAL> divflux(fOrder,0.);
  int dphr = dphi.Rows(), phr = dphi.Cols();
  int i, k, c, l;
  TPZFMatrix Beta(dphr*fOrder,fOrder);
  InverseQrt(sol,Beta);
  TPZFMatrix jacob(dphr*fOrder,fOrder);
  JacobFlux(sol,jacob);

  for(c=0;c<fOrder;c++) {
    for(l=0;l<dphr;l++)
      for(k=0;k<fOrder;k++)
        divflux[c] += jacob(l*fOrder+c,k)*dsol(l,k);
  }
  for(c=0;c<fOrder;c++) divflux[c] *= fFactor*area;

  REAL prod;
  for(i=0;i<phr;i++) {
    for(k=0;k<dphr;k++) {
      for(c=0;c<fOrder;c++) {
        prod = 0.;
        for(l=0;l<fOrder;l++)
          prod += Beta(k*fOrder+c,l)*divflux[l];
        ef(i*fOrder+c,0) -= weight*prod*dphi(k,i);
      }
    }
  }
}
*/
/***   Resolve aproximadamente lei de conservacao a partir de um arquivo   ***
       void TLinearLaw::Solver(FILE *type) {
       int i, flux, NTimes;
       double TimesChangeBC[10];
       double CFL;
       fscanf(type,"%d",&NTimes);
       if(NTimes>10) {
       printf("\n Apenas podem mudarse 10 veces as BC.");
       exit(1);
       }
       for(i=0;i<NTimes;i++)
       fscanf(type,"%lf",TimesChangeBC+i);
       // Nome do arquivo para armazenar solucao aproximada
       char Uout[20]; Uout[0]='\0';
       strncat(Uout,Name(),3);
       fscanf(type,"%lf",&CFL);
       while(CFL>0.) {
       SetCFL(CFL);
       fscanf(type,"%d",&flux);
       TLinearFlux F(*this,flux);
       while(F.FluxType()!=0) {
       for(i=3;i<20;i++) Uout[i]='\0';
       strncat(Uout,F.FluxName(),5);
       strcat(Uout,".U");

       printf("\nFLUXO ATUAL %u",F.FluxType());
       printf("\nArquivo de saida %s", Uout);

       ClearArray(fU,fNCelules*fOrder);
       ClearArray(fU1,fNCelules*fOrder);
       InitialValueU();
       GetUTimeMax(F,Uout,TimesChangeBC,NTimes);
       fscanf(type,"%d",&flux);
       F.SetFluxType(flux);
    }
    fscanf(type,"%lf",&CFL);
  }
}

***   Resolve aproximadamente a lei de conservacao(usuario)   ***
       void TLinearLaw::Solver() {
       int i,NTimes;
       double TimesChangeBC[10];
       double CFL;

       char Uout[20]; Uout[0]='\0';
       strncat(Uout,Name(),3);

       printf("\n Numero de Tempos para mudancas, Ntempos = ");
       scanf("%d",&NTimes);
       if(NTimes>10) {
       printf("\n Apenas podem mudarse 10 veces as BC.");
       exit(1);
       }
       printf("\n Ingrese os tempos :\n");
       for(i=0;i<NTimes;i++)
       scanf("%lf",TimesChangeBC+i);

       printf("\n\n (Quit 0)Valor do CFL = ");
       scanf("%lf",&CFL);
       while(CFL>0.) {
       SetCFL(CFL);
       TLinearFlux F(*this);
       while(F.FluxType()!=0) {
       for(i=3;i<20;i++) Uout[i]='\0';
       strncat(Uout,F.FluxName(),5);
       strcat(Uout,".U");

       printf("\nFLUXO ATUAL %u",F.FluxType());
       printf("\nArquivo de saida %s",Uout);

       ClearArray(fU,fNCelules*fOrder);
       ClearArray(fU1,fNCelules*fOrder);

       InitialValueU();

       GetUTimeMax(F,Uout,TimesChangeBC,NTimes);

       F.RequireFluxType();
       }
       printf("\n (Quit 0) Proximo CFL = ");
       scanf("%lf",&CFL);
       }
       }

*/
/*
void TLinearLaw::Contribute(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
     TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {
  int phc,phr,dphc,dphr,efr,efc,ekr,ekc;
  int nvar = NStateVariables();
  phc = phi.Cols();
  phr = phi.Rows();
  dphc = dphi.Cols();
  dphr = dphi.Rows();
  efr = ef.Rows();
  efc = ef.Cols();
  ekr = ek.Rows();
  ekc = ek.Cols();
  if(phc != 1 || phr != dphc ||
     ekr != phr*nvar || ekc != ekr ||
     efr != ekr || efc != 1) {
    PZError << "\nTConservationLaw. Inconsistent input data : \n" <<
      "phi.Cols() = " << phi.Cols() << " dphi.Cols() = " << dphi.Cols() <<
      " phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
      dphi.Rows() << "\nek.Rows() = " << ek.Rows() << " ek.Cols() = "
	    << ek.Cols() <<
      "\nef.Rows() = " << ef.Rows() << " ef.Cols() = "
	    << ef.Cols() << "\n";
  }

  int i,j,l,c;
  TPZVec<REAL> res(nvar,0.),flux(dphr*nvar,0.); //dphr is a spatial domain dimension
  REAL prod, soma;

  if(fForcingFunction) {
    for(i=0;i<nvar;i++) res[i] = sol[i];
    /**Observar que res contem a solucao Ul-1, e a funcao s(tk,x,sol(l-1))
       toma a tempo tk da lei de conservacao
    fForcingFunction(x,res);
  }

  Flux(sol,flux);
  for(i=0;i<phr;i++) {
    /**Calculo de keij, apenas o primeiro valor do bloco
    for(j=0;j<phr;j++) {

      //int proddphi = 0.;
      //      for(l=0;l<dphr;l++) proddphi += (dphi(l,i)*dphi(l,j));

      prod = weight*phi(i,0)*phi(j,0);
      for(c=0;c<nvar;c++) {
	//        soma = 1.;                     
	//        for(int p=0;p<nvar;p++) soma += fA[c*nvar+p];
	//        soma = phi(j,0) - (fMinDeltaX * soma * dphi(0,j));
        ek(i*nvar+c,j*nvar+c) += prod;
      }
    }

    /**Calculo de feic
    for(c=0;c<nvar;c++) {
      prod = 0.;
      for(l=0;l<dphr;l++) prod += dphi(l,i)*flux[l*nvar+c];
      ef(i*nvar+c,0) += (weight*(phi(i,0)*res[c]+prod));
    }
  }
}
*/


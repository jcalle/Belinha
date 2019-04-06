/*******       File : LinLawBi.c       *******/

#include "linlawbi.h"
#include "linflux.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include <fstream>
#include "myheader.h"
#include <math.h>

TLinearLaw2D::TLinearLaw2D(int id,int order) : TLinearLaw(id,order) {
  SetName("Linear conservation law 2D.");
}
TLinearLaw2D::TLinearLaw2D(std::istream &input) : TLinearLaw(input) {
}
TLinearLaw2D::TLinearLaw2D(int id,int order,char *name,int type) :
	TLinearLaw(id,order,name,type) {
}

/* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
int TLinearLaw2D::IdBC(double *x) {
  if(x[0]<0.5) return -1;
  else if(x[1] < 1.) return -2;
  return -1;
}
void TLinearLaw2D::Flux(TPZVec<REAL> &Ui,TPZVec<REAL> &funcao) {
  int i, j;
  for(i=0;i<fOrder;i++) {
    funcao[i]=0.;
    funcao[fOrder+i] = 0.;
    for(j=0;j<fOrder;j++) {
      funcao[i]+=(fA[i*fOrder+j]*Ui[j]);
      funcao[fOrder+i] += (fA[(fOrder+i)*fOrder+j]*Ui[j]);
    }
  }
}

void TLinearLaw2D::JacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &jacob) {
  int i,j;
  for(i=0;i<2*fOrder;i++)
    for(j=0;j<fOrder;j++)
      jacob(i,j)=fA[i*fOrder+j];
}
void TLinearLaw2D::JacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      jacob(i,j)=normal[0]*fA[i*fOrder+j]+normal[1]*fA[2*i*fOrder+j];
}

REAL TLinearLaw2D::MaxEigJacob(TPZVec<REAL> &/*U*/,TPZVec<REAL> &/*normal*/) {
  return Max(fabs(fMaxEigVal),fabs(fMaxEigValB));
}

REAL TLinearLaw2D::ValEigJacob(TPZVec<REAL> &u,int order,int dim) {
#ifndef NOTDEBUG
  if(dim!=1 && dim!=2)
    PZError << "TLinearLaw2D::ValEigJacob Bad parameter dim.\n";
#endif
  return fEigVal[(dim-1)*fOrder+order];
}

void TLinearLaw2D::Print(std::ostream &out) {
  int i,j;
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
  out << "\nInversa da matriz dos auto-vetores de A :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fInvEigVect[i] << "\t";
  }
  out << "\nMaximo autovalor = %lf\n" << fMaxEigVal;
  out << "\nMatrix |A| :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fValA[i] << "\t";
  }
  out << "\nMatrix B :";
  int index = fOrder*fOrder;
  for(i=0;i<fOrder;i++) {
    out << std::endl;
    for(j=0;j<fOrder;j++)
      out << fA[index+i] << "\t";
  }
  out << "\nAuto-Valores de B :\n";
  for(i=0;i<fOrder;i++)
    out << fEigVal[fOrder+i] << "\t";
  out << "\nMatrix dos auto-vetores direitos de B :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fEigVect[index+i] << "\t";
  }
  out << "\nInversa da matriz dos auto-vetores de B :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fInvEigVect[index+i] << "\t";
  }
  out << "\nMaximo autovalor de B = %lf\n" << fMaxEigValB;
  out << "\nMatrix |B| :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fValA[index+i] << "\t";
  }
}

void TLinearLaw2D::ValJacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &valjacob) {
  int i,j,m,index;
  double aux,aux1;
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
      aux=0.;
      for(m=0;m<fOrder;m++) {
	     aux1=fabs(fEigVal[m])*fInvEigVect[m*fOrder+j];
	     aux1*=fEigVect[index+m];
	     aux+=aux1;
      }
      valjacob(i,j)=aux;
    }
  }
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
      aux=0.;
      for(m=0;m<fOrder;m++) {
	     aux1=fabs(EigValB(m))*InvEigVectB(m*fOrder+j);
	     aux1*=EigVectB(index+m);
	     aux+=aux1;
      }
      valjacob(fOrder+i,j)=aux;
    }
  }
}
/** Revisar esta muito ruin */
int TLinearLaw2D::ValJacobFlux(TPZVec<REAL> &/*Ui*/,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) {
	return 1;
  int i,j,m,index;
  double aux,aux1;
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
      aux=0.;
      for(m=0;m<fOrder;m++) {
	     aux1=fabs(normal[0]*fEigVal[m])*fInvEigVect[m*fOrder+j];
	     aux1*=fEigVect[index+m];
	     aux+=aux1;
      }
      valjacob(i,j)=aux;
    }
  }
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
      aux=0.;
      for(m=0;m<fOrder;m++) {
	     aux1=fabs(normal[1]*EigValB(m))*InvEigVectB(m*fOrder+j);
	     aux1*=EigVectB(index+m);
	     aux+=aux1;
      }
      valjacob(fOrder+i,j)=aux;
    }
  }
	return 0;
}

void TLinearLaw2D::SetData(std::istream &input) {
  TConservationLaw::SetData(input);
  SetMatrix(input);
}

void TLinearLaw2D::SetMatrix(std::istream &input) {
  int i,index;
  double temp;
  fMaxEigVal=0.;
  fMaxEigValB = 0.;
  //Limpando para armazenar
  for(i=0;i<MAXORDER*MAXORDER;i++) fA[i]=fEigVect[i]=fValA[i]=fInvEigVect[i]=0.;
  for(i=0;i<MAXORDER;i++) fEigVal[i]=0.;
  //Ingressa dados da matriz A
  GetCommentary(input);   // Matrix A and B
  for(i=0;i<2*fOrder*fOrder;i++)
    input >> fA[i];
  //Ingresso dos autovalores, autovetores e a inversa dos autovetores
  for(i=0;i<fOrder;i++) {
    input >> temp;
    fMaxEigVal=Max(fMaxEigVal,fabs(temp));
    fEigVal[i]=temp;
  }
  for(i=0;i<fOrder;i++) {
    input >> temp;
    fMaxEigValB=Max(fMaxEigValB,fabs(temp));
    fEigVal[fOrder+i]=temp;
  }
  for(i=0;i<2*fOrder*fOrder;i++)
    input >> fEigVect[i];
  for(i=0;i<2*fOrder*fOrder;i++)
    input >> fInvEigVect[i];

  std::cout << "\nVerifique :\n Matrix A \t EigVector\t InvEigVect\t EigValue";
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(int j=0;j<fOrder;j++) {
		std::cout << "\n" << fA[index+j] << "\t" << fEigVect[index+j];
		std::cout << "\t" << fInvEigVect[index+j];
    }
	std::cout << "\t" << fEigVal[i];
  }
  std::cout << "\n";
  std::cout << "\nVerifique :\n Matrix B \t EigVector\t InvEigVect\t EigValue";
  for(i=fOrder;i<2*fOrder;i++) {
    index=i*fOrder;
    for(int j=0;j<fOrder;j++) {
		std::cout << "\n" << fA[index+j] << "\t" << fEigVect[index+j];
		std::cout << "\t" << fInvEigVect[index+j];
    }
	std::cout << "\t" << fEigVal[i];
  }
  std::cout << "\n";
}

void TLinearLaw2D::SetMatrix(TPZFMatrix<STATE> &A,TPZFMatrix<STATE> &EigVect,TPZFMatrix<STATE> &InvEigVect,
			   TPZVec<REAL> &EigVal) {
  int i,j;
  double temp;
  fMaxEigVal=0.;
  fMaxEigValB = 0.;
#ifdef DEBUG
  if(fOrder!=A.Rows() || fOrder!=EigVect.Rows() || fOrder!=InvEigVect.Rows() ||
     fOrder!=EigVal.NElements()) {
    PZError << "TLinearLaw2D1D::SetMatrix. Some matrix is out of the dimension.\n";
    return;
  }
#endif
  for(i=0;i<2*fOrder;i++)
    for(j=0;j<fOrder;j++)
      fA[i*fOrder+j] = A(i,j);

  //Ingresso dos autovalores, autovetores e a inversa dos autovetores
  for(i=0;i<fOrder;i++) {
    temp = EigVal[i];
    fMaxEigVal=Max(fMaxEigVal,fabs(temp));
    fEigVal[i]=temp;
  }
  for(i=fOrder;i<2*fOrder;i++) {
    temp = EigVal[i];
    fMaxEigValB=Max(fMaxEigValB,fabs(temp));
    fEigVal[fOrder+i]=temp;
  }
  for(i=0;i<2*fOrder;i++)
    for(j=0;j<fOrder;j++)
      fEigVect[i*fOrder+j] = EigVect(i,j);
  for(i=0;i<2*fOrder;i++)
    for(j=0;j<fOrder;j++)
      fInvEigVect[i*fOrder+j] = InvEigVect(i,j);
}

void TLinearLaw2D::RequireMatrix() {
  int i;
  std::cout << "\nUtilizara o teclado? (y/n) ";
  char c[5];
  std::cin >> c[0];
  if(c[0]=='y') {
    fMaxEigVal=0.;
    fMaxEigVal = 0.;
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
	std::cout << "5- Matrix B,\n";
    int index = fOrder*fOrder;
    for(i=index;i<2*index;i++)
		std::cin >> fA[i];
	std::cout << "6- Autovalores de B,\n";
    for(i=fOrder;i<2*fOrder;i++) {
		std::cin >> fEigVal[i];
      fMaxEigValB=Max(fMaxEigValB,fabs(fEigVal[i]));
    }
	std::cout << "7- Matriz R dos autovetores de B,\n";
    for(i=index;i<2*index;i++)
		std::cin >> fEigVect[i];
	std::cout << "8- Matriz inversa de R,\n";
    for(i=index;i<2*index;i++)
		std::cin >> fInvEigVect[i];
	std::cout << "Dados Completos.\n";
  }
  else {
	  std::cout << "\nNome do arquivo com dados relacionados a matriz A e B ";
    char Amatrix[20];
	std::cin >> Amatrix;
	std::ifstream in(Amatrix);
    SetMatrix(in);
    in.close();
  }
}

// *************************   TLinearLaw2DCircle   ***************************

TLinearLaw2DCircle::TLinearLaw2DCircle(int id) : TLinearLaw2D(id,1) {
  SetName("Linear conservation law 2D - Circular.");
}
TLinearLaw2DCircle::TLinearLaw2DCircle(std::istream &input) : TLinearLaw2D(input) {
  SetName("Linear conservation law 2D - Circular.");
  fOrder = 1;
  GetDataCommented(input,fAngle);
  SetData(fAngle);
  fNumericalFlux.SetOrder(1);
  fNumericalFlux.SetDimension(2);
}
TLinearLaw2DCircle::TLinearLaw2DCircle(int id,REAL angle,char *name,int type) :
	TLinearLaw2D(id,0,name,type) {
  fOrder = 1;
  fNumericalFlux.SetOrder(1);
  SetData(angle);
}

inline void TLinearLaw2DCircle::SetData(std::istream &input) {
  TConservationLaw::SetData(input);
  GetDataCommented(input,fAngle);
  SetData(fAngle);
}

void TLinearLaw2DCircle::SetData(REAL angle) {
  int i,index;
  fAngle = (angle*PI_Value)/180.;
  fMaxEigVal=0.;
  fMaxEigValB = 0.;
  //Limpando para armazenar
  for(i=0;i<MAXORDER*MAXORDER;i++) fA[i]=fEigVect[i]=fValA[i]=fInvEigVect[i]=0.;
  for(i=0;i<MAXORDER;i++) fEigVal[i]=0.;
  //Ingressa dados da matriz A
  fEigVal[0] = fA[0] = cos(fAngle);
  fEigVal[1] = fA[1] = sin(fAngle);
  //Ingresso dos autovalores, autovetores e a inversa dos autovetores
  fMaxEigVal = fA[0];
  fMaxEigValB = fA[1];
  fEigVect[0] = fEigVect[1] = 1.;
  fInvEigVect[0] = fInvEigVect[1] = 1.;

  std::cout << "\nVerifique :\n Matrix A \t EigVector\t InvEigVect\t EigValue";
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(int j=0;j<fOrder;j++) {
		std::cout << "\n" << fA[index+j] << "\t" << fEigVect[index+j];
		std::cout << "\t" << fInvEigVect[index+j];
    }
	std::cout << "\t" << fEigVal[i];
  }
  std::cout << "\n";
  std::cout << "\nVerifique :\n Matrix B \t EigVector\t InvEigVect\t EigValue";
  for(i=fOrder;i<2*fOrder;i++) {
    index=i*fOrder;
    for(int j=0;j<fOrder;j++) {
		std::cout << "\n" << fA[index+j] << "\t" << fEigVect[index+j];
		std::cout << "\t" << fInvEigVect[index+j];
    }
	std::cout << "\t" << fEigVal[i];
  }
  std::cout << "\n";
}

/*
void TLinearLaw2D::EigRoeMatrix(TPZVec<REAL> &Ul,TPZVec<REAL> &Ur,TPZVec<REAL> &EigRoe) {
  int i,j;
  for(i=0;i<fOrder;i++) EigRoe[i]=fEigVal[i];
  for(i=0;i<fOrder*fOrder;i++) EigRoe[fOrder+i]=fEigVect[i];
  j=fOrder*(1+fOrder);
  for(i=0;i<fOrder*fOrder;i++) EigRoe[j+i]=fInvEigVect[i];
}

void TLinearLaw2D::RoeMatrix(TPZVec<REAL> &Ul,TPZVec<REAL> &Ur,TPZFMatrix &Roe) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      Roe(i,j)=fA[i*fOrder+j];
}

void TLinearLaw2D::ValRoeMatrix(TPZVec<REAL> &Ul,TPZVec<REAL> &Ur,TPZFMatrix &ValRoe) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      ValRoe(i,j)=fValA[i*fOrder+j];
}
*/



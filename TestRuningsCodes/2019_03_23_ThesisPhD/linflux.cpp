/*******       File : linflux.c

This file contains the method definitions for class TLinearFlux.
This class is derived from TNumericalFlux and help us to calculate
a numerical flux across one edge of a cell.
In all of cases, flux lenght is fDim*fOrder, to contains flux = (A,B,C)
from   Ut + AUx + BUy + CUz = s.
To numerical flux : BackwardEuler, Godunov, LaxFriedrichs and LaxWendroff
we need  U lenght = 2*fOrder. (To Ul and Ur)
To numerical flux : BeamWarming we need  U lenght = 3*fOrder. (To Ul-1, Ul and Ur)
In the others numerical fluxs we need U lenght = 4*fOrder. (To Ul-1, Ul, Ur and Ur+1)

*******              *******/

#include "pzvec.h"
#include "linflux.h"
#include "conslaw.h"
#include "linlaw.h"
#include "pzfmatrix.h"

//Constructors and destructor
TLinearFlux::TLinearFlux(TLinearLaw *L,int dim) {
  fLaw = L;
  fLaw->SetAlfa(L->Alfa());
  fDim = dim;
  fOrder=L->NStateVariables();
  fA=L->A();
  fValA=L->ValA();
  fEigVal=L->EigVal();
  fEigVect=L->EigVect();
  fInvEigVect=L->InvEigVect();
  RequireFluxType();
}
TLinearFlux::TLinearFlux(TLinearLaw *L,int t,int dim) {
  fLaw=L;
  fLaw->SetAlfa(L->Alfa());
  fDim = dim;
  fOrder=L->NStateVariables();
  fA=L->A();
  fValA=L->ValA();
  fEigVal=L->EigVal();
  fEigVect=L->EigVect();
  fInvEigVect=L->InvEigVect();
  SetFluxType(t);
}
TLinearFlux::TLinearFlux(TConservationLaw *L,double *A,double *EigVal,double *EigVect,double *InvEigVect,int t) {
  fLaw = L;
  fLaw->SetAlfa(L->Alfa());
  fOrder=L->NStateVariables();
  fDim = L->Dimension();
  fA = A;
  fEigVal=EigVal;
  fEigVect=EigVect;
  fInvEigVect=InvEigVect;
  SetFluxType(t);
}
TLinearFlux::TLinearFlux(TConservationLaw *L,double *A,int t) {
  fOrder=L->NStateVariables();
  fDim = L->Dimension();
  fA=A;
  SetFluxType(t);
}
TLinearFlux::TLinearFlux(TConservationLaw *L,double *ValA) {
  fOrder=L->NStateVariables();
  fDim = L->Dimension();
  fValA=ValA;
}

void TLinearFlux :: ShowFluxType() {
 std::cout << "\nFLUXO NUMERICO : \n";
 std::cout << "\nAproximacao De Primeira Ordem :\n\n";
 std::cout << " 0- Backward Euler\n 1- Godunov\n 2- Lax Friedrichs\n";
 std::cout << "\nAproximacao De Segunda Ordem :\n\n 3- Lax Wendroff\n";
 std::cout << " 4- Ritchmyer\n 5- Beam Warming\n";
 std::cout << "\nAlta Resolucao : Flux Limiter\n";
 std::cout << " Godunov-LaxWendroff - Limitador :\n";
 std::cout << " 6- Superbee\t 7- Minimo Modulo\n";
 std::cout << " 8- Van Leer\t 9- Chakravarthy\n";
 std::cout << " LaxFriedrichs-LaxWendroff - Limitador :";
 std::cout << "10- Superbee\t 11- Minimo Modulo\n";
 std::cout << "12- Van Leer\t 13- Chakravarthy\n";
 std::cout << " Godunov-BeamWarming - Limitador :\n";
 std::cout << "14- Superbee\t 15- Minimo Modulo\n";
 std::cout << "16- Van Leer\t 17- Chakravarthy\n";
 std::cout << " LaxFriedrichs-BeamWarming - Limitador :\n";
 std::cout << "18- Superbee\t 19- Minimo Modulo\n";
 std::cout << "20- Van Leer\t 21- Chakravarthy\n";
 std::cout << "\nAlta resolucao : Slope-Limiter\n\n";
 std::cout << "22- Minimo modulo (MinMod)\n";
 std::cout << "\n30- Todos os metodos(exeto Backward Euler)\n";
}

void TLinearFlux::SetFluxType(int t) {
  if(t==30) fFluxType=2;
  else if(t<0 || t>22) {
   std::cout << "\nERRO na eleicao do fluxo\n";
    exit(1);
  }
  else fFluxType=t;

  if(fFluxType==0) {
    fFlux=&TLinearFlux::FluxBackwardEuler;
    SetName("Back");
  }
  else if(fFluxType==1) {
    fFlux=&TLinearFlux::FluxGodunov;
    SetName("Godu");
  }
  else if(fFluxType==2) {
    fFlux=&TLinearFlux::FluxLaxFriedrichs;
    SetName("LaxF");
  }
  else if(fFluxType==3) {
    fFlux=&TLinearFlux ::FluxLaxWendroff;
    SetName("LaxW");
  }
  else if(fFluxType==4) {
    fFlux=&TLinearFlux::FluxRitchmyer;
    SetName("Ritc");
  }
  else if(fFluxType==5) {
    fFlux=&TLinearFlux::FluxBeamWarming;
    SetName("BeaW");
  }
  else if(fFluxType<10) {
    if(fFluxType==6) {
      SetName("GdLWs");
      SetLimiter(ESuperbee);
    }
    else if(fFluxType==7) {
      SetName("GdLWm");
      SetLimiter(EMinMod);
    }
    else if(fFluxType==8) {
      SetName("GdLWv");
      SetLimiter(EVanLeer);
    }
    else if(fFluxType==9) {
      SetName("GdLWc");
      SetLimiter(EChakravarthy);
    }
    fFlux=&TLinearFlux ::FluxGodunovLaxWendroff;
  }
  else if(fFluxType<14) {
    if(fFluxType==10) {
      SetName("LFLWs");
      SetLimiter(ESuperbee);
    }
    else if(fFluxType==11) {
      SetName("LFLWm");
      SetLimiter(EMinMod);
    }
    else if(fFluxType==12) {
      SetName("LFLWv");
      SetLimiter(EVanLeer);
    }
    else if(fFluxType==13) {
      SetName("LFLWc");
      SetLimiter(EChakravarthy);
    }
    fFlux=&TLinearFlux ::FluxLaxFriedrichsLaxWendroff;
  }
  else if(fFluxType<18) {
    if(fFluxType==14) {
      SetName("GdBWs");
      SetLimiter(ESuperbee);
    }
    else if(fFluxType==15) {
      SetName("GdBWm");
      SetLimiter(EMinMod);
    }
    else if(fFluxType==16) {
      SetName("GdBWv");
      SetLimiter(EVanLeer);
    }
    else if(fFluxType==17) {
      SetName("GdBWc");
      SetLimiter(EChakravarthy);
    }
    fFlux=&TLinearFlux::FluxGodunovBeamWarming;
  }
  else if(fFluxType<22) {
    if(fFluxType==18) {
      SetName("LFBWs");
      SetLimiter(ESuperbee);
    }
    else if(fFluxType==19) {
      SetName("LFBWm");
      SetLimiter(EMinMod);
    }
    else if(fFluxType==20) {
      SetName("LFBWv");
      SetLimiter(EVanLeer);
    }
    else if(fFluxType==21) {
      SetName("LFBWc");
      SetLimiter(EChakravarthy);
    }
    fFlux=&TLinearFlux::FluxLaxFriedrichsBeamWarming;
  }
  else if(fFluxType==22) {
    fFlux=&TLinearFlux::FluxSlopeLimiter;
    SetName("Slop");
    fMaxCFL = 0.5;
  }
  else
    SetName("noflux");
}

void TLinearFlux::NumericalFlux(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  (this->*fFlux)(U,Up1,normal,flux);
}

void TLinearFlux::FluxBackwardEuler(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k,kindex,index;
  double soma;
  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      soma=0;
      for(j=0;j<fOrder;j++)
	     soma+=(((TLinearLaw *)fLaw)->A(kindex+index+j)*(U[j]+Up1[j]));
      flux[i] += (.5*soma*normal[k]);
    }
  }
}

void TLinearFlux::FluxGodunov(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k,index,kindex;
  double soma,eigval;
  double Auxiliar[MAXORDER];
  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++) {
      index=i*fOrder;
      eigval=fEigVal[k*fOrder+i];
      soma=0;
      if(eigval<0) {
        for(j=0;j<fOrder;j++)
          soma+=(eigval*Up1[j]*fInvEigVect[kindex+index+j]);
      }
      else {
	for(j=0;j<fOrder;j++)
	  soma+=(eigval*U[j]*(fInvEigVect[kindex+index+j]));
      }
      Auxiliar[i]=normal[k]*soma;
    }
    //Realizando a multiplicacao de R*Auxiliar:
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      soma=0.;
      for(j=0;j<fOrder;j++)
        soma+=((fEigVect[kindex+index+j])*Auxiliar[j]);
      flux[i] += (soma*normal[k]);
    }
  }
}

void TLinearFlux::FluxLaxFriedrichs(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k,kindex;
  TPZFMatrix<STATE> Anormal(fOrder,fOrder,0.);
  REAL alfa = fLaw->MaxEigJacob(U,normal);
  REAL aux = fLaw->MaxEigJacob(Up1,normal);
  if(alfa<aux) alfa = aux;

/*  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      eigval=normal[k]*fEigVal[k*fOrder+i];
      soma=0;
      temp=eigval-(1./(fLaw->Alfa()));
      eigval+=(1./(fLaw->Alfa()));
      for(j=0;j<fOrder;j++)
        soma+=(fInvEigVect[kindex+index+j]*((temp*Up1[j])+(eigval*U[j])));
      Auxiliar[i]=.5*soma;
    }
    //Realizando a multiplicacao de R*Auxiliar:
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      soma=0.;
      for(j=0;j<fOrder;j++)
        soma+=(fEigVect[kindex+index+j]*Auxiliar[j]);
      flux[i] += soma;
    }
  }*/
  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++)
      for(j=0;j<fOrder;j++)
        Anormal(i,j) += normal[k]*fA[kindex+i*fOrder+j];
  }
  for(i=0;i<fOrder;i++) {
    flux[i] = 0.;
    for(j=0;j<fOrder;j++)
      flux[i] += Anormal(i,j)*(U[j]+Up1[j]);
  }
  for(i=0;i<fOrder;i++) {
    flux[i] += (alfa*(U[i]-Up1[i]));
    flux[i] *= .5;
  }
  /**Flux[i] = .5 * [(A,B).n(Up1+U) - 1/alfa*(Up1-U) ] */
}

void TLinearFlux::FluxLaxWendroff(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
/*  int i,j,k,index,kindex;
  double soma,eigval,temp;
  double Auxiliar[MAXORDER];
  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      eigval=fEigVal[k*fOrder+i];
      soma=(fLaw->Alfa())*eigval*eigval;
      temp=eigval-soma;
      eigval+=soma;
      soma=0.;
      for(j=0;j<fOrder;j++)
        soma+=(fInvEigVect[kindex+index+j]*((temp*Up1[j])+(eigval*U[j])));
      Auxiliar[i]=.5*soma;
    }
    //Realizando a multiplicacao de R*Auxiliar:
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      soma=0.;
      for(j=0;j<fOrder;j++)
        soma+=(fEigVect[kindex+index+j]*Auxiliar[j]);
      flux[i] += (soma*normal[k]);
    }
  }*/
  TPZVec<REAL> upartial(fOrder,0.);
  TPZFMatrix<STATE> Anormal(fOrder,fOrder,0.);
  REAL coef = fLaw->Alfa();
  int i, j, k, kindex;
  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++)
      for(j=0;j<fOrder;j++)
        Anormal(i,j) += normal[k]*fA[kindex+i*fOrder+j];
  }
  for(i=0;i<fOrder;i++) {
    flux[i] = 0.;
    for(j=0;j<fOrder;j++)
      upartial[i] += Anormal(i,j)*(U[j]-Up1[j]);
  }
  for(i=0;i<fOrder;i++) upartial[i] *= coef;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      flux[i] += Anormal(i,j)*(Up1[i]+U[i]+upartial[i]);
  for(i=0;i<fOrder;i++) flux[i] *= .5;
  /** flux[i] = .5* [(A,B).n (Up1+U) - alfa * (A,B).n *((A,B).n (Up1-U)) ]*/
}

void TLinearFlux::FluxRitchmyer(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
}
// In this case, we need Ul-1, Ul and Ur, then U lenght is 3*fOrder
void TLinearFlux::FluxBeamWarming(TPZVec<REAL> &U,TPZVec<REAL> &Up1p2,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k,kindex,index;
  double soma,eigval,temp;
  double Auxiliar[MAXORDER];
  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++) {
      index=i*fOrder;
      eigval=ValAbs(fEigVal[k*fOrder+i]);
      soma=(fLaw->Alfa())*eigval*eigval;
      temp=3*eigval-soma;
      eigval=soma-eigval;
      soma=0.;
      //Achando V(i-sp)
      if(fEigVal[k*fOrder+i]<0.) {
	for(j=0;j<fOrder;j++)
	  soma+=(fInvEigVect[kindex+index+j]*((temp*Up1p2[j])+(eigval*Up1p2[j+fOrder])));
      }
      else {
	for(j=0;j<fOrder;j++)
	  soma+=(fInvEigVect[kindex+index+j]*((temp*Up1p2[j])+(eigval*U[j])));
      }
      Auxiliar[i]=.5*soma;
    }
    //Realizando a multiplicacao de R*Auxiliar:
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      soma=0.;
      for(j=0;j<fOrder;j++)
        soma+=(fEigVect[kindex+index+j]*Auxiliar[j]);
      flux[i] += (soma*normal[k]);
    }
  }
}
// To this numerical flux we need Ul-1,Ul, Ur and Ur+1, then U lenght is 4*fOrder
void TLinearFlux::FluxGodunovLaxWendroff(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k,kindex,index;
  double Vi,Vip1,Visgn;
  double eigval,temp;
  double Auxiliar[MAXORDER];

  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++) {
      index=i*fOrder;
      eigval=fEigVal[k*fOrder+i];
      temp=(fLaw->Alfa())*eigval*eigval;
      temp=.5*(ValAbs(eigval)-temp);
      Vi=Vip1=Visgn=0.;
      if(eigval<0) {
	//Multiplicando InvR*Ui, InvR*Uip1, InvR*Uisgn
	for(j=0;j<fOrder;j++) {
	  Vi+=(fInvEigVect[kindex+index+j]*Up1[j]);
	  Vip1+=(fInvEigVect[kindex+index+j]*Up1[j+fOrder]);
	  Visgn+=(fInvEigVect[kindex+index+j]*Up1[j+2*fOrder]);
	}
	Visgn=Limiter((Visgn-Vip1),(Vip1-Vi));
	Visgn*=temp;
	Vi*=(-1.*Visgn);
	Vip1*=(eigval+Visgn);
      }
      else {
	for(j=0;j<fOrder;j++) {
	  Vi+=(fInvEigVect[kindex+index+j]*Up1[j]);
	  Vip1+=(fInvEigVect[kindex+index+j]*Up1[j+fOrder]);
	  Visgn+=(fInvEigVect[kindex+index+j]*U[j]);
	}
	Visgn=Limiter((Vi-Visgn),(Vip1-Vi));
	Visgn*=temp;
	Vip1*=Visgn;
	Vi*=(eigval-Visgn);
      }
      Auxiliar[i]=Vi+Vip1;
    }
    //Realizando a multiplicacao de R*Auxiliar:
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      eigval=0.;
      for(j=0;j<fOrder;j++)
	eigval+=(fEigVect[kindex+index+j]*Auxiliar[j]);
      flux[i] += (eigval*normal[k]);
    }
  }
}
// As last instance, we need Ul-1, Ul, Ur and Ur+1. U lenght is 4*fOrder
void TLinearFlux::FluxLaxFriedrichsLaxWendroff(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
/*  int i,j,k,kindex,index;
  double Vi,Vip1,Visgn;
  double eigval,temp;
  double Auxiliar[MAXORDER];
  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++) {
      index=i*fOrder;
      eigval=fEigVal[k*fOrder+i];
      temp=(fLaw->Alfa())*eigval*eigval;
      temp=(1./(fLaw->Alfa()))-temp;
      Vi=Vip1=Visgn=0.;
      if(normal[k]<0) {
	//Multiplying InvR*Ui, InvR*Uip1, InvR*Uisgn
	for(j=0;j<fOrder;j++) {
//	  Vi+=(fInvEigVect[kindex+index+j]*Up1[j]);
//	  Vip1+=(fInvEigVect[kindex+index+j]*Up1[j+fOrder]);
//	  Visgn+=(fInvEigVect[kindex+index+j]*Up1[j+2*fOrder]);
	  Vi+=(fInvEigVect[kindex+index+j]*U[j]);
	  Vip1+=(fInvEigVect[kindex+index+j]*Up1[j]);
	  Visgn+=(fInvEigVect[kindex+index+j]*Up1[fOrder+j]);
	}
	Visgn=Limiter((Visgn-Vip1),(Vip1-Vi));
      }
      else {
	for(j=0;j<fOrder;j++) {
//	  Vi+=(fInvEigVect[kindex+index+j]*Up1[j]);
//	  Vip1+=(fInvEigVect[kindex+index+j]*Up1[j+fOrder]);
//	  Visgn+=(fInvEigVect[kindex+index+j]*U[j]);
	  Vi+=(fInvEigVect[kindex+index+j]*Up1[j]);
	  Vip1+=(fInvEigVect[kindex+index+j]*U[j]);
	  Visgn+=(fInvEigVect[kindex+index+j]*U[fOrder+j]);
	}
	Visgn=Limiter((Vi-Visgn),(Vip1-Vi));
      }
      Visgn*=temp;
      Vi*=(normal[k]*eigval+(1./(fLaw->Alfa()))-Visgn);
      Vip1*=(normal[k]*eigval-(1./(fLaw->Alfa()))+Visgn);
      Auxiliar[i]=.5*(Vi+Vip1);
    }
    // R*Auxiliar:
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      eigval=0.;
      flux[i] = 0.;
      for(j=0;j<fOrder;j++)
        eigval+=(normal[k]*fEigVect[kindex+index+j]*Auxiliar[j]);
      flux[i] += eigval;
    }
  }*/
  FluxLaxFriedrichs(U,Up1,normal,flux);
  TPZVec<REAL> fluxhigh(fOrder,0.);
  FluxLaxWendroff(U,Up1,normal,fluxhigh);
  int i,k;
  if(U.NElements()==fOrder) return;
  TPZVec<REAL> Visgn(fOrder,0.);
  double eigval;
  for(i=0;i<fOrder;i++) {
    Visgn[i]=0.;
    eigval = 0.;
    for(k=0;k<fDim;k++)
      eigval += normal[k]*fEigVal[k*fOrder+i];
    if(eigval<0)
      Visgn[i]=Limiter(Up1[fOrder+i]-Up1[i],U[i]-U[i+fOrder]);
    else
      Visgn[i]=Limiter(U[i]-U[fOrder+i],Up1[i+fOrder]-Up1[i]);
  }
  for(i=0;i<fOrder;i++)
    flux[i] += Visgn[i]*(fluxhigh[i]-flux[i]);
  /** flux[i] = fluxlower[i] + limiter * ( fluxhigh[i] - fluxlower[i] ) */
}
// We need  Ul-1, Ul, Ur and Ur+1. U lenght is 4*fOrder
void TLinearFlux::FluxGodunovBeamWarming(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k,kindex,index;
  double Vi,Vip1,Visgn,limitador;
  double eigval,temp;
  double Auxiliar[MAXORDER];
  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++) {
      index=i*fOrder;
      eigval=fEigVal[k*fOrder+i];
      temp=(fLaw->Alfa())*eigval*eigval;
      Vi=Vip1=Visgn=limitador=0.;
      if(eigval<0) {
	for(j=0;j<fOrder;j++) {
	  Vi+=(fInvEigVect[kindex+index+j]*Up1[j]);
	  Vip1+=(fInvEigVect[kindex+index+j]*Up1[j+fOrder]);
	  Visgn+=(fInvEigVect[kindex+index+j]*Up1[j+2*fOrder]);
	}
	limitador=Limiter((Visgn-Vip1),(Vip1-Vi));
	Visgn*=(limitador*(temp-ValAbs(eigval)));
	Vip1*=(2*eigval*(1.-limitador));
	Vi*=(limitador*((-3.*eigval)-temp));
      }
      else {
	for(j=0;j<fOrder;j++) {
	  Vi+=(fInvEigVect[kindex+index+j]*Up1[j]);
	  Vip1+=(fInvEigVect[kindex+index+j]*Up1[j+fOrder]);
	  Visgn+=(fInvEigVect[kindex+index+j]*U[j]);
	}
	limitador=Limiter((Vi-Visgn),(Vip1-Vi));
	Visgn*=(limitador*(temp-ValAbs(eigval)));
	limitador*=(eigval-temp);
	Vi*=((2*eigval)+limitador);
	Vip1=0.;
      }
      Auxiliar[i]=.5*(Vi+Vip1+Visgn);
    }
    // multiplying  R*Auxiliar:
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      limitador=0;
      for(j=0;j<fOrder;j++)
        limitador+=(fEigVect[kindex+index+j]*Auxiliar[j]);
      flux[i] += (limitador*normal[k]);
    }
  }
}

void TLinearFlux::FluxLaxFriedrichsBeamWarming(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k,kindex,index;
  double Vi,Vip1,Visgn,limitador;
  double eigval,temp;
  double inv=1./(fLaw->Alfa());
  double Auxiliar[MAXORDER];
  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++) {
      index=i*fOrder;
      eigval=fEigVal[k*fOrder+i];
      temp=(fLaw->Alfa())*eigval*eigval;
      Vi=Vip1=Visgn=limitador=0.;
      if(eigval<0) {
	for(j=0;j<fOrder;j++) {
	  Vi+=(fInvEigVect[kindex+index+j]*Up1[j]);
	  Vip1+=(fInvEigVect[kindex+index+j]*Up1[j+fOrder]);
	  Visgn+=(fInvEigVect[kindex+index+j]*Up1[j+2*fOrder]);
	}
	limitador=Limiter((Visgn-Vip1),(Vip1-Vi));
	Visgn*=(limitador*(temp+eigval));
	Vip1*=((eigval-inv)*(1.-limitador));
	limitador*=((-1)*((4*eigval)+temp+inv));
	Vi*=(limitador+eigval+inv);
      }
      else {
	for(j=0;j<fOrder;j++) {
	  Vi+=(fInvEigVect[kindex+index+j]*Up1[j]);
	  Vip1+=(fInvEigVect[kindex+index+j]*Up1[j+fOrder]);
	  Visgn+=(fInvEigVect[kindex+index+j]*U[j]);
	}
	limitador=Limiter((Vi-Visgn),(Vip1-Vi));
	Visgn*=(limitador*(temp-eigval));
	Vip1*=((eigval-inv)*(1.-limitador));
	limitador*=((2*eigval)-temp-inv);
	Vi*=(inv+eigval+limitador);
      }
      Auxiliar[i]=.5*(Vi+Vip1+Visgn);
    }
    // multiplying  R*Auxiliar:
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      limitador=0;
      for(j=0;j<fOrder;j++)
        limitador+=(fEigVect[kindex+index+j]*Auxiliar[j]);
      flux[i] += (limitador*normal[k]);
    }
  }
}
// We need Ul-1, Ul, Ur and Ur+1
void TLinearFlux::FluxSlopeLimiter(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k,kindex,index;
  double Vi,Vip1,Visgn,MinMod;
  double eigvalue;
  double Auxiliar[MAXORDER];
  for(k=0;k<fDim;k++) {
    kindex = k*fOrder*fOrder;
    for(i=0;i<fOrder;i++) {
      index=i*fOrder;
      eigvalue=fEigVal[k*fOrder+i];
      Vi=Vip1=Visgn=0.;
      if(eigvalue<0.) {
	//Determinando a compoente p-esima de Vi=InvR*Ui, Vip1 e Vip2
	for(j=0;j<fOrder;j++) {
	  Vi+=(fInvEigVect[kindex+index+j]*Up1[j]);
	  Vip1+=(fInvEigVect[kindex+index+j]*Up1[j+fOrder]);
	  Visgn+=(fInvEigVect[kindex+index+j]*Up1[j+2*fOrder]);
	}
	MinMod=minmod((Visgn-Vip1),(Vip1-Vi));
	Visgn=MinMod*(-.5)*eigvalue*(1.+eigvalue);
	Auxiliar[i]=(eigvalue*Vip1)+Visgn;
      }
      else {
	for(j=0;j<fOrder;j++) {
	  Vi+=(fInvEigVect[kindex+index+j]*Up1[j]);
	  Vip1+=(fInvEigVect[kindex+index+j]*Up1[j+fOrder]);
	  Visgn+=(fInvEigVect[kindex+index+j]*U[j]);
	}
	MinMod=minmod((Vip1-Vi),(Vi-Visgn));
	Visgn=MinMod*.5*eigvalue*(1.-eigvalue);
	Auxiliar[i]=(eigvalue*Vi)+Visgn;
      }
    }
    // multiplying R*Auxiliar:
    for(i=0;i<fOrder;i++) {
      index = i*fOrder;
      eigvalue=0;
      for(j=0;j<fOrder;j++)
	eigvalue+=(fEigVect[kindex+index+j]*Auxiliar[j]);
      flux[i] += (eigvalue*normal[k]);
    }
  }
}




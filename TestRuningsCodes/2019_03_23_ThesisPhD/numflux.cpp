/*******       File : numflux.c

This file contains the method definitions for class TNumericalFlux.

*******              *******/

#include "pzfmatrix.h"
#include "numflux.h"
#include "pzvec.h"
#include "conslaw.h"

#include <math.h>

TNumericalFlux::TNumericalFlux(TConservationLaw *L,int type) {
  fLaw=L;
  SetName("Numerical Flux to ");
  fName = L->Name();
  SetFluxType(type);
  fMaxCFL = 1.;
}
TNumericalFlux::TNumericalFlux(TConservationLaw *L,char *name,int type) {
  fLaw=L;
  SetName(name);
  SetFluxType(type);
  fMaxCFL = 1.;
}
/*TNumericalFlux::TNumericalFlux(TPZCompMesh *mesh,int t) {
  fLaw=(TConservationLaw *)mesh->MaterialVec()[0];
  fDim = 1;
  SetName("Numerical Flux to ");
  strcpy(fName,fLaw->Name());
  fOrder = fLaw->Order();
  SetFluxType(t);
  fMaxCFL = 1.;
} */
TNumericalFlux::TNumericalFlux(int order,int t) {
  fLaw=NULL;
  fDim = 1;
  SetName("No Name to Numerical Flux");
  fOrder = order;
  SetFluxType(t);
  fMaxCFL = 1.;
}

void TNumericalFlux::SetName(char *name) {
  fName = name;
}

void TNumericalFlux::ShowFluxType() {
 std::cout << "\nNUMERICAL FLUX : \n";
 std::cout << "First order approximation :\n";
 std::cout << " 1- Godunov \n 2- Lax Friedrichs\n ";
 std::cout << "Second order approximation :\n 3- Lax Wendroff\n";
 std::cout << " 4- Richtmyer-Lax Wendroff\n 5- Mac Cormack \n";
 std::cout << "High resolution : Flux Limiter\n";
 std::cout << " Godunov-LaxWendroff - Limiter :\n";
 std::cout << " 6- Superbee\t 7- Minime Module\n";
 std::cout << " 8- Van Leer\t 9- Chakravarthy\n";
 std::cout << " LaxFriedrichs-LaxWendroff - Limiter :";
 std::cout << "10- Superbee\t 11- Minime Module\n";
 std::cout << "12- Van Leer\t 13- Chakravarthy\n";
 std::cout << " Godunov-BeamWarming - Limiter :\n";
 std::cout << "14- Superbee\t 15- Minime Module\n";
 std::cout << "16- Van Leer\t 17- Chakravarthy\n";
 std::cout << " LaxFriedrichs-BeamWarming - Limiter :\n";
 std::cout << "18- Superbee\t 19- Minime Module\n";
 std::cout << "20- Van Leer\t 21- Chakravarthy\n";
 std::cout << "\nHigh resolution : Slope-Limiter\n\n";
 std::cout << "22- Minime module (MinMod)\n";
 std::cout << "\n Schemes Upstream :\n\n";
 std::cout << "23- Straightforward\n24- VanLeer\n25- Roe\n";
 std::cout << "26- ENO\n27- ENO to non-uniform mesh\n";
 std::cout << "28- Flux limiter with one minime limiter\n";
 std::cout << "29- Flux limiter with one maxime limiter\n";
 std::cout << "30- All schemes (not Godunov)\n";
}

void TNumericalFlux::RequireFluxType() {
  ShowFluxType();
 std::cout << "\n To Quit (0).\n\n  Choice = ";
  int t;
  std::cin >> t;
  if(t==0) fFluxType=0;
  else SetFluxType(t);
}

void TNumericalFlux::SetFluxType(int t) {
  if(t==30)  fFluxType=2;
  else if(t<0 || t>29) {
   std::cout << "\nERRO na eleicao do fluxo\n";
    exit(1);
  }
  else  fFluxType=t;

  if(fFluxType==1) {
    if(fOrder==1)
      fFlux = &TNumericalFlux::FluxGodunov1;
    else {
      PZError << "TNumericalFlux::SetFluxType. Flux Godunov required but order != 1.\n";
      fFlux = 0;
    }
    SetName("Godu");
  }
  else if(fFluxType==2) {
    fFlux = &TNumericalFlux::FluxLaxFriedrichs;
    SetName("LaxF");
  }
  else if(fFluxType==3) {
    fFlux = &TNumericalFlux::FluxLaxWendroff;
    SetName("LaxW");
  }
  else if(fFluxType==4) {
    fFlux = &TNumericalFlux::FluxRichtmyerLW;
    SetName("Rich");
  }
  else if(fFluxType==5) {
    fFlux = &TNumericalFlux::FluxMacCormack;
    SetName("MacC");
  }
  else if(fFluxType<22) {
    if(fFluxType==6)
      SetName("GdLWs");
    else if(fFluxType==7)
      SetName("GdLWm");
    else if(fFluxType==8)
      SetName("GdLWv");
    else if(fFluxType==9)
      SetName("GdLWc");
    else if(fFluxType==10)
      SetName("LFLWs");
    else if(fFluxType==11)
      SetName("LFLWm");
    else if(fFluxType==12)
      SetName("LFLWv");
    else if(fFluxType==13)
      SetName("LFLWc");
    else if(fFluxType==14)
      SetName("GdBWs");
    else if(fFluxType==15)
      SetName("GdBWm");
    else if(fFluxType==16)
      SetName("GdBWv");
    else if(fFluxType==17)
      SetName("GdBWc");
    else if(fFluxType==18)
      SetName("LFBWs");
    else if(fFluxType==19)
      SetName("LFBWm");
    else if(fFluxType==20)
      SetName("LFBWv");
    else if(fFluxType==21)
      SetName("LFBWc");
    fFlux = &TNumericalFlux::FluxLimiter;
  }
  else if(fFluxType==22) {
    fFlux = &TNumericalFlux::FluxSlopeLimiter;
    if(fOrder==1)
     std::cout << "CUIDADO, Verificar, f'(u) nao muda de sinal\n";
    SetName("Slop");
    fMaxCFL = .5;
  }
  else if(fFluxType==23) {
    fFlux = &TNumericalFlux::FluxUpStraightforward;
    SetName("UpSt");
  }
  else if(fFluxType==24) {
    fFlux = &TNumericalFlux::FluxUpVanLeer;
    SetName("UpVL");
  }
  else if(fFluxType==25) {
    fFlux = &TNumericalFlux::FluxUpRoe;
    SetName("UpRo");
  }
  else if(fFluxType==26) {
    fFlux = &TNumericalFlux::FluxENO;
    SetName("ENO1");
  }
  else if(fFluxType==27) {
    fFlux = &TNumericalFlux::FluxENO_NU;
    SetName("ENOnu");
  }
  else if(fFluxType==28) {
    fFlux = &TNumericalFlux::FluxLFLWMin;
    SetName("LiMin");
    SetLimiter(EVanLeer);
  }
  else if(fFluxType==29) {
    fFlux = &TNumericalFlux::FluxLFLWMax;
    SetName("LiMax");
    SetLimiter(EVanLeer);
  }
  else
    SetName("noflux");
}

void TNumericalFlux::NumericalFlux(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  /*	int j;
	double ULocal[4*MaxOrdem];
	int jmin = ((s-1) >= 0) ? -1 : 0;
	int jmax = ((*fL).NCelules()-s) > 3 ? 3 : (*fL).NCelules()-s;
	for(j=jmin;j<jmax;j++)
   	(*fL).Transform(Ui+j*fOrder,(*fL).Material(s+j),ULocal+(j+1)*fOrder,(*fL).Material(s));
	(this->*fFlux)(ULocal+fOrder,Fluxi,(*fL).Material(s));
	if((*fL).Material(s+1)!=(*fL).Material(s)) {
   	for(j=jmin;j<jmax;j++)
      	(*fL).Transform(Ui+j*fOrder,(*fL).Material(s+j),ULocal+(j+1)*fOrder,(*fL).Material(s+1));
	(this->*fFlux)(ULocal+fOrder,Fluxi+fOrder,(*fL).Material(s+1));
	}
	else {
   	for(j=0;j<fOrder;j++)
      	Fluxi[j+fOrder] = Fluxi[j];
	}
  */
  (this->*fFlux)(U,Up1,normal,flux);
}
void TNumericalFlux::FluxLaxFriedrichs(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {

  int i,j;
  TPZVec<REAL> flux1(fDim*fOrder);
  TPZVec<REAL> flux2(fDim*fOrder);

  double aux,alfa;
  fLaw->Flux(U,flux1);
  fLaw->Flux(Up1,flux2);
  for(i=0;i<fOrder;i++) {
    flux[i] = 0.;
    for(j=0;j<fDim;j++)
      flux[i] += normal[j]*(flux1[j*fOrder+i]+flux2[j*fOrder+i]);
  }

  alfa = fLaw->MaxEigJacob(U,normal);
  aux = fLaw->MaxEigJacob(Up1,normal);
  if(alfa<aux) alfa = aux;

  for(i=0;i<fOrder;i++) {
    flux[i]+= alfa*(U[i]-Up1[i]);
    flux[i]*=.5;
  }
  //Flux=.5*((*fL).MaxEigVal()*(u[i]-u[i+1])+(Funcao(u[i])+Funcao(u[i+1])))
}

void TNumericalFlux::FluxUpStraightforward(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k;
  TPZVec<REAL> flux1(fDim*fOrder);
  TPZVec<REAL> flux2(fDim*fOrder);
  TPZVec<REAL> umid(fOrder,0.);
  TPZFMatrix<STATE> jacobtemp(fOrder*fDim,fOrder,0.);
  TPZFMatrix<STATE> valjacob(fOrder,fOrder,0.);

  fLaw->Flux(U,flux1);
  fLaw->Flux(Up1,flux2);
  for(i=0;i<fOrder;i++) {
    flux[i] = 0.;
    for(j=0;j<fDim;j++)
      flux[i] += normal[j]*(flux1[j*fOrder+i]+flux2[j*fOrder+i]);
  }

  for(i=0;i<fOrder;i++) umid[i]=.5*(U[i]+Up1[i]);
  fLaw->ValJacobFlux(umid,jacobtemp);
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      for(k=0;k<fDim;k++)
        valjacob(i,j) += jacobtemp(k*fOrder+i,j);   //jacob*normal

  for(i=0;i<fOrder;i++) {   //umid armazena valjacob*(U-Up1)
    umid[i] = 0.;
    for(j=0;j<fOrder;j++)
      umid[i] += valjacob(i,j)*(U[j]-Up1[j]);
  }

  //Computando o fluxo numerico = .5((flux1+flux2)+valjacob*(U-Up1))
  for(i=0;i<fOrder;i++) {
    flux[i]+= umid[i];
    flux[i]*=.5;
  }
  // flux[j] = .5 * (Sk-> ValJacobk(.5*(U[j]+Up1[j]))*(u[j]-Up1[j]) + nk*Fk(U[j]) + nk*Fk(Up1[j]))
}

void TNumericalFlux::FluxUpVanLeer(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k;
  TPZVec<REAL> flux1(fDim*fOrder);
  TPZVec<REAL> flux2(fDim*fOrder);
  REAL temp;

  TPZVec<REAL> aux(fOrder,0.);
  TPZFMatrix<STATE> jacob1(fDim*fOrder,fOrder);
  TPZFMatrix<STATE> jacob2(fOrder*fDim,fOrder);
  TPZFMatrix<STATE> jacob(fOrder,fOrder,0.);

  fLaw->Flux(U,flux1);
  fLaw->Flux(Up1,flux2);

  fLaw->ValJacobFlux(U,jacob1);
  fLaw->ValJacobFlux(Up1,jacob2);
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      for(k=0;k<fDim;k++)
        jacob(i,j) += .5*(jacob1(i+k*fOrder,j)+jacob2(i+k*fOrder,j));

  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      aux[i] += jacob(i,j)*(U[j]-Up1[j]);

  for(i=0;i<fOrder;i++) {
    temp = 0.;
    for(k=0;k<fDim;k++)
      temp += normal[k]*(flux1[i+k*fOrder]+flux2[i+k*fOrder]);
    flux[i] = .5*(aux[i] + temp);
  }
  // flux = .5 * [Sk -> ( (u-up1)*.5*(jacob(u)+jacob(up1)) + nk*(flux(u) + flux(up1)) ) ]
}

void TNumericalFlux::FluxUpRoe(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k;
  TPZVec<REAL> flux1(fOrder*fDim);
  TPZVec<REAL> flux2(fOrder*fDim);
  REAL temp;
  TPZVec<REAL> aux(fOrder,0.);
  TPZFMatrix<STATE> ValRoe(fDim*fOrder,fOrder);
  fLaw->Flux(U,flux1);
  fLaw->Flux(Up1,flux2);

  fLaw->ValRoeMatrix(U,Up1,ValRoe);

  for(i=0;i<fOrder;i++)
    for(k=0;k<fDim;k++) {
      for(j=0;j<fOrder;j++)
        aux[i] += ValRoe(i,j)*(U[j]-Up1[j]);
    }

  for(i=0;i<fOrder;i++) {
    temp = 0.;
    for(k=0;k<fDim;k++)
      temp += normal[k]*(flux1[i+k*fOrder]+flux2[i+k*fOrder]);
    flux[i] = .5*(aux[i] + temp);
  }
}

void TNumericalFlux::FluxLaxWendroff(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int i,j,k;
  TPZVec<REAL> fluxtemp(fOrder*fDim);
  TPZVec<REAL> flux1(fOrder,0.);
  TPZVec<REAL> flux2(fOrder,0.);
  TPZVec<REAL> umid(fOrder,0.);
  TPZFMatrix<STATE> jacobtemp(fOrder*fDim,fOrder);
  TPZFMatrix<STATE> jacob(fOrder,fOrder,0.);

  /**Computando flux1.normal, flux2.normal e jacob.normal*/
  fLaw->Flux(U,fluxtemp);
  for(i=0;i<fOrder;i++)
    for(k=0;k<fDim;k++)
      flux1[i] += fluxtemp[k*fOrder+i]*normal[k];   //flux1*normal
  fLaw->Flux(Up1,fluxtemp);
  for(i=0;i<fOrder;i++)
    for(k=0;k<fDim;k++)
      flux2[i] += fluxtemp[k*fOrder+i]*normal[k];   //flux2*normal

  for(i=0;i<fOrder;i++) umid[i]=.5*(U[i]+Up1[i]);
  fLaw->JacobFlux(umid,jacobtemp);
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      for(k=0;k<fDim;k++)
        jacob(i,j) += normal[k]*jacobtemp(k*fOrder+i,j);   //jacob*normal

  REAL alfa = fLaw->Alfa();
  for(i=0;i<fOrder;i++) {   //umid armazena jacob*(flux1-flux2)
    umid[i] = 0.;
    for(j=0;j<fOrder;j++)
      umid[i] += jacob(i,j)*(flux1[j]-flux2[j]);
  }

  //Computando o fluxo numerico = .5((flux1+flux2)+alfa*jacob*(flux1-flux2))
  for(i=0;i<fOrder;i++)
    flux[i] = .5*(flux1[i]+flux2[i]+(alfa*umid[i]));
  //coef=Alfa*Jacob(.5*(u[i+1]+u[i]))
  //result=.5*((Funcao(u[i+1])+Funcao(u[i]))+coef*(Funcao(u[i])-Funcao(u[i+1]))))
}

void TNumericalFlux::FluxLFLWMin(TPZVec<REAL> &/*U*/,TPZVec<REAL> &/*Up1*/,TPZVec<REAL> &/*normal*/,TPZVec<REAL> &/*flux*/) {
 std::cout << "\nNo pode ser chamado.";
  return;
  /*	double deno,lim,argument1,argument2;
	int i;
	double FluxHigh[MaxOrdem];

	deno=(*fL).fRM[s+1]-(*fL).fRM[s];
	argument1=(*fL).fRM[s]-(*fL).fRM[s-1];
	argument2=(*fL).fRM[s+2]-(*fL).fRM[s+1];
	if(ValAbs(argument1)>ValAbs(argument2))
   	lim=Limitador(argument2,deno,s);
	else lim=Limitador(argument1,deno,s);
	FluxLaxFriedrichs(Ui,Fluxi,s);
	FluxLaxWendroff(Ui,FluxHigh,s);
	for(i=0;i<fOrder;i++)
   	Fluxi[i]+=(lim*(FluxHigh[i]-Fluxi[i]));
  for(j=0;j<fOrder;j++)
    for(k=0;k<fDim;k++) flux[j] = normal[k]*flux1[k*fOrder+j];*/
}

void TNumericalFlux::FluxLFLWMax(TPZVec<REAL> &/*U*/,TPZVec<REAL> &/*Up1*/,TPZVec<REAL> &/*normal*/,TPZVec<REAL> &/*flux*/) {
 std::cout << "\nNo pode ser chamado.";
  return;
  /*	double deno,lim,argument1,argument2;
	int i;
	double FluxHigh[MaxOrdem];

	deno=(*fL).fRM[s+1]-(*fL).fRM[s];
	argument1=(*fL).fRM[s]-(*fL).fRM[s-1];
	argument2=(*fL).fRM[s+2]-(*fL).fRM[s+1];
	if(ValAbs(argument1)>ValAbs(argument2))
   	lim=Limitador(argument1,deno,s);
	else lim=Limitador(argument2,deno,s);
	FluxLaxFriedrichs(Ui,Fluxi,s);
	FluxLaxWendroff(Ui,FluxHigh,s);
	for(i=0;i<fOrder;i++)
   	Fluxi[i]+=(lim*(FluxHigh[i]-Fluxi[i]));
  for(j=0;j<fOrder;j++)
    for(k=0;k<fDim;k++) flux[j] = normal[k]*flux1[k*fOrder+j];
*/
}

void TNumericalFlux::FluxRichtmyerLW(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int j,k,m;
  TPZVec<REAL> flux1(fOrder*fDim);
  TPZVec<REAL> flux2(fOrder*fDim);
  TPZVec<REAL> utemp(fOrder,0.);

  fLaw->Flux(U,flux1);
  fLaw->Flux(Up1,flux2);

  /** Compute utemp = .5*(u + up1 + alfa*(F(u) - F(up1))) */
  for(k=0;k<fDim;k++) {
    m = k*fOrder;
    for(j=0;j<fOrder;j++)
      utemp[j] += normal[k]*(flux1[m+j] - flux2[m+j]);
  }
  double alfa = fLaw->Alfa();
  for(j=0;j<fOrder;j++) {
    utemp[j] *= alfa;
    utemp[j] += (U[j] + Up1[j]);
    utemp[j] *= .5;
    flux[j] = 0.;
  }
  /** Compute F(utemp) and flux = normal*F */
  fLaw->Flux(utemp,flux2);
  for(k=0;k<fDim;k++) {
    m = k*fOrder;
    for(j=0;j<fOrder;j++)
      flux[j] += normal[k]*flux2[m+j];
  }
  //result=Funcao[.5*(u[i]+u[i+1]+Alfa*(Funcao(u[i])-Funcao(u[i+1])))]
}

void TNumericalFlux ::FluxMacCormack(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  int j,k,m;
  TPZVec<REAL> flux1(fOrder*fDim);
  TPZVec<REAL> flux2(fOrder*fDim);
  TPZVec<REAL> utemp(fOrder,0.);

  fLaw->Flux(U,flux1);
  fLaw->Flux(Up1,flux2);

  /** Compute utemp = .5*(u + up1 + alfa*(F(u) - F(up1))) */
  for(k=0;k<fDim;k++) {
    m = k*fOrder;
    for(j=0;j<fOrder;j++)
      utemp[j] += normal[k]*(flux1[m+j] - flux2[m+j]);
  }
  double alfa = fLaw->Alfa();
  for(j=0;j<fOrder;j++) {
    utemp[j] *= alfa;
    utemp[j] += U[j];
    flux[j] = 0.;
  }
  /** Compute F(utemp) and flux = normal*F */
  fLaw->Flux(Up1,flux1);
  fLaw->Flux(utemp,flux2);
  for(k=0;k<fDim;k++) {
    m = k*fOrder;
    for(j=0;j<fOrder;j++)
      flux[j] += .5*normal[k]*(flux1[m+j]+flux2[m+j]);
  }
  //utemp = u[i]+((*fL).Alfa())*((*fL).FuncaoFlux(u[i])-(*fL).FuncaoFlux(u[i+1]));
  //return(.5*((*fL).FuncaoFlux(u[i+1])+(*fL).FuncaoFlux(utemp)));
}

double TNumericalFlux ::Limiter(double num,double den) {
  double ratio;
  if(!IsZero(den)) {
    ratio = num/den;
    if(fLimiter==EVanLeer) {
      if(ratio>0.)
        return((2.*ratio)/(1.+ratio));
      else return(0.);
    }
    else if(fLimiter==EChakravarthy) {
      double phi;
      phi=2.;
      if(ratio>0.)
        return(Min(ratio,phi));
      else return(0.);
    }
    else {
      double phi,aux; // 1. <= phi <= 2.
      if(fLimiter==ESuperbee) phi=2.;
      else phi=1.;	//Quando fLimiter==EMinMod
      if(ratio>0.) {
        aux=Min(1.,phi*ratio);
        phi=Min(phi,ratio);
        return(Max(aux,phi));
      }
      else return(0.);
    }
  }
  return 0.;
}

#include "linflux.h"

void TNumericalFlux::FluxLimiter(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  // We construct the Roe matrix : [ Â(Ul,Ur) ]. In Roe we store
  // eigenvalues of A,B and C, Roe matrix eigenvectors RA,RB and RC and Rs inverses.
  TPZVec<REAL> Roe(120);
  fLaw->EigRoeMatrix(U,Up1,Roe);
  // Solve as a linear system, to Roe matrix data
  int k = fOrder*fOrder;
  TLinearFlux LF(fLaw,&Roe[0],&Roe[fDim*k],&Roe[fDim*(k+fOrder)],&Roe[fDim*(2*k+fOrder)],fFluxType);
  LF.NumericalFlux(U,Up1,normal,flux);
}

void TNumericalFlux::FluxSlopeLimiter(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  /*	if(fOrder==1 && fDim==1) {
	double den=U[2]-U[1];
	TFMatrix partial(1,1);
	fLaw->JacobFlux(&U[1],partial);
	if(partial(0,0)>0.) {
	fLaw->Flux(&U[1],flux);
	if(!IsZero(den)) {
	TPZVec<REAL> temp(1,U[1]-U[0]);
	partial=minmod(den,temp[0]);
	fLaw->Flux(&U[2],temp);
	temp[0]-=flux[0];
	temp[0]=temp[0]/den;
	den=(fLaw->Alfa())*temp[0];
	den=1.-den;
	den*=(.5*temp[0]);
	partial(0,0)*=den;
	flux[0]+=partial(0,0);
	}
	}
	else {
	fLaw->Flux(&U[2],flux);
	if(!IsZero(den)) {
	TPZVec<REAL> temp(1,U[3]-U[2]);
	partial=minmod(temp[0],den);
	fLaw->Flux(&U[1],temp);
	temp[0]=flux[0]-temp[0];
	temp[0]/=den;
	den=(fLaw->Alfa())*temp[0];
	den+=1.;
	den*=(.5*temp[0]);
	partial(0,0)*=den;
	flux[0]-=partial(0,0);
	}
	}
	}
	else {*/
  //Linearizamos a funcao fluxo na posicao i-esima
  //Geramos a matrix de Roe para Ul e Ur nesta celula : Â(Ul,Ur)
  TPZVec<REAL> Roe(120);     //Dados da matriz Roe
  fLaw->EigRoeMatrix(U,Up1,Roe);
  //Resolvemos em esta celula, como se fosse um sistema linear,
  //para a matrix de Roe obtida
  int k = fOrder*fOrder;
  TLinearFlux LF(fLaw,&Roe[0],&Roe[k],&Roe[k+fOrder],&Roe[2*k+fOrder],2);
  LF.NumericalFlux(U,Up1,normal,flux);
  //	}
}
//We need here  Ul-1, Ul, Ur and Ur+1
void TNumericalFlux::FluxENO(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  double soma;
  TPZVec<REAL> NewU(fOrder,0.);
  TPZVec<REAL> NewUp1(fOrder,0.);
  double MinModule[2*MAXORDER];
  TPZFMatrix<STATE> NewJacob(fOrder,fOrder);
  int j,k;
  if(fDim==1) {
    for(k=0;k<fOrder;k++) {
      NewU[k]=U[k];
      NewUp1[k]=Up1[k];
    }
    for(k=0;k<fOrder;k++) {
      MinModule[k]=minmod(U[k]-U[k+fOrder],Up1[k+fOrder]-Up1[k]);
      MinModule[k]*=.5;
      NewU[k]+=MinModule[k];
      MinModule[k+fOrder]=minmod(Up1[k+fOrder]-Up1[k],U[k]-U[k+fOrder]);  //A toa
      MinModule[k+fOrder]*=.5;													//A toa
      NewUp1[k]-=MinModule[k+fOrder];
    }
    fLaw->JacobFlux(NewU,NewJacob);
    for(k=0;k<fOrder;k++) {
      soma=0.;
      for(j=0;j<fOrder;j++)
	soma += normal[0]*(NewJacob(k,j)*MinModule[j]);
      soma*=fLaw->Alfa();
      NewU[k]-=soma;
    }
    fLaw->JacobFlux(NewUp1,NewJacob);
    for(k=0;k<fOrder;k++) {
      soma=0.;
      for(j=0;j<fOrder;j++)
        soma += normal[0]*(NewJacob(k,j)*MinModule[j+fOrder]);
      soma*=fLaw->Alfa();
      NewUp1[k]-=soma;
    }
    FluxLaxFriedrichs(NewU,NewUp1,normal,flux);
  }
  else std::cout << "TNumericalFlux::FluxENO is yet implemented.\n";
}

void TNumericalFlux::FluxENO_NU(TPZVec<REAL> &/*U*/,TPZVec<REAL> &/*Up1*/,TPZVec<REAL> &/*normal*/,TPZVec<REAL> &/*flux*/) {
 std::cout << "\nNo pode ser chamado.";
  return;
  /*	double soma,num1,num2;
	double hm1,h,hp1,hp2;
	double NewU[2*MaxOrdem];
	double NewJacob[MaxOrdem*MaxOrdem];
	double MinModule[2*MaxOrdem];
	int j,k;
	for(k=0;k<2*fOrder;k++) NewU[k]=Ui[k];
	hm1=fLaw->Area(s-1);
	h=fLaw->Area(s);
	hp1=fLaw->Area(s+1);
	hp2=fLaw->Area(s+2);
	for(k=0;k<fOrder;k++) {
	num1=Ui[k]-Ui[k-fOrder];
	num1/=(hm1+h);
	num2=Ui[k+fOrder]-Ui[k];
	num2/=(hp1+h);
	MinModule[k]=minmod(num1,num2);
	NewU[k]+=(MinModule[k]*h);
	j=k+fOrder;
	num1=Ui[j+fOrder]-Ui[j];
	num1/=(hp1+hp2);
	num2=Ui[j]-Ui[k];
	num2/=(hp1+h);
	MinModule[j]=minmod(num1,num2);
	NewU[j]-=(MinModule[j]*hp1);
	}
	fLaw->JacobFlux(NewU,NewJacob,s);
	for(k=0;k<fOrder;k++) {
	soma=0.;
	for(j=0;j<fOrder;j++)
	soma+=(NewJacob[k*fOrder+j]*MinModule[j]);
	soma*=(h*fLaw->Alfa());
	NewU[k]-=soma;
	}
	fLaw->JacobFlux(NewU+fOrder,NewJacob,s+1);
	for(k=0;k<fOrder;k++) {
	soma=0.;
	for(j=0;j<fOrder;j++)
	soma+=(NewJacob[k*fOrder+j]*MinModule[j+fOrder]);
	soma*=(h*fLaw->Alfa());
	NewU[k+fOrder]-=soma;
	}
	FluxUpRoe(NewU,Fluxi,s);
  */
}

//#include "burger1d.h"
void TNumericalFlux::FluxGodunov1(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux) {
  if(fDim!=1 || fOrder!=1) {
   std::cout << "TNumerical::FluxGodunov1 is implemented only to dim=1 and order=1.\n";
    return;
  }
  TPZFMatrix<STATE> jaco1(1,1), jaco2(1,1);

  fLaw->JacobFlux(U,jaco1);
  jaco1 *= normal[0];
  fLaw->JacobFlux(Up1,jaco2);
  jaco2 *= normal[0];
  if(jaco1(0,0)>=0. && jaco2(0,0)>=0.) {
    fLaw->Flux(U,flux);
    flux[0] *= normal[0];
  }
  else if(jaco1(0,0)<0. && jaco2(0,0)<=0.) {
    fLaw->Flux(Up1,flux);
    flux[0] *= normal[0];
  }
  else if(jaco1(0,0)>=0. && jaco2(0,0)<0.) {
    double deltau=U[0]-U[1];
    if(IsZero(deltau)) {
      fLaw->Flux(U,flux);
      flux[0] *= normal[0];
    }
    else {
      TPZVec<REAL> temp(1);
      fLaw->Flux(U,temp);
      temp[0] *= normal[0];
      fLaw->Flux(Up1,flux);
      flux[0] *= normal[0];
      if(sgn(temp[0]-flux[0])==sgn(deltau))
        flux[0]=temp[0];
    }
  }
  else {
    TPZVec<REAL> sonic(1,fLaw->SonicPoint());
    fLaw->Flux(sonic,flux);
    flux[0] *= normal[0];
  }
}

/*
  int TNumericalFlux::NCellsNeeded() {
  if(fFluxType < 6) return 2;
  else if(fFluxType < 23) return 4;
  else if(fFluxType < 26) return 2;
  else if(fFluxType < 28) return 4;
  return 2;
  }
*/


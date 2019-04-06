#include <stdio.h>
#include <math.h>

#include "euler.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "hadaptive.h"
#include "myheader.h"


/*******       Equacao de Euler Uni-dimensional       *******/

TEulerLaw1D::TEulerLaw1D(std::istream &input) : TConservationLaw(input) {
  fGamma=1.4;
  fNumericalFlux.SetOrder(NStateVariables());
  fNumericalFlux.SetDimension(Dimension());
}
TEulerLaw1D::TEulerLaw1D(int id) : TConservationLaw(id,"Euler Conservation Law") {
  fGamma=1.4;
  fNumericalFlux.SetOrder(NStateVariables());
  fNumericalFlux.SetDimension(Dimension());
}
TEulerLaw1D::TEulerLaw1D(int id,char *name,int type) : TConservationLaw(id,name) {
  fGamma=1.4;
  fNumericalFlux.SetOrder(NStateVariables());
  SetNumericalFluxType(type,type);
  fNumericalFlux.SetDimension(Dimension());
}
TEulerLaw1D::TEulerLaw1D(TEulerLaw1D &law) : TConservationLaw(law) {
  fGamma = law.Gamma();
}

int TEulerLaw1D::VariableIndex(char *name) {
  if(!strcmp(name,"density")) return 0;
  else if(!strcmp(name,"dens_velocity")) return 1;
  else if(!strcmp(name,"dens_energy")) return 2;
  else if(!strcmp(name,"velocity_x")) return 4;
  else if(!strcmp(name,"energy")) return 5;
  else if(!strcmp(name,"Pressure")) return 3;
	else if(!strcmp(name,"enthalpy")) return 6;
  return TPZMaterial::VariableIndex(name);
}

void TEulerLaw1D::VariablesName(TPZVec<std::string> &names) {
  int i, n = AllVariablesPost();
  names.Resize(n);
  for(i=0;i<n;i++) {
    switch(i) {
    case 0:
      names[0] = "density";
      break;
    case 1:
      names[1] = "velocity_x";
      break;
    case 2:
      names[2] = "energy";
      break;
    case 3:
      names[3] = "Pressure";
      break;
    case 4:
      names[4] = "dens_velocity";
      break;
    case 5:
      names[5] = "dens_energy";
      break;
    case 6:
      names[6] = "enthalpy";
      break;
    default:
      names[i] = "state_nonumber";
      break;
    }
  }
}

int TEulerLaw1D::NSolutionVariables(int index) {
  if(index<0) {
    PZError << "TEulerLaw1D::NSolutionVariables. Bad parameter index.\n";
    return -1;
  }
  if(index < AllVariablesImplemented()) return 1;
  return TPZMaterial::NSolutionVariables(index);
}

/**returns the solution associated with the var index based on the finite element approximation*/
void TEulerLaw1D::Solution(TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
                 int var,TPZVec<REAL> &Solout){
  Solout.Resize(NSolutionVariables(var));
  switch(var) {
  case 0:
    Solout[0] = Sol[0];
    return;
  case 1:
    Solout[0] = Sol[1];
    return;
  case 2:
    Solout[0] = Sol[2];
    return;
  case 3:
    Solout[0] = Pressure(Sol);
    return;
  case 4: {
	  REAL den = 1.e-7;
		if(!IsZero(Sol[0])) den = Sol[0];
    if(IsZero(Sol[1])) Solout[0] = 0.;
    else Solout[0] = Sol[1]/den;
    return;
	}
  case 5: {
	  REAL den = 1.e-7;
		if(!IsZero(Sol[0])) den = Sol[0];
    if(IsZero(Sol[2])) Solout[0] = 0.;
    else Solout[0] = Sol[2]/den;
    return;
	}
  case 6: {
    REAL den = 1.e-7;
    if(!IsZero(Sol[0])) den = Sol[0];
    Solout[0] = (Sol[2]+Pressure(Sol))/den;
    return;
  }
  default:
    PZError << "TEulerLaw1D::Solution. Bad Parameter var.\n";
  }

  TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

REAL TEulerLaw1D::Pressure(TPZVec<REAL> &U) {
  if(IsZero(U[0])) return (fGamma-1.)*U[2];
  return (fGamma-1.)*(U[2]-(.5*U[1]*(U[1]/U[0])));
}

void TEulerLaw1D::Flux(TPZVec<REAL> &U,TPZVec<REAL> &flux) {

//Variavel : u = (densidade, densidade*velocidade, energia)
//Funcao : f(u)=(densid*veloc, densid*(velocid)^2+pressao, veloc*(energia+pressao))
//Funcao : f(u)=( u2, ((u2*u2)/u1)+pressao, (u2/u1)*(u3+pressao))

  REAL pressao = Pressure(U);
  if(IsZero(U[0])) {
//    if(!IsZero(U[1])) {
//      printf("\nERRO : Densidade nula e momento nao nulo\n");
//      ReAssembling = 1;
      flux[0] = flux[2] = 0.;
      flux[1] = pressao;
      return;
  }
  REAL velocidade = U[1]/U[0];
  //Verificacao para existencia das caracteristicas
/*  if((2*U[0]*U[2])<(U[1]*U[1])) {
    printf("\n ERRO, Na determinacao das caracteristicas (f[u])\n");
    ReAssembling = 1;
    flux[0] = flux[1] = flux[2] = 0.;
    return;
  }*/
  flux[0]=U[1];
  flux[1]=(U[1]*velocidade)+pressao;
  flux[2]=velocidade*(U[2]+pressao);
}

void TEulerLaw1D::JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob) {
  if(IsZero(U[0])) {
//    if(!IsZero(U[1])) {
      printf("\nERRO : Densidade NULA. Falla jacobiano.\n");
      exit(1);
//    }
//    vel=0.;
//    if(!IsZero(U[2])) {
//      printf("\nERRO : Densidade NULA e energia nao nula.\n");
//      exit(1);
//    }
//    energy=0.;
  }
  REAL  vel = U[1]/U[0], energy = U[2]/U[0];

  jacob(0,0)=jacob(0,2)=0.;
  jacob(0,1)=1.;
  jacob(1,2)=fGamma-1.;
  jacob(1,0)=.5*(fGamma-3.)*vel*vel;
  jacob(1,1)=(3.-fGamma)*vel;
  jacob(2,0)=(-1.*fGamma*vel*energy)+((fGamma-1.)*vel*vel*vel);
  jacob(2,1)=(fGamma*energy)-(0.5*(fGamma+1.)*vel*vel);
  jacob(2,2)=fGamma*vel;
}

REAL TEulerLaw1D::MaxEigJacob(TPZVec<REAL> &U,TPZVec<REAL> &normal) {
	REAL vel = 0., velson = 0.;
	if(!IsZero(U[0])) {
    vel=normal[0]*(U[1]/U[0]);              //Que acontece se Ui[0]= 0.
    velson = fGamma*Pressure(U)/U[0];
  }
	if(velson<0.) {
//		printf("\nERRO, na determinacao de c en MaxEigValJacob.\n");
//		exit(1);
//cout << Pressure(U) << "   " << U[2] << "   " << U[1] << "   " << U[0] << endl;
    return fabs(vel);
	}
	velson=sqrt(velson);
	return Max(fabs(vel-velson),fabs(vel+velson));
}

REAL TEulerLaw1D::ValEigJacob(TPZVec<REAL> &u,int order,int dim) {
#ifndef NOTDEBUG
  if(dim!=1) PZError << "TEulerLaw1D::ValEigJacob. Bad parameter dim.\n";
#endif
  if(IsZero(u[0])) return 0.;
  REAL vel = u[1]/u[0];              //Que acontece se Ui[0]= 0.
  if(!order) return vel;
  REAL velson = fGamma*Pressure(u)/u[0];
	if(velson<0.) {
//    ReAssembling = 1;
		printf("\nERRO, na determinacao de c en ValEigJacob.\n");
    return vel;
	}
  velson=sqrt(velson);
  if(order==1) return vel+velson;
  return vel-velson;
}

void TEulerLaw1D::JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
#ifndef NOTDEBUG
  if(IsZero(U[0])) {
	  std::cout << "TEulerLaw1D::JacobFlux. Densidade NULA. Falla jacobiano.\n";
      exit(1);
  }
	if(jacob.Rows()!=3)
	  PZError << "TEulerLaw1D::JacobFlux. Matrix dimension incompatible.\n";
#endif
  REAL  vel = U[1]/U[0], energy = U[2]/U[0], norm = normal[0];

  jacob(0,0)=jacob(0,2)=0.;
  jacob(0,1)=norm;
  jacob(1,2)=norm*(fGamma-1.);
  jacob(1,0)=.5*norm*(fGamma-3.)*vel*vel;
  jacob(1,1)=(3.-fGamma)*vel*norm;
  jacob(2,0)=norm*((-1.*fGamma*vel*energy)+((fGamma-1.)*vel*vel*vel));
  jacob(2,1)=norm*((fGamma*energy)-(0.5*(fGamma+1.)*vel*vel));
  jacob(2,2)=norm*fGamma*vel;
}

int TEulerLaw1D::ValJacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) {
	valjacob.Zero();
  REAL rho = U[0];
  REAL den, temp, cval2 = (fGamma*Pressure(U))/rho;
  if(IsZero(rho) || cval2<0)
    return 1;
  REAL cval = sqrt(cval2);
	REAL vel = U[1]/rho, ener = U[2]/rho;

	/** Matriz de autovetores direitos do jacobiano*/
	TPZFMatrix<STATE> eigvecright(3,3), jacob(3,3);
	eigvecright(2,0) = eigvecright(2,1) = eigvecright(2,2) = 1.;
	eigvecright(1,0) = 2./vel;
	eigvecright(0,0) = 2./(vel*vel);
	temp = (fGamma - 1.)*(2.*ener*fGamma - (fGamma+1.)*vel*vel);
	eigvecright(1,1) = -2.*(temp+(fGamma-3.)*cval*vel);
	eigvecright(1,2) = 2.*(temp-(fGamma-3.)*vel*cval);
	den = (fGamma-1.)*vel;
	eigvecright(0,1) = 4*cval+2.*den;
	eigvecright(0,2) = 4*cval-2.*den;
	temp *= vel;
	den = 4*cval*(ener*fGamma-den*vel);
	eigvecright(1,1) /= (den-temp);
	eigvecright(0,1) /= (den-temp);
	den += temp;
	eigvecright(1,2) /= den;
	eigvecright(0,2) /= den;
	/** montando matriz valor absoluto |A| = R*|L|*InvR. Primeiro R*|L| */
	temp = fabs(normal[0]*vel);
	jacob(0,0) = temp*eigvecright(0,0);
	jacob(1,0) = temp*eigvecright(1,0);
	jacob(2,0) = temp*eigvecright(2,0);
	temp = fabs(normal[0]*(vel-cval));
	jacob(0,1) = temp*eigvecright(0,1);
	jacob(1,1) = temp*eigvecright(1,1);
	jacob(2,1) = temp*eigvecright(2,1);
	temp = fabs(normal[0]*(vel+cval));
	jacob(0,2) = temp*eigvecright(0,2);
	jacob(1,2) = temp*eigvecright(1,2);
	jacob(2,2) = temp*eigvecright(2,2);
	/** Multiplicando vezes a inversa dos autovetores */
	InverseJacob(eigvecright);
//	int i,j,k;
	Multiply(jacob,eigvecright,valjacob);
	
  return 0;
}

void TEulerLaw1D::ValRoeMatrix(TPZVec<REAL> &u,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &ValRoe) {
}
void TEulerLaw1D::RoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &Roe) {
}

//void TEulerLaw1D::JacobFluxInv_Diffussion(TPZVec<REAL> &U,TPZVec<REAL> &coef,TPZFMatrix &jacob) {
//}

void TEulerLaw1D::EigRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &Roe) {
	double velprom,entalprom,cprom,raizcprom,dif;
	double densi, densip1, vel, velp1, ental, entalp1;
	dif=fGamma-1.;
	densi=U[0];
	densip1=Up1[0];
	vel=U[1]/densi;
	velp1=Up1[1]/densip1;
	if((densi<0.) || (densip1<0.)) {
		printf("\nERRO : Densidade negativa\n");
		exit(1);
	}
	//determinando os valores de entalpia a esquerda e direita
	ental=(-.5)*dif*vel*vel;
	ental+=((fGamma*U[2])/densi);
	entalp1=(-.5)*dif*velp1*velp1;
	entalp1+=((fGamma*Up1[2])/densip1);
	densi=sqrt(densi);
	densip1=sqrt(densip1);
	entalprom=densip1*entalp1;
	entalprom+=(densi*ental);
	ental=densi+densip1;
	entalprom/=ental;
	velprom=densip1*velp1;
	velprom+=(densi*vel);
	velprom/=ental;
	cprom=(-.5)*velprom*velprom;
	cprom+=entalprom;
	if(cprom<0) {
		printf("\nERRO : Nao caracteristicas na matriz de Roe\n");
		exit(1);
	}
	cprom*=dif;
	raizcprom=sqrt(cprom);
	//armazenando os autovalores e autovetores da matriz de Roe por linha
	Roe[0]=Roe[6]=(velprom-raizcprom);
	Roe[1]=Roe[7]=velprom;
	Roe[2]=Roe[8]=(velprom+raizcprom);
	Roe[3]=Roe[4]=Roe[5]=1.;
	Roe[9]=(entalprom-(velprom*raizcprom));
	Roe[10]=(.5*velprom*velprom);
	Roe[11]=(entalprom+(velprom*raizcprom));
	//para armazenar os elementos da matriz inversa de R
	densi=1./raizcprom;
	densip1=0.5*dif*densi*densi;
	vel=densip1*velprom;
	velp1=0.5*vel*velprom;
	densi*=0.5;
	ental=densi*velprom;
	Roe[12]=velp1+ental;
	Roe[13]=(-1.)*(vel+densi);
	Roe[14]=Roe[20]=densip1;
	Roe[15]=1.-(2.*velp1);
	Roe[16]=2*vel;
	Roe[17]=(-2.)*densip1;
	Roe[18]=velp1-ental;
	Roe[19]=densi-vel;
}

void TEulerLaw1D::SetData(std::istream &input) {
  TConservationLaw::SetData(input);
}

void TEulerLaw1D::Print(std::ostream &out) {
	double dens,temp;
	fName += ".out";
	FILE * input;
	FILE * output1;
	FILE * output2;
	FILE * output3;
	input = fopen(Name().c_str(),"r");
	char *densidad="densidad.eul";
	char *velocity="velocity.eul";
	char *energy="energy.eul";
	output1 = fopen(densidad,"w");
	output2 = fopen(velocity,"w");
	output3 = fopen(energy,"w");
	while(fscanf(input,"%lf",&dens)) {
		fscanf(input,"%lf",&temp);
		fprintf(output1,"%f\t",dens);
		if(IsZero(dens)) fprintf(output2,"%f\t",0.);
		else fprintf(output2,"%f\t%",(temp/dens));
		fscanf(input,"%lf",&temp);
		fprintf(output3,"%f\t",temp);
//         }
//         fprintf(output1,"\n"); fprintf(output2,"\n"); fprintf(output3,"\n");
	}
	fclose(input);
	fclose(output1);
	fclose(output2);
	fclose(output3);
}

void TEulerLaw1D::SolRiemannProblem(double *Ul,double *Ur,double *Result) {
	double ldensidade=Ul[0],rdensidade=Ur[0];
	double pressao,velocidade,densidade2,densidade3;
	double An, Bn,gamma=fGamma;
	int ite=0;
	double lvel,rvel,lpressao,rpressao,lc,rc;
	double w,pressaom1=-1.,dif=1.;
	//Determinando os valores da pressao, velocidade, velson de Ul e Ur
	rpressao=Ur[2]; lpressao=Ul[2];
	//Verificando a compatibilidade dos dados
	if(!IsZero(ldensidade)) {
		lvel=Ul[1]/ldensidade;
		lpressao-=(.5*Ul[1]*lvel);
	}
	else {
		lvel=0.;
		if(!IsZero(Ul[1]))
			printf("\nCUIDADO, ldensidade nula, mas lmomento nao.\n");
	}
	if(!IsZero(rdensidade)) {
		rvel=Ur[1]/rdensidade;
		rpressao-=(.5*Ur[1]*rvel);
	}
	else {
		rvel=0.;
		if(!IsZero(Ul[4]))
			printf("\nCUIDADO, rdensidade nula, mas rmomento nao.\n");
	}
	if((rpressao<0.)||(lpressao<0.))
		printf("\nERRO, nao existem caracteristicas, pressao negativa.\n");

	rpressao*=(gamma-1.);
	lpressao*=(gamma-1.);
	lc=gamma*lpressao;
	rc=gamma*rpressao;
	lc=sqrt(lc/ldensidade);
	rc=sqrt(rc/rdensidade);
	//Armazenando os dados parciais obtidos
	Result[0]=lvel; Result[1]=lpressao; Result[2]=lc; Result[3]=lvel-lc;
	Result[4]=rvel; Result[5]=rpressao; Result[6]=rc; Result[7]=rvel+rc;
	//Achando os estados intermeios
	pressao=.5*(lpressao+rpressao);
	while(!IsZero(dif)&&(ite<50)) {
		ite++;
		//determinando An
		w=pressao/lpressao;
		An=Phi(w);
		w=lpressao*ldensidade;
		An*=sqrt(w);
		//determinando Bn
		w=pressao/rpressao;
		Bn=Phi(w);
		w=rpressao*rdensidade;
		Bn*=sqrt(w);
		pressao=(Bn*lpressao)+(An*rpressao);
		pressao+=(An*Bn*(lvel-rvel));
		pressao/=(An+Bn);
		dif=ValAbs(pressao-pressaom1);
		if((ite==40)&&(dif>.05)) {
			double z=rpressao+lpressao;
			z=pressaom1/z;
			dif=(gamma+1.)/(2*gamma);
			dif=pow(z,dif);
			dif-=1.;
			z=(1.-z)/dif;
			z*=((gamma-1.)/(3*gamma));
			z-=1.;
			z=Max(z,0.);
			pressao+=(z*pressaom1);
			pressao=pressao/(z+1.);
		}
		pressaom1=pressao;
	}
	velocidade=(An*lvel)+(Bn*rvel);
	velocidade+=(lpressao-rpressao);
	velocidade/=(An+Bn);
	densidade2=((gamma+1.)*pressao)+((gamma-1.)*lpressao);
	densidade3=((gamma+1.)*pressao)+((gamma-1.)*rpressao);
	densidade2/=(((gamma-1.)*pressao)+((gamma+1.)*lpressao));
	densidade3/=(((gamma-1.)*pressao)+((gamma+1.)*rpressao));
	densidade2*=ldensidade;
	densidade3*=rdensidade;
	//Armazenando os dados dos estados intermeios
	Result[8]=densidade2;
	rc=densidade2*velocidade;
	Result[9]=rc;
	lc=(1./(gamma-1.))*pressao;
	Result[10]=lc+(.5*rc*velocidade);
	Result[11]=densidade3;
	rc=densidade3*velocidade;
	Result[12]=rc;
	Result[13]=lc+(.5*rc*velocidade);
	Result[14]=velocidade;
	Result[15]=pressao;
	rc=gamma*pressao;
	lc=sqrt(rc/densidade2);
	Result[16]=velocidade-lc;
	lc=sqrt(rc/densidade3);
	Result[17]=velocidade+lc;
	//velocidades de propagacao das ondas consideradas como shocks
	Result[18]=lvel-(An/ldensidade);
	Result[19]=rvel+(Bn/rdensidade);
}

double TEulerLaw1D::Phi(double w) {
	double gamma=fGamma;
	double expoente=.5*(gamma-1.);
	if(w<1.) {
		expoente/=gamma;
		gamma=sqrt(gamma);
		gamma*=(expoente*(1.-w));
		expoente=pow(w,expoente);
		expoente=1.-expoente;
		return(gamma/expoente);
	}
	else {
		expoente+=(.5*(gamma+1.)*w);
		return(sqrt(expoente));
	}
}


/* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
int TEulerLaw1D::IdBC(double *x) {
  if(x[0]<0.5) return -1;
  return -2;
}

void TEulerLaw1D::ContributeOverInterface(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZVec<REAL> &up1,
     REAL weight,REAL area,int type,TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &phi0,TPZFMatrix<STATE> &phi1,TPZVec<REAL> &normal,
     TPZFMatrix<STATE> &/*ek*/,TPZFMatrix<STATE> &ef) {

  SetPoint(x);
//	if(Dimension()==2) area = sqrt(area);
  /**To determine maxime jacobian eigenvalue to next time step */
	if(type < 0) {
    up1[0] = u[0];
    u[1] = up1[1] = 0.;  // fazendo velocidade zero.
    up1[2] = Pressure(u)/(fGamma - 1.);   // energia igual a pressao
    type = 2;
	}
	else {
    REAL maxeigen = MaxEigJacob(u,normal);
    if(fMaxEigen<maxeigen) fMaxEigen = maxeigen;
    maxeigen = MaxEigJacob(up1,normal);
    if(fMaxEigen<maxeigen) fMaxEigen = maxeigen;
  }

  if(area!=0.) SetAlfa(fDeltaT/area);   // Preenchendo alfa para fluxos numericos de segunda ordem

  int phr0,phr1,efr;
  int nvar = NStateVariables();
  phr0 = phi0.Rows();
  phr1 = phi1.Rows();
  efr = ef.Rows();
#ifndef NOTDEBUG
  if(efr != (phr0+phr1)*nvar) {
    PZError << "\nTConservationLaw. Inconsistent input data : \n" <<
      " phi.Rows = " << (phi0.Rows()+phi1.Rows()) << "\nef.Rows() = " << ef.Rows() << "\n";
  }
#endif NOTDEBUG
  int i, c;
  TPZVec<REAL> flux(nvar,0.);

  fNumericalFlux.SetFluxType(type);
  fNumericalFlux.NumericalFlux(u,up1,normal,flux);
  for(c=0;c<nvar;c++) flux[c] *= weight;

  for(i=0;i<phr0;i++) {
    for(c=0;c<nvar;c++) {
        ef(i*nvar+c,0) += (-1.*phi0(i,0)*flux[c]);
    }
  }
  phr1 += phr0;
  /** Computing ef */
  for(i=phr0;i<phr1;i++) {
    for(c=0;c<nvar;c++) {
        ef(i*nvar+c,0) += (phi1(i-phr0,0)*flux[c]);
    }
  }
}

void TEulerLaw1D::InverseJacob(TPZFMatrix<STATE> &jac) {
	REAL a=jac.GetVal(0,0), b=jac.GetVal(0,1), c=jac.GetVal(0,2);
	REAL d=jac.GetVal(1,0), e=jac.GetVal(1,1), f=jac.GetVal(1,2);
	REAL g=jac.GetVal(2,0), h=jac.GetVal(2,1), i=jac.GetVal(2,2);
  REAL EIHF = e*i-h*f, FGDI = f*g-d*i, HDGE = h*d-g*e;
  REAL det = a*EIHF + b*FGDI + c*HDGE;
  if(IsZero(det)) {
    PZError << "TEulerLaw1D::InverseJacob. Determinante zero, matriz singular.\n";
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
}

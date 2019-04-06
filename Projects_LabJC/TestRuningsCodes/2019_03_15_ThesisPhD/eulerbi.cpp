
#include <stdio.h>
#include <math.h>

#include "eulerbi.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "hadaptive.h"
#include "myheader.h"

#include "eulerdif.h"


/*******       Equacao de Euler Bi-dimensional       *******/

TEulerLaw2D::TEulerLaw2D(std::istream &input) : TConservationLaw(input) {
  fGamma=1.4;
  fNumericalFlux.SetOrder(NStateVariables());
  fNumericalFlux.SetDimension(Dimension());
}
TEulerLaw2D::TEulerLaw2D(int id) : TConservationLaw(id,"Euler-2D Conservation Law") {
  fGamma=1.4;
  fNumericalFlux.SetOrder(NStateVariables());
  fNumericalFlux.SetDimension(Dimension());
}
TEulerLaw2D::TEulerLaw2D(int id,char *name,int type) : TConservationLaw(id,name) {
  fGamma=1.4;
  fNumericalFlux.SetOrder(NStateVariables());
  SetNumericalFluxType(type,type);
  fNumericalFlux.SetDimension(Dimension());
}
TEulerLaw2D::TEulerLaw2D(TEulerLaw2D &law) : TConservationLaw(law) {
  fGamma = law.Gamma();
}

int TEulerLaw2D::IdBC(double *x) {
	if(fName == "Euler two-dimensional") {
	  if(IsZero(x[0])) return -1;
	  else if(IsZero(x[0]-fXend)) return -2;
//	  return -3;
	}
	else if(fName == "Euler shock-obliq") {
		if(IsZero(x[0])) return -1;
		else if(IsZero(x[1]-1.)) return -2;
		else if(IsZero(x[1])) return -5;
	}
	return -3;
}

int TEulerLaw2D::VariableIndex(char *name) {
  if(!strcmp(name,"density")) return 0;
  else if(!strcmp(name,"velocity_x")) return 5;
  else if(!strcmp(name,"velocity_y")) return 6;
  else if(!strcmp(name,"energy")) return 7;
	else if(!strcmp(name,"enthalpy")) return 8;
  else if(!strcmp(name,"Pressure")) return 4;
  else if(!strcmp(name,"dens_velocity_x")) return 1;
  else if(!strcmp(name,"dens_velocity_y")) return 2;
  else if(!strcmp(name,"dens_energy")) return 3;
  return TPZMaterial::VariableIndex(name);
}

void TEulerLaw2D::VariablesName(TPZVec<std::string> &names) {
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
      names[2] = "velocity_y";
      break;
    case 3:
      names[3] = "energy";
      break;
    case 4:
      names[4] = "Pressure";
      break;
    case 5:
      names[5] = "enthalpy";
      break;
    case 6:
      names[6] = "dens_velocity_x";
      break;
    case 7:
      names[7] = "dens_velocity_y";
      break;
    case 8:
      names[8] = "dens_energy";
      break;
    default:
      names[i] = "state_nonumber";
      break;
    }
  }
}

int TEulerLaw2D::NSolutionVariables(int index) {
  if(index<0) {
    PZError << "TEulerLaw2D::NSolutionVariables. Bad parameter index.\n";
    return -1;
  }
  if(index < AllVariablesImplemented()) return 1;
  return TPZMaterial::NSolutionVariables(index);
}

/**returns the solution associated with the var index based on the finite element approximation */
void TEulerLaw2D::Solution(TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout) {
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
    Solout[0] = Sol[3];
    return;
  case 4:
    Solout[0] = Pressure(Sol);
    return;
  case 5: {
    REAL den = 1.e-7;
    if(!IsZero(Sol[0])) den = Sol[0];
    if(IsZero(Sol[1])) Solout[0] = 0.;
    else Solout[0] = Sol[1]/den;
    return;
  }
  case 6: {
    REAL den = 1.e-7;
    if(!IsZero(Sol[0])) den = Sol[0];
    if(IsZero(Sol[2])) Solout[0] = 0.;
    else Solout[0] = Sol[2]/den;
    return;
  }
  case 7: {
    REAL den = 1.e-7;
    if(!IsZero(Sol[0])) den = Sol[0];
    if(IsZero(Sol[3])) Solout[0] = 0.;
    else Solout[0] = Sol[3]/den;
    return;
  }
  case 8: {
    REAL den = 1.e-7;
    if(!IsZero(Sol[0])) den = Sol[0];
    Solout[0] = (Sol[3]+Pressure(Sol))/den;
    return;
  }
  default:
    PZError << "TEulerLaw2D::Solution. Bad Parameter var.\n";
  }
  TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

REAL TEulerLaw2D::Pressure(TPZVec<REAL> &U) {
  if(IsZero(U[0])) return(fGamma-1.)*U[3];
  REAL velocity = (U[1]*U[1]+U[2]*U[2])/U[0];
  return (fGamma-1.)*(U[3]-(.5*velocity));
}

//------------------->  Verificar  A PARTIR DAQUI   <-------------
void TEulerLaw2D::Flux(TPZVec<REAL> &U,TPZVec<REAL> &flux) {

//Variavel : u = (densidade, densidade*velocidade, energia)
//Funcao : f(u)=(densid*veloc, densid*(velocid)^2+pressao, veloc*(energia+pressao))
//Funcao : f(u)=( u2, ((u2*u2)/u1)+pressao, (u2/u1)*(u3+pressao))

  REAL pressao = Pressure(U);
  if(IsZero(U[0])) {
//    if(!IsZero(U[1]) || !IsZero(U[2])) {
//      printf("\nERRO : Densidade nula e momentos nao nulo\n");
//      ReAssembling = 1;
      flux[0] = flux[2] = flux[3] = flux[4] = flux[5] = flux[7] = 0.;
//      flux[1] = flux[6] = 0.;
	  std::cout << "EulerLaw2D::Flux find nule density.\n";
      flux[1] = flux[6] = pressao;  // Jorge 17/08
      return;
  }
  flux[0] = U[1];
  flux[1] = ((U[1]/U[0])*U[1]) + pressao;
  flux[2] = U[1]*(U[2]/U[0]);
  flux[5] = U[2]*(U[1]/U[0]);
  flux[3] = (U[1]/U[0])*(U[3]+pressao);
  flux[4] = U[2];
  flux[6] = ((U[2]/U[0])*U[2]) + pressao;
  flux[7] = (U[2]/U[0])*(U[3] + pressao);
}

void TEulerLaw2D::InverseJacob(TPZFMatrix<STATE> &jac) {
	REAL a=jac.GetVal(0,0), b=jac.GetVal(0,1), c=jac.GetVal(0,2), d=jac.GetVal(0,3);
	REAL e=jac.GetVal(1,0), f=jac.GetVal(1,1), g=jac.GetVal(1,2), h=jac.GetVal(1,3);
	REAL i=jac.GetVal(2,0), j=jac.GetVal(2,1), k=jac.GetVal(2,2), l=jac.GetVal(2,3);
	REAL m=jac.GetVal(3,0), n=jac.GetVal(3,1), r=jac.GetVal(3,2), s=jac.GetVal(3,3);
  REAL BHDF = b*h-d*f, BKCJ = b*k-c*j, BSDN = b*s-d*n, CFBG = c*f-b*g;
	REAL CLDK = c*l-d*k, CNBR = c*n-b*r, DGCH = d*g-c*h, DJBL = d*j-b*l;
	REAL DRCS = d*r-c*s, FRGN = f*r-g*n, GJFK = g*j-f*k, GSHR = g*s-h*r;
	REAL FSHN = f*s-h*n, LNJS = l*n-j*s, HJFL = h*j-f*l, HKGL = h*k-g*l;
	REAL KNJR = k*n-j*r, LRKS = l*r-k*s;
  REAL det = DJBL*(g*m-e*r)+BHDF*(k*m-i*r)+BSDN*(g*i-e*k);
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
}

void TEulerLaw2D::JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob) {

  if(IsZero(U[0])) {
//    if(!IsZero(U[1])) {
	  std::cout << "\nERRO : Densidade NULA. Falla jacobiano.\n";
      exit(1);

  }
  REAL velx, vely, ener, aux, gammam1 = fGamma-1.;
  velx = U[1]/U[0];
  vely = U[2]/U[0];
  ener = U[3]/U[0];

  jacob(0,0)=jacob(0,2)=jacob(0,3)=jacob(2,3)=jacob(4,0)=jacob(4,1)=jacob(4,3)=jacob(5,3)=0.;
  jacob(0,1)=jacob(4,2)=1.;
  jacob(1,3)=jacob(6,3)=gammam1;
  jacob(1,1)=(3.-fGamma)*velx;
  jacob(6,2)=(3.-fGamma)*vely;
  jacob(1,2)=-1.*gammam1*vely;
  jacob(6,1)=-1.*gammam1*velx;
  jacob(3,2)=jacob(7,1)=jacob(1,2)*velx;
  jacob(2,0)=jacob(5,0)=-1.*velx*vely;
  jacob(2,1)=jacob(5,1)=vely;
  jacob(2,2)=jacob(5,2)=velx;
  jacob(3,3)=fGamma*velx;
  jacob(7,3)=fGamma*vely;
  aux = gammam1*(velx*velx+vely*vely);
  jacob(1,0) = 0.5*aux - (velx*velx);
  jacob(3,1) = (fGamma*ener) - (.5*gammam1*(3*velx*velx+vely*vely));
  jacob(7,2) = (fGamma*ener) - (.5*gammam1*(velx*velx+3*vely*vely));
  jacob(6,0) = 0.5*aux - (vely*vely);
  jacob(3,0) = (aux-(fGamma*ener))*velx;
  jacob(7,0) = (aux-(fGamma*ener))*vely;
}

void TEulerLaw2D::JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
#ifndef NOTDEBUG
  if(IsZero(U[0])) {
	  std::cout << "TEulerLaw2D::JacobFlux. Densidade NULA. Falla jacobiano.\n";
      exit(1);
  }
	if(jacob.Rows()!=4)
	  PZError << "TEulerLaw2D::JacobFlux. Matrix dimension incompatible.\n";
#endif

  REAL velx, vely, ener, aux, gammam1 = fGamma-1.;
  velx = U[1]/U[0];
  vely = U[2]/U[0];
  ener = U[3]/U[0];
	REAL alfa = normal[0], beta= normal[1];
  aux = gammam1*(velx*velx+vely*vely);

  jacob(0,0)=jacob(0,3)=0.;
  jacob(0,1)=alfa; jacob(0,2)=beta;
  jacob(1,3)= alfa*gammam1;
  jacob(2,3)= beta*gammam1;
  jacob(1,1)= alfa*((3.-fGamma)*velx)+beta*vely;
  jacob(1,2)=-alfa*gammam1*vely+beta*velx;
  jacob(3,2)=-alfa*gammam1*velx*vely+beta*((fGamma*ener) - .5*gammam1*(3*vely*vely+velx*velx));
  jacob(2,0)=-alfa*velx*vely+beta*(0.5*aux - (vely*vely));
  jacob(2,1)=alfa*vely-beta*gammam1*velx;
  jacob(2,2)=alfa*velx+beta*(3.-fGamma)*vely;
  jacob(3,3)=fGamma*(alfa*velx+beta*vely);
  jacob(1,0) = alfa*(0.5*aux - (velx*velx))-beta*(velx*vely);
  jacob(3,1) = alfa*((fGamma*ener) - .5*gammam1*(3*velx*velx+vely*vely))-beta*gammam1*vely*velx;
  jacob(3,0) = (alfa*velx+beta*vely)*(aux - fGamma*ener);
}

REAL TEulerLaw2D::MaxEigJacob(TPZVec<REAL> &U,TPZVec<REAL> &normal) {
  REAL rho = U[0];
	if(IsZero(rho)) return 0.;
  REAL lambda, c, d;
  REAL alfa = normal[0], beta = normal[1];
  lambda = (alfa*U[1] + beta*U[2])/rho;
  c = (fGamma*Pressure(U))/rho;
  if(c < 0.) return lambda;
  c = sqrt(c);
  d = sqrt(alfa*alfa + beta*beta);
	return Max(fabs(lambda - d*c),fabs(lambda + d*c));
}
REAL TEulerLaw2D::ValEigJacob(TPZVec<REAL> &U,int order,int dim) {
  if(IsZero(U[0])) return 0.;
  REAL velx = U[1]/U[0];
  if(dim==1) {
    if(!order || order==1) return velx;
    REAL velson = (fGamma*Pressure(U))/U[0];
		if(velson<0.) return velx;
		velson = sqrt(velson);
    if(order==2) return velx+velson;
    return velx-velson;
  }
  REAL vely = U[2]/U[0];
  if(!order || order==1) return vely;
  REAL velson = (fGamma*Pressure(U))/U[0];
  if(velson<0.) return vely;
  if(order==2) return vely+velson;
  return vely-velson;
}

void TEulerLaw2D::ValRoeMatrix(TPZVec<REAL> &u,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &ValRoe) {
}
void TEulerLaw2D::RoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &Roe) {
}

void TEulerLaw2D::EigRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &Roe) {
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

void TEulerLaw2D::SetData(std::istream &input) {
  TConservationLaw::SetData(input);
  GetDataCommented(input,fXend);
}

void TEulerLaw2D::Print(std::ostream &out) {
}

void TEulerLaw2D::ContributeOverInterface(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZVec<REAL> &up1,
     REAL weight,REAL area,int type,TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &phi0,TPZFMatrix<STATE> &phi1,TPZVec<REAL> &normal,
     TPZFMatrix<STATE> &/*ek*/,TPZFMatrix<STATE> &ef) {

// type = -8 condicao fronteira na presenca de uma parede. Acertando a pressao
// type = -2 caso de saida livre
// type =  0 caso Dirichlet
  SetPoint(x);
  /**To determine maxime jacobian eigenvalue to next time step */
  if(type == -8) {
    up1[0] = u[0];
		REAL Vtangent, VNormal;
    Vtangent = ((-normal[1])*u[1]) + normal[0]*u[2];   // fazendo velocidade normal igual a zero.
		VNormal = u[1]*normal[0]+u[2]*normal[1];
    up1[3] = Pressure(u)/(fGamma-1.);   // energia igual a pressao
    up1[1] = Vtangent*(-normal[1])-VNormal*normal[0];
    up1[2] = Vtangent*normal[0]-VNormal*normal[1];
		if(!IsZero(u[0]))
      up1[3] += (.5*(u[1]*u[1]+u[2]*u[2])/u[0]);
  }
  else if(type == -2) {
    up1[0] = u[0];
    up1[1] = u[1];
    up1[2] = u[2];
    up1[3] = u[3];
  }
//  else {
    REAL maxeigen = MaxEigJacob(u,normal);
    if(fMaxEigen<maxeigen) fMaxEigen = maxeigen;
    maxeigen = MaxEigJacob(up1,normal);
    if(fMaxEigen<maxeigen) fMaxEigen = maxeigen;
//  }

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

  fNumericalFlux.SetFluxType(FluxType());
  fNumericalFlux.NumericalFlux(u,up1,normal,flux);
  for(c=0;c<nvar;c++) flux[c] *= weight;

  for(i=0;i<phr0;i++) {
    for(c=0;c<nvar;c++) {
      ef(i*nvar+c,0) += (-1.*phi0(i,0)*flux[c]);
    }
  }
  phr1 += phr0;
  // Computing ef
  for(i=phr0;i<phr1;i++) {
    for(c=0;c<nvar;c++) {
      ef(i*nvar+c,0) += (phi1(i-phr0,0)*flux[c]);
    }
  }
}

int TEulerLaw2D::ValJacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) {
	valjacob.Zero();
  REAL rho = U[0];
  REAL vx, vy, lambda, cval, cval2, u, kt, st, aa, bb;
  REAL gM1 = fGamma-1.;
  vx = U[1]/rho;
  vy = U[2]/rho;
  cval2 = (fGamma*Pressure(U))/rho;
  if(cval2<0) {
    return 1;
  }
  cval = sqrt(cval2);
 
  u = sqrt(vx*vx + vy*vy);
  REAL Mach = u/cval;
  REAL Mach2 = u*u/cval2;
  if(IsZero(u)) {
    kt=1.;
    st=0.;
  } else {
    kt = vx/u;
    st = vy/u;
  }
  aa = normal[0]*kt + normal[1]*st;
  bb = normal[1]*kt - normal[0]*st;

	REAL aabbnorm2 = aa*aa+bb*bb;
  REAL aabbnorm = sqrt(aabbnorm2);
  //Para o computo da matrix direita
  TPZFMatrix<STATE> right(4,4);
  TPZVec<REAL> diag(4);
  REAL aascal = aa/aabbnorm;
  REAL bbscal = bb/aabbnorm;
  right(0,0) = -(bbscal*cval*Mach);
  right(0,1) = bbscal*kt + aascal*st;
  right(0,2) = -(aascal*kt) + bbscal*st;
  right(0,3) = 0.;

  right(1,0) = 1 - gM1*Mach2/2.;
  right(1,1) = (gM1*kt*Mach)/cval;
  right(1,2) = (gM1*Mach*st)/cval;
  right(1,3) = -gM1/cval2;

  right(2,0) = (cval2*Mach*(2*aascal + gM1*Mach))/4.;
  right(2,1) = -(cval*(aascal*kt + gM1*kt*Mach - bbscal*st))/2.;
  right(2,2) = -(cval*(bbscal*kt + (aascal + gM1*Mach)*st))/2.;
  right(2,3) = gM1/2.;

  right(3,0) = (cval2*Mach*(-2*aascal + gM1*Mach))/4.;
  right(3,1) = (cval*(aascal*kt - gM1*kt*Mach - bbscal*st))/2.;
  right(3,2) = (cval*(bbscal*kt + (aascal - gM1*Mach)*st))/2.;
  right(3,3) = gM1/2.;

  //matriz esquerda RtQX*val(autovalores)
  TPZFMatrix<STATE> left(4,4);
  left(0,0) = 0.;
  left(0,1) = 1.;
  left(0,2) = 1/cval2;
  left(0,3) = 1/cval2;
  
  left(1,0) = bbscal*kt + aascal*st;
  left(1,1) = cval*kt*Mach;
  left(1,2) = (-(aascal*kt) + kt*Mach + bbscal*st)/cval;
  left(1,3) = (aascal*kt + kt*Mach - bbscal*st)/cval;

  left(2,0) = -(aascal*kt) + bbscal*st;
  left(2,1) = cval*Mach*st;
  left(2,2) = -((bbscal*kt + aascal*st - Mach*st)/cval);
  left(2,3) = (bbscal*kt + (aascal + Mach)*st)/cval;


  left(3,0) = bbscal*cval*Mach;
  left(3,1) = (cval2*Mach2)/2.;
  left(3,2) = 1/gM1 + (Mach*(-2*aascal + Mach))/2.;
  left(3,3) = 1/gM1 + (Mach*(2*aascal + Mach))/2.;

  //right.Print("right eigenvalues");
  //left.Print("left eigenvalues");
  diag[0] = fabs(aa*cval*Mach);
  diag[1] = fabs(aa*cval*Mach);
  diag[2] = fabs(-(aabbnorm*cval) + aa*cval*Mach);
  diag[3] = fabs(cval*(aabbnorm + aa*Mach));

  int i,j,k;
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      valjacob(i,j)=0.;
      for(k=0;k<4;k++) {
        valjacob(i,j) += left(i,k)*diag[k]*right(k,j);
      }
    }
  }
  return 0;
}

/******************************   **************************************/
void TEulerLaw2D::JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &Ajacob,TPZFMatrix<STATE> &Bjacob) {
  if(IsZero(U[0])) {
	  std::cout << "\nERRO : Densidade NULA. Falla jacobiano.\n";
      exit(1);
  }

  REAL velx, vely, ener, aux, gammam1 = fGamma-1.;
  velx = U[1]/U[0];
  vely = U[2]/U[0];
  ener = U[3]/U[0];

  Ajacob(0,0)=Ajacob(0,2)=Ajacob(0,3)=Ajacob(2,3)=0.;
  Bjacob(0,0)=Bjacob(0,1)=Bjacob(0,3)=Bjacob(1,3)=0.;
  Ajacob(0,1)=Bjacob(0,2)=1.;
  Ajacob(1,3)=Bjacob(2,3)=gammam1;
  Ajacob(1,1)=(3.-fGamma)*velx;
  Bjacob(2,2)=(3.-fGamma)*vely;
  Ajacob(1,2)=-1.*gammam1*vely;
  Bjacob(2,1)=-1.*gammam1*velx;
  Ajacob(3,2)=Bjacob(3,1)=Ajacob(1,2)*velx;
  Ajacob(2,0)=Bjacob(1,0)=-1.*velx*vely;
  Ajacob(2,1)=Bjacob(1,1)=vely;
  Ajacob(2,2)=Bjacob(1,2)=velx;
  Ajacob(3,3)=fGamma*velx;
  Bjacob(3,3)=fGamma*vely;
  aux = gammam1*(velx*velx+vely*vely);
  Ajacob(1,0) = 0.5*aux - (velx*velx);
  Ajacob(3,1) = (fGamma*ener) - (.5*gammam1*(3*velx*velx+vely*vely));
  Bjacob(3,2) = (fGamma*ener) - (.5*gammam1*(velx*velx+3*vely*vely));
  Bjacob(2,0) = 0.5*aux - (vely*vely);
  Ajacob(3,0) = (aux-(fGamma*ener))*velx;
  Bjacob(3,0) = (aux-(fGamma*ener))*vely;
}

void TEulerLaw2D::InvJacob2d(TPZFMatrix<STATE> &axes,TPZFMatrix<STATE> &jacinv) {

  REAL tmp[2][2];
  // dksi/dx
  tmp[0][0] = jacinv(0,0)*axes(0,0)+jacinv(0,1)*axes(1,0);
  //d ksi/d y
  tmp[0][1] = jacinv(0,0)*axes(0,1)+jacinv(0,1)*axes(1,1);
  //d eta / d x
  tmp[1][0] = jacinv(1,0)*axes(0,0)+jacinv(1,1)*axes(1,0);
  // d eta / d y
  tmp[1][1] = jacinv(1,0)*axes(0,1)+jacinv(1,1)*axes(1,1);
  
  jacinv(0,0) = tmp[0][0];
  jacinv(0,1) = tmp[0][1];
  jacinv(1,0) = tmp[1][0];
  jacinv(1,1) = tmp[1][1];
}

void TEulerLaw2D::MatrixDiff(TPZVec<REAL> &sol,TPZFMatrix<REAL> &axes, TPZFMatrix<STATE> &jacinv,TPZFMatrix<STATE>
&ATauA,TPZFMatrix<STATE> &ATauB,TPZFMatrix<STATE> &BTauA,TPZFMatrix<STATE> &BTauB) {

  //Computando o jacobiano da transformacao do elemento mestre ao elemento triangular
  InvJacob2d(axes,jacinv);

  // Vetor de coeficientes dos jacobianos (alfa,beta) onde alfa*B+beta*B
  TPZVec<REAL> Beta(2);

  TPZFMatrix<STATE> Ajacob(4,4);
  TPZFMatrix<STATE> Bjacob(4,4);
  TPZFMatrix<STATE> valjacob(4,4);
  TPZFMatrix<STATE> invvaljacob(4,4,0.);

  JacobFlux(sol,Ajacob,Bjacob);
	
  // Calculando Tau
  // Matrix |	jacinv(0,0)*A + jacinv(0,1)*B|
	Beta[0] = jacinv(0,0); Beta[1] = jacinv(0,1);
	ValJacobFlux(sol,valjacob,Beta);
	//valjacob.Print("Jacob ValAbs dksi");
	
	invvaljacob = valjacob;
	Beta[0] = jacinv(1,0); Beta[1] = jacinv(1,1);
	ValJacobFlux(sol,valjacob,Beta);
	//valjacob.Print("Jacob ValAbs deta");
	
	invvaljacob += valjacob;
	InverseJacob(invvaljacob);
	//invvaljacob.Print("InvValJacob");

	ATauA.Zero();
	ATauB.Zero();
	BTauA.Zero();
	BTauB.Zero();
	int i,j,k,l;
	for(i=0;i<4;i++) {
	  for(j=0;j<4;j++) {
	    for(k=0;k<4;k++) {
	      for(l=0;l<4;l++) {
					ATauA(i,j) += Ajacob(i,k)*invvaljacob(k,l)*Ajacob(l,j);
					ATauB(i,j) += Ajacob(i,k)*invvaljacob(k,l)*Bjacob(l,j);
					BTauA(i,j) += Bjacob(i,k)*invvaljacob(k,l)*Ajacob(l,j);
					BTauB(i,j) += Bjacob(i,k)*invvaljacob(k,l)*Bjacob(l,j);
	      }
	    }
	  }
	}
}

/*void TEulerLaw2D::IncrementDiffusion(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,
        REAL weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,
        TPZFMatrix &ek,TPZFMatrix &ef,REAL cfldiff) {
  SetPoint(x);
  int phr = phi.Rows();
  int i,j,k,p,c;

  TPZFMatrix dphixy(2,phr,0.);
  for(i=0;i<phr;i++)
    for(j=0;j<2;j++)
      for(k=0;k<2;k++) dphixy(j,i) += dphi(k,i)*axes(k,j);

  //Computando o jacobiano da transformacao do elemento mestre ao elemento triangular
  TPZFMatrix jacinvxy(2,2);
	jacinvxy(0,0) = jacinv(0,0)*axes(0,0)+jacinv(0,1)*axes(1,0);
  jacinvxy(0,1) = jacinv(0,0)*axes(0,1)+jacinv(0,1)*axes(1,1);
  jacinvxy(1,0) = jacinv(1,0)*axes(0,0)+jacinv(1,1)*axes(1,0);
  jacinvxy(1,1) = jacinv(1,0)*axes(0,1)+jacinv(1,1)*axes(1,1);

  int dsolr = dsol.Rows(), dsolc = dsol.Cols();
  REAL Prod, coeff;
  TPZVec<REAL> Beta(2);
  TPZFMatrix vjacob(4,4);
  TPZFMatrix matprod(4,4);
  TPZFMatrix tau(4,4,0.);
	TPZFMatrix temp(4,4,0.);

  coeff = fCoef*cfldiff*weight*fDeltaT;
	
  //Computa a soma das matrices valores absolutos e a inversa da soma
  Beta[0] = jacinvxy(0,0); Beta[1] = jacinvxy(0,1);
  p = ValJacobFlux(sol,vjacob,Beta);
  tau = vjacob;
  Beta[0] = jacinvxy(1,0); Beta[1] = jacinvxy(1,1);
  p += ValJacobFlux(sol,vjacob,Beta);
  tau += vjacob;
  if(!p)
    InverseJacob(tau);
  else
    tau.Identity();

	tau *= coeff;

 if(fExplicit<1) {
  for(i=0;i<phr;i++) {
		Beta[0] = dphixy(0,i); Beta[1] = dphixy(1,i);
		JacobFlux(sol,vjacob,Beta);
		Multiply(vjacob,tau,matprod);*/
/*	TEulerDiffusivity eulerdif;
		TPZFMatrix ATA(4,4,0.);
		TPZFMatrix ATB(4,4,0.), BTA(4,4,0.), BTB(4,4,0.);
		TPZFMatrix K(temp);

		eulerdif.MatrixDiff(sol,axes,jacinv,ATA,ATB,BTA,BTB);
		for(c=0;c<4;c++) {
			for(k=0;k<4;k++) {
				K(c,k) -= cfldiff*fDeltaT*weight*(
					dphixy(0,i)*ATA(c,k)*dphixy(0,j)
					+dphixy(0,i)*ATB(c,k)*dphixy(1,j)
					+dphixy(1,i)*BTA(c,k)*dphixy(0,j)
					+dphixy(1,i)*BTB(c,k)*dphixy(1,j)
					);
      }
    }
//		ATA.Print("ATA1");
//		ATB.Print("ATB1");
//		BTA.Print("BTA1");
//		BTB.Print("BTB1");
		for(c=0;c<4;c++) {
			for(k=0;k<4;k++) {
				if(fabs(K(c,k)) > 1.e-18) cout << "K"<< c << k << " =  " << K(c,k) << "    time = " << Time() << endl;
			}
		}

//		K.Print("Versao Philippe1");
//		pause = Pause(pause);
*/
/*		for(j=0;j<phr;j++) {
		  Beta[0] = dphixy(0,j); Beta[1] = dphixy(1,j);
		  JacobFlux(sol,vjacob,Beta);
			Multiply(matprod,vjacob,temp);

      for(c=0;c<4;c++) {
			  for(k=0;k<4;k++) {
          ek(4*i+c,4*j+k) += temp(c,k);
				}
			}
    }
  }

 }
 else if(fExplicit==1) {
  TPZVec<REAL> divflux(4,0.);
  for(i=0;i<phr;i++) {
    p = 4*i;
    for(j=0;j<4;j++) {
//      module = 0.;
//      for(c=0;c<dim;c++) module += Beta[c*nvar+j]*Beta[c*nvar+j];
//      module = sqrt(module);
//      if(IsZero(module)) break;
      for(c=0;c<phr;c++) {
        ek(p+j,4*c+j) += weight*cfldiff*phi(c,0);
      }
      ef(p+j,0) -= weight*cfldiff*divflux[j];
    }
  }
 }
 else {
  TPZVec<REAL> divflux(4,0.);
  for(i=0;i<phr;i++) {
    p = 4*i;
    for(j=0;j<4;j++) {
//      module = 0.;
//      for(c=0;c<2;c++) module += Beta[c*nvar+j]*Beta[c*nvar+j];
//      module = sqrt(module);
//      if(IsZero(module)) break;
      ef(p+j,0) -= weight*cfldiff*divflux[j];
    }
  }
 }
}
*/
/*
int TEulerDiffusivity::main() {

  TEulerDiffusivity eul;
  TPZVec<REAL> U(4);
  U[0] = 5.4;
  U[1] = -.5;
  U[2] = -.25;
  U[3] = 50.8;
  TPZVec<REAL> flux(8);
  eul.Flux(U,flux);

	TPZVec<REAL> normal(2);
	normal[0] = 0.5;
	normal[1] = sqrt(.75);
	TPZFMatrix A(4,4);
	TPZFMatrix B(4,4);
	TPZFMatrix jacob(4,4);
	eul.JacobFlux(U,A,B);
	eul.JacobFlux(U,jacob,normal);
	A *= normal[0];
	B *= normal[1];
	int i,j;
	for(j=0;j<4;j++)
	   for(i=0;i<4;i++)
		    A(i,j) += B(i,j);
	
	eul.ValJacobFlux(U,jacob,normal);
	jacob.Print("ValJacob");
	TPZFMatrix axes(3,3,0.),jacinv(2,2,0.), ATA(4,4), ATB(4,4), BTA(4,4),BTB(4,4);
	axes(0,0)=1.5; axes(0,1)=2.3; axes(1,0)=.1; axes(1,1)=-1.9;
	jacinv(0,0) = 1.3;jacinv(0,1) = -2.4;jacinv(1,0)=3.5;jacinv(1,1)=5.2;
	eul.MatrixDiff(U,axes,jacinv,ATA,ATB,BTA,BTB);
	ATA.Print("ATA ");
	ATB.Print("ATB ");
	BTA.Print("BTA ");
	BTB.Print("BTB ");
  return 0;
}

*/

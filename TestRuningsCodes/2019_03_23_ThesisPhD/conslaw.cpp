/*******       FILE :   conslaw.c

Contains the definitions of the methods for the class TConservationLaw.

*******                               *******/

#include "conslaw.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "filedat1d.h"
#include "filedat2d.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzbndcond.h"
#include "bndcond.h"
#include "myheader.h"

TConservationLaw::TConservationLaw(int id,char *name) : TPZMaterial(id),
         fNumericalFlux(this,2) {
  SetName(name);
  fFluxType = 2;
  fMaxEigen = 0.;
  fExplicit = 0;
}
TConservationLaw::TConservationLaw(std::istream &input) : TPZMaterial(0), fNumericalFlux(this,1) {
  fFluxType = 2;
  SetData(input);
  MakeMesh(input);
  fMaxEigen = 0.;
  fExplicit = 0;
}
TConservationLaw::TConservationLaw(TConservationLaw &law) : TPZMaterial(law),
     fNumericalFlux(&law,2) {
  fFluxType = 2;
  fName = law.Name();
  fMaxEigen = 0.;
  fExplicit = 0;
}

/**Create an object TPZBndCond derived of TPZMaterial*/
TPZBndCond *TConservationLaw::CreateBC(int id,int typ,TPZFMatrix<REAL> &val1,TPZFMatrix<REAL> &val2) {
  return new TBC(this,id,typ,val1,val2);
}

REAL TConservationLaw::MinDeltaX() {
  return fMinDeltaX;
}

void TConservationLaw::SetName(char *name) {
  if(!name) {
    if(Dimension()==1) fName = "One-dimensional Conservation Law";
    else if(Dimension()==2) fName = "Two-dimensional Conservation Law";
    else if(Dimension()==3) fName = "Three-dimensional Conservation Law";
    else fName = "Conservation Law undefined";
    return;
  }
  fName = name;
}

void TConservationLaw::Print(std::ostream &out) {
  out << "CONSERVATION LAW" << std::endl << fName << std::endl;
  out << Id() << "\tOrder = " << NStateVariables() << std::endl;
}

void TConservationLaw::Clean() {
  fMaxEigen = 0.;
  fMinDeltaX = 0.;
  fAlfa = fDeltaT = fCurrentTime = 0.;
  fNumericalFlux.SetFluxType(fFluxTypeLower);
  fNumericalFlux.SetOrder(NStateVariables());
  fNumericalFlux.SetDimension(Dimension());
  fCoef = 0.;
  Point[0] = Point[1] = Point[2] = 0.;
}

void TConservationLaw::MakeMesh(std::istream &input) {
  char lbuf[256];
  int reg,dim = Dimension();
  input.get(lbuf[0]);
  input.getline(lbuf,256); //Regular mesh(1) or irregular mesh(0)
  input >> reg;
  input.get(lbuf[0]);
  input.getline(lbuf,256); //Data file name.

  TPZGeoMesh *gmesh = new TPZGeoMesh();
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

  if(dim==1) {
    //Inicializando malha computacional e geometrica
	  std::cout << std::endl << "Inicializing a geometrical and computational mesh 1d.";
    TDatafile1d RMesh(input);
    input.get(lbuf[0]);
    RMesh.InitializeMesh(cmesh,this);
	std::cout << "\nGeometrical and computational mesh 1d are done.\n";
  }
  else if(dim==2) {
	  std::cout << "\nInicializing a geometrical and computational mesh 2d.";
    TDatafile2d Mod(input);
    input.get(lbuf[0]);
    Mod.InitializeMesh(cmesh,this);
	std::cout << "\nComputational mesh 2d is done.\n";
  }
}

void TConservationLaw::SetData(std::istream &input) {
  GetDataCommented(input,Name(),256);
}

void TConservationLaw::SetNumericalFluxType(int typehigh,int typelower) {
  fFluxTypeHigher = typehigh;
  fFluxTypeLower = typelower;
  fNumericalFlux.SetFluxType(typelower);
  fNumericalFlux.SetOrder(NStateVariables());
  fNumericalFlux.SetDimension(Dimension());
	fFluxType = fFluxTypeLower;
}

void TConservationLaw::SetFluxTypeLower() {
  fFluxType = fFluxTypeLower;
}
void TConservationLaw::SetFluxTypeHigher() {
  fFluxType = fFluxTypeHigher;
}
/**
* @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
* @param data [in] stores all input data
* @param weight [in] is the weight of the integration rule
* @param ek [out] is the stiffness matrix
* @param ef [out] is the load vector
* @since April 16, 2007
*/
void TConservationLaw::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
	TPZVec<REAL> x = data.x; 
	TPZFMatrix<REAL> jacinv = data.jacinv; 
	TPZVec<REAL> sol = data.sol[0];
	TPZFMatrix<STATE> dsol = data.dsol[0];
	TPZFMatrix<REAL> axes = data.axes;
	TPZFMatrix<STATE> phi = data.phi; 
	TPZFMatrix<STATE> dphi = data.dphi;
	int phc, phr, dphc, dphr, efr, efc, ekr, ekc;
	int nvar = NStateVariables();
	phc = phi.Cols();
	phr = phi.Rows();
	dphc = dphi.Cols();
	dphr = dphi.Rows();   // It is dimension.
	if (!phr || dphr != Dimension()) return;

	/**To determine maxime jacobian eigenvalue to next time step */
	TPZVec<REAL> normal(2, 0.);
	normal[0] = 1.;
	REAL maxeigen = MaxEigJacob(sol, normal);
	if (fMaxEigen < maxeigen) fMaxEigen = maxeigen;

	SetPoint(x);

	efr = ef.Rows();
	efc = ef.Cols();
	ekr = ek.Rows();
	ekc = ek.Cols();
	if (phc != 1 || phr != dphc ||
		ekr != phr * nvar || ekc != ekr ||
		efr != ekr || efc != 1) {
		PZError << "\nTConservationLaw. Inconsistent input data : \n" <<
			"phi.Cols() = " << phi.Cols() << " dphi.Cols() = " << dphi.Cols() <<
			" phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
			dphi.Rows() << "\nek.Rows() = " << ek.Rows() << " ek.Cols() = "
			<< ek.Cols() <<
			"\nef.Rows() = " << ef.Rows() << " ef.Cols() = "
			<< ef.Cols() << "\n";
	}

	TPZVec<REAL> res(nvar, 0.), flux(dphr*nvar, 0.); //dphr is a spatial domain dimension
	REAL prod;

	/**To apply axes to dphi*/
	int i, j, l, c;
	TPZFMatrix<STATE> dphixy(dphr, phr, 0.);
	for (i = 0; i < phr; i++)
		for (c = 0; c < dphr; c++)
			for (l = 0; l < dphr; l++)
				dphixy(c, i) += dphi(l, i)*axes(l, c);

	if (fForcingFunction) {
		for (i = 0; i < nvar; i++) res[i] = sol[i];
		/**Observar que res contem a solucao Ul-1, e a funcao s(tk,x,sol(l-1))
		   toma a tempo tk da lei de conservacao*/
		fForcingFunction->Execute(x, res);
	}

	Flux(sol, flux);
	for (i = 0; i < phr; i++) {
		/**Calculo de keij, apenas o primeiro valor do bloco*/
		for (j = 0; j < phr; j++) {
			for (c = 0; c < nvar; c++) {
				ek(i*nvar + c, j*nvar + c) += (weight*phi(j, 0)*phi(i, 0));
			}
		}

		/**Calculo de feic*/
		for (c = 0; c < nvar; c++) {
			prod = 0.;
			for (l = 0; l < dphr; l++) prod += dphixy(l, i)*flux[l*nvar + c];
			ef(i*nvar + c, 0) += weight * (phi(i, 0)*res[c] + prod);
		}
	}
	/** Incrementando o termo difusivo segundo seja o esquema explicito ou implicito */
	/** Se a ordem de interpolação é zero (volumes finitos) não faz nada */
	if (nvar == phr)
		return;
	int k, p;

	//Computando o jacobiano da transformacao do elemento mestre ao elemento triangular
	TPZFMatrix<STATE> jacinvxy(dphr, dphr, 0.);
	for (i = 0; i < dphr; i++)
		for (j = 0; j < dphr; j++)
			for (k = 0; k < dphr; k++)
				jacinvxy(i, j) += (jacinv(i, k)*axes(k, j));

	int dsolr = dsol.Rows(), dsolc = dsol.Cols();
	REAL coeff = fCoef * fCFLDifusion*weight*fDeltaT;
	TPZFMatrix<STATE> vjacob(nvar, nvar);
	TPZFMatrix<STATE> matprod(nvar, nvar);
	TPZFMatrix<STATE> tau(nvar, nvar, 0.);
	TPZFMatrix<STATE> temp(nvar, nvar, 0.);
	TPZVec<REAL> Beta(dphr);

	/** Computa a matrix tau */
	if (dphr == 1) tau.Identity();
	else Tau(jacinvxy, sol, tau);
	tau *= coeff;

	if (fExplicit < 1) {
		for (i = 0; i < phr; i++) {
			for (k = 0; k < dphr; k++) Beta[k] = dphixy(k, i);
			JacobFlux(sol, vjacob, Beta);
			Multiply(vjacob, tau, matprod);
			for (j = 0; j < phr; j++) {
				for (k = 0; k < dphr; k++) Beta[k] = dphixy(k, j);
				JacobFlux(sol, vjacob, Beta);
				Multiply(matprod, vjacob, temp);
				for (k = 0; k < nvar; k++)
					for (p = 0; p < nvar; p++)
						ek(nvar*i + k, nvar*j + p) -= temp(k, p);   ///?? Mudanca de sinal 07/03/2003
			}
		}
	}
	else if (fExplicit > 1) {
		TPZVec<REAL> divflux(nvar, 0.);
		for (i = 0; i < phr; i++) {
			p = nvar * i;
			for (j = 0; j < nvar; j++) {
				//      module = 0.;
				//      for(c=0;c<pdhr;c++) module += Beta[c*nvar+j]*Beta[c*nvar+j];
				//      module = sqrt(module);
				//      if(IsZero(module)) break;
				for (c = 0; c < phr; c++) {
					ek(p + j, nvar*c + j) += weight * fCFLDifusion*phi(c, 0);
				}
				ef(p + j, 0) -= weight * fCFLDifusion*divflux[j];
			}
		}
	}
	else {
		TPZVec<REAL> divflux(nvar, 0.);
		for (i = 0; i < phr; i++) {
			p = nvar * i;
			for (j = 0; j < nvar; j++) {
				//      module = 0.;
				//      for(c=0;c<pdhr;c++) module += Beta[c*nvar+j]*Beta[c*nvar+j];
				//      module = sqrt(module);
				//      if(IsZero(module)) break;
				ef(p + j, 0) -= weight * fCFLDifusion*divflux[j];
			}
		}
	}
}

void TConservationLaw::Contribute(TPZVec<REAL> &x,TPZFMatrix<REAL> &jacinv,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,REAL weight,
     TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
	/** Evita trabalho desnecessario */
  int phc,phr,dphc,dphr,efr,efc,ekr,ekc;
  int nvar = NStateVariables();
  phc = phi.Cols();
  phr = phi.Rows();
  dphc = dphi.Cols();
  dphr = dphi.Rows();   // It is dimension.
	if(!phr || dphr!=Dimension()) return;
	
  /**To determine maxime jacobian eigenvalue to next time step */
  TPZVec<REAL> normal(2,0.);
  normal[0] = 1.;
  REAL maxeigen = MaxEigJacob(sol,normal);
  if(fMaxEigen<maxeigen) fMaxEigen = maxeigen;

  SetPoint(x);

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

  TPZVec<REAL> res(nvar,0.),flux(dphr*nvar,0.); //dphr is a spatial domain dimension
  REAL prod;

  /**To apply axes to dphi*/
  int i, j, l, c;
  TPZFMatrix<STATE> dphixy(dphr,phr,0.);
  for(i=0;i<phr;i++)
    for(c=0;c<dphr;c++)
      for(l=0;l<dphr;l++) 
        dphixy(c,i) += dphi(l,i)*axes(l,c);

  if(fForcingFunction) {
    for(i=0;i<nvar;i++) res[i] = sol[i];
    /**Observar que res contem a solucao Ul-1, e a funcao s(tk,x,sol(l-1))
       toma a tempo tk da lei de conservacao*/
    fForcingFunction->Execute(x,res);
  }

  Flux(sol,flux);
  for(i=0;i<phr;i++) {
    /**Calculo de keij, apenas o primeiro valor do bloco*/
    for(j=0;j<phr;j++) {
      for(c=0;c<nvar;c++) {
        ek(i*nvar+c,j*nvar+c) += (weight*phi(j,0)*phi(i,0));
      }
    }

    /**Calculo de feic*/
    for(c=0;c<nvar;c++) {
      prod = 0.;
      for(l=0;l<dphr;l++) prod += dphixy(l,i)*flux[l*nvar+c];
        ef(i*nvar+c,0) += weight*(phi(i,0)*res[c]+prod);
    }
  }
	/** Incrementando o termo difusivo segundo seja o esquema explicito ou implicito */
	/** Se a ordem de interpolação é zero (volumes finitos) não faz nada */
	if(nvar == phr) 
		return;
  int k, p;

  //Computando o jacobiano da transformacao do elemento mestre ao elemento triangular
  TPZFMatrix<STATE> jacinvxy(dphr,dphr,0.);
	for(i=0;i<dphr;i++)
		for(j=0;j<dphr;j++)
			for(k=0;k<dphr;k++)
				jacinvxy(i,j) += (jacinv(i,k)*axes(k,j));

  int dsolr = dsol.Rows(), dsolc = dsol.Cols();
  REAL coeff = fCoef*fCFLDifusion*weight*fDeltaT;
  TPZFMatrix<STATE> vjacob(nvar,nvar);
  TPZFMatrix<STATE> matprod(nvar,nvar);
  TPZFMatrix<STATE> tau(nvar,nvar,0.);
	TPZFMatrix<STATE> temp(nvar,nvar,0.);
	TPZVec<REAL> Beta(dphr);
	
	/** Computa a matrix tau */
	if(dphr==1) tau.Identity();
	else Tau(jacinvxy,sol,tau);
	tau *= coeff;

 if(fExplicit<1) {
 	for(i=0;i<phr;i++) {
		for(k=0;k<dphr;k++) Beta[k] = dphixy(k,i);
		JacobFlux(sol,vjacob,Beta);
		Multiply(vjacob,tau,matprod);
		for(j=0;j<phr;j++) {
			for(k=0;k<dphr;k++) Beta[k] = dphixy(k,j);
			JacobFlux(sol,vjacob,Beta);
			Multiply(matprod,vjacob,temp);
      for(k=0;k<nvar;k++)
        for(p=0;p<nvar;p++)
					ek(nvar*i+k,nvar*j+p) -= temp(k,p);   ///?? Mudanca de sinal 07/03/2003
		}
	}
 }
 else if(fExplicit>1) {
  TPZVec<REAL> divflux(nvar,0.);
  for(i=0;i<phr;i++) {
    p = nvar*i;
    for(j=0;j<nvar;j++) {
//      module = 0.;
//      for(c=0;c<pdhr;c++) module += Beta[c*nvar+j]*Beta[c*nvar+j];
//      module = sqrt(module);
//      if(IsZero(module)) break;
      for(c=0;c<phr;c++) {
        ek(p+j,nvar*c+j) += weight*fCFLDifusion*phi(c,0);
      }
      ef(p+j,0) -= weight*fCFLDifusion*divflux[j];
    }
  }
 }
 else {
  TPZVec<REAL> divflux(nvar,0.);
  for(i=0;i<phr;i++) {
    p = nvar*i;
    for(j=0;j<nvar;j++) {
//      module = 0.;
//      for(c=0;c<pdhr;c++) module += Beta[c*nvar+j]*Beta[c*nvar+j];
//      module = sqrt(module);
//      if(IsZero(module)) break;
      ef(p+j,0) -= weight*fCFLDifusion*divflux[j];
    }
  }
 }
}

void TConservationLaw::ContributeOverInterface(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZVec<REAL> &up1,
     REAL weight,REAL area,int type,TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &phi0,TPZFMatrix<STATE> &phi1,TPZVec<REAL> &normal,
     TPZFMatrix<STATE> &/*ek*/,TPZFMatrix<STATE> &ef) {

  SetPoint(x);
//	if(Dimension()==2) area = sqrt(area);
  /**To determine maxime jacobian eigenvalue to next time step */
  REAL maxeigen = MaxEigJacob(u,normal);
  if(fMaxEigen<maxeigen) fMaxEigen = maxeigen;
  maxeigen = MaxEigJacob(up1,normal);
  if(fMaxEigen<maxeigen) fMaxEigen = maxeigen;

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
        ef(i*nvar+c,0) = (-1.*phi0(i,0)*flux[c]*fCoef);
    }
  }
  phr1 += phr0;
  /** Computing ef */
  for(i=phr0;i<phr1;i++) {
    for(c=0;c<nvar;c++) {
        ef(i*nvar+c,0) = (phi1(i-phr0,0)*flux[c]*fCoef);
    }
  }
}

void TConservationLaw::SetPoint(TPZVec<REAL> &x) {
  int i, dim = x.NElements();
  for(i=0;i<dim;i++) Point[i] = x[i];
}
/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
 * @param data [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @param bc [in] is the boundary condition material
 * @since October 07, 2011
 */
void TConservationLaw::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {

}

/*Compute contribution to the stiffness matrix and right hand side at the integration point of a boundary*/
void TConservationLaw::ContributeBC(TPZVec<REAL> &/*x*/,TPZVec<REAL> &/*sol*/,REAL /*weight*/,TPZFMatrix<REAL> &/*axes*/,
	       TPZFMatrix<STATE> &/*phi*/,TPZFMatrix<STATE> &/*ek*/,TPZFMatrix<STATE> &/*ef*/,TPZBndCond &/*bc*/) {
}

/** To compute the tau matrix to diffusive term */
void TConservationLaw::Tau(TPZFMatrix<REAL> &jacinv,TPZVec<REAL> &sol,TPZFMatrix<STATE> &tau) {
	tau.Zero();
	int i, j, p, dim=jacinv.Rows(), nvar=NStateVariables();
	TPZVec<REAL> Beta(Dimension(),0.);
	TPZFMatrix<STATE> vjacob(nvar,nvar,0.);
	for(i=0;i<dim;i++) {
		for(j=0;j<dim;j++)
			Beta[j] = jacinv(j,i);
		p = ValJacobFlux(sol,vjacob,Beta);
		tau += vjacob;
	}
	if(!p) InverseJacob(tau);
	else tau.Identity();
}

/** To compute the inverse jacobian matrix */
void TConservationLaw::InverseJacob(TPZFMatrix<STATE> &mat) {
	mat(0,0) = 1./mat(0,0);
}

void TConservationLaw::VariablesName(TPZVec<std::string> &names) {
  names[0] = "state";
}

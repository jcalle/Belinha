/*******       FILE : ConsLaw.h   ConsLaw -> Conservation Law

Cont\'em a declara\c c\~ao da classe base "TConservationLaw". 
Esta \'e uma classe abstrata que permite resolver uma lei de conserva\c c\~ao
do tipo hiperb\'olico.
Para esta resolu\c c\~ao emplea-se o m\'etodo do volume finito e/ou
o m\'etodo de elementos finitos.
Aplica tamb\'em um algoritmo de multiresolu\c c\~ao sobre as medias celulares 
para detectar as singularidades da fun\c c\~ao solu\c c\~ao.
A lei de conserva\c c\~ao pode ser uni- ou bi-dimensional e linear ou n\~ao linear.
*******                                                *******/

#ifndef CONSERVATIONLAWH
#define CONSERVATIONLAWH

#include <stdlib.h>
#include <string.h>

#include "tpzmaterial.h"
#include "numflux.h"
#include "pzerror.h"

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZInterpolatedElement;

class TConservationLaw : public TPZMaterial {

  /** Alfa used into numerical flux computations */
  REAL fAlfa;
  /** Current Time = Tk*/
  REAL fCurrentTime;
 protected:
  REAL fDeltaT;
	REAL fCFLDifusion;

  /** To store conservation law name*/
  std::string fName;
  /** Minime of the element measures -> min measure(Cell) */
  REAL fMinDeltaX;
  /** Store the maxime jacobian eigenvalue over all solution values */
  REAL  fMaxEigen;

  /** Pointer to numerical flux to be will used*/
  TNumericalFlux fNumericalFlux;

 public:
  /** Flux types high and lower resolution to adapt over numerical flux boundary cells */
  int fFluxType;
  int fFluxTypeHigher;
  int fFluxTypeLower;

    /**CONSTRUCTORS and DESTRUTOR */
  TConservationLaw(int id,char *name=0);
  TConservationLaw(std::istream &input);
  TConservationLaw(TConservationLaw &law);

  /**Create an object TPZBndCond derived of TPZMaterial*/
  virtual TPZBndCond *CreateBC(int id,int typ, TPZFMatrix<REAL> &val1, TPZFMatrix<REAL> &val2);

  /** To set and put private data*/
  virtual double SonicPoint();
  std::string Name() { return fName; }
  void SetName(char *name);
  REAL Time() { return fCurrentTime; }
  REAL MinDeltaX();
  REAL MaxEigen() { return fMaxEigen; }
    
  REAL Alfa() { return fAlfa; }
  void SetAlfa(REAL alfa) { fAlfa = alfa; }
  void SetTime(double time,double deltat=0.);
  void SetMinDeltaX(double deltax) { fMinDeltaX = deltax; }
  void SetMaxEigen(REAL maxeigen) { fMaxEigen = maxeigen; }
  
  /** To set the numerical flux function choosed*/
  void SetNumericalFluxType(int typehigh,int typelower);
  void SetFluxTypeLower();
  void SetFluxTypeHigher();
	void SetFactor(REAL cfldif) { fCFLDifusion = cfldif; }
  TNumericalFlux &NumericalFlux() { return fNumericalFlux; }

  /** return the order of the system conservation law system*/
  int NFluxes() { return NStateVariables(); }

  //u tem m componentes, e flux(jacob, ...) tem N*m componentes
  virtual void Flux(TPZVec<REAL> &u,TPZVec<REAL> &flux)=0;
  virtual void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &jacob)=0;
	virtual void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal)=0;

  virtual void ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob);
  virtual int ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta);
  /**Implement the maxime eigenvalue of the jacobian matrix : n*jacob(u) */
  virtual REAL MaxEigJacob(TPZVec<REAL> &u,TPZVec<REAL> &normal)=0;
  virtual REAL ValEigJacob(TPZVec<REAL> &u,int order,int dim)=0;
	
  virtual void RoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &Roe);
  // EigRoe must to have Roe matrix, Roe eigenvalues, Roe matrix R and R inverse
  virtual void EigRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &EigRoe);
  virtual void ValRoeMatrix(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZFMatrix<STATE> &ValRoe);

  /**
 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);

  /**Compute contribution to the stiffness matrix and right hand side at an integration point*/
  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix<REAL> &,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,REAL weight,
	       TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
  /**Compute contribution to the stiffness matrix and right hand side at an integration point
  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol,
	REAL weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix
	&ek,TPZFMatrix &ef) { }*/
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 * @since October 07, 2011
	 */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
  /*Compute contribution to the stiffness matrix and right hand side at the integration point of a boundary*/
  void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,TPZFMatrix<REAL> &axes,
         TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

	/** To compute the tau matrix to diffusive term */
	virtual void Tau(TPZFMatrix<STATE> &jacinv,TPZVec<REAL> &sol,TPZFMatrix<STATE> &tau);
	/** To compute the inverse jacobian matrix */
	virtual void InverseJacob(TPZFMatrix<STATE> &mat);

  /**To construct appropiate computational mesh with interface elements*/
  void MakeMesh(std::istream &input);
  /**To evaluate true-error, L2-error and estimate error*/
  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<REAL> &flux,
        TPZVec<REAL> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val);

  /**To clean the conservation law data*/
  virtual void Clean();

  /**Print conservation law information*/
  void Print(std::ostream &out);

  /**To compute the contribution over edges of the elements*/
  void ContributeOverInterface(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZVec<REAL> &up1,
    REAL weight,REAL area,int type,TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &phileft,TPZFMatrix<STATE> &phiright,TPZVec<REAL> &normal,
    TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

  /**Read data of the material from a istream (file data)*/
  virtual void SetData(std::istream &data);

  int FluxType() { return fFluxType; }
	
  void SetExplicit(int explicito) { fExplicit = explicito; }
  void SetCoef(REAL coef) { fCoef = coef; }
  virtual void ApplyPeriodic(TPZCompMesh &cmesh) { }
  void SetPoint(TPZVec<REAL> &x);

  /** Return number of variables to print in postprocessing file 	*/
  virtual int AllVariablesPost() { return NStateVariables(); }
  virtual int AllVariablesImplemented() { return NStateVariables(); }
  /**To return the variables name of the material*/
  virtual void VariablesName(TPZVec<std::string> &names);

  /** Factor to diffussive term*/
  virtual int IdBC(double *x) { return 5; }

 protected:
  int        fExplicit;
  /** Coefficient to RungeKutta method : cl*fDeltaT */
  REAL       fCoef;
  /** Point to application*/
  double Point[3];
};

inline int TConservationLaw::ValJacobFlux(TPZVec<REAL> &/*u*/,TPZFMatrix<STATE> &/*valjacob*/,TPZVec<REAL> &/*beta*/) {
  PZError << "TConservationLaw::ValJacobFlux is called." << std::endl;
	return 0;
}
inline void TConservationLaw::ValJacobFlux(TPZVec<REAL> &/*u*/,TPZFMatrix<STATE> &/*valjacob*/) {
  PZError << "TConservationLaw::ValJacobFlux (void) is called." << std::endl;
}
inline void TConservationLaw::RoeMatrix(TPZVec<REAL> &/*u*/,TPZVec<REAL> &/*up1*/,TPZFMatrix<STATE> &/*Roe*/) {
  PZError << "TConservationLaw::RoeMatrix is called." << std::endl;
}
inline void TConservationLaw::EigRoeMatrix(TPZVec<REAL> &/*u*/,TPZVec<REAL> &/*up1*/,TPZVec<REAL> &/*EigRoe*/) {
  PZError << "TConservationLaw::EigRoeMatrix is called." << std::endl;
}
inline void TConservationLaw::ValRoeMatrix(TPZVec<REAL> &/*u*/,TPZVec<REAL> &/*up1*/,TPZFMatrix<STATE> &/*ValRoe*/) {
  PZError << "TConservationLaw::ValRoeMatrix is called." << std::endl;
}

inline void TConservationLaw::Errors(TPZVec<REAL> &,TPZVec<REAL> &,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,
				     TPZVec<REAL> &,TPZVec<REAL> &,TPZFMatrix<STATE> &,TPZVec<REAL> &) {
  PZError << "TConservationLaw::Errors is called." << std::endl;
}
inline double TConservationLaw::SonicPoint() {
  PZError << "TConservationLaw::SonicPoint is called.\n";
  return 0.;
}
inline void TConservationLaw::SetTime(double time,double deltat) {
  fCurrentTime = time;
	fDeltaT = deltat;
}

#endif

  /**Set initial function at time = 0
  void SetInitialFunction(void (*UZero)(TPZVec<REAL> &loc,TPZVec<REAL> &u)) {
    fUZero = UZero;
  }

//	void SetAlfa(REAL alfa) { fAlfa = alfa; }
//   double Alfa() { return fAlfa; }

//When fRegularMesh==1 then the first point of the mesh is the first node of
//the first element of the fElementMap
//	virtual void InitialPRiemann(double *U,int n);
//	virtual void InitialValueU();
//	virtual void InitialValueU(char * ) { }
//	virtual void EigRoeMatrixX(double *Ul,double *Ur,double *Roe) { }
//	virtual void ValRoeMatrixX(double *Ul,double *Ur,double *ValRoe) { }
// virtual int Material(int s) {return 0; }
//	void SetNonHomogeneous() { fNonHomogeneous=1; }
//	void SetU(double *U,double *U1) { fU=U; fU1=U1; }
//	void SetUi(double *U,int i);
//	void SetUij(double U,int i,int j) { fU[i*fOrdem+j]=U; }
//	int fPRiemann;		//Identificador do problema de Riemann
//	int fNonHomogeneous;	//Identifica equacao diferencial nao homogenea
//	int fDefinedRoeMatrix;	//Identifica se a matrix de Roe foi definida para a lei
//	int fDefinedGodunov;	//Identifica se definido o esquema de Godunov para a lei
//	double fAlfa;		  //fAlfa = dt/dx = k/h ou = k/Area(Celula)
//	double fCFL;		  //fCFL = |fAlfa*Max(f'(u))|
//	double fMaxCFL;     //Constante verifica condicao estabilidade(CFL)
//	int fHaveExact;		//Se diferente de zero, problema tem sol. exata
//	double fPontoSonico;	//Ponto onde f'(u)=0 (caso escalar nao linear)
//	TCompGrid *c;	  //Ponteiro para a malha computacional, o dominio espacial
//	void SetAlfa(double alfa) { fAlfa = alfa; }
//	void SetCFL(double CFL);
//	double Alfa() { return fAlfa; }
//	double CFL() { return fCFL; }

  /** Maxime eigenvalue of the jacobian at current time      -> 16/10/99
      This value is zero when is putting the current time 
  REAL fMaxEigen;
  void SetMaxEigen(REAL maxeig) { fMaxEigen = maxeig; }

  /**Return the maxime jacobian eigenvalue over all the solution values
  virtual REAL MaxEigen();
  
   */

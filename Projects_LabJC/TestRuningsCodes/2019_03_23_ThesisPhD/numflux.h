/*******   FILE :  numflux.h : NumFlux -> Numerical Flux

Header file for class TNumericalFlux. This class let us to calculate
of numerical flux across a edge in the discretized mesh, when we use
a conservative method for a conservation law (linear or non linear).

*******/

#ifndef NUMFLUXH
#define NUMFLUXH

#define ValAbs(x) ( (x) < 0 ? -(x) : (x) )
#define Max(x,y)  ( (x) > (y) ? (x) : (y) )
#define Min(x,y)  ( (x) > (y) ? (y) : (x) )
#define minmod(x,y) ( (x*y) > 0 ? ( ValAbs(x) < ValAbs(y) ? (x) : (y) ) : (0.) )
#define sgn(x) ( (x) > 0 ? (1) : ( (x) < 0 ? (-1) : (0) ) )
#ifndef ISZERO
#define ISZERO
#define IsZero( x )    ( (x) < 1.e-10 && (x) > -1.e-10 )
#endif
#define MAXORDER 5

#include <string.h>

template<class T>
class TPZVec;
class TConservationLaw;

class TNumericalFlux {

/*******DADOS MEMBRO*******/
 protected:
  TConservationLaw *fLaw;	 // pointer to conservation law

  std::string fName;          // Flux name
  int fFluxType;
  int fOrder;              // order of the conservation law
  double fMaxCFL;
  int fDim;                // geometrical dimension mesh
  enum MLimiter { ESuperbee, EMinMod, EVanLeer, EChakravarthy };
  MLimiter fLimiter;      // Identify type of limiter
  
  
/*******FUNCOES MEMBRO*******/
// ACCESS TO DATA MEMBERS :
 public:
  void SetName(char *Name);
  void SetOrder(int order) { fOrder = order; }
  void SetDimension(int dim) { fDim = dim; }
  virtual void ShowFluxType();
  virtual void SetFluxType(int t);
  void RequireFluxType();

  void SetLimiter(MLimiter En) { fLimiter = En; }
  double Limiter(double num,double den);
  
  std::string Name() { return fName; }
  int FluxType() { return fFluxType; }
  TConservationLaw *Equation() { return fLaw; }
  MLimiter Limiter() { return fLimiter; }

//CONSTRUTORES E DESTRUTOR

  TNumericalFlux(TConservationLaw *L,int type);
  TNumericalFlux(TConservationLaw *L,char *name,int type);
//  TNumericalFlux(TPZCompMesh *mesh,int t);
  TNumericalFlux(int order,int t);
  TNumericalFlux() { }
  virtual ~TNumericalFlux() {
  }

//Calcula o fluxo da celula Ui para a celula U1i
//No caso uni-dimensional o segundo ponteiro guarda o fluxo de Ui - Uip1
  virtual void NumericalFlux(TPZVec<REAL> &U0,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  void (TNumericalFlux::*fFlux)(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);

  virtual void FluxLaxFriedrichs(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxGodunov1(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxLaxWendroff(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxRichtmyerLW(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxMacCormack(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxLimiter(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxSlopeLimiter(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxUpStraightforward(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxUpVanLeer(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxUpRoe(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxENO(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxENO_NU(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxLFLWMin(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);
  virtual void FluxLFLWMax(TPZVec<REAL> &U,TPZVec<REAL> &Up1,TPZVec<REAL> &normal,TPZVec<REAL> &flux);

};


#endif


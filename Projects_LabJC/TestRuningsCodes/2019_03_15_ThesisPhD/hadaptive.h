#ifndef MYHPADAPTIVEHH
#define MYHPADAPTIVEHH

#include <iostream>
#include "pzreal.h"
#include "pzmanvector.h"

class TWavelet;

#define IsZeroToCoarsing(x) ( ((x) < 1.e-2 && (x) > -1.e-2) ?  1 : 0 )

class TPZCompMesh;
class TPZCompEl;
class TPZGeoEl;
class TPZGraphMesh;
class TPZInterpolatedElement;
class TTimeRungeKutta;
template<class T> class TPZFMatrix;
template<class T> class TPZMatrix;

//template<int> class TPZManVector;
template<class T> class TPZVec;

class TAdaptive {

  int hpadaptive;
  int MaxLevelRefinements;
	TWavelet *fWavelet;
	TWavelet *fWaveletSec;

  int p_decrement;
	int p_increment;
	int segurityzone; 

  TPZManVector<int> elfathers;
  TPZManVector<int> elsons;
	
	// When the current analysis is over a stationary problem
	int fSteadyState;

 public:
  TAdaptive(std::istream &input,int dim,int hadap = 0);
	~TAdaptive();

  /** To adapt numerical flux for high or lower resolution.
    elsons is a vector of computational elements over continuity region */
  void AdaptNumericalFlux(TPZCompMesh &cmesh,int fluxtype,TPZManVector<int> &geoel);
  /** Chama ao objeto da classe TWavelet para calcular os coeficientes wavelet e
  determinar se eh necessario refinar elementos da malha*/
  int IsNecessaryToRefine(TPZCompMesh &cmesh,int level,int onlyfluxadaptive);

  void WaveletHPAdaptive(TPZCompMesh &mesh);

  void WaveletDecomposition(TPZCompMesh &cmesh,TPZManVector<int> &geoel,int var,
	       int level);

  int DecrementOrder(TPZCompMesh &cmesh,TPZManVector<int> &elvector);
  int IncrementOrder(TPZCompMesh &cmesh);

  void Refinements(TPZCompMesh &cmesh,int ncycles,int minelem,int maxelem);
  void Refinements(TPZCompMesh &cmesh,int ncycles);
  void Refinements(TPZCompMesh &cmesh,TPZManVector<TPZGeoEl *> &elvec,int ncycles);

  void Coarsing(TPZCompMesh &cmesh);
	
	void CleanVectors();
	int MaxLevelRef() { return MaxLevelRefinements; }

  void CheckVectorsToAdaptivity(TPZCompMesh &cmesh,TPZManVector<int> &smalls);

  int SegurityZone(TPZCompMesh &cmesh,TPZManVector<int> &elvector);

  void AllIncrementInterpolation(TPZCompMesh &cmesh);
	
  // Minime level of refinements user informed.
  int MinLevel;
	
	void SteadyState(int steady) { fSteadyState = steady; }
};

#endif

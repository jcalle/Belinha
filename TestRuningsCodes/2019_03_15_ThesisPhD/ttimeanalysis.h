/*******       File : tarungek.h

Header file for class TTimeAnalysis. This class has a weighted residual analysis
to solve partial differential equations with first order time derivated. An object
of this class apply the L-order Runge-Kutta method.

*******              *******/

#ifndef TIMEANALYSISHH
#define TIMEANALYSISHH

#include <stdlib.h>
#include <string.h>

#include "pzanalysis.h"
#include "pzbndmat.h"
#include "pzfmatrix.h"

#define BIGMAXEIG 3.e2
#define BIGNUMBER 1.e5
#define MAXITERATIONS 3000000000

#ifdef PARALLELVERSION
#include "mpi.h"
#endif

class TConservationLaw;
class TPZInterpolatedElement;
class TAdaptive;
template<class T> class TPZVec;

extern int ReAssembling;

class TTimeAnalysis : public TPZAnalysis {

 protected:
//  int          fOrder;        //Order of the Runge-Kutta method
//  REAL         fCoef[4][2];   //Coefficients to operator H(u)
  REAL         fStartTime;    //Initial time of the time domain
	int          fNEndTimes;    //Number of the partial end times to write post-processing files
  REAL*        fEndTime;      //End values of the time domain
  char         fName[128];    //Associated name to time analysis
  REAL         fTimeStep;     //Time step intervale :  DeltaT
  REAL         fCFL;          //The spatial CFL to convergence
//  REAL         fCFLDiffussion; //The time CFL to increment diffusion
//	int          fExplicit;     // 0 (implicito) =1(explicito com difusividade) >1(explicito Runge-Kutta)
	int          fSteadyState;  // 0 se o processo é nao estacionario, 1 se estacionario

	/** Pointer to object does adaptivity */
	TAdaptive *fAdaptive;

  double MINIMETIMESTEP;    //Minimo passo de tempo permitido
  int MAXROWS;              //Maximo numero de linhas ao particionar as matrices

  TPZFMatrix<STATE>   *fUk;          //Solution vector to store solution in current time
  TConservationLaw  *fLaw;    //Conservation law to solve. Improve to several equations
  /** Pointer to initial function at T=0 */
  void (*fUZero)(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel);

  /** Parameter to use or not use limiter operator after computing solution */
//  int fUseLimiter;
  /** Norm of the residual vector */
  REAL fRhsNorm;
  /** Threashold value for aplying adaptivity on the norm */
  REAL fRhsAdaptNorm;
	TPZFMatrix<STATE> fPreviousSol;

 public:
  /** To indicate the level until the comp. mesh was refined */
  int fLevel;
	int fNPlots;

  TTimeAnalysis(std::istream &input, std::ostream &out,TPZCompMesh *cmesh,REAL TEnd=0.);
  TTimeAnalysis(std::istream &input, std::ostream &out,TPZCompMesh *mesh,int level);
  TTimeAnalysis(std::ostream &out);
  ~TTimeAnalysis();

  void SetDeltaT(double dt) { fTimeStep = dt; }
  REAL DeltaT() { return fTimeStep; }
  void SetName(char *name) { strcpy(fName,name); }
  void SetLaw(TConservationLaw *law) { fLaw = law; }
  void SetTimeDomain(double StartT,double EndT);
  double CurrentTime();
	void SetAdaptive(TAdaptive *adaptive) { fAdaptive = adaptive; }

  char *Name() { return fName; }
  double StartTime() { return fStartTime; }
  double EndTime(int n) { return fEndTime[n]; }
  double LastTimeReached() { return CurrentTime(); }
  int OrderLaw();

	/**Read complementary data to analysis */
	virtual void ReadData(std::ifstream &input);
  /**To clean TTimeAnalysis data members */
  virtual void CleanToStartRun();
	
	/**Put cfldiff as diffusion coefficient into the computational elements of the vector */
//	void SetFactor(TPZAdmChunkVector<TPZCompEl *> &elvec,REAL cfldiff);

  /**Compute time step from CFL condition and compute next time adjusting it at end time*/
  void ComputeNextTime(double currenttime,REAL EndTime);
  double StabilityCondition();

  virtual void AdjustBeforeAssemble(REAL ctime);

  /**To assembler only the vector on the right hand*/
  void AssembleRhs();
  void ReAssemble(REAL time);
  virtual void Run(std::istream &input, std::ostream &out);
  void Solve();
  void SolveToImplicit();
  /**Metodos adequados para resolver o sistema em banda em paralelo*/
  void SolveToImplicitChunks();
  void SolveToImplicitParallel();
  int MatrixTruncated(int n);

  /**Return maxime area and return minime area using pointer parameter*/
  REAL ComputeAreas(TPZCompMesh *cmesh,REAL *areamin);

  /**To put real dimensions in stiffness, rhs and solution matrices*/
  void AdequateMatrix();

  /**Applying adaptive scheme*/
  virtual int Adapting(int &step,int onlyflux);

  void Print(char *namevar, std::ostream &out= std::cout,int all=0);
  void PrintSolution(std::ostream &out= std::cout);
  void AppendMethodName(char *name);

  /**Compute initial function into fSolution of the mesh*/
  void ApplyUZero();
  void SetUZero(void (*uzero)(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel));

  /**To apply boundary condition*/
  void ApplyBC();
  void ApplyRhsBC();
  void ApplyBCDirichlet();
//  void ApplyBCNonFlux();
	
  void Integral(TPZInterpolatedElement *el,
         void (*fp)(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel),TPZVec<REAL> &integral);

  /**Operator Cockburn limiter - Generalized slope limiter*/
  /**It is used when the interpolated order is greater than zero : TPZCompEl::gOrder>0
  void CockburnLimiter(TPZAdmChunkVector<TPZCompEl *> &elvec,int var);
  void CockburnLimiter1d(TPZAdmChunkVector<TPZCompEl *> &elvec,int var);
  void CockburnLimiter2d(TPZAdmChunkVector<TPZCompEl *> &elvec,int var);

  /**Compute the minime module of the three real values
  REAL MinMod(REAL first,REAL second,REAL third);*/

  /** To Redefine the graphmeshes after a hp-adaptivity */
  void ReDefineGraphMesh(int dim);
  void ReDraw(int stepgraph,int ending=0);
  virtual void GetSchemeType(std::string &filename);
  /**Imprime estado atual da malha geometrica, computacional, solucao, etc*/
  void PrintCurrentState();

  /**Tentando melhorar o desempenho na resolucao da matriz banda fStiffness
  TPZFBMatrix *PartialStiff;
  TPZFMatrix *PartialRhs;
  TPZFMatrix *PartialSol;
*/	
#ifdef PARALLELVERSION
  int    fRank;   // process id
  int    fSize;   // number of process
  MPI_Status fStatus;
  int    fTypePar;
#endif
};

inline void TTimeAnalysis::AppendMethodName(char *name) {
  strcat(name," Runge-Kutta Time analysis.");
}

inline void TTimeAnalysis::SetUZero(void (*uzero)(TPZVec<REAL> &x,TPZVec<REAL> &u,
       TPZInterpolatedElement *cel)) {
  fUZero = uzero;
}

/** defines the postprocessing parameters for the graphical grid 
void DefineGraphMesh(int dimension, TPZVec<char *> &scalnames,TPZVec<char *> &vecnames,
		   TPZGraphMesh *graph);
*/
#endif

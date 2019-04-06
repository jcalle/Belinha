#ifndef MYHEADERFILEHPP
#define MYHEADERFILEHPP

#define PI_Value 3.1415926535897932

#include <iostream>
#include "pzreal.h"

#include "pzvec.h"
#include "pzmanvector.h"

class TPZCompMesh;
class TPZCompEl;
class TPZGeoEl;
class TPZGraphMesh;
class TPZInterpolatedElement;
class TTimeAnalysis;
template<class T>
class TPZFMatrix;
template<class T>
class TPZMatrix;


void GetCommentary(std::istream &input,int nlines = 1);
void GetDataCommented(std::istream &input,int &value);
void GetDataCommented(std::istream &input,int const nlines,int &value);
void GetDataCommented(std::istream &input,REAL &value);
void GetDataCommented(std::istream &input,TPZVec<REAL> &vector);
void GetDataCommented(std::istream &input,char *string,int size);
void GetDataCommented(std::istream &input, std::string &str, int size);

void SetCreateFunctionToElement(int dim,int type);
int BandWidth(TPZCompMesh *mesh);

char Itoa(int num);

extern int npausas;

void Pause(char *data,int nelem=0);
int Pause(int data);
void Pause(float data);
void Pause(double data);

void Refinements(TPZCompMesh &cmesh,int ncycles,int MaxLevel);
void Refinements(TPZCompMesh &cmesh,int ncycles,TPZVec<int> &typecel,
       int interpol=1,int MaxLevel=5);

void Multiply(TPZMatrix<STATE> &A,TPZMatrix<STATE> &B,TPZMatrix<STATE> &result);
void InvJacob2d(TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &invjac);

void CountElements(TPZCompMesh &cmesh,std::ostream &out= std::cout);
REAL MeanSolutionFromSubElements(TPZGeoEl *gel,int var);
REAL LinearApproximation(TPZVec<REAL> &qsi,TPZInterpolatedElement *el,int var);

/**Print information of the interfaces into comp. elements and the
   computational elements related with each interface*/
void PrintInterfaces(TPZCompMesh *cmesh,std::ostream &out);
/**Print vector elements*/
void PrintVector(TPZVec<int> &vec, std::ostream &out,char *vecname="Vector");
/**Print erros of several runnings in "mathematica" format*/
void PrintErrors(TPZVec<REAL> &err,int ngrouped,char *name);
/**Print erros of several runnings in "dx" format*/
void PrintErrorsDX(TPZVec<REAL> &err,int ngrouped,std::string &name);
/** To asign name to file for post-processing - name contains of law name */
void PlotFileName(char *name,int num,char *plotfile);
/** To asign name to file for post-processing - name contains of law name */
void PlotFileName(std::string name, int num, std::string plotfile);
int PrintDiagonal(TPZMatrix<STATE> *matrix, int mask = 0);
/** Recupera erros computados em rodadas anteriores*/
void GetOldErros(TPZVec<REAL> &errovec,int &erroini,char *ex_out);

/** To determine lateral values on connect with discontinuity
    n represents the number of connects into the computational element */
void PutDiscontinuousHyperplane(REAL term_c,REAL x=1.,REAL y=0.,REAL z=0.);
void PutSecondDiscontinuousHyperplane(REAL term_c=-1.e12,REAL x=1.,REAL y=0.,REAL z=0.);
int DiscontinuousHyperplane(TPZVec<REAL> &x);
int SecondDiscontinuousHyperplane(TPZVec<REAL> &x);

/** Return 1 if UZero is not right */
int PutUZero(int nuzero,TTimeAnalysis *an, std::istream &input,int dim=1);
/** Lateral intial functions when exist discontinuous hyperplanes */
void UZeroLeft(TPZVec<REAL> &x,TPZVec<REAL> &u);
void UZeroRight(TPZVec<REAL> &x,TPZVec<REAL> &u);
void UZeroRight_2(TPZVec<REAL> &x,TPZVec<REAL> &u);
/** Initial function with discontinuous hyperplanes */
void UZero(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel);
void UZero(TPZVec<REAL> &x,TPZVec<REAL> &u);
/** Initial continuous functions. Sinoidal functions*/
void UZeroSin(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel);
void UZeroBolha(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel);
void UZeroTwo(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel);
void UZeroSin2d(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel);
void UZeroCone(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel);
void UZeroNull(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel);
void UZeroUnity(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel);

/** Initial continuous function to circular displacement around the origen */
void UZeroLinearCircular(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel);

/** Global variables to Initial functions and Dirichlet conditions */
/** Determine whether the x point is into, at left or at right of the hyperplane
    a*x + b*y + c*z + d = 0   */
extern REAL Coef[];
extern REAL Term_c;
extern REAL SecondCoef[];
extern REAL SecondTerm_c;
/** Global variables to Dirichlet conditions */
extern REAL u0[];   // To Dirichlet values

/**Detect the Dirichlet values from boundary conditions into the cmesh */
REAL DetectDirichletValue(TPZCompMesh *cmesh,int nvar);
/**Count which different Dirichlet conditions has the boundary conditions */
int NumberOfDirichletConditions(TPZCompMesh &cmesh);

/**Make continuous all connects into the cmesh except boundary connects*/
int AllConnectContinuous(TPZCompMesh *cmesh);
/**Make discontinuous all connects into the cmesh except boundary connects*/
int AllConnectDiscontinuous(TPZCompMesh *cmesh);
/**Create a *.xls file with excell format, to data visualization */
void CreateXLS(char *name,int nobjects=1,int ncoords=1);

int NoZeroNumbers(TPZVec<int> &vec);
int DifferentNumbers(TPZVec<int> &vec);

extern std::string plotfile;   // nome do arquivo para pos-processamento

#endif

#ifndef TECHNIQUEWAVELET
#define TECHNIQUEWAVELET

#include <iostream>
#include "pzreal.h"
#include "pzmanvector.h"

template<class T> class TPZVec;

class TPZCompMesh;
class TTimeAnalysis;
class TPZGeoEl;
template<class T>
class TPZFMatrix;

class TWavelet {
 protected:
  double THREADHOLDER;
	double MAXWAVELETCOEFFICIENT;
	
 public:
	TWavelet(double ,double );
	~TWavelet() { }
	
  virtual int Decomposition(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
     TPZManVector<int> &elsons,int minlevel) = 0;
  virtual void DecompositionOnlyOrderZero(TPZCompMesh &cmesh,int var,
	     TPZManVector<int> &elfathers,TPZManVector<int> &elsons,int minlevel);

  virtual REAL OneLevelDetails(TPZVec<REAL> &values) = 0;
  virtual REAL OneLevelEven(TPZVec<REAL> &value);
  virtual REAL OneLevelEven(REAL value);

  void PrintDecomposition(int var,TPZVec<REAL> &coefs,int nsubel, std::ostream &out);

  virtual void DrawWavelets(int var,TTimeAnalysis &an,int dim,int &step,REAL &time,
		      TPZVec<std::string> &scal,TPZVec<std::string> &vec,int last);

};

inline void TWavelet::DrawWavelets(int var,TTimeAnalysis &an,int dim,int &step,REAL &time,
		      TPZVec<std::string> &scal,TPZVec<std::string> &vec,int last) {

	std::cout << "TWavelet::DrawWavelets is called.\n";
}
inline REAL TWavelet::OneLevelEven(REAL value) {
	std::cout << "TWavelet::OneLevelEven is called.\n";
	return 0.;
}
inline REAL TWavelet::OneLevelEven(TPZVec<REAL> &value) {
	std::cout << "TWavelet::OneLevelEven (vector) is called.\n";
	return 0.;
}
inline void TWavelet::DecompositionOnlyOrderZero(TPZCompMesh &cmesh,int var,
       TPZManVector<int> &elfathers,TPZManVector<int> &elsons,int minlevel) {
	std::cout << "TWavelet::DecompositionOnlyOrderZero is called.\n";
}


class THaar1D : public TWavelet {

 public:
  THaar1D(double ,double);
	
  int Decomposition(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
     TPZManVector<int> &elsons,int minlevel);

  void DecompositionOnlyOrderZero(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
     TPZManVector<int> &elsons,int minlevel);
  REAL OneLevelEven(TPZVec<REAL> &values);   // at even is the mean value
  REAL OneLevelDetails(TPZVec<REAL> &values);

  void Decomposition(TPZCompMesh &cmesh,int var,TPZManVector<REAL> &coeff,
			TPZManVector<TPZGeoEl *> &els,int level);

  void DrawWavelets(int var,TTimeAnalysis &an,int dim,int &step,REAL &time,
		      TPZVec<std::string> &scal,TPZVec<std::string> &vec,int last);

};

class THaar2D : public TWavelet {

 public:
  THaar2D(double ,double);

  int Decomposition(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
       TPZManVector<int> &elsons,int minlevel);
  void DecompositionOnlyOrderZero(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
       TPZManVector<int> &elsons,int minlevel);
  REAL OneLevelEven(TPZVec<REAL> &values);
  REAL OneLevelDetails(TPZVec<REAL> &values);

};

class TSchauder1D : public TWavelet {

 public:
  TSchauder1D(double ,double);
	
  int Decomposition(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
     TPZManVector<int> &elsons,int minlevel);

  void Preconditioning(int var,int dim,TPZFMatrix<STATE> &pre,int nlevels);
  REAL OneLevelEven(REAL value);
  REAL OneLevelDetails(TPZVec<REAL> &values);

  void DrawWavelets(int var,TTimeAnalysis &an,int posdim,int &step,
			  REAL &time,TPZVec<std::string> &scal,TPZVec<std::string> &vec,int last);

};

class TSchauder2D : public TWavelet {

 public:
  TSchauder2D(double ,double);
	
  int Decomposition(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
       TPZManVector<int> &elsons,int minlevel);

  REAL OneLevelEven(REAL value);
  REAL OneLevelDetails(TPZVec<REAL> &values);

  void DrawWavelets(int var,TTimeAnalysis &an,int posdim,int &step,
			  REAL &time,TPZVec<std::string> &scal,TPZVec<std::string> &vec,int last);

};

#endif

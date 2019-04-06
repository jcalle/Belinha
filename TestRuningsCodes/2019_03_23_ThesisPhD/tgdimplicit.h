#ifndef TGALERKINDESCONTINUOUSIMPLICITHH
#define TGALERKINDESCONTINUOUSIMPLICITHH

#include "ttimeanalysis.h"

class TGDImplicit : public TTimeAnalysis {

  REAL  fCFLDiffussion; //The time CFL to increment diffusion

 public:
  TGDImplicit(std::istream &input, std::ostream &out,TPZCompMesh *mesh,int level);
  ~TGDImplicit() { }

  void Run(std::istream &input, std::ostream &out);

	void AdjustBeforeAssemble(REAL ctime);
	void CleanDiffusion();
	void CleanToStartRun();

	/**Read complementary data to analysis */
	virtual void ReadData(std::ifstream &input);
  /**Applying adaptive scheme*/
  virtual int Adapting(int &step,int onlyflux);

  virtual void GetSchemeType(char *filename);
	
	/** Acelerando integracao sobre os elementos finitos */
	TPZFMatrix<STATE> fSolutionNext;
  void AssembleOnlyIntegral(TPZFMatrix<STATE> &rhs);
};

#endif

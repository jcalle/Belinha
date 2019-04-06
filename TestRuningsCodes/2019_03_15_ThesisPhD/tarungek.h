#ifndef TIMERUNGEKUTTAHH
#define TIMERUNGEKUTTAHH

#include "ttimeanalysis.h"

class TTimeRungeKutta : public TTimeAnalysis {

  int          fOrder;        //Order of the Runge-Kutta method
  REAL         fCoef[4][2];   //Coefficients to operator H(u)

  /** Parameter to use or not use limiter operator after computing solution */
  int fUseLimiter;

 public:
  TTimeRungeKutta(std::istream &input, std::ostream &out,TPZCompMesh *mesh,REAL TEnd=0.);
  TTimeRungeKutta(std::istream &input, std::ostream &out,TPZCompMesh *mesh,int level);
  TTimeRungeKutta(std::ostream &out);
  ~TTimeRungeKutta() { }
  void SetUseLimiter(int use) { fUseLimiter = use; }

  void Run(std::istream &input, std::ostream &out);
//  void Solve();

  /**Operator Cockburn limiter - Generalized slope limiter*/
  /**It is used when the interpolated order is greater than zero : TPZCompEl::gOrder>0*/
  void CockburnLimiter(TPZAdmChunkVector<TPZCompEl *> &elvec,int var);
  void CockburnLimiter1d(TPZAdmChunkVector<TPZCompEl *> &elvec,int var);
  void CockburnLimiter2d(TPZAdmChunkVector<TPZCompEl *> &elvec,int var);

  /**Compute the minime module of the three real values*/
  REAL MinMod(REAL first,REAL second,REAL third);

	/**Read complementary data to analysis */
	virtual void ReadData(std::ifstream &input);

  virtual void GetSchemeType(char *filename);
};

#endif

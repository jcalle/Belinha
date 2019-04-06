#ifndef EULERDIFFUSIONHH
#define EULERDIFFUSIONHH


#include "pzvec.h" 
#include "pzfmatrix.h"

class TEulerDiffusivity {

static  REAL  fGamma;
	
 public:

static  REAL Pressure(TPZVec<REAL> &U);

static  void Flux(TPZVec<REAL> &u,TPZVec<REAL> &flux);
static  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &Ajacob,TPZFMatrix<STATE> &Bjacob);
static	void JacobFlux(TPZVec<REAL> &U,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal);

static  void ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal);

static  void MatrixDiff(TPZVec<REAL> &sol, TPZFMatrix<REAL> &axes, TPZFMatrix<STATE> &jacinv,TPZFMatrix<STATE>
&ATauA,TPZFMatrix<STATE> &ATauB,TPZFMatrix<STATE> &BTauA,TPZFMatrix<STATE> &BTauB);

static  void InvJacob2d(TPZFMatrix<STATE> &axes,TPZFMatrix<STATE> &jacinv);

static  void InverseJacob(TPZFMatrix<STATE> &jac);

static int main();
	
};

#endif

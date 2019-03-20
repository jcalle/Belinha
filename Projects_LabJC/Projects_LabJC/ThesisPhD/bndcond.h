
#ifndef BOUNDCONDITIONHPP
#define BOUNDCONDITIONHPP

#include "pzbndcond.h"
#include "tpzmaterial.h"

class TBC : public TPZBndCond {

//  TBC(TBC &bccopy) : TPZBndCond(bcopy) { } 
 public:
  ~TBC() { }

  TBC(TPZMaterial *material,int id,int type,TPZFMatrix<REAL> &val1,TPZFMatrix<REAL> &val2) :
    TPZBndCond(material,id,type,val1,val2) { }

  /**returns the variable index associated with the name*/
  int VariableIndex(char *name) { return Material()->VariableIndex(name); }

  /** returns the number of variables associated with the variable indexed by var.
      var is obtained by calling VariableIndex*/
  int NSolutionVariables(int var) {
    return fMaterial->NSolutionVariables(var); }

  /**To return the variables name of the material*/
//  void VariablesName(TPZVec<char *> &names) { Material()->VariablesName(names); }

  /**returns the solution associated with the var index based on the finite element approximation*/
  virtual void Solution(TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, int var, TPZVec<STATE> &Solout)
  {
	  TPZMaterialData data;
	  data.sol[0] = Sol;
	  data.dsol[0] = DSol;
	  data.axes = axes;
    fMaterial->Solution(data,var,Solout); }

  /**compute the value of the flux function to be used by ZZ error estimator*/
  void Flux(TPZVec<REAL> &x,TPZVec<REAL> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,TPZVec<REAL> &flux) {
    Material()->Flux(x,Sol,DSol,axes,flux); }

  /**Compute contribution to the stiffness matrix and right hand side at an integration point*/
  void Contribute(TPZVec<REAL> &x,TPZFMatrix<REAL> &jac,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol,REAL weight,TPZFMatrix<REAL> &axes,
      TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    ((TConservationLaw *)Material())->Contribute(x,jac,sol,dsol,weight,axes,phi,dphi,ek,ef); }

  /*Compute contribution to the stiffness matrix and right hand side at the integration point of a boundary*/
  void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,TPZFMatrix<REAL> &axes,
			    TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	  TPZMaterialData data;
	  data.x = x;
	  data.sol[0] = sol;
	  data.axes = axes;
	  data.phi = phi;
    Material()->ContributeBC(data,weight,ek,ef,bc); }

  /**Compute the error due to the difference between the interpolated flux and
  the flux computed based on the derivative of the solution*/
  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
		      TPZVec<REAL> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val) {
    Material()->Errors(x,sol,dsol,axes,flux,uexact,duexact,val); }

};

#endif


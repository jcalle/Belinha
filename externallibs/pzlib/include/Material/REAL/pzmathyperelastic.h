/**
 * @file
 * @brief Contains the TPZMatHyperElastic class which implements a hyper elasticity material.
 */

#ifndef MATHYPERELASTICHPP
#define MATHYPERELASTICHPP

#include "TPZMaterial.h"
#include "pzfmatrix.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

/**
 * @ingroup material
 * @brief Implements a hyper elasticity material.
 */
class TPZMatHyperElastic : public TPZMaterial {
	
    STATE fXf[3];
#ifndef _AUTODIFF
	STATE fK2[3][3], fK3[3][3], fK4[3][3], fK6[3][3], fK7[3][3], fK8[3][3];
	STATE fL1[3][3],fL2[3][3],fL3[3][3],fL4[3][3],fL5[3][3],fL6[3][3],fL7[3][3],fL8[3][3],fL9[3][3],fGradtrC[3][3];
    STATE fGradDetF[3][3];
#endif
	STATE fE1[3],fE5[3],fE9[3];

    
	STATE fLambda,fNu,fE,fMu;
	STATE fCoef1,fCoef2,fCoef3;
	
	public :
	
	TPZMatHyperElastic(int nummat,STATE e,STATE mu,STATE nu=-1.,STATE lambda=-1.,STATE coef1=-1.,STATE coef2=-1.,STATE coef3=-1.);
	
	virtual ~TPZMatHyperElastic();
	
	void SetMaterial(TPZFMatrix<STATE> &xfin){
		fXf[0] = xfin(0,0);
		fXf[1] = xfin(1,0);
		fXf[2] = xfin(2,0);
	}
	
	int Dimension() const { return 3;}
	
	int NStateVariables();
	
	virtual void Print(std::ostream & out);
	
	std::string Name() { return "TPZMatHyperElastic"; }
	
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
	
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
#ifdef _AUTODIFF
	
	/** @brief Computes contribution to the energy at an integration point*/
	virtual void ContributeEnergy(TPZVec<REAL> &x,
								  TPZVec<FADFADREAL> &sol,
								  TPZVec<FADFADREAL> &dsol,
								  FADFADREAL &U,
								  REAL weight);
	
	static void ComputeEnergy(STATE lambda, STATE mu,  TPZFMatrix<STATE> &dsol, TFad<9,TFad<9,STATE> > &energy);
	
	/** @brief Computes contribution of BC to the Energy*/
	virtual void ContributeBCEnergy(TPZVec<REAL> & x,
									TPZVec<FADFADREAL> & sol, FADFADREAL &U,
									REAL weight, TPZBndCond &bc);
	
#endif
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 9;}
	
protected:
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	/** @brief Computes the value of the flux function to be used by ZZ error estimator */
	virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
		        TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);//Cedric
    public:
virtual int ClassId() const;
 
};

#endif

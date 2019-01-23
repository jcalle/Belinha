/**
 * \file
 * @brief Contains the TPZCoupledTransportDarcy class which implements two equations to transport problem.(Jorge?)
 */

#ifndef MATCOUPLEDTRANSPDARCY
#define MATCOUPLEDTRANSPDARCY

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include "pzpoisson3d.h"

class TPZCoupledTransportDarcyBC;
/**
 * @ingroup material
 * @brief Implements two equations where the second one requires the solution of the first.
 */
/** 
 * First equation is  \f$ -div[ K Grad[p] ] = 0 \f$ \n
 * Second equation is \f$ -div[ A Grad[u] ] + beta . Grad[u] = 0 \f$ \n
 * where \f$ K \f$ and \f$ A \f$ are constants and \f$ beta = alpha * (-K Grad[p] ) \f$ \n
 * Both equations are implemented by TPZMatPoisson3d material. \n
 * This material has two instances of TPZMatPoisson3d. In a first moment
 * it behaves like first equation. \n Once the problem is solved, it behaves
 * like the second equation. \n In the latter case beta value is set when Contribute
 * is called by processing the solution of the previous equation.
 */
class TPZCoupledTransportDarcy : public TPZDiscontinuousGalerkin {
	
	protected :
	
	/** @brief In second equation: \f$ beta = alpha * (-K Grad[p] ) \f$. Here alpha is stored. */
	REAL fAlpha;
	
	/** @name Two instances of TPZMatPoisson3d */
	// @{
	TPZMaterial * fMaterialRefs[2];
	TPZMatPoisson3d * fMaterials[2];
	// @}
	
	static int gCurrentEq;
	
public:
	
	virtual TPZBndCond *CreateBC(TPZMaterial * mat, int id, int typ, TPZFMatrix<STATE> &val1,TPZFMatrix<STATE> &val2){
		PZError << "Error! - This method should not be called - " << __PRETTY_FUNCTION__ << std::endl;
		return 0;
	}
	
	TPZCoupledTransportDarcyBC *CreateBC2(int id);
	
	static void SetCurrentMaterial(const int i);
	
	static int CurrentEquation();
	
	void UpdateConvectionDir(TPZFMatrix<STATE> &dsol);
	
	void UpdateConvectionDirInterface(TPZFMatrix<STATE> &dsolL, TPZFMatrix<STATE> &dsolR);
	
	virtual int HasForcingFunction() {return this->GetCurrentMaterial()->HasForcingFunction();}
	
	TPZMatPoisson3d * GetCurrentMaterial() const {
#ifdef PZDEBUG
		if (!this->fMaterials[0] || !this->fMaterials[1]){
			PZError << "Error! - " << __PRETTY_FUNCTION__ << std::endl;
			exit (-1);
		}
#endif
		if (TPZCoupledTransportDarcy::gCurrentEq == 0 ||
			TPZCoupledTransportDarcy::gCurrentEq == 1) return this->fMaterials[TPZCoupledTransportDarcy::gCurrentEq];
		else {
			PZError << "Error! - " << __PRETTY_FUNCTION__ << std::endl;
			exit (-1);
		}
		// the code will never reach this point
		return 0;
	}
	
	TPZMatPoisson3d *GetMaterial(int eq){
#ifdef PZDEBUG
		if (!this->fMaterials[0] || !this->fMaterials[1]){
			PZError << "Error! - " << __PRETTY_FUNCTION__ << std::endl;
			exit (-1);
		}
#endif
		if (eq == 0 || eq ==1) return this->fMaterials[eq];
		else {
			PZError << " Error - " << __PRETTY_FUNCTION__ << std::endl;
			exit (-1);
		}
		// the code will never reach this point
		return 0;
	}
	
	TPZCoupledTransportDarcy(int nummat, int nummat0, int nummat1, int dim);
	
	virtual ~TPZCoupledTransportDarcy();
	
	TPZCoupledTransportDarcy(TPZCoupledTransportDarcy &copy) : TPZRegisterClassId(&TPZCoupledTransportDarcy::ClassId),TPZDiscontinuousGalerkin(copy),
    fAlpha(copy.fAlpha){
		this->fMaterialRefs[0] = copy.fMaterialRefs[0];
		this->fMaterialRefs[1] = copy.fMaterialRefs[1];
		fMaterials[0] = copy.fMaterials[0];
		fMaterials[1] = copy.fMaterials[1];
	}
	
	virtual TPZMaterial * NewMaterial(){
		return new TPZCoupledTransportDarcy(*this);
	}
	
	int Dimension() const {
		return this->GetCurrentMaterial()->Dimension();
	}
	
	int NStateVariables();
	
	void SetAlpha(REAL alpha);
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZCoupledTransportDarcy"; }
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc);
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ef)
	{
		TPZDiscontinuousGalerkin::Contribute(data,weight,ef);
	}
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc)
	{
		TPZDiscontinuousGalerkin::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 3;}
	
protected:
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
public:
	/**
	 * @brief Returns the solution associated with the var index based on
	 * the finite element approximation
	 */
	virtual void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<STATE> &Solout)
	{
		TPZDiscontinuousGalerkin::SolutionDisc(data,dataleft,dataright,var,Solout);
	}
	
	
	/** @brief Compute the value of the flux function to be used by ZZ error estimator */
	virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);//Cedric
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix<STATE> &ek,
									 TPZFMatrix<STATE> &ef);
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ek,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix<STATE> &ef)
	{
		TPZDiscontinuousGalerkin::ContributeInterface(data,dataleft,dataright,weight,ef);
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc)
	{
		TPZDiscontinuousGalerkin::ContributeBCInterface(data,dataleft,weight,ef,bc);
	}
	
	void InterfaceErrors(TPZVec<REAL> &/*x*/,
						 TPZVec<STATE> &leftu, TPZFMatrix<STATE> &leftdudx, /* TPZFMatrix<REAL> &leftaxes,*/
						 TPZVec<STATE> &rightu, TPZFMatrix<STATE> &rightdudx, /* TPZFMatrix<REAL> &rightaxes,*/
						 TPZVec<STATE> &/*flux*/,
						 TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values,
						 TPZVec<REAL> normal, REAL elsize);
	
	virtual int IsInterfaceConservative(){ return 1;}
        
        virtual int ClassId() const;

};

#endif

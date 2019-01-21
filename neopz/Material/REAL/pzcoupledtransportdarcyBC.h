/**
 * \file
 * @brief Contains the TPZCoupledTransportDarcyBC class.
 */

#ifndef MATCOUPLEDTRANSPDARCYBC
#define MATCOUPLEDTRANSPDARCYBC

#include <iostream>

#include "pzreal.h"
#include "pzbndcond.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzcoupledtransportdarcy.h"

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZCoupledTransportDarcyBC : public TPZBndCond{
	
// #warning THIS CLASS IS NOT THREADSAFE!!!
protected:
	
	TPZBndCond * fMaterials[2];
	
	TPZBndCond * GetNonNullMaterial() const {
		if (this->fMaterials[0]) return this->fMaterials[0];
		if (this->fMaterials[1]) return this->fMaterials[1];
		PZError << "Error! - "  << __PRETTY_FUNCTION__ << std::endl;
		exit (-1);
		// the code will never reach this point
		return 0;
	}
	
	void UpdateConvectionDir(TPZFMatrix<STATE> &dsol);
	void UpdateConvectionDirInterface(TPZFMatrix<STATE> &dsolL, TPZFMatrix<STATE> &dsolR, TPZFMatrix<REAL> &phiL, TPZFMatrix<REAL> &phiR);
	
	public :
	
	TPZCoupledTransportDarcyBC(TPZCoupledTransportDarcy * material, int id);
	
	~TPZCoupledTransportDarcyBC();
	
	TPZBndCond * GetCurrentMaterial(){
		const int eq = TPZCoupledTransportDarcy::CurrentEquation();
		if (eq == 0 || eq == 1) return this->fMaterials[eq];
		else {
			PZError << "Error! - " << __PRETTY_FUNCTION__ << std::endl;
			exit (-1);
		}
		// the code will never reach this point
		return 0;
	}
	
	virtual int HasForcingFunction() {
		TPZBndCond * bc = this->GetCurrentMaterial();
		if (bc) return bc->HasForcingFunction();
		return 0;
	}
	
	void SetMaterial(int eq, TPZBndCond * mat){
		if (eq == 0 || eq == 1) this->fMaterials[eq] = mat;
		else {
			PZError << "Error! - " << __PRETTY_FUNCTION__ << std::endl;
			exit (-1);
		}
	}
	
	/** @brief Returns the integrable dimension of the material */
	int Dimension() const {
		return this->GetNonNullMaterial()->Dimension();
	}
	
	virtual int NFluxes(){ return this->GetNonNullMaterial()->NFluxes(); }
	
	int NStateVariables() { return this->GetNonNullMaterial()->NStateVariables(); }
	
	/** @brief Returns the number of norm errors. Default is 3: energy, L2 and H1. */
	virtual int NEvalErrors() {return this->GetNonNullMaterial()->NEvalErrors();}
	
	/** @brief Computes the value of the flux function to be used by ZZ error estimator */
	void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux){
		flux.Fill(0.);
	}
	
	void Print(std::ostream & out = std::cout) {
		out << " Boundary condition number = " << Id() << "\n";
	}
	
	void Contribute(TPZMaterialData &data,
					REAL weight,
					TPZFMatrix<STATE> &ek,
					TPZFMatrix<STATE> &ef);
	
	void Contribute(TPZMaterialData &data,
					REAL weight,
					TPZFMatrix<STATE> &ef)
	{
		TPZBndCond::Contribute(data,weight,ef);
	}
	
	void ContributeBC(TPZMaterialData &data,
					  REAL weight,
					  TPZFMatrix<STATE> &ek,
					  TPZFMatrix<STATE> &ef,
					  TPZBndCond &bc) {  }
	
    void ContributeBC(TPZMaterialData &data,
					  REAL weight,
					  TPZFMatrix<STATE> &ef,
					  TPZBndCond &bc)
	{
		TPZBndCond::ContributeBC(data,weight,ef,bc);
	}
	
	
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
				TPZVec<STATE> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val){
		val.Fill(0.);
	}
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                     REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                     REAL weight,
                                     TPZFMatrix<STATE> &ef);
	
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ek,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc) {
		//NOTHING TO BE DONE HERE
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc)
	{
		TPZBndCond::ContributeBCInterface(data,dataleft,weight,ef,bc);
	}
public:
virtual int ClassId() const;

};

#endif

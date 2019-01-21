/**
 * \file
 * @brief Contains the TPZIncNavierStokesKEps class which implements an imcompressible Navier-Stokes formulation with modified KEpsilon turbulence model.
 */

#ifndef PZINCNSKEPS
#define PZINCNSKEPS

#include "TPZMaterial.h"

#include <iostream>
#include <string>

#ifdef _AUTODIFF
#include "fadType.h"
#endif

class TPZBndCond;

/**
 * @ingroup material
 * @brief This class implements an imcompressible Navier-Stokes formulation with modified KEpsilon turbulence model.
 * @author Professor Paulo Vatavuk
 * @author Professor Philippe Devloo
 * @author Edimar Cesar Rylo
 * @author Luis Fernando
 * @author Roberto H. Heiderich 
 * @author Tiago Forti
 * @since June 29, 2005
 */
/** 
 * Variables are: {K, Eps, Pressure, Vx, Vy, Vz}.
 * This class is homework:
 */
class TPZIncNavierStokesKEps : public TPZMaterial {
	
private:
	
    int fDimension;
    
    STATE fMU, fRHO, fCmu, fSigmaK, fSigmaEps, fCepsilon1, fCepsilon2;
    
    TPZVec<STATE> fBodyForce; //fc = {0,0,-fc}
    
    /** @brief Dot for matrices with same dimensions. \f$ Tr[A B]. \f$ No consistence test is made. */
    STATE Dot(TPZFMatrix<STATE> &A, TPZFMatrix<STATE> &B);
    
    /** @brief Dot of vector A with row BRow of matrix B. */
    STATE Dot(TPZVec<STATE> &A, TPZFMatrix<STATE> &B, int BRow);
	
    /** @brief Dot for vectors with same dimensions. No consistence test is made. */        
    STATE Dot(TPZVec<STATE> &A, TPZVec<STATE> &B);
	
public:
	
    enum VARIABLES {ENone = -1, EK = 0, EEpsilon, EPressure, EVx, EVy, EVz, EVvector};
	
    TPZIncNavierStokesKEps(int id, int dimension);
    
    void SetParameters(STATE MU, STATE RHO, STATE Cmu, STATE SigmaK, STATE SigmaEps, STATE Cepsilon1, STATE Cepsilon2, TPZVec<STATE> &BodyForce );
    
    void GetParameters(STATE &MU, STATE &RHO, STATE &Cmu, STATE &SigmaK, STATE &SigmaEps, STATE &Cepsilon1, STATE &Cepsilon2, TPZVec<STATE> &BodyForce );    
	
    virtual ~TPZIncNavierStokesKEps();
	
    virtual int Dimension() const;
	
    /** @brief Returns the number of state variables associated with the material*/
    virtual int NStateVariables();
	
    /** @brief Print out the data associated with the material*/
    virtual void Print(std::ostream &out = std::cout);
	
    /** @brief Returns the number of variables associated with the variable
     *  indexed by var. \n var is obtained by calling VariableIndex*/
    virtual int NSolutionVariables(int var);
	
protected:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol,
						  TPZFMatrix<REAL> &axes, int var, TPZVec<STATE> &Solout);
public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	
    /** @brief Computes contribution to the tangent matrix and residual at an integration point */
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef);
	
    /** @brief Computes contribution to the residual at an integration point */
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ef);
	
	
	/** @brief Computes contribution to the stiffness matrix and right hand side at the integration point of a boundary */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc);
	
	/** @brief Computes contribution to the stiffness matrix and right hand side at the integration point of a boundary */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc)
	{
        TPZMaterial::ContributeBC(data,weight,ef,bc);
    }
	
	
    /** 
	 * @brief Computes the error due to the difference between the interpolated flux and the flux computed \n
	 * based on the derivative of the solution
     */
    virtual void Errors(TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,
                        TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
                        TPZVec<STATE> &uexact, TPZFMatrix<STATE> &duexact,
                        TPZVec<REAL> &val){
        PZError << __PRETTY_FUNCTION__ << std::endl;
        PZError << "Method not implemented! Error comparison not available. Please, implement it." << std::endl;
    }
    public:
virtual int ClassId() const;

};


#endif


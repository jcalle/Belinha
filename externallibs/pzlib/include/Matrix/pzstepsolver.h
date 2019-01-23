/**
 * @file
 * @brief Contains TPZStepSolver class which defines step solvers class.
 */

#ifndef TPZSTEPSOLVER_H
#define TPZSTEPSOLVER_H
#include "pzsolve.h"

#include "pzfmatrix.h"

#include "TPZStream.h"

#include <list>

/**
 * @brief Defines step solvers class. \ref solver "Solver"
 * @ingroup solver
 */
template<class TVar>
class TPZStepSolver: public TPZMatrixSolver<TVar>
{
public:
	TPZStepSolver(TPZAutoPointer<TPZMatrix<TVar> > refmat = 0);
	
	TPZStepSolver(const TPZStepSolver<TVar> & copy);
	
	virtual ~TPZStepSolver();
	
	void SetSOR(const int64_t numiterations, const REAL overrelax, const REAL tol, const int64_t FromCurrent);
	
	void SetSSOR(const int64_t numiterations, const REAL overrelax, const REAL tol, const int64_t FromCurrent);
	
	void SetJacobi(const int64_t numiterations, const REAL tol, const int64_t FromCurrent);
	
	void SetCG(const int64_t numiterations, const TPZMatrixSolver<TVar> &pre, const REAL tol, const int64_t FromCurrent);
	
	void SetGMRES(const int64_t numiterations, const int numvectors, const TPZMatrixSolver<TVar> &pre, const REAL tol, const int64_t FromCurrent);
	
	void SetBiCGStab(const int64_t numiterations, const TPZMatrixSolver<TVar> &pre, const REAL tol, const int64_t FromCurrent);
	
	void SetDirect(const DecomposeType decomp);
	
	void SetMultiply();
	
	virtual TPZSolver<TVar> *Clone() const
	{
		return new TPZStepSolver<TVar>(*this);
	}
	
	void SetTolerance(REAL tol)
	{
		fTol = tol;
	}
    
    /** @brief return the value of tolerance from the solver */
    REAL GetTolerance() const
    {
        return fTol;
    }
	
    /** @brief reset the data structure of the solver object */
	void ResetSolver();
    
    virtual typename TPZMatrixSolver<TVar>::MSolver Solver()
    {
        return fSolver;
    }
	
	/** @brief returns the equations for which the equations had zero pivot */
	std::list<int64_t> &Singular()
	{
		return fSingular;
	}
	
	/** @brief This method will reset the matrix associated with the solver */
	/** This is useful when the matrix needs to be recomputed in a non linear problem */
	virtual void ResetMatrix();

	/** @brief Updates the values of the current matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > matrix)
	{
		if (fPrecond)
			fPrecond->UpdateFrom(matrix);
		TPZMatrixSolver<TVar>::UpdateFrom(matrix);
	}
	
    
	void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual = 0);
    
    /** @brief Decompose the system of equations if a direct solver is used */
    virtual void Decompose();
    
    /** @brief Define the preconditioner as a solver object */
	void SetPreconditioner(TPZSolver<TVar> &solve);
    
    /** @brief Number of iterations of last solve */
    int NumIterations()
    {
        return fNumIterations;
    }
    
    /** @brief access method to the preconditioner */
    TPZSolver<TVar> *PreConditioner()
    {
        return fPrecond;
    }
	
	/** @brief Serialization methods */
	public:
virtual int ClassId() const;

	virtual void Write(TPZStream &buf, int withclassid) const;
	virtual void Read(TPZStream &buf, void *context);
	
	
private:
	typename TPZMatrixSolver<TVar>::MSolver fSolver;
	DecomposeType fDecompose;
    
    /// Maximum number of iterations
    int64_t fMaxIterations;
    
    /// Number of iterations of last solve
	int64_t fNumIterations;
	int fNumVectors;
	REAL fTol;
	REAL fOverRelax;
	
	/** @brief Solver using preconditioner matrix */
	TPZSolver<TVar> *fPrecond;
	int64_t fFromCurrent;
	
	std::list<int64_t> fSingular;
};

template<class TVar>
int TPZStepSolver<TVar>::ClassId() const{
    return Hash("TPZStepSolver") ^ TPZMatrixSolver<TVar>::ClassId() << 1;
}
#endif

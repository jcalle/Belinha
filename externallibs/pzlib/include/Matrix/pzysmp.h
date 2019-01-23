/**
 * @file
 * @brief Contains the TPZFYsmpMatrix class which implements a non symmetric sparse matrix.
 */

#ifndef YSMPMATH
#define YSMPMATH

#include "pz_config.h"

#ifdef USING_BLAS
#ifdef USING_MKL
#include <mkl.h>
#elif MACOSX
#include <Accelerate/Accelerate.h>
#else
#ifdef MACOSX
#include <Accelerate/Accelerate.h>
#else
extern "C"{
     #include "cblas.h"
     };
#endif
#endif
#endif

template<class TVar>
class TPZVerySparseMatrix;

#include "pzmatrix.h"
#include "pzfmatrix.h"

#ifdef USING_MKL
#include "TPZPardisoControl.h"
#endif

/**
 * @brief Implements a non symmetric sparse matrix (Yale Sparse Matrix Storage). \ref matrix "Matrix"
 * @ingroup matrix
 */
/**
 * Defines operations on general sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 */
template<class TVar>
class TPZFYsmpMatrix : public TPZMatrix<TVar> {
	
#ifdef USING_MKL
    friend class TPZPardisoControl<TVar>;
#endif
    
	public :
	
	/** @brief An auxiliary structure to hold the data of the subset \n of equations used to multiply in a multi-threaded environment */
	/**
	 In future versions this structure should be defined in a derived class
	 */
	struct TPZMThread {
		const TPZFYsmpMatrix<TVar> *target;
		int64_t fFirsteq;
		int64_t fLasteq;
		const TPZFMatrix<TVar> *fX;
		TPZFMatrix<TVar> *fZ;
		TVar fAlpha;
		int fOpt;
	};
	
private:
	
	static void * ExecuteMT(void *entrydata);
	
public:
    
    TPZFYsmpMatrix();
    
    TPZFYsmpMatrix(const int64_t rows,const int64_t cols );
	
	TPZFYsmpMatrix(const TPZVerySparseMatrix<TVar> &cp);
	
	TPZFYsmpMatrix &operator=(const TPZFYsmpMatrix<TVar> &copy);
	
	TPZFYsmpMatrix &operator=(const TPZVerySparseMatrix<TVar> &cp);
    
	CLONEDEF(TPZFYsmpMatrix)
	
	virtual ~TPZFYsmpMatrix();	
	
    /** @brief Fill matrix storage with randomic values */
    /** This method use GetVal and PutVal which are implemented by each type matrices */
    void AutoFill(int64_t nrow, int64_t ncol, int symmetric);
    

    
	/** @brief Get the matrix entry at (row,col) without bound checking */
	virtual const TVar &GetVal(const int64_t row,const int64_t col ) const;
	
	int64_t NumTerms()
	{
		return fIA[this->Rows()];
	}
	
	int PutVal(const int64_t row, const int64_t col, const TVar &Value);
	
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha=1.,const TVar beta = 0., const int opt = 0) const;
	
	virtual void MultAddMT(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						   const TVar alpha=1.,const TVar beta = 0., const int opt = 0);
	
	virtual int GetSub(const int64_t sRow,const int64_t sCol,const int64_t rowSize,
					   const int64_t colSize, TPZFMatrix<TVar> & A ) const;
	
	void GetSub(const TPZVec<int64_t> &indices,TPZFMatrix<TVar> &block) const;
	
	/** @brief Pass the data to the class. */
	virtual void SetData( int64_t *IA, int64_t *JA, TVar *A );
    
    /** @brief Pass the data to the class. */
    virtual void SetData( TPZVec<int64_t> &IA, TPZVec<int64_t> &JA, TPZVec<TVar> &A );
	
	/** @brief Print the matrix along with a identification title */
	virtual void Print(const char *title, std::ostream &out = std::cout , const MatrixOutputFormat form = EFormatted) const;
	
	/**
	 * @name Solvers
	 * @brief Linear system solvers. \n
	 */
	 /** For symmetric decompositions lower triangular matrix is used. \n
	 * Solves a system A*X = B returning X in B
	 */  
	//@{
	/**
	 * @brief Solves the linear system using Jacobi method. \n
	 * @param numiterations The number of interations for the process.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param scratch Available manipulation area on memory.
	 * @param tol The tolerance value.
	 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
	 */
	virtual void SolveJacobi(int64_t & numiterations, const TPZFMatrix<TVar> & F, TPZFMatrix<TVar> & result,
							 TPZFMatrix<TVar> * residual, TPZFMatrix<TVar> & scratch, REAL & tol, const int FromCurrent = 0) ;
	
	void SolveSOR(int64_t &numiterations, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &x,
				  TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &scratch,
				  const REAL overrelax, REAL &tol,
				  const int FromCurrent = 0,const int direction = 1 ) ;    
	// @}
	
	/**
	 * @brief Add a contribution of a stiffness matrix
	 * putting it on destination indexes position
	 */
	virtual void AddKelOld(
						   TPZFMatrix<TVar> & elmat //! Member stiffness matrix beeing added
						   , TPZVec < int > & destinationindex //! Positioning of such members on global stiffness matrix
						   );    
	
	virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & destinationindex);
	
	virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex);
	
	void MultiplyDummy(TPZFYsmpMatrix<TVar> & B, TPZFYsmpMatrix<TVar> & Res);
	
	virtual int Zero();
	
	/**
	 * @name Factorization
	 * @brief Those member functions are related to matrices factorization
	 */
	//@{
	/**
	 * @brief Decomposes the current matrix using LU decomposition.
	 */
	virtual int Decompose_LU(std::list<int64_t> &singular);
	virtual int Decompose_LU();
	
	//@}
	
	/**
	 * @name Substitutions
	 * @brief Substitutions forward and backward
	 */
	//@{  
	/**
	 * @brief Computes Forward and Backward substitution for a "LU" decomposed matrix.
	 * @param B right hand side and result after all
	 */
	virtual int Substitution( TPZFMatrix<TVar> * B ) const;
	
	//@}
	
    public:
virtual int ClassId() const;

private:
	
	void ComputeDiagonal();
	
	/*
	 * @brief Perform row update of the sparse matrix
	 */
	void RowLUUpdate(int64_t sourcerow, int64_t destrow);
	
protected:
	TPZVec<int64_t>  fIA;
	TPZVec<int64_t>  fJA;
	TPZVec<TVar> fA;
	
	TPZVec<TVar> fDiag;
	
	int   fSymmetric;
	
#ifdef USING_MKL    
    TPZPardisoControl<TVar> fPardisoControl;
#endif
protected:
	
	/**
	 * @brief Implements a initialization method for the sparse structure. It sets the initial value for the fIA and fJA.
	 */ 
	/**
	 * -fIA will contain the initial positions for all the equations
	 * -fJA will contain (-1) on all its positions
	 * -fA will contain 0 on all its value 
	 */
	void InitializeData();
};


template<class TVar>
inline void TPZFYsmpMatrix<TVar>::SetData( int64_t *IA, int64_t *JA, TVar *A ) {
	// Pass the data to the class.
    int nel = this->Rows()+1;
    fIA.resize(nel);
    memccpy(&fIA[0], IA, nel, sizeof(int64_t));
    int64_t nval = fIA[nel-1];
    fJA.resize(nval);
    memccpy(&fJA[0], JA, nval, sizeof(int64_t));
    fA.resize(nval);
    memccpy(&fA[0], A, nval, sizeof(TVar));
	ComputeDiagonal();
}

/** @brief Pass the data to the class. */
template<class TVar>
inline void TPZFYsmpMatrix<TVar>::SetData( TPZVec<int64_t> &IA, TPZVec<int64_t> &JA, TPZVec<TVar> &A ){
    
    if (IA.size() != this->Rows() + 1 ) {
        DebugStop();
    }
    
    if (JA.size() != IA[this->Rows()]) {
        DebugStop();
    }
    
    fIA = IA;
    fJA = JA;
    fA = A;
}

#endif

/**
 * @file
 * @brief Contains TPZSBMatrix class which implements symmetric band matrices(hermitian, for the complex case. assumed to be
 * upper triangular).
 */

#ifndef TSBNDMATH
#define TSBNDMATH

#include "pzmatrix.h"
#include "pzfmatrix.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

/**
 * @brief Implements symmetric band matrices. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSBMatrix : public TPZMatrix<TVar>
{
public:
    TPZSBMatrix() : TPZRegisterClassId(&TPZSBMatrix::ClassId),
    TPZMatrix<TVar>() , fDiag() { fBand = 0; }
    TPZSBMatrix(const int64_t dim,const int64_t band );
    TPZSBMatrix(const TPZSBMatrix<TVar> &A ) : TPZRegisterClassId(&TPZSBMatrix::ClassId),
    TPZMatrix<TVar>(A)  { Copy(A); }
    
    CLONEDEF(TPZSBMatrix)
    
    ~TPZSBMatrix() { Clear(); }
    
    int    PutVal(const int64_t row,const int64_t col,const TVar& element );
    const TVar &GetVal(const int64_t row,const int64_t col ) const;
    
    TVar &operator()(int64_t row, int64_t col);
    
    /** @brief Checks if the current matrix is symmetric */
    virtual int IsSimetric() const
    {
        return 1;
    }

    friend class TPZSBMatrix<float>;
    friend class TPZSBMatrix<double>;
    
    /// copy the values from a matrix with a different precision
    template<class TVar2>
    void CopyFrom(TPZSBMatrix<TVar2> &orig)
    {
        TPZMatrix<TVar>::CopyFrom(orig);
        fDiag.resize(orig.fDiag.size());
        int64_t nel = fDiag.size();
        for (int64_t el=0; el<nel; el++) {
            fDiag[el] = orig.fDiag[el];
        }
    }
    
    /** @brief Computes z = beta * y + alpha * opt(this)*x */
    /** z and x cannot overlap in memory */
    void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const;
    
    void Print(const char *name = NULL, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const;
    //friend std::ostream & operator<< <>(std::ostream& out,const TPZSBMatrix<TVar>  &A); Leonardo removendo o '<>' antes do (std...
    template<class TT>friend std::ostream & operator<< (std::ostream& out,const TPZSBMatrix<TT>  &A);
    
    /** Fill the matrix with random values (non singular matrix) */
    void AutoFill(int64_t nrow, int64_t ncol, int symmetric);
    
    
    /// Operadores com matrizes SKY LINE.
    // @{
    TPZSBMatrix &operator= (const TPZSBMatrix<TVar> &A );
    TPZSBMatrix operator+  (const TPZSBMatrix<TVar> &A ) const;
    TPZSBMatrix operator-  (const TPZSBMatrix<TVar> &A ) const;
    TPZSBMatrix &operator+=(const TPZSBMatrix<TVar> &A );
    TPZSBMatrix &operator-=(const TPZSBMatrix<TVar> &A );
    // @}
    TPZSBMatrix<TVar> operator*  (const TVar v ) const;
    TPZSBMatrix<TVar> &operator*=(const TVar v );
    
    TPZSBMatrix<TVar> operator-() const { return operator*(-1.0); }
    
    /// Redimension the matrix keeping original elements.
    int Resize(const int64_t newDim ,const int64_t);
    
    /// Redimension the matrix and zeroes its elements
    int Redim(const int64_t newDim) {return Redim(newDim,newDim);}
    int Redim(const int64_t newRows ,const int64_t newCols);
    
    /// Zeroes the elements of the matrix
    int Zero();
    
    int64_t GetBand() const { return fBand; }
    int   SetBand(const int64_t newBand );
    
    /// To solve linear systems
    // @{
#ifdef USING_LAPACK
    int Decompose_Cholesky();  // Faz A = GGt.
    int Decompose_Cholesky(std::list<int64_t> &singular);
#endif
    
    int Subst_Forward( TPZFMatrix<TVar>*B ) const;
    int Subst_Backward ( TPZFMatrix<TVar> *b ) const;

    int Decompose_LDLt(std::list<int64_t> &singular);
    int Decompose_LDLt();
    int Subst_LForward( TPZFMatrix<TVar> *B ) const;
    int Subst_LBackward( TPZFMatrix<TVar> *B ) const;
    int Subst_Diag( TPZFMatrix<TVar> *B ) const;
//    int Subst_Forward( TPZFMatrix<TVar>*B ) const;
//    int Subst_Backward( TPZFMatrix<TVar> *B ) const;
    
    // @}
    
#ifdef USING_LAPACK
    /*** @name Solve eigenvalues ***/
    /** @{ */
    
    /// Computes the eigenvalues and eigenvectors of the symmetric matrix
    // on exit the matrix contains the eigenvectors
    /** @brief Solves the Ax=w*x eigenvalue problem and calculates the eigenvectors
     * @param w Stores the eigenvalues
     * @param Stores the correspondent eigenvectors
     */
    int SolveEigenProblem(TPZVec < std::complex<double> > &w, TPZFMatrix < std::complex<double> > &eigenVectors);
    /** @brief Solves the Ax=w*x eigenvalue problem and does NOT calculates the eigenvectors
     * @param w Stores the eigenvalues
     */
    int SolveEigenProblem(TPZVec < std::complex<double> > &w);
    /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
     * @param w Stores the eigenvalues
     * @param Stores the correspondent eigenvectors
     */
    int SolveGeneralisedEigenProblem(TPZSBMatrix< TVar > &B , TPZVec < std::complex<double> > &w, TPZFMatrix < std::complex<double> > &eigenVectors);
    /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and does NOT calculates the eigenvectors
     * @param w Stores the eigenvalues
     */
    int SolveGeneralisedEigenProblem(TPZSBMatrix< TVar > &B , TPZVec < std::complex<double> > &w);
    
    /** @} */
#endif
    public:
virtual int ClassId() const;

private:
    
    int64_t  Size() const
    {
        return( this->Dim() * (fBand + 1) );
    }
//    int  PutZero();
    //static int  Error(const char *msg1,const char* msg2="" ) ;
    int  Clear();
    void Copy (const TPZSBMatrix<TVar> & );
    
    int64_t Index(int64_t i, int64_t j) const
    {
#ifdef PZDEBUG
        if (i>j) {
            DebugStop();
        }
#endif
        return fBand+i-j+(fBand+1)*j;
    }
    TPZVec<TVar> fDiag;
    int64_t  fBand;
};

#endif

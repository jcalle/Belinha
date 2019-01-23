/**
 * @file
 * @brief Contains TPZMatrixclass which implements full matrix (using column major representation).
 */


#ifndef _TFULLMATRIXH_
#define _TFULLMATRIXH_

#include <math.h>            // for sqrt
#include <string.h>          // for NULL, memset
#include <complex>           // for complex
#include <iostream>          // for operator<<, ostream, cout
#include <list>              // for list
#include <sstream>           // for basic_stringbuf<>::int_type, basic_strin...
#include "pzmatrix.h"		 // for TPZMatrix
#include "pzerror.h"         // for DebugStop
#include "pzmanvector.h"     // for TPZManVector
#include "pzreal.h"          // for TPZFlopCounter, IsZero, REAL, sqrt, fabs
#include "pzvec.h"           // for TPZVec
#include "tpzautopointer.h"  // for TPZAutoPointer

#ifdef _AUTODIFF
#include "tfad.h"
#include "fad.h"
#include "pzextractval.h"
#endif

class TPZSavable;
class TPZStream;
template <class TVar> class TPZFMatrix;
template <class TVar> class TPZVerySparseMatrix;


template <class T>
class TPZVec;

/**
 * @addtogroup matrix
 * @{
 */

/** @brief MACRO to get MAT(row,col) entry */
#define GETVAL(MAT,rows,row,col) MAT->fElem[((unsigned)col)*rows+row]
/** @brief MACRO to put value val into MAT(row,col) entry */
#define PUTVAL(MAT,rows,row,col,val) MAT->fElem[((unsigned)col)*rows+row]=val
/** @brief MACRO to get the entry of the vector (ptr[col*rows+row]) as matrix ( ptr(row,col) )*/
#define SELECTEL(ptr,rows,row,col) ptr[col*rows+row]


/** @brief Returns a dot product to matrices */
template<class TVar>
TVar Dot(const TPZFMatrix<TVar> &A,const TPZFMatrix<TVar> &B);

/** @brief Returns the norm of the matrix A */
template<class TVar>
TVar Norm(const TPZFMatrix<TVar> &A);


/**
 * @brief Full matrix class. \ref matrix "Matrix"
 * @note The full matrix class is special in that the data is stored column wise
 * @author Misael Mandujano
 * @since 04/1996
 */

template<class TVar=REAL>
class TPZFMatrix: public TPZMatrix<TVar> {
    
public:
    
    /** @brief Simple constructor */
    TPZFMatrix() : TPZRegisterClassId(&TPZFMatrix<TVar>::ClassId),
    TPZMatrix<TVar>( 0, 0 ), fElem(0),fGiven(0),fSize(0) {}
    /**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
     @param buf Preallocated memory area which can be used by the matrix object
     @param size Size of the area pointed to by buf
     */
    TPZFMatrix(const int64_t rows ,const int64_t columns, TVar* buf,const int64_t size);
    /**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
     @param val Inital value fill all elements
     */
    TPZFMatrix(const int64_t rows ,const int64_t columns,const TVar & val );
    /**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
     */
    inline  TPZFMatrix(const int64_t rows ,const int64_t columns = 1) : 
    TPZRegisterClassId(&TPZFMatrix<TVar>::ClassId),
    TPZMatrix<TVar>(rows,columns), fElem(0),fGiven(0),fSize(0) {
        if(rows*columns) fElem = new TVar[rows*columns];
    }
    
    /**
     * @brief Copy constructor specialized form TPZVerySparseMatrix
     * @param refmat Used as a model for current object
     */
    TPZFMatrix(TPZVerySparseMatrix<TVar> const & A);
    
    /**
     * @brief Copy constructor
     * @param refmat Used as a model for current object
     */
    TPZFMatrix(const TPZFMatrix<TVar> & refmat);
    
    
    
    CLONEDEF(TPZFMatrix<TVar>)
    TPZFMatrix(const TPZMatrix<TVar> & refmat);
    
    /** @brief Simple destructor */
    virtual  ~TPZFMatrix();
    
    int64_t MemoryFootprint() const
    {
        return (sizeof(TVar)*this->Rows()*this->Cols());
    }
    
    TVar *Adress()
    {
        return fElem;
    }
    
    friend class TPZFMatrix<float>;
    friend class TPZFMatrix<double>;

    /// copy the values from a matrix with a different precision
    template<class TVar2>
    void CopyFrom(TPZFMatrix<TVar2> &orig)
    {
        Resize(orig.Rows(), orig.Cols());
        TPZMatrix<TVar>::CopyFrom(orig);
        int64_t nel = orig.Rows()*orig.Cols();
        for (int64_t el=0; el<nel; el++) {
            fElem[el] = orig.fElem[el];
        }
    }

    /** @brief Updates the values of the matrix based on the values of the matrix */
    virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> >  mat)
    {
        TPZMatrix<TVar> *matptr = mat.operator->();
        TPZFMatrix<TVar> *from = dynamic_cast<TPZFMatrix<TVar> *>(matptr);
        if (from) {
            *this = *from;
        }
        else
        {
            std::cout << "TPZMatrix<TVar>::UdateFrom is not implemented\n";
            DebugStop();
        }
    }

    int PutVal(const int64_t row,const int64_t col,const TVar & value );
    const TVar &GetVal(const int64_t row,const int64_t col ) const;
    
    virtual TVar &s(const int64_t row, const int64_t col);
    
    TVar &g(const int64_t row, const int64_t col) const;
    /**
     * @brief Performs a right hand side assemblage
     * @param rhs Load vector
     * @param destination Destine index on current matrix
     */
    void AddFel(TPZFMatrix<TVar> &rhs,TPZVec<int64_t> &destination);
    /**
     * @brief Performs a right hand side assemblage
     * @param rhs Load vector
     * @param source Source index on rhs
     * @param destination Destine index on current matrix
     */
    void AddFel(TPZFMatrix<TVar> &rhs,TPZVec<int64_t> &source, TPZVec<int64_t> &destination);
    
    
    /**
     * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
     * @param x Is x on the above operation
     * @param y Is y on the above operation
     * @param z Is z on the above operation
     * @param alpha Is alpha on the above operation
     * @param beta Is beta on the above operation
     * @param opt Indicates if is Transpose or not
     */
    virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                         const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const ;
    
    
    static void MultAdd(const TVar *ptr, int64_t rows, int64_t cols, const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                        const TVar alpha=1.,const TVar beta = 0.,const int opt = 0);
    
    /**
     * @name Generic operator with TVar type
     * @{
     */
    TVar &operator()(const int64_t row,const int64_t col);
    TVar &operator()(const int64_t row);
    /** @} */
    
    /**
     * @name Operations with FULL matrices
     * @{
     */
    
    /** @brief Generic operator with FULL matrices */
    virtual TPZFMatrix&operator= (const TPZFMatrix<TVar> &A );
    TPZFMatrix<TVar> operator+  (const TPZFMatrix<TVar> &A ) const;
    TPZFMatrix<TVar> operator-  (const TPZFMatrix<TVar> &A ) const;
    TPZFMatrix<TVar> operator*  ( TPZFMatrix<TVar> A ) const ;
    TPZFMatrix<TVar> &operator+=(const TPZFMatrix<TVar> &A );
    TPZFMatrix<TVar> &operator-=(const TPZFMatrix<TVar> &A );
    
    /** @} */
    
    /**
     * @brief Performs an ZAXPY operation being *this += alpha * p
     * @param alpha Being alpha on above opereation
     * @param p Being p on above operation
     */
    void ZAXPY(const TVar alpha, const TPZFMatrix<TVar> &p);
    /**
     * @brief Performs an operation *this = this * beta + z
     * @param beta Being beta on above opereation
     * @param z Being z on above operation
     */
    void TimesBetaPlusZ(const TVar beta, const TPZFMatrix<TVar> &z);
    
    /** @brief Generic operator with matrices */
    TPZFMatrix<TVar> &operator= (const TPZMatrix<TVar> &A );
    
    /**
     * @name Numerics
     * @brief Numeric operations with matrices
     * @{
     */
    
    /** @brief Numeric operator with matrices */
    TPZFMatrix<TVar> &operator= (const TVar val );
    TPZFMatrix<TVar> operator+  (const TVar val ) const;
    TPZFMatrix<TVar> operator-  (const TVar val ) const;
    TPZFMatrix<TVar> operator*  (const TVar val ) const;
    TPZFMatrix<TVar> &operator+=(const TVar val );
    TPZFMatrix<TVar> &operator-=(const TVar val )  { return operator+=( -val ); }
    TPZFMatrix<TVar> &operator*=(const TVar val );
    
    //	TPZFMatrix<TVar> operator-() const;// { return operator*( -1.0 ); }
    
    /** @} */
    
    /** @brief Redimension a matrix, but maintain your elements. */
    int Resize(const int64_t newRows,const int64_t wCols );
    
    /** @brief Redimension the matrix doing nothing with the elements */
    int SetSize(int64_t newRows, int64_t newCols);
    
    /** @brief Remodel the shape of the  matrix, but keeping the same dimension. */
    int Remodel(const int64_t newRows,const int64_t wCols );
    
    /** @brief Redimension a matrix and ZERO your elements. */
    int Redim(const int64_t newRows,const int64_t newCols );
    
    /** @brief Makes Zero all the elements */
    int Zero();
    
#ifdef USING_LAPACK
    /** @brief Initialize pivot with i = i  */
    void InitializePivot();
#endif

    
    /**
     * @brief This method implements a Gram Schimidt method. \n this = Orthog.TransfToOrthog
     * @param Orthog [out] each column represents a vector orthogonalized with respect to the first vector (first column of *this). Vectors are normalized
     * @param TransfToOrthog [out] is the basis change from *this to Orthog
     * @author Caju
     * @since 2007
     */
    void GramSchmidt(TPZFMatrix<TVar> &Orthog, TPZFMatrix<TVar> &TransfToOrthog);
    
    void DeterminantInverse(TVar &determinant, TPZFMatrix<TVar> &inverse);
    
    void Transpose(TPZMatrix<TVar> *const T) const;
    
    /** @see TPZMatrix<TVar>::Transpose */
    
    void Transpose();
    
    /*** @name Solve linear system of equations ***/
    /** @{ */
    
    /** @brief Cholesky Decomposition Optmized. for walks in the direction of the vector that composes the matrix */
    virtual int Decompose_Cholesky();
    virtual int Decompose_Cholesky(std::list<int64_t> &singular);
    
    /** @brief LU Decomposition. Stores L and U matrices at the storage of the same matrix */
    virtual int Decompose_LU(std::list<int64_t> &singular);
    virtual int Decompose_LU();
    
    /**
     * @brief Decomposes the current matrix using LDLt. \n
     * The current matrix has to be symmetric.
     * "L" is lower triangular with 1.0 in its diagonal and "D" is a Diagonal matrix.
     */
    virtual int Decompose_LDLt();
    
    static int Substitution(const TVar *ptr, int64_t rows, TPZFMatrix<TVar> *B);
    
    virtual int Substitution( TPZFMatrix<TVar> *B ) const;
    
    /** @brief LU Decomposition using pivot */
    virtual int Decompose_LU(TPZVec<int> &index);
    
    /** @brief LU substitution using pivot. */
    virtual int Substitution( TPZFMatrix<TVar> *B, const TPZVec<int> &index ) const;
    
    /** @brief LU substitution using pivot. Static version. */
    static int Substitution(const TVar *ptr, int64_t rows,  TPZFMatrix<TVar> *B, const TPZVec<int> &index );
    
#ifdef USING_LAPACK
    /**
     * @brief Computes B = Y, where A*Y = B, A is lower triangular.
     * @param b right hand side and result after all
     */
    virtual int Subst_Forward( TPZFMatrix<TVar>* b ) const;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is upper triangular.
     * @param b right hand side and result after all
     */
    virtual int Subst_Backward( TPZFMatrix<TVar>* b ) const;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
     * @param b right hand side and result after all
     */
    virtual int Subst_LForward( TPZFMatrix<TVar>* b ) const;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
     * @param b right hand side and result after all
     */
    virtual int Subst_LBackward( TPZFMatrix<TVar>* b ) const;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
     * @param b right hand side and result after all
     */
    virtual int Subst_Diag( TPZFMatrix<TVar>* b ) const;
#endif
    
    /** @} */
    
#ifdef USING_LAPACK
    /*** @name Solve eigenvalues ***/
    /** @{ */
    /** @brief Solves the Ax=w*x eigenvalue problem and calculates the eigenvectors
     * @param w Stores the eigenvalues
     * @param Stores the correspondent eigenvectors
     */
    virtual int SolveEigenProblem(TPZVec < std::complex<double> > &w, TPZFMatrix < std::complex<double> > &eigenVectors);
    /** @brief Solves the Ax=w*x eigenvalue problem and does NOT calculates the eigenvectors
     * @param w Stores the eigenvalues
     */
    virtual int SolveEigenProblem(TPZVec < std::complex<double> > &w);

    /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
     * @param w Stores the eigenvalues
     * @param Stores the correspondent eigenvectors
     */
    virtual int SolveGeneralisedEigenProblem(TPZFMatrix< TVar > &B , TPZVec < std::complex<double> > &w, TPZFMatrix < std::complex<double> > &eigenVectors);
    /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and does NOT calculates the eigenvectors
     * @param w Stores the eigenvalues
     */
    virtual int SolveGeneralisedEigenProblem(TPZFMatrix< TVar > &B , TPZVec < std::complex<double> > &w);    
    /** @} */
#endif
    
    /** @brief Routines to send and receive messages */
    public:
virtual int ClassId() const;

    
    virtual void Read( TPZStream &buf, void *context );
    virtual void Write(TPZStream &buf, int withclassid ) const;
    
    /** @brief Compare the object for identity with the object pointed to, eventually copy the object */
    /**
     * compare both objects bitwise for identity. Put an entry in the log file if different
     * overwrite the calling object if the override flag is true
     */
    virtual bool Compare(TPZSavable *copy, bool override = false);
    
    /** @brief Compare the object for identity with the object pointed to, eventually copy the object */
    /**
     * compare both objects bitwise for identity. Put an entry in the log file if different
     * generate an interupt if the override flag is true and the objects are different
     */
    virtual bool Compare(TPZSavable *copy, bool override = false) const;
    
    operator const TVar*() const { return fElem; }
    
    static void PrintStatic(const TVar *ptr, int64_t rows, int64_t cols, const char *name, std::ostream& out,const MatrixOutputFormat form);
    
private:
    
    static int Error(const char *msg1,const char *msg2=0 );
    int Clear();
    
    TVar *fElem;
    TVar *fGiven;
    int64_t fSize;
    
#ifdef USING_LAPACK
    TPZManVector<int,5> fPivot;
    
    TPZVec<TVar> fWork;
#endif
};

/** @} */

template<class TVar>
inline TPZFMatrix<TVar>::TPZFMatrix(const int64_t rows,const int64_t cols,TVar * buf,const int64_t sz)
: TPZRegisterClassId(&TPZFMatrix<TVar>::ClassId),TPZMatrix<TVar>( rows, cols ), 
fElem(buf),fGiven(buf),fSize(sz) {
    int64_t size = rows * cols;
    if(size == 0)
    {
        fElem = NULL;
    }
    else if(size > sz)
    {
        fElem=new TVar[size];
#ifndef NODEBUG
        if ( fElem == NULL && size) Error( "Constructor <memory allocation error>." );
#endif
    }
}
template<class TVar>
inline TPZFMatrix<TVar>::TPZFMatrix(const int64_t rows,const int64_t cols,const TVar & val )
: TPZRegisterClassId(&TPZFMatrix<TVar>::ClassId), TPZMatrix<TVar>( rows, cols ), 
 fElem(0), fGiven(0), fSize(0) {
    int64_t size = rows * cols;
    if(!size) return;
    fElem=new TVar[size];
#ifdef PZDEBUG
    if ( fElem == NULL && size) Error( "Constructor <memory allocation error>." );
#endif
    for(int64_t i=0;i<size;i++) fElem[i] = val;
}

template<class TVar>
/** @brief Implements a scalar product val*A */
inline TPZFMatrix<TVar> operator*(TVar val, const TPZFMatrix<TVar> &A)
{
    return A*val;
}


/*******************************/
/*** Operator*( TPZMatrix<TVar> & ) ***/
template<class TVar>
inline TPZFMatrix<TVar> TPZFMatrix<TVar>::operator*( TPZFMatrix<TVar> A ) const {
    if ( this->Cols() != A.Rows() )
        Error( "Operator* <matrixs with incompatible dimensions>" );
    
    TPZFMatrix<TVar> res;
    res.Redim( this->Rows(), A.Cols() );
    MultAdd(A,A,res,1.,0.,0);
    return( res );
}

/**************/
/*** PutVal ***/
template<class TVar>
inline int TPZFMatrix<TVar>::PutVal(const int64_t row, const int64_t col,const TVar & value ) {
    fElem[ ((unsigned)col) * this->Rows() + row ] = value;
    return( 1 );
}

/******************/
/*** Destructor ***/
template<class TVar>
inline TPZFMatrix<TVar>::~TPZFMatrix() {
    if(fElem && fElem != fGiven) delete[]( fElem );
    fElem = 0;
    fSize = 0;
}


/**************/
/*** GetVal ***/
template<class TVar>
inline const TVar &TPZFMatrix<TVar>::GetVal( const int64_t row, const int64_t col ) const {
#ifdef PZDEBUG
    if(row >=  this->Rows() || row<0 || col >=  this->Cols() || col<0) {
        Error("TPZFMatrix::operator() "," Index out of bounds");
        DebugStop();
    }
#endif
    return( fElem[ col*this->fRow + row ] );
}

template<class TVar>
inline TVar &TPZFMatrix<TVar>::operator()( const int64_t row, const int64_t col) {
#ifndef NODEBUG
    if(row >=  this->Rows() || row<0 || col >=  this->Cols() || col<0) {
        Error("TPZFMatrix<TVar>::operator() "," Index out of bounds");
        DebugStop();
    }
#endif
    return *(this->fElem+col*this->fRow+row);
}

template<class TVar>
inline TVar &TPZFMatrix<TVar>::s(const int64_t row, const int64_t col) {
    // verificando se o elemento a inserir esta dentro da matriz
    return operator()(row,col);
}

template<class TVar>
inline TVar &TPZFMatrix<TVar>::g( const int64_t row, const int64_t col) const {
#ifdef PZDEBUG
    if(row >=  this->Rows() || row<0 || col >=  this->Cols() || col<0) {
        Error("TPZFMatrix<TVar>::operator() "," Index out of bounds");
        DebugStop();
    }
#endif
    return *(this->fElem+col*this->fRow+row);
}

template<class TVar>
inline TVar &TPZFMatrix<TVar>::operator()(const int64_t row) {
#ifdef PZDEBUG
    if(row >=  this->Rows() || row<0) {
        Error("TPZFMatrix<TVar>::operator() "," Index out of bounds");
        DebugStop();
    }
#endif
    return *(this->fElem+row);
}

template<class TVar>
inline int TPZFMatrix<TVar>::Redim(const int64_t newRows,const int64_t newCols) {
    int64_t newsize = newRows*newCols;
    int64_t size = this->fRow*this->fCol;
    if ( newsize == size) {
        this->fRow = newRows;
        this->fCol = newCols;
        Zero();
        return( 1 );
    }
    if(this->fElem && this->fElem != this->fGiven) delete []this->fElem;
    
    if(this->fGiven && newsize <= this->fSize) {
        this->fElem = this->fGiven;
    } else if(newsize == 0) {
        this->fElem = NULL;
    } else {
        this->fElem = new TVar[ newsize ] ;
    }
#ifndef NODEBUG
    if (newsize && this->fElem == NULL )
        Error( "Resize <memory allocation error>." );
#endif
    
    this->fRow  = newRows;
    this->fCol  = newCols;
    
    Zero();
    
    return( 1 );
}

/***************/
/****Zero*******/

template<class TVar>
inline int TPZFMatrix<TVar>::Zero() {
    int64_t size = this->fRow * this->fCol * sizeof(TVar);
    memset(((void*)this->fElem),'\0',size);
    this->fDecomposed = 0;
    return( 1 );
}

/**************************/
/*** Operations Global ****/

inline int64_t Norm(const TPZFMatrix<int64_t> &A) {
    return (int64_t)sqrt((REAL)Dot(A,A));
}

inline int Norm(const TPZFMatrix<int> &A) {
    return (int)sqrt((REAL)Dot(A,A));
}

inline float Norm(const TPZFMatrix<float> &A) {
    return sqrt(Dot(A,A));
}

inline double Norm(const TPZFMatrix<double> &A) {
    return sqrt(Dot(A,A));
}

inline long double Norm(const TPZFMatrix<long double> &A) {
    return sqrt(Dot(A,A));
}

inline float Norm(const TPZFMatrix< std::complex <float> > &A) {
    return sqrt(Dot(A,A).real());
}

inline double Norm(const TPZFMatrix< std::complex <double> > &A) {
    return sqrt(Dot(A,A).real());
}

inline long double Norm(const TPZFMatrix< std::complex <long double> > &A) {
    return sqrt(Dot(A,A).real());
}

#ifdef _AUTODIFF
inline float Norm(const TPZFMatrix< Fad <float> > &A) {
    return (float)TPZExtractVal::val(sqrt(Dot(A,A)));
}

inline double Norm(const TPZFMatrix< Fad <double> > &A) {
    return TPZExtractVal::val(sqrt(Dot(A,A)));
}

inline long double Norm(const TPZFMatrix< Fad <long double> > &A) {
    return TPZExtractVal::val(sqrt(Dot(A,A)));
}
#endif

inline TPZFlopCounter Norm(const TPZFMatrix<TPZFlopCounter> &A)
{
    return sqrt(Dot(A, A));
}
/**
 * @brief Non abstract class which implements full matrices with preallocated storage with (N+1) entries. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<int N, class TVar=REAL>
class TPZFNMatrix : public TPZFMatrix<TVar> {
    
    TVar fBuf[N+1];
    
public:
    /*
     * @brief Constructor which does not initialize the data. \n
     * WARNING : this class will dynamically allocate memory if the template parameter N is smaller than row*col
     * @param row Number of rows
     * @param col Number of cols
     */
    inline TPZFNMatrix(int64_t row, int64_t col) : TPZRegisterClassId(&TPZFNMatrix<N,TVar>::ClassId), 
    TPZFMatrix<TVar>(row,col,fBuf,N) {}
    
    inline TPZFNMatrix() : TPZRegisterClassId(&TPZFNMatrix<N,TVar>::ClassId), TPZFMatrix<TVar>(0,0,fBuf,N)
    {
    }
    
    inline TPZFNMatrix(const TPZFMatrix<TVar> &copy) : TPZRegisterClassId(&TPZFNMatrix<N,TVar>::ClassId), 
    TPZFMatrix<TVar>(0,0,fBuf,N)
    {
        *this = copy;
    }
    
    inline TPZFNMatrix(const TPZFNMatrix<N,TVar> &copy) : TPZRegisterClassId(&TPZFNMatrix<N,TVar>::ClassId), 
    TPZFMatrix<TVar>(0,0,fBuf,N)
    {
        *this = copy;
    }
    
    virtual ~TPZFNMatrix()
    {
    }
    
    CLONEDEF(TPZFNMatrix)
    /*
     * @brief Constructor which initializes the data. \n
     * WARNING : this class will dynamically allocate memory if the template parameter N is smaller than row*col
     * @param row Number of rows
     * @param col Number of cols
     * @param val initial value of the matrix elements
     */
    inline  TPZFNMatrix(int64_t row, int64_t col, const TVar &val) : TPZRegisterClassId(&TPZFNMatrix<N,TVar>::ClassId), 
    TPZFMatrix<TVar>(row,col,fBuf,N) {
        TPZFMatrix<TVar>::operator=(val);
    }
    
    inline  TPZFMatrix<TVar> &operator=(const TPZFMatrix<TVar> &copy) {
        return TPZFMatrix<TVar>::operator=(copy);
    }
    inline  TPZFNMatrix<N, TVar> &operator=(const TPZFNMatrix<N, TVar> &copy) {
        TPZFMatrix<TVar>::operator=(copy);
        return *this;
    }
    
    virtual int ClassId() const;
    void Read(TPZStream& buf, void* context);
    void Write(TPZStream& buf, int withclassid) const;
    
};

template<int N, class TVar>
int TPZFNMatrix<N, TVar>::ClassId() const{
    return Hash("TPZFNMatrix") ^ TPZFMatrix<TVar>::ClassId() << 1 ^ (N << 2);
}

template<int N, class TVar>
void TPZFNMatrix<N, TVar>::Read(TPZStream& buf, void* context) {
    TPZFMatrix<TVar>::Read(buf, context);
    buf.Read(fBuf, N+1);
}

template<int N, class TVar>
void TPZFNMatrix<N, TVar>::Write(TPZStream& buf, int withclassid) const {
    TPZFMatrix<TVar>::Write(buf, withclassid);
    buf.Write(fBuf, N+1);
}

#endif

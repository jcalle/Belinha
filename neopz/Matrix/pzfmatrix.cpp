/**
 * @file
 * @brief Contains the implementation of the TPZFMatrix<>methods.
 */


#include "pzfmatrix.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <complex>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include "pzerror.h"
#include "pzaxestools.h"
#include "pzextractval.h"
#include "pzlog.h"
#include "pzmatrix.h"
#include "TPZSavable.h"
#include "pzvec.h"
#include "tpzverysparsematrix.h"

#ifdef _AUTODIFF
#include "tfad.h"
#include "fad.h"
#endif

class TPZStream;

#ifdef PZDEBUG
#define DEBUG2
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzfmatrix"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.checkconsistency"));
#endif

#ifdef USING_LAPACK
/** CBlas Math Library */
#include "TPZLapack.h"
#define BLAS_MULT
#endif


//#define IsZero( a )  ( fabs(a) < 1.e-20)


using namespace std;

/********************/
/*** Constructors ***/

template <class TVar>
TPZFMatrix<TVar>::TPZFMatrix(const TPZMatrix<TVar> &mat) : 
TPZRegisterClassId(&TPZFMatrix::ClassId),
TPZMatrix<TVar>(mat), fElem(0),fGiven(0),fSize(0) {
    if(this->fRow*this->fCol) {
        
        fElem = new TVar[this->fRow*this->fCol];
        TVar * p = fElem;
        int64_t i,j;
        for(j=0; j<this->fCol; j++) {
            for(i=0; i<this->fRow; i++) {
                *p++ = mat.GetVal(i,j);
            }
        }
    }
}


/********************************/
/*** Constructor( TPZFMatrix<TVar> & ) ***/

template<class TVar>
TPZFMatrix<TVar>::TPZFMatrix(const TPZFMatrix<TVar> &A)
: TPZRegisterClassId(&TPZFMatrix::ClassId),
TPZMatrix<TVar>( A.fRow, A.fCol ), fElem(0), fGiven(0), fSize(0) {
    int64_t size = this->fRow * this->fCol;
    if(!size) return;
    fElem = new TVar[ size ] ;
#ifdef PZDEBUG2
    if ( size && fElem == NULL ) Error( "Constructor <memory allocation error>." );
#endif
    // Copia a matriz
    TVar * src = A.fElem;
    TVar * p = fElem;
    memcpy((void *)(p),(void *)(src),(size_t)size*sizeof(TVar));
}

/********************************/
/*** Constructor( TPZVerySparseMatrix<TVar> & ) ***/

template<class TVar>
TPZFMatrix<TVar>::TPZFMatrix(TPZVerySparseMatrix <TVar> const & A)
: TPZRegisterClassId(&TPZFMatrix::ClassId),
TPZMatrix<TVar>( A.Rows(), A.Cols() ), fElem(0), fGiven(0), fSize(0) {
    
    int64_t size = this->fRow * this->fCol;
    if(!size) return;
    fElem = new TVar[ size ] ;
    
#ifdef PZDEBUG2
    if ( size && fElem == NULL ) Error( "Constructor <memory allocation error>." );
#endif
    
    typename std::map <std::pair<int64_t, int64_t>, TVar>::const_iterator it = A.MapBegin();
    typename std::map <std::pair<int64_t, int64_t>, TVar>::const_iterator end = A.MapEnd();
    
    TVar * p = fElem;
    memset((void *)p, 0, (size_t)size*sizeof(TVar));
    
    for (; it != end; it++) {
        const std::pair<int64_t, int64_t>& key = it->first;
        PutVal(key.first, key.second, it->second);
    }
    
}



/******** Operacoes com matrizes FULL  ********/

/******************/
/*** Operator = ***/
template<class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator=(const TPZFMatrix<TVar> &A ) {
    if(this == &A) return *this;
    int64_t size = A.fRow * A.fCol;
    
    TVar * newElem = fElem;
    if(fSize < size && size != this->fRow*this->fCol) {
        newElem = new TVar
        [size] ;
    } else if (fSize >= size) {
        newElem = fGiven;
    }
    
    if ( newElem == NULL && size > 0) Error( "Operator= <memory allocation error>." );
    if (fElem && fElem != newElem && fElem != fGiven) delete[]( fElem );
    this->fRow  = A.fRow;
    this->fCol  = A.fCol;
    fElem = newElem;
    
    // Copia a matriz
    memcpy((void *)(fElem),(void *)(A.fElem),(size_t)size*sizeof(TVar));
    
    TPZMatrix<TVar>::operator=(A);
    
    
    return *this;
}

template <class TVar>
void TPZFMatrix<TVar>::AddFel(TPZFMatrix<TVar> &rhs,TPZVec<int64_t> &destination) {
    if(rhs.Cols() != this->Cols()) {
        PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
        DebugStop();
        return;
    }
    int64_t ncol = this->Cols();
    int64_t nrow = rhs.Rows();
    int64_t i,j;
    for(j=0; j<ncol; j++) {
        for(i=0; i<nrow; i++) {
            operator()(destination[i],j) += rhs(i,j);
        }
    }
}

template<class TVar>
void TPZFMatrix<TVar>::AddFel(TPZFMatrix<TVar> &rhs,TPZVec<int64_t> &source, TPZVec<int64_t> &destination) {
    if(rhs.Cols() != this->Cols() && source.NElements()) {
        PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
        DebugStop();
        return;
    }
    int64_t ncol = this->Cols();
    int64_t nrow = source.NElements();
    int64_t i,j;
    for(j=0; j<ncol; j++) {
        for(i=0; i<nrow; i++) {
            operator()(destination[i],j) += rhs(source[i],j);
        }
    }
}

template<>
void TPZFMatrix<double>::AddFel(TPZFMatrix<double> &rhs,TPZVec<int64_t> &source, TPZVec<int64_t> &destination) {
    if(rhs.Cols() != this->Cols() && source.NElements()) {
        PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
        DebugStop();
        return;
    }
    int64_t ncol = this->Cols();
    int64_t nrow = source.NElements();
    int64_t i,j;
    for(j=0; j<ncol; j++) {
        for(i=0; i<nrow; i++) {
#pragma omp atomic
            operator()(destination[i],j) += rhs(source[i],j);
        }
    }
}

template<>
void TPZFMatrix<float>::AddFel(TPZFMatrix<float> &rhs,TPZVec<int64_t> &source, TPZVec<int64_t> &destination) {
    if(rhs.Cols() != this->Cols() && source.NElements()) {
        PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
        DebugStop();
        return;
    }
    int64_t ncol = this->Cols();
    int64_t nrow = source.NElements();
    int64_t i,j;
    for(j=0; j<ncol; j++) {
        for(i=0; i<nrow; i++) {
#pragma omp atomic
            operator()(destination[i],j) += rhs(source[i],j);
        }
    }
}



/*******************************/
/*** Operator+( TPZFMatrix>& ) ***/

template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator+(const TPZFMatrix<TVar> &A ) const {
    if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
        Error( "Operator+ <matrixs with different dimensions>" );
    
    TPZFMatrix<TVar> res;
    res.Redim( this->Rows(), this->Cols() );
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    TVar * pm = fElem, *plast = fElem+size;
    TVar * pa = A.fElem;
    TVar * pr = res.fElem;
    
    while(pm < plast) *pr++ = (*pm++) + (*pa++);
    
    return( res );
}

/*******************************/
/*** Operator-( TPZFMatrix<>& ) ***/
template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator-(const TPZFMatrix<TVar> &A ) const {
    if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
        Error( "Operator- <matrixs with different dimensions>" );
    
    TPZFMatrix<TVar> res;
    res.Redim( this->Rows(), this->Cols() );
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    TVar * pm = fElem;
    TVar * pa = A.fElem;
    TVar * pr = res.fElem, *prlast =pr+size;
    
    while(pr < prlast) *pr++ = (*pm++) - (*pa++);
    return( res );
}

template<>
void TPZFMatrix<int>::GramSchmidt(TPZFMatrix<int> &Orthog, TPZFMatrix<int> &TransfToOrthog)
{
    std::cout << "Nothing to do\n";
    DebugStop();
}

#ifdef _AUTODIFF
template <>
void TPZFMatrix<TFad<6,REAL> >::GramSchmidt(TPZFMatrix<TFad<6,REAL> > &Orthog, TPZFMatrix<TFad<6, REAL> > &TransfToOrthog)
{
    DebugStop();
}

template <>
void TPZFMatrix<Fad<double> >::GramSchmidt(TPZFMatrix<Fad<double> > &Orthog, TPZFMatrix<Fad<double> > &TransfToOrthog)
{
    DebugStop();
}

#endif

template <class TVar>
void TPZFMatrix<TVar>::GramSchmidt(TPZFMatrix<TVar> &Orthog, TPZFMatrix<TVar> &TransfToOrthog)
{
#ifdef LOG4CXX2
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        Print("GrSchmidt Entrada",sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    double scale = 1.;
    for(int64_t j = 0; j < this->Cols(); j++)
    {
        double norm = 0.;
        for(int64_t i = 0; i < this->Rows(); i++)
        {
            norm += fabs(TPZExtractVal::val(this->GetVal(i,j)*this->GetVal(i,j)));
        }
        norm = sqrt(norm);
        if(norm > 1.e-10)
        {
            if((1.)/norm > scale) scale = (1.)/norm;
        }
    }
    
    this->operator *=( scale );
    
    int64_t QTDcomp = this->Rows();
    int64_t QTDvec = this->Cols();
    Orthog.Resize(QTDcomp,QTDvec);
    Orthog.Zero();
    /// Making a copy of *this (Ortog = *this)
    for(int64_t r = 0; r < QTDcomp; r++)
    {
        for(int64_t c = 0; c < QTDvec; c++)
        {
            Orthog(r,c) = GetVal(r,c);
        }
    }
    
#ifdef PZDEBUG
    int check = 0;
    for(int64_t c = 0; c < QTDvec; c++)
    {
        TVar summ = 0.;
        for(int64_t r = 0; r < QTDcomp; r++)
        {
            summ += fabs(GetVal(r,c));
        }
        if(fabs(summ) < 0.00001)
        {
            std::stringstream sout;
            sout << "Null Vector on Gram-Schmidt Method! Col = " << c << "\n";
            LOGPZ_ERROR(logger,sout.str())
            check = 1;
        }
    }
#endif
    
    TVar dotUp, dotDown;
    for(int64_t c = 1; c < QTDvec; c++)
    {
        for(int64_t stop = 0; stop < c; stop++)
        {
            dotUp = 0.;
            dotDown = 0.;
            for(int64_t r = 0; r < QTDcomp; r++)
            {
                dotUp += GetVal(r,c)*Orthog(r,stop);
                dotDown += Orthog(r,stop)*Orthog(r,stop);
            }
            if(fabs(dotDown) < 1.E-8)
            {
#ifdef PZDEBUG
                if(check == 0)
                {
                    std::stringstream sout;
                    sout << "Parallel Vectors on Gram-Schmidt Method! Col = " << stop << "\n";
                    LOGPZ_ERROR(logger,sout.str())
                }
#endif
                
                for(int64_t r = 0; r < QTDcomp; r++)
                {
                    Orthog(r,stop) = 0.;
                }
            }
            else
            {
#ifdef LOG4CXX2
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "dotdown = " << dotDown << " dotup = " << dotUp;
                    LOGPZ_DEBUG(logger,sout.str())
                }
#endif
                
                for(int64_t r = 0; r < QTDcomp; r++)
                {
                    Orthog(r,c) -= dotUp*Orthog(r,stop)/dotDown;
                }
            }
        }
    }
    for(int64_t c = 0; c < QTDvec; c++)
    {
        dotUp = 0.;
        for(int64_t r = 0; r < QTDcomp; r++)
        {
            dotUp += Orthog(r,c)*Orthog(r,c);
        }
        if(fabs(dotUp) > 1.e-8)
        {
            for(int64_t r = 0; r < QTDcomp; r++)
            {
                Orthog(r,c) = Orthog(r,c)/sqrt(dotUp);
            }
        }
        else {
#ifdef LOG4CXX
            std::stringstream sout;
            sout << "Linearly dependent columns dotUp = " << dotUp;
            LOGPZ_ERROR(logger,sout.str())
#endif
            
            for(int64_t r = 0; r < QTDcomp; r++)
            {
                Orthog(r,c) = 0.;
            }
        }
        
    }
    Orthog.Multiply(*this,TransfToOrthog,1);
    
    this->operator*= ( 1./scale );
    TransfToOrthog.operator*= ( 1./scale );
    
#ifdef LOG4CXX2
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << endl;
        sout << "Output GS" << endl;
        Orthog.Print("Orthog matrix",sout);
        TransfToOrthog.Print("TransfToOrthog matrix",sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
#ifdef PZDEBUG
    TPZFNMatrix<9, TVar> OrthogT;
    Orthog.Transpose(&OrthogT);
    TPZAxesTools<TVar>::VerifyAxes(OrthogT);
#endif
}

template <>
void TPZFMatrix<TPZFlopCounter>::GramSchmidt(TPZFMatrix<TPZFlopCounter> &Orthog, TPZFMatrix<TPZFlopCounter> &TransfToOrthog)
{
    std::cout << __PRETTY_FUNCTION__ << " please implement me\n";
    DebugStop();
}

template <class TVar>
void TPZFMatrix<TVar>::DeterminantInverse(TVar &determinant, TPZFMatrix<TVar> &inverse)
{
    TPZFNMatrix<100, TVar> copy(*this);
    inverse.Redim(this->Rows(),this->Rows());
    int64_t r;
    for(r=0; r<this->Rows(); r++) inverse(r,r) = 1.;
    copy.Solve_LU(&inverse);
    determinant = 1.;
    for(r=0; r<this->Rows(); r++) determinant *= copy(r,r);
}

#ifdef USING_LAPACK

template <class TVar>
/** @brief Initialize pivot with i = i  */
void TPZFMatrix<TVar>::InitializePivot()
{
    fPivot.Resize(this->Rows());
    for(int64_t i = 0; i < this->Rows(); i++){
        fPivot[i] = i+1; // Fortran based indexing
    }
}

#endif

template <class TVar>
void TPZFMatrix<TVar>::MultAdd(const TVar *ptr, int64_t rows, int64_t cols, const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                               const TVar alpha,const TVar beta ,const int opt)
{
    
    
    if ((!opt && cols != x.Rows()) || (opt && rows != x.Rows())) {
        Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
        return;
    }
    if(beta != (TVar)0. && ((!opt && rows != y.Rows()) || (opt && cols != y.Rows()) || y.Cols() != x.Cols())) {
        Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
        return;
    }
    if(!opt) {
        if(z.Cols() != x.Cols() || z.Rows() != rows) {
            z.Redim(rows,x.Cols());
        }
    } else {
        if(z.Cols() != x.Cols() || z.Rows() != cols) {
            z.Redim(cols,x.Cols());
        }
    }
    unsigned numeq = opt ? cols : rows;
    int64_t xcols = x.Cols();
    int64_t ic, c;
    if(!(rows*cols)) return;
    for (ic = 0; ic < xcols; ic++) {
        TVar *zp = &z(0,ic), *zlast = zp+numeq;
        if(beta != (TVar)0.) {
            const TVar *yp = &y.g(0,ic);
            if(&z != &y) {
                memcpy((void *)zp,(void *)yp,numeq*sizeof(TVar));
            }
            for(int64_t i=0; i< numeq; i++) for(int64_t c=0; c<xcols; c++) z(i,c) *= beta;
        } else {
            while(zp != zlast) {
                *zp = 0.;
                zp ++;
            }
        }
    }
    
    
    for (ic = 0; ic < xcols; ic++) {
        if(!opt) {
            for ( c = 0; c<cols; c++) {
                TVar * zp = &z(0,ic), *zlast = zp+rows;
                const TVar * fp = ptr +rows*c;
                const TVar * xp = &x.g(c,ic);
                while(zp < zlast) {
                    *zp += alpha* *fp++ * *xp;
                    zp ++;
                }
            }
        } else {
            const TVar * fp = ptr;
            TVar *zp = &z(0,ic);
            for (c = 0; c<cols; c++) {
                TVar val = 0.;
                // bug correction philippe 5/2/97
                //					 REAL * xp = &x(0,ic), xlast = xp + numeq;
                const TVar *xp = &x.g(0,ic);
                const TVar *xlast = xp + rows;
                while(xp < xlast) {
                    val += *fp++ * *xp;
                    xp ++;
                }
                *zp += alpha *val;
                zp ++;
            }
        }
    }
    
    
    
    
    
}

#ifdef USING_LAPACK
template<>
void TPZFMatrix<double>::MultAdd(const TPZFMatrix<double> &x,const TPZFMatrix<double> &y, TPZFMatrix<double> &z,
                                 const double alpha,const double beta,const int opt) const {
    
#ifdef PZDEBUG
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows())) {
        Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
        return;
    }
    if(beta != (double)0. && ((!opt && this->Rows() != y.Rows()) || (opt && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
        Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
        return;
    }
#endif
    if(!opt) {
        if(z.Cols() != x.Cols() || z.Rows() != this->Rows()) {
            z.Redim(this->Rows(),x.Cols());
        }
    } else {
        if(z.Cols() != x.Cols() || z.Rows() != this->Cols()) {
            z.Redim(this->Cols(),x.Cols());
        }
    }
    if(this->Cols() == 0) {
        z.Zero();
        if (beta != 0) {
            z = y;
            z *= beta;
        }
        return;
    }
    if (beta != (double)0.) {
        z = y;
    }
    if (Rows() == 0 || Cols() == 0 || x.Rows() == 0 || x.Cols() == 0) {
        return;
    }
    if (!opt) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, this->Rows(), x.Cols(), this->Cols(),
                    alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), beta, z.fElem, z.Rows());
    } else {
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, this->Cols(), x.Cols(), this->Rows(),
                    alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), beta, z.fElem, z.Rows());
    }
    
}
template<>
void TPZFMatrix<float>::MultAdd(const TPZFMatrix<float> &x,const TPZFMatrix<float> &y, TPZFMatrix<float> &z,
                                const float alpha,const float beta,const int opt) const {
    
#ifdef PZDEBUG
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows())) {
        Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
        return;
    }
    if(beta != (float)0. && ((!opt && this->Rows() != y.Rows()) || (opt && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
        Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
        return;
    }
#endif
    if(!opt) {
        if(z.Cols() != x.Cols() || z.Rows() != this->Rows()) {
            z.Redim(this->Rows(),x.Cols());
        }
    } else {
        if(z.Cols() != x.Cols() || z.Rows() != this->Cols()) {
            z.Redim(this->Cols(),x.Cols());
        }
    }
    if(this->Cols() == 0) {
        z.Zero();
    }
    if (beta != (float)0.) {
        z = y;
    }
    if (!opt) {
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, this->Rows(), x.Cols(), this->Cols(),
                    alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), beta, z.fElem, z.Rows());
    } else {
        cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, this->Cols(), x.Cols(), this->Rows(),
                    alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), beta, z.fElem, z.Rows());
    }
    
}

template<>
void TPZFMatrix<std::complex<double> >::MultAdd(const TPZFMatrix<std::complex<double> > &x,const TPZFMatrix<std::complex<double> > &y, TPZFMatrix<std::complex<double> > &z,
                                                const std::complex<double> alpha,const std::complex<double> beta,const int opt) const {
    
#ifdef PZDEBUG
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows())) {
        Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
        return;
    }
    if(beta.real() != 0. && ((!opt && this->Rows() != y.Rows()) || (opt && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
        Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
        return;
    }
#endif
    if(!opt) {
        if(z.Cols() != x.Cols() || z.Rows() != this->Rows()) {
            z.Redim(this->Rows(),x.Cols());
        }
    } else {
        if(z.Cols() != x.Cols() || z.Rows() != this->Cols()) {
            z.Redim(this->Cols(),x.Cols());
        }
    }
    if(this->Cols() == 0) {
        z.Zero();
    }
    if (beta.real() != 0.) {
        z = y;
    }
    if (!opt) {
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, this->Rows(), x.Cols(), this->Cols(),
                    &alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), &beta, z.fElem, z.Rows());
    } else {
        cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans, this->Cols(), x.Cols(), this->Rows(),
                    &alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), &beta, z.fElem, z.Rows());
    }
    
}

#endif // USING_LAPACK

/**
 * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
 * @param x Is x on the above operation
 * @param y Is y on the above operation
 * @param z Is z on the above operation
 * @param alpha Is alpha on the above operation
 * @param beta Is beta on the above operation
 * @param opt Indicates if is Transpose or not
 */
template <class TVar>
void TPZFMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                               const TVar alpha,const TVar beta,const int opt) const {
    
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows())) {
        Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
        return;
    }
    if(beta != (TVar)0. && ((!opt && this->Rows() != y.Rows()) || (opt && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
        Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
        return;
    }
    if(!opt) {
        if(z.Cols() != x.Cols() || z.Rows() != this->Rows()) {
            z.Redim(this->Rows(),x.Cols());
        }
    } else {
        if(z.Cols() != x.Cols() || z.Rows() != this->Cols()) {
            z.Redim(this->Cols(),x.Cols());
        }
    }
    if(this->Cols() == 0)
    {
        z.Zero();
    }
    unsigned numeq = opt ? this->Cols() : this->Rows();
    int64_t rows = this->Rows();
    int64_t cols = this->Cols();
    int64_t xcols = x.Cols();
    int64_t ic, c;
    if (numeq)
    {
        for (ic = 0; ic < xcols; ic++) {
            TVar *zp = &z(0,ic), *zlast = zp+numeq;
            if(beta != (TVar)0.) {
                const TVar *yp = &y.g(0,ic);
                if(&z != &y) {
                    memcpy((void *)zp,(void *)yp,numeq*sizeof(TVar));
                }
                for(int64_t i=0; i< numeq; i++) z(i,ic) *= beta;
                
            } else {
                while(zp != zlast) {
                    *zp = 0.;
                    zp ++;
                }
            }
        }
    }
    
    if(!(rows*cols)) return;
    
    for (ic = 0; ic < xcols; ic++) {
        if(!opt) {
            for ( c = 0; c<cols; c++) {
                TVar * zp = &z(0,ic), *zlast = zp+rows;
                TVar * fp = fElem +rows*c;
                const TVar * xp = &x.g(c,ic);
                while(zp < zlast) {
                    *zp += alpha* *fp++ * *xp;
                    zp ++;
                }
            }
        } else {
            TVar * fp = fElem,  *zp = &z(0,ic);
            for (c = 0; c<cols; c++) {
                TVar val = 0.;
                // bug correction philippe 5/2/97
                //					 REAL * xp = &x(0,ic), xlast = xp + numeq;
                const TVar *xp = &x.g(0,ic);
                const TVar *xlast = xp + rows;
                while(xp < xlast) {
                    val += *fp++ * *xp;
                    xp ++;
                }
                *zp += alpha *val;
                zp ++;
            }
        }
    }
    
}

/********************************/
/*** Operator+=( TPZFMatrix<>& ) ***/
template <class TVar>
TPZFMatrix<TVar> & TPZFMatrix<TVar>::operator+=(const TPZFMatrix<TVar> &A ) {
    if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
        Error( "Operator+= <matrixs with different dimensions>" );
    
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    TVar * pm = fElem, *pmlast=pm+size;
    TVar * pa = A.fElem;
    while(pm < pmlast) (*pm++) += (*pa++);
    return( *this );
}

/*******************************/
/*** Operator-=( TPZFMatrix<>& ) ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator-=(const TPZFMatrix<TVar> &A ) {
    if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
        Error( "Operator-= <matrixs with different dimensions>" );
    
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    TVar * pm = fElem;
    TVar * pa = A.fElem;
    
    for ( int64_t i = 0; i < size; i++ ) *pm++ -= *pa++;
    
    return( *this );
}

#ifdef USING_LAPACK
template <>
void TPZFMatrix<double>::ZAXPY(const double alpha,const TPZFMatrix<double> &p) {
    // 	//Como definir o tamanho dos vetores
    int size  = (fRow*fCol) ;
    cblas_daxpy(size, alpha, &p.fElem[0], 1, &fElem[0], 1);
}
#endif

template <class TVar>
void TPZFMatrix<TVar>::ZAXPY(const TVar alpha,const TPZFMatrix<TVar> &p) {
    
    TVar * pt = fElem;
    TVar * pp = p.fElem;
    TVar * ptlast = fElem + this->fRow*this->fCol;
    while(pt < ptlast) *pt++ += alpha * *pp++;
}

#ifdef USING_LAPACK
template<>
void TPZFMatrix<double>::TimesBetaPlusZ(const double beta,const TPZFMatrix<double> &z)
{
    int size = fRow*fCol;
    cblas_dscal(size,beta,fElem,1);
    cblas_daxpy(size,1.,z.fElem,1,fElem,1);
}
#endif

template<class TVar>
void TPZFMatrix<TVar>::TimesBetaPlusZ(const TVar beta,const TPZFMatrix<TVar> &z) {
    // #ifdef USING_ATLAS
    // 	int size = fRow*fCol;
    // 	cblas_dscal(size,beta,fElem,1);
    // 	cblas_daxpy(size,1.,z.fElem,1,fElem,1);
    // #else
    
    TVar * pt = fElem,  *ptlast = fElem + this->fRow*this->fCol;
    TVar * pz = z.fElem;
    //	while(pt < ptlast) *pt++ *= (beta) + *pz++;
    while(pt < ptlast) {
        *pt *= (beta);
        *pt++ += *pz++;
    }
    // #endif
}

/******** Operacoes com MATRIZES GENERICAS ********/

/******************/
/*** Operator = ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator=(const TPZMatrix<TVar> &A ) {
    int64_t arows  = A.Rows();
    int64_t acols  = A.Cols();
    int64_t size = arows * acols;
    if(fElem != fGiven) {
        delete []fElem;
        fElem = 0;
    }
    this->fRow  =  arows;
    this->fCol  = acols;
    if(fSize < size) {
        fElem = new TVar[ arows * acols ] ;
    } else {
        fElem = fGiven;
    }
    TVar * dst = fElem;
    for ( int64_t c = 0; c < this->fCol; c++ )
        for ( int64_t r = 0; r < this->fRow; r++ )
            *dst++ = A.Get( r, c );
    return( *this );
}


/******** Operacoes com valores NUMERICOS ********/

/******************/
/*** Operator = ***/
template <class TVar>
TPZFMatrix<TVar>& TPZFMatrix<TVar>::operator=(const TVar value ) {
    int64_t size = ((int64_t)this->fRow) * this->fCol;
    TVar * dst   = fElem;
    for ( int64_t i = 0; i < size; i++ )
        *dst++ = value;
    this->fDecomposed = 0;
    return *this;
}



/***************************/
/*** Operator+=( value ) ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator+=(const TVar value ) {
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    
    TVar * dst = fElem, *dstlast = dst+size;
    while ( dst < dstlast ) *dst++ += value;
    return( *this );
}



/**************************/
/*** Operator+( value ) ***/
template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator+(const TVar value ) const {
    TPZFMatrix<TVar> res( *this );
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    
    TVar * dst = res.fElem,  *dstlast = dst+size;
    while ( dst < dstlast )
        *dst++ += value;
    
    return( res );
}

template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator-  (const TVar val ) const {
    return operator+( -val );
}


/**************************/
/*** Operator*( value ) ***/

template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator*(const TVar value ) const
{
    TPZFMatrix<TVar> res( *this );
    res *= value;
    return( res );
}

/***************************/
/*** Operator*=( value ) ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator*=( const TVar value ) {
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    TVar * dst = fElem, *dstlast = dst+size;
    while ( dst < dstlast ) *dst++ *= value;
    return( *this );
}

/**************/
/*** Resize ***/
template <class TVar>
int TPZFMatrix<TVar>::Resize(const int64_t newRows,const int64_t newCols) {
    if ( newRows == this->Rows() && newCols == this->Cols() ) return( 1 );
    int64_t newsize = ((int64_t)newRows)*newCols;
    TVar * newElem;
    if(fGiven && fElem != fGiven && newsize <= fSize)
    {
        newElem = fGiven;
    } else
    {
        newElem = new TVar[ newRows * newCols ] ;
    }
    if ( newElem == NULL )
        Error( "Resize <memory allocation error>." );
    
    int64_t minRow  = ( this->fRow < newRows ? this->fRow : newRows );
    int64_t minCol  = ( this->fCol < newCols ? this->fCol : newCols );
    TVar * src;
    TVar * dst;
    int64_t r, c;
    
    for ( c = 0; c < minCol; c++ ) {
        // Copia as linhas da matriz antiga para a nova.
        // Copia os elementos de uma linha.
        dst = newElem + c*newRows;
        src = fElem + c*this->fRow;
        for ( r = 0; r < minRow; r++ ) *dst++ = *src++;
        
        // Se a nova linha for maior (mais colunas), preenche o
        //  resto da linha com ZEROS.
        for ( ; r < newRows; r++ ) *dst++ = 0.0;
    }
    
    // Preenche as linha que sobrarem (se sobrarem) com ZEROS.
    for ( ;c < newCols; c++ )
    {
        dst = newElem + c*newRows;
        for (r = 0 ; r < newRows; r++ ) *dst++ = 0.0;
    }
    
    if (fElem && fElem != fGiven )delete[]( fElem );
    fElem = newElem;
    this->fRow  = newRows;
    this->fCol  = newCols;
    return( 1 );
}

template <class TVar>
int TPZFMatrix<TVar>::Remodel(const int64_t newRows,const int64_t newCols) {
    if(newRows*newCols != this->fRow*this->fCol) return -1;
    this->fRow = newRows;
    this->fCol = newCols;
    return 1;
}

/********************/
/*** Transpose () ***/
template <class TVar>
void TPZFMatrix<TVar>::Transpose(TPZMatrix<TVar> *const T) const{
    T->Resize( this->Cols(), this->Rows() );
    //Transposta por filas
    TVar * p = fElem;
    for ( int64_t c = 0; c < this->Cols(); c++ ) {
        for ( int64_t r = 0; r < this->Rows(); r++ ) {
            T->PutVal( c, r, *p++ );
            //            cout<<"(r,c)= "<<r<<"  "<<c<<"\n";
        }
    }
}

template <class TVar>
void TPZFMatrix<TVar>::Transpose() {
    TPZFMatrix<TVar> temp;
    Transpose(&temp);
    *this = temp;
}

#ifdef USING_LAPACK
template <>
int TPZFMatrix<float>::Decompose_LU(TPZVec<int> &index) {
    
    
    if (this->fDecomposed != ENoDecompose && this->fDecomposed != ELUPivot) DebugStop();
    
    if (this->fDecomposed != ENoDecompose) {
        return ELUPivot;
    }
    
    if ( this->Rows() != this->Cols() ) {
        cout << "TPZFPivotMatrix::DecomposeLU ERRO : A Matriz não é quadrada" << endl;
        return 0;
    }
    
    
    int nRows = this->Rows();
    int zero = 0;
    float b;int info;
    
    fPivot.Resize(nRows);
    
    //    int sgesv_(__CLPK_integer *__n, __CLPK_integer *__nrhs, __CLPK_real *__a,
    //               __CLPK_integer *__lda, __CLPK_integer *__ipiv, __CLPK_real *__b,
    //               __CLPK_integer *__ldb,
    //               __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
    //                                                                __IPHONE_4_0);
    
    
    sgesv_(&nRows,&zero,fElem,&nRows,&fPivot[0],&b,&nRows,&info);
    index = fPivot;
    this->fDecomposed = ELUPivot;
    return 1;
}
template <>
int TPZFMatrix<double>::Decompose_LU(TPZVec<int> &index) {
    
    
    if (this->fDecomposed != ENoDecompose && this->fDecomposed != ELUPivot) DebugStop();
    
    if (this->fDecomposed != ENoDecompose) {
        return ELUPivot;
    }
    
    if ( this->Rows() != this->Cols() ) {
        cout << "TPZFPivotMatrix::DecomposeLU ERRO : A Matriz não é quadrada" << endl;
        return 0;
    }
    
    
    int nRows = this->Rows();
    int zero = 0;
    double b;int info;
    
    fPivot.Resize(nRows);
    
    //    int sgesv_(__CLPK_integer *__n, __CLPK_integer *__nrhs, __CLPK_real *__a,
    //               __CLPK_integer *__lda, __CLPK_integer *__ipiv, __CLPK_real *__b,
    //               __CLPK_integer *__ldb,
    //               __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
    //                                                                __IPHONE_4_0);
    
    
    dgesv_(&nRows,&zero,fElem,&nRows,&fPivot[0],&b,&nRows,&info);
    index = fPivot;
    this->fDecomposed = ELUPivot;
    return 1;
}
#endif //USING_LAPACK

template <class TVar>
int TPZFMatrix<TVar>::Decompose_LU(TPZVec<int> &index) {
    
    if (this->fDecomposed) return 0;
    
    if ( this->Rows() != this->Cols() ) {
        cout << "TPZFPivotMatrix::DecomposeLU ERRO : A Matriz não é quadrada" << endl;
        return 0;
    }
    
    int64_t i,j,k;
    TVar sum = 0.;
    int64_t nRows = this->Rows();
    int64_t nCols = this->Cols();
    
    index.Resize(nRows);
    //inicializo o vetor de índices para o caso de pivotamento
    for (i=0;i<nRows;i++) index[i] = i;
    
    //LU
    for (j=0;j<nCols;j++){
        // cout << "line..." << j << endl;
        for (i=0;i<=j;i++){
            sum = 0.;
            for (k=0;k<i;k++){
                sum += this->Get(i,k) * this->Get(k,j);
            }
            TVar aux = this->Get(i,j);
            PutVal(i,j,(aux - sum));
            //cout << "0_A[" << i << "," << j << "]= " << Get(i,j) << endl;
        }
        //Print(cout);
        TVar piv = this->Get(j,j);
        //  cout << "Pivo 1 =" << piv << endl;
        int64_t row = j;
        for (i=j+1;i<nRows;i++){
            sum = 0.;
            for (k=0;k<(j);k++){
                sum += this->Get(i,k) * this->Get(k,j);
            }
            TVar aux = this->Get(i,j);
            PutVal(i,j,(aux - sum));
            //cout << "1_A[" << i << "," << j << "]= " << Get(i,j) << endl;
            
            if (fabs(this->Get(i,j)) > fabs(piv)){
                piv = this->Get(i,j);
                //  cout << "Pivo 2 =" << piv << endl;
                row = i;
            }
        }
        //    Print(cout);
        if (row > j){
            for (k=0;k<nCols;k++){
                TVar aux = this->Get(j,k);
                PutVal(j,k,this->Get(row,k));
                //cout << "2_A[" << j << "," << k << "]= " << Get(j,k) << endl;
                PutVal(row,k,aux);
                //cout << "3_A[" << row << "," << k << "]= " << Get(row,k) << endl;
            }
            k = index[j];
            index[j] = index[row];
            index[row] = k;
        }
        //    cout << "Pivo = " << piv << endl;
        for (i=j+1;i<nRows;i++){
            if (fabs(piv) < fabs((TVar)1e-12)) {
                cout << "Pivot < 1e-12. Probably matrix is singular." << endl;
                DebugStop();
            }
            TVar aux = this->Get(i,j) / piv;
            PutVal(i,j,aux);
            //cout << "4_A[" << i << "," << j << "]= " << Get(i,j) << endl;
        }
        //Print(cout);
    }
    this->fDecomposed = ELUPivot;
    return 1;
}


/*****************/
/*** DecomposeLU ***/
template <class TVar>
int TPZFMatrix<TVar>::Decompose_LU(std::list<int64_t> &singular) {
    //return Decompose_LU();
#ifndef USING_LAPACK
    if (  this->fDecomposed && this->fDecomposed != ELU)  Error( "Decompose_LU <Matrix already Decomposed with other scheme>" );
#else
    if (  this->fDecomposed && this->fDecomposed != ELUPivot)  Error( "Decompose_LU <Matrix already Decomposed with other scheme>" );
#endif
    if (this->fDecomposed) return 1;
    
    const int  min = ( this->Cols() < (this->Rows()) ) ? this->Cols() : this->Rows();
    const int nrows = this->Rows();
    const int ncols = this->Cols();
    
    for ( int k = 0; k < min ; k++ ) {
        TVar pivot = GetVal(k, k);
        if (IsZero( pivot )){
            singular.push_back(k);
            PutVal(k,k,1.);
            pivot = 1.;
        }
        for ( int i = k+1; i < nrows; i++ ) {
            this->operator()(i,k) /= pivot;
        }
        for ( int j = k+1; j < ncols; j++ ){
            int i = k+1;
            TVar * elemPtr = &(this->operator()(i,j));
            TVar * ikPtr = &(this->operator()(i,k));
            const TVar kjVal = GetVal(k,j);
            for ( ; i < nrows; i++, elemPtr++, ikPtr++ ) {
                (*elemPtr) -= (*ikPtr)*kjVal;
                //          this->operator()(i,j) -= GetVal( i, k )*GetVal(k,j);
            }
        }
    }
    this->fDecomposed=ELU;
#ifdef USING_LAPACK
    fPivot.resize(nrows);
    for (int i=0; i<nrows; i++) {
        fPivot[i] = i+1;
    }
    this->fDecomposed = ELUPivot;
#endif //USING_LAPACK
    return 1;
}

#ifdef USING_LAPACK
template <>
int TPZFMatrix<float>::Decompose_LU() {
    
    
    return this->Decompose_LU(fPivot);
}
template <>
int TPZFMatrix<double>::Decompose_LU() {
    
    
    return this->Decompose_LU(fPivot);
}
#endif //USING_LAPACK

template <class TVar>
int TPZFMatrix<TVar>::Decompose_LU() {
    
    std::list<int64_t> fake;
    return this->Decompose_LU(fake);
}


#ifdef _AUTODIFF
template <class TVar>
int TPZFMatrix<TVar>::Substitution(const TVar *ptr, int64_t rows, TPZFMatrix<TVar> *B) {
    std::cout << __PRETTY_FUNCTION__ << " bailing out\n";
    DebugStop();
    return 1;
}
#else
template <class TVar>
int TPZFMatrix<TVar>::Substitution(const TVar *ptr, int64_t rows, TPZFMatrix<TVar> *B)
{
    int64_t rowb = B->Rows();
    int64_t colb = B->Cols();
    if ( rowb != rows ) Error( "static::SubstitutionLU <incompatible dimensions>" );
    int64_t i,j;
    for ( i = 0; i < rowb; i++ ) {
        for ( int64_t col = 0; col < colb; col++ )
            for (j = 0; j < i; j++ )
                //B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
                PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - SELECTEL(ptr, rows, i, j) * GETVAL(B, rowb, j, col));
    }
    
    for (int64_t col=0; col<colb; col++){
        for ( i = rowb-1; i >= 0; i-- ) {
            for (j = i+1; j < rowb ; j++ )
                //B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
                PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - SELECTEL(ptr, rows, i, j) * GETVAL(B, rowb, j, col));
            if ( IsZero( SELECTEL(ptr, rows, i, i)/*GetVal(i, i)*/ ) ) {
                if (fabs(SELECTEL(ptr, rows, i, i)/*GetVal(i, i)*/) > fabs((TVar)0.)) {
                    if (fabs(GETVAL(B, rowb, i, col) - SELECTEL(ptr, rows, i, i)/*B->GetVal(i, col) - GetVal(i, i)*/) > fabs(((TVar)1e-12))) {
                        Error( "static::BackSub(SubstitutionLU) <Matrix is singular even after Power Plus..." );
                    }
                }else  Error( "static::BackSub(SubstitutionLU) <Matrix is singular" );
            }
            PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col)/SELECTEL(ptr, rows, i, i));
            //B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
        }
    }
    return( 1 );
}
#endif

#ifdef _AUTODIFF
#include "fadType.h"
#endif

/****************/
/*** Substitution ***/
template <class TVar>
int TPZFMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B ) const {
#ifdef USING_LAPACK    
	if (this->fDecomposed != ELUPivot) {
		Error("TPZFMatrix::Decompose_LU substitution called for a wrongly decomposed matrix");
	}
#else
	if(this->fDecomposed != ELU) {
        Error("TPZFMatrix::Decompose_LU substitution called for a wrongly decomposed matrix");
    }
#endif
    int64_t rowb = B->Rows();
    int64_t colb = B->Cols();
    int64_t row = this->Rows();
    if ( rowb != this->Rows() ) Error( "SubstitutionLU <incompatible dimensions>" );
    
    
    int64_t i,j;
    for ( i = 0; i < rowb; i++ ) {
        for ( int64_t col = 0; col < colb; col++ )
            for (j = 0; j < i; j++ )
                //B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
                PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - GETVAL(this, row, i, j) * GETVAL(B, rowb, j, col));
    }
    
    for (int64_t col=0; col<colb; col++){
        for ( i = rowb-1; i >= 0; i-- ) {
            for (j = i+1; j < rowb ; j++ )
                //B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
                PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - GETVAL(this, row, i, j) * GETVAL(B, rowb, j, col));
            if ( IsZero( GETVAL(this, row, i, i)/*GetVal(i, i)*/ ) ) {
                if (fabs(GETVAL(this, row, i, i)/*GetVal(i, i)*/) > 0.){
                    TVar diff = GETVAL(B, rowb, i, col) - GETVAL(this, row, i, i)/*B->GetVal(i, col) - GetVal(i, i)*/;
                    if (fabs(diff) > 1e-12){
                        Error( "BackSub(SubstitutionLU) <Matrix is singular even after Power Plus..." );
                    }
                }else  Error( "BackSub(SubstitutionLU) <Matrix is singular" );
            }
            PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col)/GETVAL(this, row, i, i));
            //B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
        }
    }
    return( 1 );
    
}


#ifdef USING_LAPACK
template<>
int TPZFMatrix<float>::Substitution( TPZFMatrix<float> *B, const TPZVec<int> &index ) const{
    
    if(!B){
        PZError << __PRETTY_FUNCTION__ << "TPZFMatrix<>*B eh nulo" << endl;
        return 0;
    }
    
    TPZFMatrix<float> &b = *B;
    
    if (!this->fDecomposed){
        PZError <<  __PRETTY_FUNCTION__ << "Matriz não decomposta" << endl;
        return 0;
    }
    
    if (this->fDecomposed != ELUPivot){
        PZError << __PRETTY_FUNCTION__ << "\nfDecomposed != ELUPivot" << endl;
    }
    
    //    int sgetrs_(char *__trans, __CLPK_integer *__n, __CLPK_integer *__nrhs,
    //                __CLPK_real *__a, __CLPK_integer *__lda, __CLPK_integer *__ipiv,
    //                __CLPK_real *__b, __CLPK_integer *__ldb,
    //                __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
    //                                                                 __IPHONE_4_0);
    int nRows = this->Rows();
    char notrans = 'N';
    int BCols = B->Cols();
    int info = 0;
    
    sgetrs_(&notrans,&nRows,&BCols,fElem,&nRows,&fPivot[0],B->fElem,&nRows,&info);
    
#ifdef PZDEBUG
    if(info != 0)
    {
        DebugStop();
    }
#endif
    
    return 1;
}

template<>
int TPZFMatrix<double>::Substitution( TPZFMatrix<double> *B, const TPZVec<int> &index ) const{
    
    if(!B){
        PZError << __PRETTY_FUNCTION__ << "TPZFMatrix<>*B eh nulo" << endl;
        return 0;
    }
    
    
    if (!this->fDecomposed){
        PZError <<  __PRETTY_FUNCTION__ << "Matriz não decomposta" << endl;
        return 0;
    }
    
    if (this->fDecomposed != ELUPivot){
        PZError << __PRETTY_FUNCTION__ << "\nfDecomposed != ELUPivot" << endl;
    }
    
    //    int sgetrs_(char *__trans, __CLPK_integer *__n, __CLPK_integer *__nrhs,
    //                __CLPK_real *__a, __CLPK_integer *__lda, __CLPK_integer *__ipiv,
    //                __CLPK_real *__b, __CLPK_integer *__ldb,
    //                __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
    //                                                                 __IPHONE_4_0);
    int nRows = this->Rows();
    char notrans = 'N';
    int BCols = B->Cols();
    int info = 0;
    
    dgetrs_(&notrans,&nRows,&BCols,fElem,&nRows,&fPivot[0],B->fElem,&nRows,&info);
    
#ifdef PZDEBUG
    if(info != 0)
    {
        DebugStop();
    }
#endif
    
    return 1;
}

#endif //USING_LAPACK


template<class TVar>
int TPZFMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B, const TPZVec<int> &index ) const{
    
    if(!B){
        PZError << __PRETTY_FUNCTION__ << "TPZFMatrix<>*B eh nulo" << endl;
        return 0;
    }
    
    TPZFMatrix<TVar> &b = *B;
    
    if (!this->fDecomposed){
        PZError <<  __PRETTY_FUNCTION__ << "Matriz não decomposta" << endl;
        return 0;
    }
    
    if (this->fDecomposed != ELUPivot){
        PZError << __PRETTY_FUNCTION__ << "\nfDecomposed != ELUPivot" << endl;
    }
    
    int nRows = this->Rows();
    
    if (index.NElements() != nRows || b.Rows() != nRows)
    {
        cout << "TMatrix::Substituicao ERRO : vetores com dimensões incompativeis:\n"
        << "this->fIndex = " << index.NElements() << "  b = " << b.Rows() << endl;
        return 0;
    }
    
    int64_t ncols = B->Cols();
    for(int64_t ic = 0; ic<ncols; ic++)
    {
        int64_t i,j;
        TVar sum = 0;
        
        TPZVec<TVar> v(nRows);
        
        
        for (i=0;i<nRows;i++)
        {
            v[i] = b(index[i],ic);
        }
        
        //Ly=b
        for (i=0;i<nRows;i++)
        {
            sum = 0.;
            for (j=0;j<(i);j++) sum +=this->Get(i,j) * v[j];
            v[i] -= sum;
        }
        
        //Ux=y
        for (i=(nRows-1);i>-1;i--)
        {
            sum = 0.;
            for (j=(i+1);j<nRows;j++) sum += this->Get(i,j) * v[j];
            v[i] = (v[i] - sum) / this->Get(i,i);
        }
        
        for (i=0;i<nRows;i++) b(i,ic) = v[i];
    }
    return 1;
}

#ifdef USING_LAPACK
template <>
int TPZFMatrix<float>::Substitution( TPZFMatrix<float> *B ) const {
    
    return this->Substitution(B,fPivot);
}
template <>
int TPZFMatrix<double>::Substitution( TPZFMatrix<double> *B ) const {
    
    return this->Substitution(B,fPivot);
}
#endif //USING_LAPACK


//NAO TESTADO
template <class TVar>
int TPZFMatrix<TVar>::Decompose_Cholesky(){
    std::list<int64_t> fake;
    int res = this->Decompose_Cholesky(fake);
    if(fake.size()){
        DebugStop();
    }
    return res;
}

#ifdef USING_LAPACK
template <>
int TPZFMatrix<float>::Decompose_Cholesky(std::list<int64_t> &singular) {
    if (  this->fDecomposed && this->fDecomposed != ECholesky) Error( "Decompose_Cholesky <Matrix already Decomposed>" );
    if (  this->fDecomposed ) return ECholesky;
    if ( this->Rows() != this->Cols() ) Error( "Decompose_Cholesky <Matrix must be square>" );
    int dim=this->Dim();
    
    TPZFMatrix<float> B(*this);
    int nrhs = 0;
    float *A = fElem;
    char uplo = 'U';
    int info;
    //    sposv_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
    spotrf_(&uplo, &dim, A, &dim, &info);
    this->fDecomposed = ECholesky;
    
    if (info != 0) {
        DebugStop();
    }
    return 1;
}
template <>
int TPZFMatrix<double>::Decompose_Cholesky(std::list<int64_t> &singular) {
    if (  this->fDecomposed && this->fDecomposed != ECholesky) Error( "Decompose_Cholesky <Matrix already Decomposed>" );
    if (  this->fDecomposed ) return ECholesky;
    if ( this->Rows() != this->Cols() ) Error( "Decompose_Cholesky <Matrix must be square>" );
    int dim=this->Dim();
    
    double B;
    int nrhs = 0;
    double *A = fElem;
    char uplo = 'U';
    int info;
    dpotrf_(&uplo, &dim, A, &dim, &info);
    this->fDecomposed = ECholesky;
    
    if (info != 0) {
        DebugStop();
    }
    return 1;
}
#endif //USING_LAPACK

template <class TVar>
int TPZFMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular) {
    
    if (  this->fDecomposed && this->fDecomposed != ECholesky) Error( "Decompose_Cholesky <Matrix already Decomposed>" );
    if (  this->fDecomposed ) return ECholesky;
    if ( this->Rows() != this->Cols() ) Error( "Decompose_Cholesky <Matrix must be square>" );
    //return 0;
    
    int dim=this->Dim();
    
    for (int i=0 ; i<dim; i++) {
        
        TVar * diagPtr = &(this->operator()(i,i));
        for(int k=0; k<i; k++) {             //elementos da diagonal
            (*diagPtr) -= this->operator()(i,k)*this->operator()(i,k);
        }
        
        
        if( IsZero(*diagPtr) ){
            singular.push_back(i);
            (*diagPtr) = 1.;
        }
        
        (*diagPtr) = sqrt(*diagPtr);
        
        for (int j=i+1;j<dim; j++) {           //elementos fora da diagonal
            TVar sum = 0.;
            { ///escopo
                int k=0;
                TVar * ikPtr = &(this->operator()(k,i));//&(this->operator()(i,k));///(k,i) = (i,k) pela simetria da matriz, mas o alinhamento acelera a execucao
                TVar * kjPtr = &(this->operator()(k,j));
                for(; k<i; k++, kjPtr++, ikPtr++) {
                    sum += (*ikPtr)*(*kjPtr);
                }
            }
            TVar *ijPtr = &(this->operator()( i,j ));
            (*ijPtr) -= sum;
            
            (*ijPtr) /= (*diagPtr);
            this->operator()(j,i) = (*ijPtr);
        }
    }
    
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    this->fDecomposed = ECholesky;
    return ECholesky;
}

template <class TVar>
int TPZFMatrix<TVar>::Substitution(const TVar *ptr, int64_t rows, TPZFMatrix<TVar> *B, const TPZVec<int> &index )
{
    
    if(!B){
        PZError << __PRETTY_FUNCTION__ << "TPZFMatrix<>*B eh nulo" << endl;
        return 0;
    }
    
    TPZFMatrix<TVar> &b = *B;
    
    
    
    if (index.NElements() != rows || b.Rows() != rows)
    {
        cout << "TMatrix::Substituicao ERRO : vetores com dimensões incompativeis:\n"
        << "this->fIndex = " << index.NElements() << "  b = " << b.Rows() << endl;
        return 0;
    }
    
    int64_t i,j;
    TVar sum = 0;
    
    TPZVec<TVar> v(rows);
    
    
    for (i=0;i<rows;i++)
    {
        v[i] = b(index[i]);
    }
    
    //Ly=b
    for (i=0;i<rows;i++)
    {
        sum = 0.;
        for (j=0;j<(i);j++) sum += SELECTEL(ptr,rows,i,j) * v[j];
        v[i] -= sum;
    }
    
    //Ux=y
    for (i=(rows-1);i>-1;i--)
    {
        sum = 0.;
        for (j=(i+1);j<rows;j++) sum += SELECTEL(ptr,rows,i,j) * v[j];
        v[i] = (v[i] - sum) / SELECTEL(ptr,rows,i,i);
    }
    
    for (i=0;i<rows;i++) b(i) = v[i];
    return 1;
}

#ifdef USING_LAPACK
template <>
int TPZFMatrix<float>::Decompose_LDLt() {
    
    if (  this->fDecomposed && this->fDecomposed != ELDLt) {
        Error( "Decompose_LDLt <Matrix already Decomposed with other scheme> " );
    } else if(this->fDecomposed ) {
        return ELDLt;
    }
    if ( this->Rows()!=this->Cols() ) Error( "Decompose_LDLt <Matrix must be square>" );
    char uplo = 'U';
    int dim = Rows();
    int nrhs = 0;
    fPivot.Resize(dim,0);
    float B  = 0.;
    int worksize = 3*dim;
    fWork.Resize(worksize);
    int info;
    
    if (dim == 0) {
        this->fDecomposed  = ELDLt;
        this->fDefPositive = 0;
        return( 1 );
    }
    
    //    ssysv_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_real *__work#>, <#__CLPK_integer *__lwork#>, <#__CLPK_integer *__info#>)
    
    ssysv_(&uplo, &dim, &nrhs, fElem, &dim, &fPivot[0], &B, &dim, &fWork[0], &worksize, &info);
    fDecomposed = ELDLt;
    return 1;
}

template <>
int TPZFMatrix<double>::Decompose_LDLt() {
    
    if (  this->fDecomposed && this->fDecomposed != ELDLt) {
        Error( "Decompose_LDLt <Matrix already Decomposed with other scheme> " );
    } else if(this->fDecomposed ) {
        return ELDLt;
    }
    if ( this->Rows()!=this->Cols() ) Error( "Decompose_LDLt <Matrix must be square>" );
    char uplo = 'L';
    int dim = Rows();
    int nrhs = 0;
    fPivot.Resize(dim,0);
    double B  = 0.;
    int worksize = 3*dim;
    fWork.Resize(worksize);
    int info;
    
    if (dim == 0) {
        this->fDecomposed  = ELDLt;
        this->fDefPositive = 0;
        return( 1 );
    }
    
    //    ssysv_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_real *__work#>, <#__CLPK_integer *__lwork#>, <#__CLPK_integer *__info#>)
    
    dsysv_(&uplo, &dim, &nrhs, fElem, &dim, &fPivot[0], &B, &dim, &fWork[0], &worksize, &info);
    fDecomposed = ELDLt;
    return 1;
}

#endif //USING_LAPACK

template <class TVar>
int TPZFMatrix<TVar>::Decompose_LDLt() {
    
    if (  this->fDecomposed && this->fDecomposed != ELDLt) {
        Error( "Decompose_LDLt <Matrix already Decomposed with other scheme> " );
    } else if(this->fDecomposed ) {
        return ELDLt;
    }
    if ( this->Rows()!=this->Cols() ) Error( "Decompose_LDLt <Matrix must be square>" );
    
    int64_t j,k,l,dim=this->Rows();
    
    for ( j = 0; j < dim; j++ ) {
        for ( k=0; k<j; k++) {
            PutVal( j,j,GetVal(j,j) - GetVal(k,k)*GetVal(k,j)*GetVal(k,j) );
        }
        for ( k=0; k<j; k++) {
            for( l=j+1; l<dim;l++) {
                PutVal(l,j, GetVal(l,j)-GetVal(k,k)*GetVal(j,k)*GetVal(l,k) );
                PutVal(j,l,GetVal(l,j) );
            }
        }
        TVar tmp = GetVal(j,j);
        if ( IsZero(tmp) ) Error( "Decompose_LDLt <Zero on diagonal>" );
        for( l=j+1; l<dim;l++) {
            PutVal(l,j, GetVal(l,j)/GetVal(j,j) ) ;
            PutVal(j,l, GetVal(l,j) );
        }
    }
    this->fDecomposed  = ELDLt;
    this->fDefPositive = 0;
    return( 1 );
}


#ifdef USING_LAPACK
/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<float>::Subst_Forward( TPZFMatrix<float>* b ) const
{
    if (fDecomposed == ECholesky) {
        //        CALL strsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs,
        //                   179      $               one, a, lda, b, ldb )
        char left[]= "Left", upper[] = "Upper", transpose[] = "Transpose", non_unit[] = "Non-unit";
        int ncol = b->Cols();
        int dim = Rows();
        //        b->Print("before =",std::cout,EMathematicaInput);
        cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, dim, ncol, 1., fElem, dim, b->fElem, dim);
        //        TPZMatrix<float>::Subst_Forward(b);
        //        b->Print("after =",std::cout,EMathematicaInput);
        return 1;
    }
    else
    {
        return TPZMatrix<float>::Subst_Forward(b);
    }
}

/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<double>::Subst_Forward( TPZFMatrix<double>* b ) const
{
    if (fDecomposed == ECholesky) {
        //        CALL strsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs,
        //                   179      $               one, a, lda, b, ldb )
        char left[]= "Left", upper[] = "Upper", transpose[] = "Transpose", non_unit[] = "Non-unit";
        int ncol = b->Cols();
        int dim = Rows();
        //        b->Print("before =",std::cout,EMathematicaInput);
        cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, dim, ncol, 1., fElem, dim, b->fElem, dim);
        //        TPZMatrix<float>::Subst_Forward(b);
        //        b->Print("after =",std::cout,EMathematicaInput);
        return 1;
    }
    else
    {
        return TPZMatrix<double>::Subst_Forward(b);
    }
}

/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZFMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar>* b ) const
{
    return TPZMatrix<TVar>::Subst_Forward(b);
}

/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZFMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar>* b ) const
{
    return TPZMatrix<TVar>::Subst_Backward(b);
}


template<>
int TPZFMatrix<float>::Subst_Backward( TPZFMatrix<float>* b ) const
{
    if (fDecomposed == ECholesky) {
        
        //        spotrs_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
        //        CALL strsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n,
        //                   184      $               nrhs, one, a, lda, b, ldb )
        char left[]= "Left", upper[] = "Upper", transpose[] = "Transpose", non_unit[] = "Non-unit";
        int ncol = b->Cols();
        int dim = Rows();
        int info = 0;
        //        b->Print("before =",std::cout,EMathematicaInput);
        //        TPZMatrix<float>::Subst_Backward(b);
        cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, dim, ncol, 1., fElem, dim, b->fElem, dim);
        //        spotrs_(upper, &dim, &ncol, fElem, &dim, b->fElem, &ncol, &info);
        //        b->Print("after =",std::cout,EMathematicaInput);
        
        return 1;
        
        if (info != 0) {
            DebugStop();
        }
    }
    else
    {
        return TPZMatrix<float>::Subst_Backward(b);
    }
}

template<>
int TPZFMatrix<double>::Subst_Backward( TPZFMatrix<double>* b ) const
{
    if (fDecomposed == ECholesky) {
        
        //        spotrs_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
        //        CALL strsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n,
        //                   184      $               nrhs, one, a, lda, b, ldb )
        char left[]= "Left", upper[] = "Upper", transpose[] = "Transpose", non_unit[] = "Non-unit";
        int ncol = b->Cols();
        int dim = Rows();
        int info = 0;
        //        b->Print("before =",std::cout,EMathematicaInput);
        //        TPZMatrix<float>::Subst_Backward(b);
        cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, dim, ncol, 1., fElem, dim, b->fElem, dim);
        //        spotrs_(upper, &dim, &ncol, fElem, &dim, b->fElem, &ncol, &info);
        //        b->Print("after =",std::cout,EMathematicaInput);
        
        return 1;
        
        if (info != 0) {
            DebugStop();
        }
    }
    else
    {
        return TPZMatrix<double>::Subst_Backward(b);
    }
}

/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<float>::Subst_LForward( TPZFMatrix<float>* b ) const
{
    //    ssytrs_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
    
    if (fDecomposed != ELDLt) {
        DebugStop();
    }
    
    char uplo = 'U';
    int dim = Rows();
    int nrhs = b->Cols();
    float B  = 0.;
    int info;
    
    //    ssytrs_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
    ssytrs_(&uplo, &dim, &nrhs, fElem, &dim, &fPivot[0], b->fElem, &dim, &info);
    return 1;
    //    return TPZMatrix<TVar>::Subst_LForward(b);
}

/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<double>::Subst_LForward( TPZFMatrix<double>* b ) const
{
   
    if (fDecomposed != ELDLt) {
        DebugStop();
    }
    
    char uplo = 'L';
    int dim = Rows();
    int nrhs = b->Cols();
    double B  = 0.;
    int info;
    if (dim == 0 || nrhs == 0) {
        return 0;
    }
    dsytrs_(&uplo, &dim, &nrhs, fElem, &dim, &fPivot[0], b->fElem, &dim, &info);
    return 1;
    //    return TPZMatrix<TVar>::Subst_LForward(b);
}


/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZFMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar>* b ) const
{
    //    ssytrs2
    return TPZMatrix<TVar>::Subst_LForward(b);
}

/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<float>::Subst_LBackward( TPZFMatrix<float>* b ) const
{
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<double>::Subst_LBackward( TPZFMatrix<double>* b ) const
{
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZFMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar>* b ) const
{
    return TPZMatrix<TVar>::Subst_LBackward(b);
}

/**
 * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<float>::Subst_Diag( TPZFMatrix<float>* b ) const
{
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<double>::Subst_Diag( TPZFMatrix<double>* b ) const
{
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZFMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar>* b ) const
{
    return TPZMatrix<TVar>::Subst_Diag(b);
}
#endif //USING_LAPACK

/** @brief Implement dot product for matrices */
template<class TVar>
TVar Dot(const TPZFMatrix<TVar> &A, const TPZFMatrix<TVar> &B) {
    int64_t size = (A.Rows())*A.Cols();
    TVar result = 0.;
    if(!size) return result;
    // #ifdef USING_ATLAS
    // 	result = cblas_ddot(size, &A.g(0,0), 1, &B.g(0,0), 1);
    // 	return result;
    //
    // #elif USING_BLAS
    // 	result = cblas_ddot(size, &A.g(0,0), 1, &B.g(0,0), 1);
    // 	return result;
    //
    // #else
    const TVar *fpA = &A.g(0,0), *fpB = &B.g(0,0);
    const TVar *fpLast = fpA+size;
    while(fpA < fpLast)
    {
        result += (*fpA++ * *fpB++);
    }
    return result;
    // #endif
}

template
std::complex<float> Dot(const TPZFMatrix< std::complex<float> > &A, const TPZFMatrix< std::complex<float> > &B);

template
std::complex<double> Dot(const TPZFMatrix< std::complex<double> > &A, const TPZFMatrix< std::complex<double> > &B);

template
std::complex<long double> Dot(const TPZFMatrix< std::complex<long double> > &A, const TPZFMatrix< std::complex<long double> > &B);

template
long double Dot(const TPZFMatrix<long double> &A, const TPZFMatrix<long double> &B);

template
double Dot(const TPZFMatrix<double> &A, const TPZFMatrix<double> &B);

template
float Dot(const TPZFMatrix<float> &A, const TPZFMatrix<float> &B);

template
int64_t Dot(const TPZFMatrix<int64_t> &A, const TPZFMatrix<int64_t> &B);

template
int Dot(const TPZFMatrix<int> &A, const TPZFMatrix<int> &B);

#ifdef _AUTODIFF
template
Fad<float> Dot(const TPZFMatrix<Fad<float> > &A, const TPZFMatrix<Fad<float> > &B);

template
Fad<double> Dot(const TPZFMatrix<Fad<double> > &A, const TPZFMatrix<Fad<double> > &B);

template
Fad<long double> Dot(const TPZFMatrix<Fad<long double> > &A, const TPZFMatrix<Fad<long double> > &B);
#endif

template
TPZFlopCounter Dot(const TPZFMatrix<TPZFlopCounter> &A, const TPZFMatrix<TPZFlopCounter> &B);

/** @brief Increments value over all entries of the matrix A. */
template <class TVar>
TPZFMatrix<TVar> operator+(const TVar value, const TPZFMatrix<TVar> &A ) {
    return( A + value );
}

/** @brief Decrements value over all entries of the matrix A. */
template <class TVar>
TPZFMatrix<TVar> operator-(const TVar value, const TPZFMatrix<TVar> &A ) {
    return( A - value );
}

/************************** Private **************************/

/*************/
/*** Error ***/
template<class TVar>
int TPZFMatrix<TVar>::Error(const char *msg1,const char *msg2 ) {
    ostringstream out;
    out << "TPZFMatrix::" << msg1;
    if(msg2) out << msg2;
    out << ".\n";
    LOGPZ_ERROR (logger, out.str().c_str());
    DebugStop();
    return 0;
}


/*************/
/*** Clear ***/
template <class TVar>
int TPZFMatrix<TVar>::Clear() {
    if(fElem && fElem != fGiven) delete[]( fElem );
    fElem = NULL;
    this->fRow  = this->fCol = 0;
    return( 1 );
}

#ifdef _AUTODIFF
template <>
void TPZFMatrix<TFad<6,REAL> >::Read( TPZStream &buf, void *context ){
    DebugStop();
}

template <>
void TPZFMatrix<TFad<6,REAL> >::Write( TPZStream &buf, int withclassid ) const {
    DebugStop();
}

template <>
void TPZFMatrix<Fad<REAL> >::Read( TPZStream &buf, void *context ){
    DebugStop();
}

template <>
void TPZFMatrix<Fad<REAL> >::Write( TPZStream &buf, int withclassid ) const {
    DebugStop();
}
#endif

template <class TVar>
void TPZFMatrix<TVar>::Read( TPZStream &buf, void *context ){ //ok
    TPZMatrix<TVar>::Read(buf,context);
    int64_t row = this->fRow;
    int64_t col = this->fCol;
    //this is odd, but necessary.
    this->fRow = this->fCol = 0;
    Resize(row,col);
    buf.Read(fElem,this->fRow*this->fCol);
}

//template <class TVar>
//void TPZFMatrix<TVar>::Write( TPZStream &buf, int withclassid ) {
//    const TPZFMatrix<TVar> *cp = this;
//    cp->Write(buf,withclassid);
//    //    const Write(buf, withclassid);
//    //	TPZMatrix<TVar>::Write(buf,withclassid);
//    //	buf.Write(fElem,this->fRow*this->fCol);
//}

template <class TVar>
void TPZFMatrix<TVar>::Write( TPZStream &buf, int withclassid ) const { //ok
    TPZMatrix<TVar>::Write(buf,withclassid);
    buf.Write(fElem,this->fRow*this->fCol);
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
template<class TVar>
bool TPZFMatrix<TVar>::Compare(TPZSavable *copy, bool override)
{
    TPZFMatrix<TVar> *fmat = dynamic_cast<TPZFMatrix<TVar> *> (copy);
    if(!fmat) return false;
    
    bool matresult = TPZMatrix<TVar>::Compare(copy,false);
    int64_t nel = this->fRow*this->fCol;
    TVar diff=0.;
    int64_t numdif = 0;
    int64_t iel;
    for(iel=0; iel<nel; iel++)
    {
        if(fElem[iel] != fmat->fElem[iel])
        {
            matresult = false;
            numdif++;
        }
        TVar exp = fElem[iel]-fmat->fElem[iel];
        diff += fabs(exp);
    }
    if(!matresult)
    {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << " did not compare ";
        sout << " number different terms " << numdif << " number terms " << this->fRow*this->fCol;
        sout << " difference in norm L1 " << diff;
        LOGPZ_ERROR(loggerCheck,sout.str())
    }
    if(!matresult && override)
    {
        this->operator=(*fmat);
    }
    return matresult;
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
template<class TVar>
bool TPZFMatrix<TVar>::Compare(TPZSavable *copy, bool override) const
{
    TPZFMatrix<TVar> *fmat = dynamic_cast<TPZFMatrix<TVar> *> (copy);
    if(!fmat) return false;
    
    bool matresult = TPZMatrix<TVar>::Compare(copy,false);
    int64_t nel = this->fRow*this->fCol;
    TVar diff=0.;
    int64_t numdif = 0;
    int64_t iel;
    for(iel=0; iel<nel; iel++)
    {
        if(fElem[iel] != fmat->fElem[iel])
        {
            matresult = false;
            numdif++;
        }
        TVar exp = fElem[iel]-fmat->fElem[iel];
        diff += fabs(exp);
    }
    if(!matresult)
    {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << " did not compare ";
        sout << " number different terms " << numdif << " number terms " << this->fRow*this->fCol;
        sout << " difference in norm L1 " << diff;
        LOGPZ_ERROR(loggerCheck,sout.str())
    }
    if(!matresult && override)
    {
        DebugStop();
    }
    return matresult;
}

#ifdef _AUTODIFF
template <>
void TPZFMatrix<TFad<6,REAL> >::PrintStatic(const TFad<6,REAL> *ptr, int64_t rows, int64_t cols, const char *name, std::ostream& out,const MatrixOutputFormat form)
{
    DebugStop();
}
#endif

template <class TVar>
void TPZFMatrix<TVar>::PrintStatic(const TVar *ptr, int64_t rows, int64_t cols, const char *name, std::ostream& out,const MatrixOutputFormat form){
    
    if(form == EFormatted) {
        out << "Writing matrix '";
        if(name) out << name;
        out << "' (" << rows << " x " << cols << "):\n";
        
        for ( int64_t row = 0; row < rows; row++) {
            out << "\t";
            for ( int64_t col = 0; col < cols; col++ ) {
                out << SELECTEL(ptr,rows, row, col) << "  ";
            }
            out << "\n";
        }
        out << "\n";
    } else if (form == EInputFormat) {
        out << rows << " " << cols << endl;
        for ( int64_t row = 0; row < rows; row++) {
            for ( int64_t col = 0; col < cols; col++ ) {
                TVar val = SELECTEL(ptr,rows,row, col);
                if(val != (TVar)0.) out << row << ' ' << col << ' ' << val << std::endl;
            }
        }
        out << "-1 -1 0.\n";
    } else if( form == EMathematicaInput)
    {
        char number[32];
        out << name << "\n{ ";
        for ( int64_t row = 0; row < rows; row++) {
            out << "\n{ ";
            for ( int64_t col = 0; col < cols; col++ ) {
                TVar val = SELECTEL(ptr,rows,row, col);
                sprintf(number, "%16.16lf", (double)fabs(val));
                out << number;
                if(col < cols-1)
                    out << ", ";
                if((col+1) % 6 == 0)out << std::endl;
            }
            out << " }";
            if(row < rows-1)
                out << ",";
        }
        
        out << " }\n";
        
    }
    
}

template<class TVar>
int TPZFMatrix<TVar>::ClassId() const{
    return Hash("TPZFMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}

template <class TVar>
int TPZFMatrix<TVar>::SetSize(const int64_t newRows,const int64_t newCols) {
    int64_t newsize = ((int64_t)newRows)*newCols;
    int64_t oldsize = this->fRow*this->fCol;
    if(newsize == oldsize) return 1;
    if(fElem && fElem != fGiven)
    {
        delete []fElem;
        fElem = 0;
    }
    if(fGiven && newsize <= fSize)
    {
        fElem = fGiven;
    } else
    {
        fElem = new TVar[ newRows * newCols ] ;
    }
    if (newsize && fElem == NULL )
        Error( "Resize <memory allocation error>." );
    
    return( 1 );
}

#ifdef USING_LAPACK

template <class TVar>
int TPZFMatrix<TVar>::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues)
{
    DebugStop();
	return -1;
}

template <class TVar>
int TPZFMatrix<TVar>::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues, TPZFMatrix < std::complex<double> > &eigenvectors)
{
    DebugStop();
	return -1;
}

template <>
int TPZFMatrix<float>::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues)
{
    if (Rows() != Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< float > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<float> realeigen(dim,0.);
    TPZVec<float> imageigen(dim,0.);
    
    TPZFMatrix<float> temp(*this);
    TPZVec<float> work(lwork);
    sgeev_(jobvl, jobvr, &dim, temp.fElem, &dim, &realeigen[0], &imageigen[0], VL.fElem, &dim, VR.fElem, &dim, &work[0], &lwork, &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        
        eigenvalues[i] = realeigen[i] + I*imageigen[i];
    }
    return 1;
}

template <>
int TPZFMatrix<float>::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues, TPZFMatrix < std::complex<double> > &eigenvectors)
{
    if (Rows() != Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< float > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<float> realeigen(dim,0.);
    TPZVec<float> imageigen(dim,0.);
    
    TPZFMatrix<float> temp(*this);
    TPZVec<float> work(lwork);
    sgeev_(jobvl, jobvr, &dim, temp.fElem, &dim, &realeigen[0], &imageigen[0], VL.fElem, &dim, VR.fElem, &dim, &work[0], &lwork, &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = realeigen[i] + I*imageigen[i];
    }
    for(int i = 0 ; i < dim ; i ++){
        if(imageigen[i] == 0){
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i);
            }
        }
        else{
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i) + I * VR(iV,i+1) ;
                eigenvectors(iV,i + 1) = VR(iV,i) - I * VR(iV,i+1) ;
            }
            i++;
        }
    }
    
    return 1;
}

template <>
int TPZFMatrix<double>::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues)
{
    if (Rows() != Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< double > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<double> I(0,1.);
    TPZVec<double> realeigen(dim,0.);
    TPZVec<double> imageigen(dim,0.);
    
    TPZFMatrix<double> temp(*this);
    TPZVec<double> work(lwork);
    dgeev_(jobvl, jobvr, &dim, temp.fElem, &dim, &realeigen[0], &imageigen[0], VL.fElem, &dim, VR.fElem, &dim, &work[0], &lwork, &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = realeigen[i] + I*imageigen[i];
    }
    return 1;
}

template <>
int TPZFMatrix<double>::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues, TPZFMatrix < std::complex<double> > &eigenvectors)
{
    if (Rows() != Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< double > VL(Rows(),Cols(),0.),VR(Rows(),Cols(),0.);
    int dim = Rows();
//    double testwork;
    int lwork = 10+50*dim;
    int info = 0;
    std::complex<double> I(0,1.);
    TPZVec<double> realeigen(dim,0.);
    TPZVec<double> imageigen(dim,0.);
    
    TPZFMatrix<double> temp(*this);
    TPZVec<double> work(lwork,0.);
    dgeev_(jobvl, jobvr, &dim, temp.fElem, &dim, &realeigen[0], &imageigen[0], VL.fElem, &dim, VR.fElem, &dim, &work[0], &lwork, &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = realeigen[i] + I*imageigen[i];
    }
    for(int i = 0 ; i < dim ; i ++){
        if(imageigen[i] == 0){
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i);
            }
        }
        else{
            double *realVRptr = VR.fElem;
            double *imagVRptr = VR.fElem + dim;
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i) + I * VR(iV,i+1) ;
                eigenvectors(iV,i + 1) = VR(iV,i) - I * VR(iV,i+1) ;
            }
            i++;
        }
    }
    
    return 1;
}

template <>
int TPZFMatrix<complex<double> >::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues)
{
    if (Rows() != Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< complex<double> > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<double> I(0,1.);
    TPZVec<complex<double> > eigen(dim,0.);
    
    TPZFMatrix<complex<double> > temp(*this);
    TPZVec<complex<double> > work(lwork);
    TPZVec< double > rwork( 2 * dim);
    

    
    zgeev_(jobvl, jobvr, &dim, (vardoublecomplex *)temp.fElem, &dim, (vardoublecomplex *)&eigen[0], (vardoublecomplex *)VL.fElem, &dim, (vardoublecomplex *)VR.fElem, &dim, (vardoublecomplex *)&work[0], &lwork, &rwork[0], &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
        eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = eigen[i];
    }
    
    return 1;
}

template <>
int TPZFMatrix<complex<double> >::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues, TPZFMatrix < std::complex<double> > &eigenvectors)
{
    if (Rows() != Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< complex<double> > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<double> I(0,1.);
    TPZVec<complex<double> > eigen(dim,0.);
    
    TPZFMatrix<complex<double> > temp(*this);
    TPZVec<complex<double> > work(lwork);
    TPZVec< double > rwork( 2 * dim);
   
#ifdef USING_MKL
    typedef MKL_Complex16 vardoublecomplex;
#elif MACOSX
    typedef __CLPK_doublecomplex vardoublecomplex ;
#endif


    zgeev_(jobvl, jobvr, &dim, (vardoublecomplex *)temp.fElem, &dim, (vardoublecomplex *)&eigen[0], (vardoublecomplex *)VL.fElem, &dim, (vardoublecomplex *)VR.fElem, &dim, (vardoublecomplex *)&work[0], &lwork, &rwork[0], &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = eigen[i];
    }
    for(int i = 0 ; i < dim ; i ++){
    
        for( int iV = 0 ; iV < dim ; iV++ ){
            eigenvectors(iV,i) = VR(iV,i);
        }
    
    }
    
    return 1;
}

template <>
int TPZFMatrix<complex<float> >::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues)
{
    if (Rows() != Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< complex<float> > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<complex<float> > eigen(dim,0.);
    
    TPZFMatrix<complex<float> > temp(*this);
    TPZVec<complex<float> > work(lwork);
    TPZVec< float > rwork( 2 * dim);
    
    cgeev_(jobvl, jobvr, &dim, (varfloatcomplex *)temp.fElem, &dim, (varfloatcomplex *)&eigen[0], (varfloatcomplex *)VL.fElem, &dim, (varfloatcomplex *)VR.fElem, &dim, (varfloatcomplex *)&work[0], &lwork, &rwork[0], &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = eigen[i];
    }
	return 0;
}

template <>
int TPZFMatrix<complex< float> >::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues, TPZFMatrix < std::complex<double> > &eigenvectors)
{
    if (Rows() != Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< complex<float> > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<complex<float> > eigen(dim,0.);
    
    TPZFMatrix<complex<float> > temp(*this);
    TPZVec<complex<float> > work(lwork);
    TPZVec< float > rwork( 2 * dim);
    
    cgeev_(jobvl, jobvr, &dim, (varfloatcomplex *)temp.fElem, &dim, (varfloatcomplex *)&eigen[0], (varfloatcomplex *)VL.fElem, &dim, (varfloatcomplex *)VR.fElem, &dim, (varfloatcomplex *)&work[0], &lwork, &rwork[0], &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = eigen[i];
    }
    for(int i = 0 ; i < dim ; i ++){
        
        for( int iV = 0 ; iV < dim ; iV++ ){
            eigenvectors(iV,i) = VR(iV,i);
        }
        
    }
    
    return 1;
}


template< class TVar>
int
TPZFMatrix<TVar>::SolveGeneralisedEigenProblem(TPZFMatrix<TVar> &B , TPZVec < complex<double> > &w, TPZFMatrix < complex<double> > &eigenVectors)
{
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <LAPACK does not support this specific data type>" );
    return( 0 );
}
template< class TVar>
int
TPZFMatrix<TVar>::SolveGeneralisedEigenProblem(TPZFMatrix<TVar> &B , TPZVec < complex<double> > &w)
{
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <LAPACK does not support this specific data type>" );
    return( 0 );
}

template<>
int
TPZFMatrix<float>::SolveGeneralisedEigenProblem(TPZFMatrix<float> &B , TPZVec <complex<double> > &eigenvalues, TPZFMatrix < complex<double> > &eigenvectors)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< float > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<float> realeigen(dim,0.);
    TPZVec<float> imageigen(dim,0.);
    
    TPZVec<float> beta(dim);
    
    TPZFMatrix<float> temp(*this), tempB(B);
    TPZVec<float> work(lwork);

    sggev_(jobvl, jobvr, &dim, temp.fElem, &dim , tempB.fElem, &dim , &realeigen[0], &imageigen[0], &beta[0]  , VL.fElem, &dim , VR.fElem, &dim, &work[0], &lwork, &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = (realeigen[i] + I*imageigen[i]) / beta[i];
        }
    }
    for(int i = 0 ; i < dim ; i ++){
        if(imageigen[i] == 0){
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i);
            }
        }
        else{
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i) + I * VR(iV,i+1) ;
                eigenvectors(iV,i + 1) = VR(iV,i) - I * VR(iV,i+1) ;
            }
            i++;
        }
    }
    
    return 1;
}


template<>
int
TPZFMatrix<float>::SolveGeneralisedEigenProblem(TPZFMatrix<float> &B , TPZVec <complex<double> > &eigenvalues)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< float > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<float> realeigen(dim,0.);
    TPZVec<float> imageigen(dim,0.);
    
    TPZVec<float> beta(dim);
    
    TPZFMatrix<float> temp(*this), tempB(B);
    TPZVec<float> work(lwork);
    
    sggev_(jobvl, jobvr, &dim, temp.fElem, &dim , tempB.fElem, &dim , &realeigen[0], &imageigen[0], &beta[0]  , VL.fElem, &dim , VR.fElem, &dim, &work[0], &lwork, &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = (realeigen[i] + I*imageigen[i]) / beta[i];
        }
    }
    
    return 1;
}
template<>
int
TPZFMatrix<double>::SolveGeneralisedEigenProblem(TPZFMatrix<double> &B , TPZVec <complex<double> > &eigenvalues, TPZFMatrix < complex<double> > &eigenvectors)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< double > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<double> I(0,1.);
    TPZVec<double> realeigen(dim,0.);
    TPZVec<double> imageigen(dim,0.);
    
    TPZVec<double> beta(dim);
    
    TPZFMatrix<double> temp(*this), tempB(B);
    TPZVec<double> work(lwork);
    
    dggev_(jobvl, jobvr, &dim, temp.fElem, &dim , tempB.fElem, &dim , &realeigen[0], &imageigen[0], &beta[0]  , VL.fElem, &dim , VR.fElem, &dim, &work[0], &lwork, &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = (realeigen[i] + I*imageigen[i]) / beta[i];
        }
    }
    for(int i = 0 ; i < dim ; i ++){
        if(imageigen[i] == 0){
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i);
            }
        }
        else{
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i) + I * VR(iV,i+1) ;
                eigenvectors(iV,i + 1) = VR(iV,i) - I * VR(iV,i+1) ;
            }
            i++;
        }
    }
    
    return 1;
}


template<>
int
TPZFMatrix<double>::SolveGeneralisedEigenProblem(TPZFMatrix<double> &B , TPZVec <complex<double> > &eigenvalues)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< double > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<double> I(0,1.);
    TPZVec<double> realeigen(dim,0.);
    TPZVec<double> imageigen(dim,0.);
    
    TPZVec<double> beta(dim);
    
    TPZFMatrix<double> temp(*this), tempB(B);
    TPZVec<double> work(lwork);
    
    dggev_(jobvl, jobvr, &dim, temp.fElem, &dim , tempB.fElem, &dim , &realeigen[0], &imageigen[0], &beta[0]  , VL.fElem, &dim , VR.fElem, &dim, &work[0], &lwork, &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = (realeigen[i] + I*imageigen[i]) / beta[i];
        }
    }
    
    return 1;
}

template<>
int
TPZFMatrix<complex<float> >::SolveGeneralisedEigenProblem(TPZFMatrix<complex<float> > &B , TPZVec <complex<double> > &eigenvalues, TPZFMatrix < complex<double> > &eigenvectors)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< complex<float> > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    TPZVec<complex<float> > eigen(dim,0.);
    
    TPZVec<complex<float> > beta(dim);
    
    TPZFMatrix<complex<float> > temp(*this), tempB(B);
    TPZVec<complex<float> > work(lwork);
    TPZVec<float> rwork( 8 * dim );

    cggev_(jobvl, jobvr, &dim, (varfloatcomplex *)temp.fElem, &dim , (varfloatcomplex *)tempB.fElem, &dim , (varfloatcomplex *)&eigen[0], (varfloatcomplex *)&beta[0]  , (varfloatcomplex *)VL.fElem, &dim , (varfloatcomplex *)VR.fElem, &dim, (varfloatcomplex *)&work[0], &lwork, &rwork[0], &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = eigen[i] / beta[i];
        }
    }
    for(int i = 0 ; i < dim ; i ++){
        for( int iV = 0 ; iV < dim ; iV++ ){
            eigenvectors(iV,i) = VR(iV,i);
        }
    }
    
    return 1;
}


template<>
int
TPZFMatrix<complex<float> >::SolveGeneralisedEigenProblem(TPZFMatrix<complex<float> > &B , TPZVec <complex<double> > &eigenvalues)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< complex<float> > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    TPZVec<complex<float> > eigen(dim,0.);
    
    TPZVec<complex<float> > beta(dim);
    
    TPZFMatrix<complex<float> > temp(*this), tempB(B);
    TPZVec<complex<float> > work(lwork);
    TPZVec<float> rwork( 8 * dim );
    
    cggev_(jobvl, jobvr, &dim, (varfloatcomplex *)temp.fElem, &dim , (varfloatcomplex *)tempB.fElem, &dim , (varfloatcomplex *)&eigen[0], (varfloatcomplex *)&beta[0]  , (varfloatcomplex *)VL.fElem, &dim , (varfloatcomplex *)VR.fElem, &dim, (varfloatcomplex *)&work[0], &lwork, &rwork[0], &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = eigen[i] / beta[i];
        }
    }
    return 1;

}

template<>
int
TPZFMatrix<complex<double> >::SolveGeneralisedEigenProblem(TPZFMatrix<complex<double> > &B , TPZVec <complex<double> > &eigenvalues, TPZFMatrix < complex<double> > &eigenvectors)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< complex<double> > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    TPZVec<complex<double> > eigen(dim,0.);
    
    TPZVec<complex<double> > beta(dim);
    
    TPZFMatrix<complex<double> > temp(*this), tempB(B);
    TPZVec<complex<double> > work(lwork);
    TPZVec<double> rwork( 8 * dim );

    zggev_(jobvl, jobvr, &dim, (vardoublecomplex *)temp.fElem, &dim , (vardoublecomplex *)tempB.fElem, &dim , (vardoublecomplex *)&eigen[0], (vardoublecomplex *)&beta[0]  , (vardoublecomplex *)VL.fElem, &dim , (vardoublecomplex *)VR.fElem, &dim, (vardoublecomplex *)&work[0], &lwork, &rwork[0], &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = eigen[i] / beta[i];
        }
    }
    for(int i = 0 ; i < dim ; i ++){
        for( int iV = 0 ; iV < dim ; iV++ ){
            eigenvectors(iV,i) = VR(iV,i);
        }
    }
    
    return 1;
}


template<>
int
TPZFMatrix<complex<double> >::SolveGeneralisedEigenProblem(TPZFMatrix<complex<double> > &B , TPZVec <complex<double> > &eigenvalues)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< complex<double> > VL(Rows(),Cols()),VR(Rows(),Cols());
    int dim = Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    TPZVec<complex<double> > eigen(dim,0.);
    
    TPZVec<complex<double> > beta(dim);
    
    TPZFMatrix<complex<double> > temp(*this), tempB(B);
    TPZVec<complex<double> > work(lwork);
    TPZVec<double> rwork( 8 * dim );
    
    zggev_(jobvl, jobvr, &dim, (vardoublecomplex *)temp.fElem, &dim , (vardoublecomplex *)tempB.fElem, &dim , (vardoublecomplex *)&eigen[0], (vardoublecomplex *)&beta[0]  , (vardoublecomplex *)VL.fElem, &dim , (vardoublecomplex *)VR.fElem, &dim, (vardoublecomplex *)&work[0], &lwork, &rwork[0], &info);
    
    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = eigen[i] / beta[i];
        }
    }
    
    return 1;
    
}

#endif // USING_LAPACK


#ifdef _AUTODIFF
/** @brief Returns the norm of the matrix A */
template<>
TFad<6,REAL> Norm(const TPZFMatrix<TFad<6,REAL> > &A)
{
    DebugStop();
    TFad<6,REAL> res;
    return res;
}
template<>
Fad<REAL> Norm(const TPZFMatrix<Fad<REAL> > &A)
{
    DebugStop();
    Fad<REAL> res;
    return res;
}



#endif

#include <complex>

template class TPZFMatrix<int >;
template class TPZFMatrix<int64_t >;
template class TPZFMatrix<float >;
template class TPZFMatrix<double >;
template class TPZFMatrix<long double>;

template class TPZFMatrix< std::complex<float> >;
template class TPZFMatrix< std::complex<double> >;
template class TPZFMatrix< std::complex<long double> >;

template class TPZFMatrix<TPZFlopCounter>;

template class TPZRestoreClass< TPZFMatrix<int> >;
template class TPZRestoreClass< TPZFMatrix<int64_t> >;
template class TPZRestoreClass< TPZFMatrix<double> >;
template class TPZRestoreClass< TPZFMatrix<float> >;
template class TPZRestoreClass< TPZFMatrix<long double> >;

template class TPZRestoreClass< TPZFMatrix<std::complex<float> > >;
template class TPZRestoreClass< TPZFMatrix<std::complex<double> > >;
template class TPZRestoreClass< TPZFMatrix<std::complex<long double> > >;
template class TPZRestoreClass< TPZFMatrix<TPZFlopCounter > >;

#ifdef _AUTODIFF
#include "fad.h"
template class TPZFMatrix<TFad<6,REAL> >;
template class TPZFMatrix<Fad<double> >;
template class TPZFMatrix<Fad<float> >;
template class TPZFMatrix<Fad<long double> >;

template class TPZRestoreClass<TPZFMatrix<TFad<6,REAL> >>;
template class TPZRestoreClass<TPZFMatrix<Fad<double> >>;
template class TPZRestoreClass<TPZFMatrix<Fad<float> >>;
template class TPZRestoreClass<TPZFMatrix<Fad<long double> >>;
#endif

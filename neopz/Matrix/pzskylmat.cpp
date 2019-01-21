//EBORIN: Uncomment the following lines to enable matrices do be dumped before decompose/subst
//#define DUMP_BEFORE_DECOMPOSE
//#define DUMP_BEFORE_SUBST

#include "arglib.h"


#ifdef USING_NEW_SKYLMAT

/**
 @file
 @brief Contains the implementation of the TPZSkylMatrix methods.
 */

#include <math.h>
#include <stdlib.h>

#ifdef BLAS
extern "C" {
#include <cblas.h>
}
#endif

#include "pzfmatrix.h"
#include "pzskylmat.h"

#include <sstream>
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzskylmatrix"));
#endif

using namespace std;

/**************************** PUBLIC ****************************/
template<class TVar>
TPZSkylMatrix<TVar>::TPZSkylMatrix(const int64_t dim ) :
TPZRegisterClassId(&TPZSkylMatrix::ClassId),TPZMatrix<TVar>( dim, dim ), fElem(dim+1), fStorage(0)
{
    // Initializes the diagonal with zeros.
    fElem.Fill(0);
}

template<class TVar>
TPZSkylMatrix<TVar>::TPZSkylMatrix(const int64_t dim, const TPZVec<int64_t> &skyline ) :
TPZRegisterClassId(&TPZSkylMatrix::ClassId),TPZMatrix<TVar>( dim, dim ), fElem(dim+1), fStorage(0)
{
    fElem.Fill(0);
    InitializeElem(skyline,fStorage,fElem);
}

template<class TVar>
void TPZSkylMatrix<TVar>::SetSkyline(const TPZVec<int64_t> &skyline)
{
    fElem.Fill(0);
    InitializeElem(skyline,fStorage,fElem);
}

template<class TVar>
void TPZSkylMatrix<TVar>::AddSameStruct(TPZSkylMatrix<TVar> &B, double k){
#ifdef PZDEBUG
    {
        int64_t size = this->fElem.NElements();
        if(size != B.fElem.NElements()){
            PZError << "\nFATAL ERROR at " << __PRETTY_FUNCTION__ << "\n";
            PZError.flush();
            DebugStop();
        }
        for(int64_t i = 0; i < size; i++){
            if((this->fElem[i]-this->fElem[0]) != (B.fElem[i]-B.fElem[0])){
                PZError << "\nFATAL ERROR at " << __PRETTY_FUNCTION__ << "\n";
                PZError.flush();
                DebugStop();
            }
        }
    }
#endif
    
    const int64_t n = this->fStorage.NElements();
    for(int64_t i = 0; i < n; i++)
        this->fStorage[i] += TVar(k) * B.fStorage[i];
}

template<class TVar>
void TPZSkylMatrix<TVar>::UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat)
{
    TPZMatrix<TVar> *matrix = mat.operator->();
    TPZSkylMatrix<TVar> *skylmat = dynamic_cast<TPZSkylMatrix<TVar> *>(matrix);
    if (!skylmat) {
        DebugStop();
    }
    if (fStorage.NElements() != skylmat->fStorage.NElements()) {
        DebugStop();
    }
    memcpy(&fStorage[0], &(skylmat->fStorage[0]) , fStorage.NElements()*sizeof(TVar));
    this->fDecomposed = skylmat->fDecomposed;
    this->fDefPositive = skylmat->fDefPositive;
}

template<class TVar>
int
TPZSkylMatrix<TVar>::PutVal(const int64_t r,const int64_t c,const TVar & value )
{
    // Adjusting row and col to work with superior triangular
    if (r > c) 
        return PutVal(c, r, value);
    
    // Column array index
    int64_t sz = Size(c);
    if ( (c-r) >= sz) {
        if (IsZero(value)) {
            return 1; // Return OK
        }
        else {
            cerr << "TPZSkylMatrix::PutVal Size" << sz;
            TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Index out of range");
            return 0;
        }
    } 
    
    // See pzskylmat.h for more info on how to compute the index.
    int64_t index = r + sz - 1 - c;
    fElem[c][index] = value;
    return 1;
}

template<class TVar>
const TVar &
TPZSkylMatrix<TVar>::GetVal(const int64_t r,const int64_t c ) const
{
    if (r > c) 
        return GetVal(c,r);
    
#ifdef PZDEBUG
    unsigned dim = this->Dim();
    
    if(r >= dim || c >= dim  || r < 0 || c < 0) {
        cerr << "TPZSkylMatrix::GetVal index out of range row = " << r
        << " col = " << c << endl;
        return this->gZero;
    }
#endif
    
    int64_t sz = Size(c);
    if ( (c-r) < sz ) {
        int64_t index = r + sz - 1 - c;
        return (fElem[c][index]);
    }
    else {
#ifdef PZDEBUG
        if(this->gZero != TVar(0.)) {
            cerr << "TPZSkylMatrix gZero = " << this->gZero << endl;
            DebugStop();
        }
#endif
        return(this->gZero );
    }
}

template<class TVar>
TVar &
TPZSkylMatrix<TVar>::operator()(int64_t ri, int64_t ci)
{
    int64_t r = ri;
    int64_t c = ci;
    if (ri > ci) {
        r = ci;
        c = ri;
    }
    
    int64_t sz = Size(c);
    
    if ( (c-r)  >= sz ) {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Index out of range");
        DebugStop();
    }
    
    int64_t index = r + sz - 1 - c;
    
    return fElem[c][index];
}

template<class TVar>
void TPZSkylMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                                  const TVar alpha,const TVar beta ,const int opt) const
{
    // Computes z = beta * y + alpha * opt(this)*x
    //          z and x cannot overlap in memory
	
    if (this->fDecomposed != ENoDecompose) {
        //TODO: EBORIN: Figure out why the next line was commented. 
        //		DebugStop();
    }
    
    if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows()) {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__," <matrix with incompatible dimensions>" );
    }
    
    if(z.Rows() != x.Rows() || z.Cols() != x.Cols()) 
        z.Redim(x.Rows(),x.Cols());
    
    if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
        cerr << "x.Cols = " << x.Cols() << " y.Cols()"<< y.Cols() << " z.Cols() " << z.Cols() 
        << " x.Rows() " << x.Rows() << " y.Rows() "<< y.Rows() << " z.Rows() "<< z.Rows() << endl;
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__," incompatible dimensions\n");
    }
    
    this->PrepareZ(y,z,beta,opt);
    
    int64_t rows = this->Rows();
    int64_t xcols = x.Cols();
    int64_t ic, r;
    for (ic = 0; ic < xcols; ic++) {
        for( r = 0 ; r < rows ; r++ ) {
            int64_t offset = Size(r);
            TVar val = 0.;
            const TVar *p = &x.g((r-offset+1),ic);
            TVar *diag = fElem[r];
            TVar *diaglast = fElem[r+1]-1;
            while( diag < diaglast ) {
                val += *diag++ * *p;
                p ++;
            }
            if( diag == diaglast ) 
                val += *diag * *p;
            
            z(r,ic) += val*alpha;
            
            TVar *zp = &z((r-offset+1),ic);
            val = x.g(r,ic);
            diag = fElem[r];
            while( diag < diaglast ) {
                *zp += alpha * *diag++ * val;
                zp ++;
            }
        }
    }
}

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator=(const TPZSkylMatrix<TVar> &A )
{
    Clear();
    Copy(A);
    return(*this);
}

template<class TVar>
TPZSkylMatrix<TVar>
TPZSkylMatrix<TVar>::operator+(const TPZSkylMatrix<TVar> &A) const
{
    if ( this->Dim() != A.Dim() )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"<incompatible dimensions>" );
    
    TPZVec<int64_t> skylinesize(this->Dim());
    ComputeMaxSkyline(*this,A,skylinesize);
    TPZSkylMatrix res( this->fRow, skylinesize );
    
    TVar *elemMax;
    TVar *elemMin;
    int64_t  sizeMax;
    int64_t  sizeMin;
    
    for (int64_t col = 0; col < this->Dim(); col++) {
        // Select the size and the elements of the smallest and highest column.
        if ( Size(col) > A.Size(col) ) {
            sizeMax = Size(col);
            elemMax = fElem[col];
            sizeMin = A.Size(col);
            elemMin = A.fElem[col];
        }
        else {
            sizeMax = A.Size(col);
            elemMax = A.fElem[col];
            sizeMin = Size(col);
            elemMin = fElem[col];
        }
        
        TVar *dest = res.fElem[col];
        int64_t i = 0;
        for ( ; i < sizeMin; i++ )
            *dest++ = (*elemMax++) + (*elemMin++);
        for ( ; i < sizeMax; i++ )
            *dest++ = *elemMax++;
    }
    
    return( res );
}

template<class TVar>
TPZSkylMatrix<TVar>
TPZSkylMatrix<TVar>::operator-(const TPZSkylMatrix<TVar> &A ) const
{
    if ( this->Dim() != A.Dim() )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator-( TPZSkylMatrix ) <incompatible dimensions>" );
    
    TPZVec<int64_t> skylinesize(this->Dim());
    ComputeMaxSkyline(*this,A,skylinesize);
    TPZSkylMatrix<TVar> res( this->fRow, skylinesize );
    
    for ( int64_t col = 0; col < this->fRow; col++ ) {
        int64_t  sizeThis  = Size(col);
        TVar *elemThis = fElem[col];
        int64_t  sizeA     = A.Size(col);
        TVar *elemA    = A.fElem[col];
        
        TVar *dest = res.fElem[col];
        int64_t i;
        for ( i = 0; (i < sizeThis) && (i < sizeA); i++ ) 
            *dest++ = (*elemThis++) - (*elemA++);
        
        if ( i == sizeA ) 
            for ( ; i < sizeThis; i++ ) *dest++ = *elemThis++;
        else 
            for ( ; i < sizeA; i++ ) *dest++ = -(*elemA++);
    }
    
    return( res );
}

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator+=(const TPZSkylMatrix<TVar> &A )
{
    if ( this->Dim() != A.Dim() )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator+=( TPZSkylMatrix ) <incompatible dimensions>" );
	
    TPZSkylMatrix res((*this)+A);
    *this = res;
    return *this;
}

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator-=(const TPZSkylMatrix<TVar> &A )
{
    if ( this->Dim() != A.Dim() )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator-=( TPZSkylMatrix ) <incompatible dimensions>" );
    
    TPZSkylMatrix res(*this-A);
    *this = res;
    return *this;
}

template<class TVar>
TPZSkylMatrix<TVar>
TPZSkylMatrix<TVar>::operator*(const TVar value ) const
{
    TPZSkylMatrix res( *this );
    int64_t nelem = res.fStorage.NElements();
    TVar *elemRes = res.fElem[0];
    
    for (int64_t i=0; i<nelem; i++) {
        *elemRes++ *= value;
    }
    
    return( res );
}

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator*=(const TVar value )
{
    if ( IsZero( value ) ) {
        Zero();
        return( *this );
    }
    
    int64_t nelem = fStorage.NElements();
    TVar *elemRes = fElem[0];
    
    for (int64_t i=0; i<nelem; i++) {
        *elemRes++ *= value;
    }
	
    this->fDecomposed = 0;
    return( *this );
}

/*** Resize ***/
// Change the dimensions of the matrix, but keep the old values. New elements are zeroed.
template<class TVar>
int TPZSkylMatrix<TVar>::Resize( int64_t newDim ,int64_t )
{
    if ( newDim == this->Dim() )
        return( 1 );
    
    fElem.Resize(newDim+1);
    
    // Copy the elements into the new matrix.
    int64_t min = MIN( newDim, this->Dim() );
    for (int64_t i = min+1; i <= newDim; i++ )
        fElem[i] = fElem[i-1];
    
    // Set the remaining positions with zero.
    fStorage.Resize(fElem[newDim]-fElem[0]);
    this->fRow = this->fCol = newDim;
    this->fDecomposed = 0;
    return( 1 );
}

template<class TVar>
int TPZSkylMatrix<TVar>::Redim( int64_t newDim , int64_t)
{
    if ( newDim == this->Dim() ) {
        Zero();
        return( 1 );
    }
	
    Clear();
    fElem.Resize(newDim);
    fElem.Fill(0);
    this->fRow = this->fCol = newDim;
    this->fDecomposed = 0;
    return( 1 );
}

template<class TVar>
int TPZSkylMatrix<TVar>::Zero()
{
    fStorage.Fill(0.);
    this->fDecomposed = 0;
    this->fDefPositive = 0;
    return( 1 );
}

template<class TVar>
int64_t TPZSkylMatrix<TVar>::NumElements(const TPZVec<int64_t> &skyline)
{
    int64_t dim = skyline.NElements();
    
    int64_t nelem=0;
    for(int64_t i=0; i<dim; i++) {
        nelem += i-skyline[i]+1;
    }
    
    return nelem;
}

template<class TVar>
void TPZSkylMatrix<TVar>::InitializeElem(const TPZVec<int64_t> &skyline, TPZVec<TVar> &storage, TPZVec<TVar *> &point)
{
    int64_t dim = skyline.NElements();
    int64_t nel = NumElements(skyline);
    storage.Resize(nel);
    storage.Fill(0.);
    point.Resize(dim+1);
    if(dim) {
        point[0] = &storage[0];
        point[dim] = &storage[0]+nel;
    } 
    else {
        point[0] = 0;
    }
    
    for(int64_t i=1; i<dim+1; i++) { 
        point[i] = point[i-1] + (i-1)-skyline[i-1]+1;
    }
}

template<class TVar>
void TPZSkylMatrix<TVar>::ComputeMaxSkyline(const TPZSkylMatrix<TVar> &first, 
                                            const TPZSkylMatrix<TVar> &second, TPZVec<int64_t> &res)
{
    if (first.Rows() != second.Rows()) {
        cout<<"ComputeMaxSkyline : incompatible dimension";
        return;
    }
    int64_t dim = first.Rows();
    res.Resize(dim);
    
    for(int64_t i=0; i<dim; i++) {
        int64_t sz_first = first.Size(i);
        int64_t sz_secon = second.Size(i);
        if (sz_first > sz_secon)
            res[i] = i-(sz_first-1);
        else
            res[i] = i-(sz_secon-1);
    }
}

template<class TVar>
void TPZSkylMatrix<TVar>::SolveSOR(int64_t & numiterations,const TPZFMatrix<TVar> &F,
                                   TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &scratch,const REAL overrelax,
                                   REAL &tol,const int FromCurrent,const int direction)  
{  
    if(residual == &F) {
        cout << "TPZMatrix::SolveSOR called with residual and F equal, no solution\n";
        return;
    }
    REAL res = 2.*tol+1.;
    if(residual) res = Norm(*residual);
    if(!FromCurrent) {
        result.Zero();
    }
    TVar over = overrelax;
    int64_t r = this->Dim();
    int64_t c = F.Cols();
    int64_t i,ifirst = 0, ilast = r, iinc = 1;
    if(direction == -1) {
        ifirst = r-1;
        ilast = 0;
        iinc = -1;
    }
    int64_t it;
    for(it=0; it<numiterations && res > tol; it++) {
        res = 0.;
        scratch = F;
        for(int64_t ic=0; ic<c; ic++) {
            if(direction == 1) {
                //
                // compute the upper triangular part first and put into the scractch vector
                //
                for(i=ifirst; i!=ilast; i+= iinc) {
                    //TPZColuna *mydiag = &fDiag[i];
                    int64_t offset = Size(i);
                    TVar val;
                    TVar *diaglast = (fElem[i+1]-1);
                    TVar *scratchp = &scratch(i-offset+1,ic);
                    val = result(i,ic);
                    TVar *p = fElem[i];
                    int64_t lastid = diaglast-p;
                    int64_t id;
                    for(id=0; id<=lastid; id++) 
                        *scratchp++ -= *p++ * val;
                    /* codeguard fix
                     while( diag >= diaglast ) *scratchp++ -= *diag-- * val;
                     */
                }
                //
                // perform the SOR operation
                //
                for(i=ifirst; i!=ilast; i+= iinc) {
                    //TPZColuna *mydiag = &fDiag[i];
                    int64_t offset = Size(i);
                    TVar val = scratch(i,ic);
                    TVar *p = &result(i-offset+1,ic);
                    TVar *diag = fElem[i];
                    TVar *diaglast = (fElem[i+1]-1);
                    while( diag < diaglast ) 
                        val -= *diag++ * *p++;
                    res += abs(val*val);
                    result(i,ic) += val*over/ (*diag);
                }
            } else {
                //
                // the direction is upward
                //
                // put the lower triangular part of the multiplication into the scratch vector
                //
                for(i=ifirst; i!=ilast; i+= iinc) {
                    //TPZColuna *mydiag = &fDiag[i];
                    int64_t offset = Size(i);
                    TVar val = scratch(i,ic);
                    TVar *p = &result(i-offset+1,ic);
                    TVar *diag = fElem[i];
                    TVar *diaglast = (fElem[i+1]-1);
                    while( diag < diaglast ) 
                        val -= *diag++ * *p++;
                    //					res += val*val;
                    scratch(i,ic) = val;
                }
                //
                // perform the SOR operation
                //
                for(i=ifirst; i!=ilast; i+= iinc) {
                    //TPZColuna *mydiag = &fDiag[i];
                    int64_t offset = Size(i);
                    //	TVar val = scratch(i,ic);
                    TVar *diaglast = (fElem[i+1]-1);
                    TVar *scratchp = &scratch(i-offset+1,ic);
                    //val= result(i,ic);
                    TVar val = scratch(i,ic);
                    val -= *diaglast * result(i,ic);
                    res += abs(val*val);
                    val = over * val / *diaglast;
                    result(i,ic) += val;
                    val = result(i,ic);
                    TVar *diag = fElem[i];
                    while( diag < diaglast ) 
                        *scratchp++ -= *diag++ * val;
                }
            }
        }
        res = sqrt(res);
    }
    if(residual) {
        this->Residual(result,F,*residual);
    }
    numiterations = it;
    tol = res;
}

template<class TVar>
void TPZSkylMatrix<TVar>::AddKel(TPZFMatrix<TVar>&elmat,
                                 TPZVec<int64_t> &source,
                                 TPZVec<int64_t> &destination)
{
    
    std::cout << "Passei por aqui\n";
    int64_t nelem = source.NElements();
    int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
    for(icoef=0; icoef<nelem; icoef++) {
        ieq = destination[icoef];
        ieqs = source[icoef];
        for(jcoef=icoef; jcoef<nelem; jcoef++) {
            jeq = destination[jcoef];
            jeqs = source[jcoef];
            int64_t row(ieq), col(jeq);
            // invertendo linha-coluna para triangular superior
            if (row > col)
                this->Swap(&row, &col);
#ifdef PZDEBUG
            // checando limites
            if(row >= this->Dim() || col >= this->Dim()) {
                cout << "TPZSkylMatrix::GetVal index out of range row = " <<
                row << " col = " << col << endl;
                DebugStop();
            }
#endif
            int64_t sz = Size(col);
            
            int64_t index = row + sz - 1 - col;
#pragma omp atomic
            fElem[col][index] += elmat(ieqs,jeqs);

            
        }
    }
}

template<>
void TPZSkylMatrix<double>::AddKel(TPZFMatrix<double>&elmat,
                                   TPZVec<int64_t> &source,
                                   TPZVec<int64_t> &destination)
{
    int64_t nelem = source.NElements();
    int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
    for(icoef=0; icoef<nelem; icoef++) {
        ieq = destination[icoef];
        ieqs = source[icoef];
        for(jcoef=icoef; jcoef<nelem; jcoef++) {
            jeq = destination[jcoef];
            jeqs = source[jcoef];
            int64_t row(ieq), col(jeq);
            // invertendo linha-coluna para triangular superior
            if (row > col)
                this->Swap(&row, &col);
#ifdef PZDEBUG
            // checando limites
            if(row >= this->Dim() || col >= this->Dim()) {
                cout << "TPZSkylMatrix::GetVal index out of range row = " <<
                row << " col = " << col << endl;
                DebugStop();
            }
#endif
            int64_t sz = Size(col);

            int64_t index = row + sz - 1 - col;
#pragma omp atomic
            fElem[col][index] += elmat(ieqs,jeqs);
        }
    }
}


template<>
void TPZSkylMatrix<float>::AddKel(TPZFMatrix<float>&elmat,
                                  TPZVec<int64_t> &source,
                                  TPZVec<int64_t> &destination)
{
    int64_t nelem = source.NElements();
    int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
    for(icoef=0; icoef<nelem; icoef++) {
        ieq = destination[icoef];
        ieqs = source[icoef];
        for(jcoef=icoef; jcoef<nelem; jcoef++) {
            jeq = destination[jcoef];
            jeqs = source[jcoef];
            int64_t row(ieq), col(jeq);
            // invertendo linha-coluna para triangular superior
            if (row > col)
                this->Swap(&row, &col);
#ifdef PZDEBUG
            // checando limites
            if(row >= this->Dim() || col >= this->Dim()) {
                cout << "TPZSkylMatrix::GetVal index out of range row = " <<
                row << " col = " << col << endl;
                DebugStop();
            }
#endif
            // executando contribuição
            int64_t sz = Size(col);

            int64_t index = row + sz - 1 - col;
#pragma omp atomic
            fElem[col][index] += elmat(ieqs,jeqs);

            
        }
    }
}

template<>
int TPZSkylMatrix<std::complex<float> >::Decompose_Cholesky(std::list<int64_t> &singular)
{
    DebugStop();
    return -1;
}

template<>
int TPZSkylMatrix<std::complex<double> >::Decompose_Cholesky(std::list<int64_t> &singular)
{
    DebugStop();
    return -1;
}

template<>
int TPZSkylMatrix<std::complex<long double> >::Decompose_Cholesky(std::list<int64_t> &singular)
{
    DebugStop();
    return -1;
}

//EBORIN: Define these if you want to use the experimental version.
//#define DECOMPOSE_CHOLESKY_OPT2
//#define SKYLMATRIX_GETVAL_OPT1

/**************************/
/*** Decompose Cholesky ***/
template<class TVar>
int TPZSkylMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular)
{
    if(this->fDecomposed == ECholesky) 
        return 1;
    
    if (this->fDecomposed )  
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
    
#ifdef DUMP_BEFORE_DECOMPOSE
    dump_matrix(this, "TPZSkylMatrix::Decompose_Cholesky()");
#endif
    
    singular.clear();
    
    TVar pivot;
    int64_t dimension = this->Dim();
    
    for ( int64_t k = 0; k < dimension; k++ ) {
        
        if ( Size(k) == 0 )	
            return 0;
        
        // sum = Sum( A(k,p) * A(k,p) ), p = 1, ..., k-1.
        TVar sum = 0.0;
        TVar *elem_k = fElem[k];
        TVar *end_k  = fElem[k+1]-1;
        
#pragma clang loop vectorize_width(2)
        for ( ; elem_k < end_k; elem_k++ )
            sum += (*elem_k) * (*elem_k);
        
        elem_k = fElem[k+1]-1;    // Diagonal element.
        TVar* first_k = fElem[k]; // First element at row k
        TVar* last_k = elem_k;    // Last element at row k (diagonal element)
        int64_t k_sz = last_k - first_k;
        
        // Faz A(k,k) = sqrt( A(k,k) - sum ).
        pivot = *elem_k - sum;
        
        //EBORIN: FIXME: Shouldn't this be IsZero(pivot)???
        if (pivot < TVar(1.e-10)) {
            singular.push_back(k);
            pivot = 1.;
        }
        
        pivot = sqrt( pivot );
        *elem_k = pivot; 
        
        // Loop para i = k+1 ... Dim().
        int64_t i=k+1;
        
        for ( int64_t j = 2; i<dimension; j++,i++ ) {
            
            TVar* upd_elem = fElem[i+1] - j; // Element to update (row i, column k)
            TVar* first_i = fElem[i];        // First element at row i (also the first element at rowi)
            TVar* last_i = upd_elem;         
            
            if (first_i < last_i) {
                
                // Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1,..,k-1.
                sum = 0.0;
                
                int64_t min_sz = last_i - first_i;
                if (min_sz > k_sz) 
                    min_sz = k_sz;
                
                TVar* ip = last_i - min_sz;
                TVar* kp = last_k - min_sz;
                
#pragma clang loop vectorize_width(2)
                for(unsigned l=0; l<min_sz; l++)
                    sum += (*ip++) * (*kp++);
                
                // A(i,k) = (A(i,k) - sum) / A(k,k)
                *upd_elem = (*upd_elem -sum) / pivot;
            }
            else if (first_i == last_i) {
                // sum == 0.0
                *upd_elem = *upd_elem / pivot;
            }
            // else { no elements to update }
        }
    }
    
    if(this->Rows() && (GetVal(this->Rows()-1,this->Rows()-1)) < TVar(1.e-15)) {
        singular.push_back(this->Rows()-1);
        PutVal(this->Rows()-1,this->Rows()-1,1.);
    }
    
    this->fDecomposed  = ECholesky;
    this->fDefPositive = 1;
    return( 1 );
}

template<>
int TPZSkylMatrix<std::complex<float> >::Decompose_Cholesky()
{
    DebugStop();
    return -1;
}
template<>
int TPZSkylMatrix<std::complex<double> >::Decompose_Cholesky()
{
    DebugStop();
    return -1;
}
template<>
int TPZSkylMatrix<std::complex<long double> >::Decompose_Cholesky()
{
    DebugStop();
    return -1;
}

clarg::argBool clk_mig("-skl_chk_pm", "Migrate Skyline pages before Cholesky Decomposition");
clarg::argBool clk_rea("-skl_chk_rea", "Reallocate Skyline data before Cholesky Decomposition");

template<class TVar>
int TPZSkylMatrix<TVar>::Decompose_Cholesky()
{
    if(this->fDecomposed == ECholesky) 
        return 1;
    
    if (this->fDecomposed )  
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	
#ifdef DUMP_BEFORE_DECOMPOSE
    dump_matrix(this, "TPZSkylMatrix::Decompose_Cholesky()");
#endif
    
    if (clk_mig.was_set())
        MigratePages();
    if (clk_rea.was_set())
        ReallocForNuma();
    
    TVar pivot;
    TVar minpivot = 10000.;
    int64_t dimension = this->Dim();
    
    for ( int64_t k = 0; k < dimension; k++ ) {
        
        if ( Size(k) == 0 )	
            return 0;
        
        // sum = Sum( A(k,p) * A(k,p) ), p = 1, ..., k-1.
        TVar sum = 0.0;
        TVar *elem_k = fElem[k];
        TVar *end_k  = fElem[k+1]-1; 
#pragma clang loop vectorize_width(2)
        for ( ; elem_k < end_k; elem_k++ )
            sum += (*elem_k) * (*elem_k);
        
        elem_k = fElem[k+1]-1;    // Diagonal element.
        TVar* first_k = fElem[k]; // First element at row k
        TVar* last_k = elem_k;    // Last element at row k (diagonal element)
        int64_t k_sz = last_k - first_k;
        
        // Faz A(k,k) = sqrt( A(k,k) - sum ).
        pivot = *elem_k - sum;
        minpivot = minpivot < pivot ? minpivot : pivot;
        if ( pivot < TVar(0.) || IsZero(pivot) ) {
            cout << "TPZSkylMatrix::DecomposeCholesky! Matrix is not positive definite" << pivot << endl;
            return 0;
        }
        
        pivot = sqrt( pivot );
        *elem_k = pivot; 
        
        // Loop para i = k+1 ... Dim().
        int64_t i=k+1;
        
        for ( int64_t j = 2; i<dimension; j++,i++ ) {
            
            TVar* upd_elem = fElem[i+1] - j; // Element to update (row i, column k)
            TVar* first_i = fElem[i];        // First element at row i (also the first element at rowi)
            TVar* last_i = upd_elem;         
            
            if (first_i < last_i) {
                
                // Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1,..,k-1.
                sum = 0.0;
                
                int64_t min_sz = last_i - first_i;
                if (min_sz > k_sz) 
                    min_sz = k_sz;
                
                TVar* ip = last_i - min_sz;
                TVar* kp = last_k - min_sz;
                
#pragma clang loop vectorize_width(2)
                for(unsigned l=0; l<min_sz; l++)
                    sum += (*ip++) * (*kp++);
                
                // A(i,k) = (A(i,k) - sum) / A(k,k)
                *upd_elem = (*upd_elem -sum) / pivot;
            }
            else if (first_i == last_i) {
                // sum == 0.0
                *upd_elem = *upd_elem / pivot;
            }
            // else { no elements to update }
        }
    }
    
    this->fDecomposed  = ECholesky;
    this->fDefPositive = 1;
#ifdef PZDEBUG
    std::cout << __PRETTY_FUNCTION__ << " minpivot " << minpivot << std::endl;
#endif
    return( 1 );
}

template<>
int TPZSkylMatrix<std::complex<float> >::Decompose_Cholesky_blk(int64_t blk_sz)
{
    DebugStop();
    return -1;
}

template<>
int TPZSkylMatrix<std::complex<double> >::Decompose_Cholesky_blk(int64_t blk_sz)
{
    DebugStop();
    return -1;
}

template<>
int TPZSkylMatrix<std::complex<long double> >::Decompose_Cholesky_blk(int64_t blk_sz)
{
    DebugStop();
    return -1;
}

template<class TVar>
int TPZSkylMatrix<TVar>::Decompose_Cholesky_blk(int64_t blk_sz)
{
    DebugStop();
    return -1;
}

template<class TVar>
int
TPZSkylMatrix<TVar>::Decompose_LDLt(std::list<int64_t> &singular)
{
    if( this->fDecomposed == ELDLt) 
        return 1;
    
    if( this->fDecomposed ) {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, " Decompose_LDLt <Matrix already Decomposed with different decomposition>" );
    }
    
#ifdef DUMP_BEFORE_DECOMPOSE
    dump_matrix(this, "TPZSkylMatrix::Decompose_LDLt(singular)");
#endif
    
    singular.clear();
    
    // Third try
    TVar *elj,*ell;
    int64_t j,l,minj,minl,minrow,dimension = this->Dim();
    TVar sum;
    j = 1;
    while(j < dimension) {
        
        minj = j-Size(j)+1;
        l = minj;
        
        while(l <= j) {
            minl = l-Size(l)+1;
            minrow = (minj<minl)? minl:minj;
            int64_t k = minrow;
            
            // DiagkPtr = fDiag+minrow;
            elj = fElem[j+1] - (j+1) + minrow ; // fElem[j]+j-minrow;
            ell = fElem[l+1] - (l+1) + minrow ; // fElem[l]+l-minrow;
            sum = 0.;
            
            while(k < l) {
                sum += *elj++  *  *ell++ * *(fElem[k+1]-1);
                k++;
            }
            
            *elj -= sum;
            
            if(ell != elj) {
                *elj /= *ell;
            }
            else if(IsZero(*elj)) {
                singular.push_back(l);
                *elj = 1.;
            }
            l++;
        }
        j++;
    }
    
    if(this->Rows() && IsZero(GetVal(this->Rows()-1,this->Rows()-1))) {
        singular.push_back(this->Rows()-1);
        PutVal(this->Rows()-1,this->Rows()-1,1.);
    }
    
    this->fDecomposed  = ELDLt;
    this->fDefPositive = 0;
    //if(Dim() > 100) cout << endl;
    return( 1 );
}

template<class TVar>
int TPZSkylMatrix<TVar>::Decompose_LDLt()
{
    if( this->fDecomposed == ELDLt) 
        return 1;
    if (  this->fDecomposed )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, " Decompose_LDLt <Matrix already Decomposed with different decomposition>" );
    
#ifdef DUMP_BEFORE_DECOMPOSE
    dump_matrix(this, "TPZSkylMatrix::Decompose_LDLt()");
#endif
    
    // Third try
    TVar *elj,*ell;
    int64_t j,l,minj,minl,minrow,dimension = this->Dim();
    TPZVec<TVar> diag(dimension);
    for(j=0; j<dimension; j++)
        diag[j] = *(fElem[j+1]-1);
    
    //std::cout << "TPZSkylMatrix<TVar>::Decompose_LDLt: dimension = " << dimension  << std::endl;
    
    TVar sum;
    j = 1;
    while(j < dimension) {
        /*    if(!(j%100) && Dim() > 100) {
         cout <<  j << ' ';
         cout.flush();
         }
         if(!(j%1000)) cout << endl;*/
        minj = j-Size(j)+1;
        l = minj;
        while(l <= j) {
            minl = l-Size(l)+1;
            minrow = (minj<minl)? minl:minj;
            int64_t k = minrow;
            //			DiagkPtr = fDiag+minrow;
            elj = fElem[j]+minrow-minj; // fElem[j]+minrow-(j-Size(j)+1); // Pointer to A[minrow,j]
            ell = fElem[l]+minrow-minl; // fElem[l]+minrow-(l-Size(l)+1); // Pointer to A[minrow,l]
            TVar *diagptr = &diag[k];
            sum = 0.;
            while(k < l) {
                sum += *elj++ * *ell++ * *diagptr++;
                k++;
            }
            *elj -= sum;
            if(ell != elj) *elj /= *ell;
            else if(IsZero(*elj)) {
#ifdef LOG4CXX
                std::stringstream sout;
                sout << "col = " << j << " diagonal " << *elj;
                LOGPZ_DEBUG(logger,sout.str())
#endif
                
                *diagptr = *elj;
                cout << "TPZSkylMatrix pivot = " << *elj << endl;
                cout << "TPZSkylMatrix::DecomposeLDLt zero pivot\n";
                cout << "j = " << j << " l = " << l << endl;
            }
            else {
                *diagptr = *elj;
            }
            l++;
        }
        j++;
    }
    this->fDecomposed  = ELDLt;
    this->fDefPositive = 0;
    //if(Dim() > 100) cout << endl;
    return( 1 );
}

#ifdef USING_MKL
/** Intel Math Kernel Library */
#include <mkl.h> 

/*********************/
/*** Subst Forward ***/
//
//  Faz Ax = b, onde A e' triangular superior.
//	Utilizando MKL
template<>
int
TPZSkylMatrix<double>::Subst_Forward( TPZFMatrix<double> *B ) const
{
	if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
		TPZMatrix<double>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_Forward not decomposed with cholesky");
	
	int n = this->Dim();
	TPZVec<int> pntr(n+1);
    
	for (int i=0; i<n+1; i++)
		pntr[i] = fElem[i] - &fStorage[0] + 1;
    
	char desc[4] = { 'T', 'L', 'N', 'F' };
	char trans = 'N';
	double alfa = 1.0;
    
	for(int j=0; j<B->Cols(); j++) 
		mkl_dskysv(&trans, &n, &alfa, desc, &fStorage[0], &pntr[0], &(*B)(0,j), &(*B)(0,j));
    
	return 1;
}

/*** Subst Backward ***/
//  Perform Ax = b, where A is triangular inferior.
//  Utilizando MKL
template<>
int TPZSkylMatrix<double>::Subst_Backward( TPZFMatrix<double> *B ) const
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_Backward not decomposed with cholesky");
	
	int n = this->Dim();
	TPZVec<int> pntr(n+1);
    
	for (int i=0; i<n+1; i++)
		pntr[i] = fElem[i] - &fStorage[0] + 1;
    
	char desc[4] = { 'T', 'L', 'N', 'F' };
	char trans = 'T';
	double alfa = 1.0;
    
	for(int j=0; j<B->Cols(); j++) 
		mkl_dskysv(&trans, &n, &alfa, desc, &fStorage[0], &pntr[0], &(*B)(0,j), &(*B)(0,j));
    
	return 1;
}
template<>
int
TPZSkylMatrix<float>::Subst_Forward( TPZFMatrix<float> *B ) const
{
	if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
		TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_Forward not decomposed with cholesky");
	
	int n = this->Dim();
	TPZVec<int> pntr(n+1);
    
	for (int i=0; i<n+1; i++)
		pntr[i] = fElem[i] - &fStorage[0] + 1;
    
	char desc[4] = { 'T', 'L', 'N', 'F' };
	char trans = 'N';
	float alfa = 1.0;
    
	for(int j=0; j<B->Cols(); j++) 
		mkl_sskysv(&trans, &n, &alfa, desc, &fStorage[0], &pntr[0], &(*B)(0,j), &(*B)(0,j));
    
	return 1;
}

/*** Subst Backward ***/
//  Perform Ax = b, where A is triangular inferior.
//  Utilizando MKL
template<>
int TPZSkylMatrix<float>::Subst_Backward( TPZFMatrix<float> *B ) const
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_Backward not decomposed with cholesky");
	
	int n = this->Dim();
	TPZVec<int> pntr(n+1);
    
	for (int i=0; i<n+1; i++)
		pntr[i] = fElem[i] - &fStorage[0] + 1;
    
	char desc[4] = { 'T', 'L', 'N', 'F' };
	char trans = 'T';
	float alfa = 1.0;
    
	for(int j=0; j<B->Cols(); j++) 
		mkl_sskysv(&trans, &n, &alfa, desc, &fStorage[0], &pntr[0], &(*B)(0,j), &(*B)(0,j));
    
	return 1;
}
#endif

/*** Subst Forward ***/
//  Perform Ax = b, where A is triangular inferior.
template<class TVar>
int TPZSkylMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__," TPZSkylMatrix::Subst_Forward not decomposed with cholesky");
    
#ifdef DUMP_BEFORE_SUBST
	dump_matrices(this, B, "TPZSkylMatrix::Subst_Forward(B)");
#endif
    
    //	std::cout << "SubstForward this " << (void *) this << " neq " << Dim() << " normb " << Norm(*B) << std::endl;
    int64_t dimension=this->Dim();
    for ( int64_t j = 0; j < B->Cols(); j++ ) {
        int64_t k=0;
        
        // Find the first non-zero element at column j;
        while (k<dimension && (*B)(k,j) == TVar(0)) k++;
        
        //		std::cout << "kstart " << k << std::endl;
        for (; k < dimension; k++ ) {
            
            // sum = SUM( A[k,i] * B[i,j] ), for i = 1, ..., k-1.
            TVar sum = 0.0;
            TVar *elem_ki = fElem[k];
            TVar *end_ki  = fElem[k+1]-1;
            int64_t  array_sz = end_ki - elem_ki;
            TVar* BPtr = &(*B)(k,j) - array_sz; 
            
            for (unsigned l=0; l<array_sz; l++)
                sum += (*elem_ki++) * (*BPtr++);
            
            // B[k,j] = (B[k,j] - sum) / A[k,k].
            BPtr  = &(*B)(k,j);
            *BPtr = (*BPtr - sum) / fElem[k+1][-1];
        }
    }
    
    return 1;
}

/*** Subst Backward ***/
//  Perform Ax = b, where A is triangular inferior.
template<class TVar>
int TPZSkylMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_Backward not decomposed with cholesky");
    
#ifdef DUMP_BEFORE_SUBST
	dump_matrices(this, B, "TPZSkylMatrix::Subst_Backward(B)");
#endif
    
    int64_t Dimension = this->Dim();
    if(!Dimension) 
        return 1;	// nothing to do
    
    for (int64_t j = 0; j < B->Cols(); j++ ) {
        
        int64_t k = Dimension-1;
        
        // Find the last non-zero element at column j of B.
        while (k>0 && (*B)(k,j) == TVar(0.)) k--;
        
        //		std::cout << "kstart " << k << std::endl;
        for ( ; k > 0; k--) {
            // sum = SUM( A[k,i] * B[i,j] ), for i = 1, ..., k-1.
            
            TVar *elem_ki = fElem[k];
            TVar *end_ki  = fElem[k+1]-1; // diagonal(k)
            int64_t  array_sz = end_ki - elem_ki;
            
            TVar *BPtr = &(*B)(k,j);
            *BPtr = *BPtr / *end_ki;
            TVar val = *BPtr;
            
            BPtr = BPtr - array_sz;
            
            for (unsigned l=0; l<array_sz; l++) {
                *BPtr++ -= ((*elem_ki++) * val);
            }
        }
    }
    
    for(int64_t j = 0; j < B->Cols(); j++)
        (*B)(0,j) /= fElem[0][0];
    
    return 1;
}

/*** Subst L Forward ***/
// Forward Substitution suposing diagonal elements are equal to 1.
template<class TVar>
int
TPZSkylMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar> *B ) const 
{
    if ( (B->Rows() != this->Dim()) || (this->fDecomposed != ELDLt && this->fDecomposed != ELU) )
        return( 0 );
    
    int64_t dimension =this->Dim();
    for ( int64_t k = 0; k < dimension; k++ ) {
        for ( int64_t j = 0; j < B->Cols(); j++ ) {
            
            // sum = SUM( A[k,i] * B[i,j] ), for i = 1, ..., k-1.
            TVar sum = 0.0;
            TVar *elem_ki = fElem[k];
            TVar *end_ki  = fElem[k+1]-1;
            int64_t  array_sz = end_ki - elem_ki;
            TVar* BPtr = &(*B)(k,j) - array_sz; 
            for (unsigned l=0; l<array_sz; l++)
                sum += (*elem_ki++) * (*BPtr++);
            
            //EBORIN: Not sure if the code or the commentary is incorrect.
            // Faz B[k,j] = (B[k,j] - sum) / A[k,k].
            BPtr = &(*B)(k,j);
            *BPtr-= sum;
        }
    }
    
    return 1 ;
}

template<class TVar>
int TPZSkylMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar> *B ) const 
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ELDLt) 
        return  0;
    
    int64_t dimension = this->Dim();
    
    for (int64_t j = 0; j < B->Cols(); j++ ) {
        TVar *BPtr = &(*B)(0,j);
        int64_t k=0;
        
        while(k < dimension) {
            //*BPtr = *BPtr / diagonal(k)
            *BPtr++ /= *(fElem[k+1]-1);
            k++;
        }
    }
    
    return 1;
}

template<class TVar>
int TPZSkylMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed || this->fDecomposed == ECholesky) {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_LBackward not decomposed properly");
    }
    
    int64_t Dimension = this->Dim();
    
    for ( int64_t k = Dimension-1; k > 0; k-- ) {
        for ( int64_t j = 0; j < B->Cols(); j++ ) {
            
            TVar *elem_ki = fElem[k];
            TVar *end_ki  = fElem[k+1]-1; // diagonal(k)
            int64_t  array_sz = end_ki - elem_ki;
            TVar* BPtr = &(*B)(k,j) - array_sz; 
            TVar val = (*B)(k,j);
            
            for (unsigned l=0; l<array_sz; l++) {
                *BPtr++ -= (*elem_ki++) * val;
            }
        }
    }
    
    return 1;
}

/**************************** PRIVATE ****************************/
template<class TVar>
int
TPZSkylMatrix<TVar>::Clear()
{
    this->fStorage.Resize(0);
    this->fElem.Resize(0);
    this->fRow = this->fCol = 0;
    this->fDecomposed = 0;
    return( 1 );
}

template<class TVar>
void TPZSkylMatrix<TVar>::Copy(const TPZSkylMatrix<TVar> &A )
{
    int64_t dimension = A.Dim();
    this->fRow = this->fCol = dimension;
    fElem.Resize(A.fElem.NElements());
    fStorage = A.fStorage;
    TVar *firstp = 0;
    if(fStorage.NElements()) 
        firstp = &fStorage[0];
    
    for(int64_t i=0; i<=dimension; i++)
        fElem[i]=firstp+(A.fElem[i]-A.fElem[0]);
    
    this->fDecomposed  = A.fDecomposed;
    this->fDefPositive = A.fDefPositive;
}

template <class TVar>
void TPZSkylMatrix<TVar>::Read(TPZStream &buf, void *context )
{
    TPZMatrix<TVar>::Read(buf, context);
    buf.Read( fStorage);
    TPZVec<int64_t> skyl(this->Rows()+1,0);
    buf.Read( skyl);
    TVar *ptr = 0;
    if (this->Rows()) {
        ptr = &fStorage[0];
    }
    fElem.Resize(this->Rows()+1);
    for (int64_t i=0; i<this->Rows()+1; i++) {
        fElem[i] = skyl[i] + ptr;
    }
}

template <class TVar>
void TPZSkylMatrix<TVar>::Write( TPZStream &buf, int withclassid ) const
{
    TPZMatrix<TVar>::Write(buf,withclassid);
    buf.Write( fStorage);
    TPZVec<int64_t> skyl(this->Rows()+1,0);
    TVar *ptr = 0;
    if (this->Rows()) {
        ptr = &fStorage[0];
    }
    for (int64_t i=0; i<this->Rows()+1; i++) {
        skyl[i] = fElem[i] - ptr;
    }
    buf.Write( skyl);
}

template <class TVar>
void TPZSkylMatrix<TVar>::Write( TPZStream &buf, int withclassid )
{
    TPZMatrix<TVar>::Write(buf,withclassid);
    buf.Write( fStorage);
    TPZVec<int64_t> skyl(this->Rows()+1,0);
    TVar *ptr = 0;
    if (this->Rows()) {
        ptr = &fStorage[0];
    }
    for (int64_t i=0; i<this->Rows()+1; i++) {
        skyl[i] = fElem[i] - ptr;
    }
    buf.Write( skyl);
}

template<class TVar>
void TPZSkylMatrix<TVar>::DecomposeColumn(int64_t col, int64_t prevcol)
{
    TVar *ptrprev;     //Pointer to prev column
    TVar *ptrcol;      //Pointer to col column
    int64_t skprev, skcol; //prev and col Skyline height respectively
    int64_t minline;
    
    skprev = SkyHeight(prevcol);
    skcol = SkyHeight(col);
    
    ptrprev = Diag(prevcol);
    ptrcol = Diag(col);
    
    if((prevcol-skprev) > (col-skcol)){
        minline = prevcol - skprev;
    }
    else {
        minline = col - skcol;
    }
    if(minline > prevcol) {
        cout << "error condition\n";
        cout.flush();
        return;
    }
    TVar *run1 = ptrprev - (prevcol-minline);
    TVar *run2 = ptrcol - (col-minline);
    TVar sum = 0;
    
    while(run1 != ptrprev) 
        sum += (*run1++)*(*run2++);
    *run2-=sum;
    if(run1 != run2){
        *run2 /= *run1;
    }
    else{
        *run2=sqrt(*run2);
    }
}

template<>
void TPZSkylMatrix<std::complex<float> >::DecomposeColumn(int64_t col, int64_t prevcol,std::list<int64_t> &singular)
{
    DebugStop();
}

template<>
void TPZSkylMatrix<std::complex<double> >::DecomposeColumn(int64_t col, int64_t prevcol,std::list<int64_t> &singular)
{
    DebugStop();
}

template<>
void TPZSkylMatrix<std::complex<long double> >::DecomposeColumn(int64_t col, int64_t prevcol,std::list<int64_t> &singular)
{
    DebugStop();
}

template<class TVar>
void TPZSkylMatrix<TVar>::DecomposeColumn(int64_t col, int64_t prevcol,std::list<int64_t> &singular)
{
    TVar *ptrprev;     //Pointer to prev column
    TVar *ptrcol;      //Pointer to col column
    int64_t skprev, skcol; //prev and col Skyline height respectively
    int64_t minline;
    
    skprev = SkyHeight(prevcol);
    skcol = SkyHeight(col);
    
    ptrprev = Diag(prevcol);
    ptrcol = Diag(col);
    
    if((prevcol-skprev) > (col-skcol)){
        minline = prevcol - skprev;
    }
    else {
        minline = col - skcol;
    }
    if(minline > prevcol) {
        cout << "error condition\n";
        cout.flush();
        return;
    }
    TVar *run1 = ptrprev - (prevcol-minline);
    TVar *run2 = ptrcol - (col-minline);
    TVar sum = 0;
    
    while(run1 != ptrprev) 
        sum += (*run1++)*(*run2++);
    *run2-=sum;
    if(run1 != run2){
        *run2 /= *run1;
    }
    else{
        TVar pivot = *run2;
        if ( pivot < TVar(1.e-10) ) {
#ifdef LOG4CXX
            std::stringstream sout;
            sout << "equation " << col << " is singular pivot " << pivot;
            LOGPZ_WARN(logger,sout.str())
#endif
            singular.push_back(col);
            pivot = 1.;
        }    
        *run2=sqrt(pivot);
    }
}

template<>
void TPZSkylMatrix<std::complex<float> >::DecomposeColumn2(int64_t col, int64_t prevcol)
{
    DebugStop();
}

template<>
void TPZSkylMatrix<std::complex<double> >::DecomposeColumn2(int64_t col, int64_t prevcol)
{
    DebugStop();
}

template<>
void TPZSkylMatrix<std::complex<long double> >::DecomposeColumn2(int64_t col, int64_t prevcol)
{
    DebugStop();
}

template<class TVar>
void TPZSkylMatrix<TVar>::DecomposeColumn2(int64_t col, int64_t prevcol)
{
    //Cholesky Decomposition
    TVar *ptrprev;     //Pointer to prev column
    TVar *ptrcol;      //Pointer to col column
    int64_t skprev, skcol; //prev and col Skyline height respectively
    int64_t minline;
    
    skprev = SkyHeight(prevcol);
    skcol = SkyHeight(col);
    
    ptrprev = Diag(prevcol);
    ptrcol = Diag(col);
    
    if((prevcol-skprev) > (col-skcol)){
        minline = prevcol - skprev;
    }
    else {
        minline = col - skcol;
    }
    if(minline > prevcol) {
        cout << "error condition\n";
        cout.flush();
        return;
    }
    //EBORIN: TODO: Improve this code (change run1 and run2 so that dot product increment both)
    TVar *run1 = ptrprev - 1;
    TVar *run2 = ptrcol - ((col-prevcol)+1);
    TVar *lastptr = ptrprev - (prevcol-minline+1);
    TVar sum = 0;
    TVar *modify = ptrcol-(col-prevcol);
#ifndef BLAS
    while(run1 > lastptr) 
        sum += (*run1--)*(*run2--);
#else
    int64_t n=lastptr-run1;
    sum = cblas_ddot(n,run1-(n-1),1,run2-(n-1),1);
#endif
    *modify-=sum;
    if(col != prevcol){
        *modify /= *ptrprev;
    }
    else{
        if ( *modify < TVar(1.e-25) ) {
            cout << "TPZSkylMatrix::DecomposeCholesky a matrix nao e positiva definida" << *modify << endl;
            *modify = 1.e-10;
        }
        *modify=sqrt(*modify);
    }
}

template <class TVar>
void TPZSkylMatrix<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric) {
    if (nrow != ncol || !symmetric)
    {
        DebugStop();
    }
    TPZMatrix<TVar>::Redim(nrow,ncol);
    TPZVec<int64_t> skyline(nrow);
    fElem.resize(nrow+1);
    fElem.Fill(0);
    for (int64_t i=0; i<nrow; i++) {
        skyline[i]=(i*(rand()+RAND_MAX/2))/RAND_MAX;
    }
    InitializeElem(skyline,fStorage,fElem);

    for (int64_t i=0; i<nrow; i++) {
        TVar sum = 0.;
        for (int64_t j=skyline[i]; j<i; j++) {
            TVar val = ((TVar)rand())/RAND_MAX;
            PutVal(j,i,val);
            sum += fabs(val);
        }
        PutVal(i,i,sum+1.);
    }
}

template class TPZSkylMatrix<float>;
template class TPZSkylMatrix<std::complex<float> >;

template class TPZSkylMatrix<double>;
template class TPZSkylMatrix<std::complex<double> >;

#ifndef BORLAND
template class TPZRestoreClass<TPZSkylMatrix<float>>;
template class TPZRestoreClass<TPZSkylMatrix<double>>;
template class TPZRestoreClass<TPZSkylMatrix<long double>>;

template class TPZRestoreClass<TPZSkylMatrix<std::complex<float>>>;
template class TPZRestoreClass<TPZSkylMatrix<std::complex<double>>>;
template class TPZRestoreClass<TPZSkylMatrix<std::complex<long double>>>;
#endif

template class TPZSkylMatrix<long double>;
template class TPZSkylMatrix<std::complex<long double> >;

template class TPZSkylMatrix<TPZFlopCounter>;

#else // USING_NEW_SKYLMAT

/**
 * @file
 * @brief Contains the implementation of the TPZSkylMatrix methods.
 */

#include <math.h>
#include <stdlib.h>

#ifdef BLAS
extern "C" {
#include <cblas.h>
}
#endif

#include "pzfmatrix.h"
#include "pzskylmat.h"

#include <sstream>
#include "pzlog.h"
#include "tpzverysparsematrix.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzskylmatrix"));
#endif

using namespace std;

/*******************/
/*** TPZSkylMatrix ***/


/**************************** PUBLIC ****************************/

/*****************************/
/*** Construtor (int) ***/

template<class TVar>
TPZSkylMatrix<TVar>::TPZSkylMatrix(const int64_t dim )
: TPZRegisterClassId(&TPZSkylMatrix::ClassId),TPZMatrix<TVar>( dim, dim ), fElem(dim+1), fStorage(0)
{
	
	// Inicializa a diagonal (vazia).
	fElem.Fill(0);
}
template<class TVar>
TPZSkylMatrix<TVar>::TPZSkylMatrix(const int64_t dim, const TPZVec<int64_t> &skyline )
: TPZRegisterClassId(&TPZSkylMatrix::ClassId),TPZMatrix<TVar>( dim, dim ), fElem(dim+1), fStorage(0)
{
	
	// Inicializa a diagonal (vazia).
	fElem.Fill(0);
	InitializeElem(skyline,fStorage,fElem);
}

template<class TVar>
void TPZSkylMatrix<TVar>::AddSameStruct(TPZSkylMatrix<TVar> &B, double k){
#ifdef PZDEBUG
	{
		int64_t size = this->fElem.NElements();
		if(size != B.fElem.NElements()){
			PZError << "\nFATAL ERROR at " << __PRETTY_FUNCTION__ << "\n";
			PZError.flush();
			DebugStop();
		}
		for(int64_t i = 0; i < size; i++){
			if((this->fElem[i]-this->fElem[0]) != (B.fElem[i]-B.fElem[0])){
				PZError << "\nFATAL ERROR at " << __PRETTY_FUNCTION__ << "\n";
				PZError.flush();
				DebugStop();
			}
		}
	}
#endif
	
	const int64_t n = this->fStorage.NElements();
	for(int64_t i = 0L; i < n; i++) this->fStorage[i] += TVar(k)*B.fStorage[i];
	
}

/**
 * @brief Updates the values of the matrix based on the values of the matrix
 */
template<class TVar>
void TPZSkylMatrix<TVar>::UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat)
{
    TPZMatrix<TVar> *matrix = mat.operator->();
    TPZSkylMatrix<TVar> *skylmat = dynamic_cast<TPZSkylMatrix<TVar> *>(matrix);
    if (!skylmat) {
        DebugStop();
    }
    if (fStorage.NElements() != skylmat->fStorage.NElements()) {
        DebugStop();
    }
    memcpy(&fStorage[0], &(skylmat->fStorage[0]) , fStorage.NElements()*sizeof(TVar));
    this->fDecomposed = skylmat->fDecomposed;
    this->fDefPositive = skylmat->fDefPositive;
}


template<class TVar>
void TPZSkylMatrix<TVar>::SetSkyline(const TPZVec<int64_t> &skyline)
{
#ifdef PZDEBUG
	for (int64_t i = 0 ; i < this->Rows() ; i++){
		if (skyline[i] < 0 || skyline[i] > i) DebugStop();
	}
#endif
	fElem.Fill(0);
	InitializeElem(skyline,fStorage,fElem);
}
template<class TVar>
int64_t TPZSkylMatrix<TVar>::NumElements(const TPZVec<int64_t> &skyline) {
	int64_t dim = skyline.NElements();
	int64_t i;
	int64_t nelem=0;
	for(i=0; i<dim; i++) {
		nelem += i-skyline[i]+1;
	}
	return nelem;
}

template<class TVar>
void TPZSkylMatrix<TVar>::InitializeElem(const TPZVec<int64_t> &skyline, TPZVec<TVar> &storage, TPZVec<TVar *> &point) {   // JORGE 2013 OUTUBRO ???
	int64_t dim = skyline.NElements();
	int64_t nel = NumElements(skyline);
#ifdef PZDEBUG
    //	std::cout << "Skyline Matrix, Number of elements : " << nel << " in floating point " << nel*sizeof(TVar) << std::endl;
#endif
	storage.Resize(nel);
	storage.Fill(0.);
	int64_t i;
	point.Resize(dim+1);
	if(dim) {
		point[0] = &storage[0];
		point[dim] = &storage[0]+nel;
	} else {
		point[0] = 0;
	}
	for(i=1; i<dim+1; i++)
		point[i] = point[i-1]+(i-1)-skyline[i-1]+1;
}

/**
 Computes the highest skyline of both objects
 */
template<class TVar>
void TPZSkylMatrix<TVar>::ComputeMaxSkyline(const TPZSkylMatrix<TVar> &first, const TPZSkylMatrix<TVar> &second, TPZVec<int64_t> &res) {
	
	if (first.Rows() != second.Rows()) {
		cout<<"ComputeMaxSkyline : incompatible dimension";
		return;
	}
	int64_t i, dim = first.Rows();
	res.Resize(dim+1);
	
	for(i=1; i<dim+1; i++) {
		
		int64_t aux = ( first.Size(i) > second.Size(i) ) ? first.Size(i) : second.Size(i);
		res[i] = i-aux-1;
	}
}

template<class TVar>
TVar &
TPZSkylMatrix<TVar>::operator()(const int64_t r, const int64_t c) {
	int64_t row(r),col(c);
	if ( row > col ) this->Swap( &row, &col );
	
	// Indice do vetor coluna.
	int64_t index = col - row;
	if ( index >= Size(col) ) {
		//Error("TPZSkylMatrix::operator()","Index out of range");
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Index out of range");
		DebugStop();
	}
	return fElem[col][index];
}

template<class TVar>
TVar &TPZSkylMatrix<TVar>::s(const int64_t row, const int64_t col) {
	return operator()(row,col);
}

template<class TVar>
TVar &
TPZSkylMatrix<TVar>::operator()(const int64_t r) {
	return operator()(r,r);
}

//EBORIN: Define these if you want to use the experimental version.
#define DECOMPOSE_CHOLESKY_OPT2
//#define SKYLMATRIX_PUTVAL_OPT1
//#define SKYLMATRIX_GETVAL_OPT1

#ifdef SKYLMATRIX_PUTVAL_OPT1
#pragma message ( "warning: Using experimental version of TPZSkylMatrix<TVAr>::PutVal(...)" )
/**************/
/*** PutVal ***/
template<class TVar>
int
TPZSkylMatrix<TVar>::PutVal(const int64_t r,const int64_t c,const TVar & value )
{
	// inicializando row e col para trabalhar com a triangular superior
    if (r > c) return PutVal(c, r, value);
    
    // Indice do vetor coluna.
    int64_t index = c - r;
    // Se precisar redimensionar o vetor.
    //EBORIN: Do we really need to check this?
    if (index >= Size(c)) {
        if (!IsZero(value)) {
            cout << "TPZSkylMatrix::PutVal Size" << Size(c);
            cout.flush();
            TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Index out of range");
        } 
        else 
            return 1;
    }
    
    fElem[c][index] = value;
    this->fDecomposed = 0;
    return (1);
}
#else
/**************/
/*** PutVal ***/
template<class TVar>
int
TPZSkylMatrix<TVar>::PutVal(const int64_t r,const int64_t c,const TVar & value )
{
	// inicializando row e col para trabalhar com a triangular superior
	int64_t row(r),col(c);
	if ( row > col )
		this->Swap( &row, &col );
	
	// Indice do vetor coluna.
	int64_t index = col - row;
	// Se precisar redimensionar o vetor.
	//EBORIN: Do we really need to check this?
	if ( index >= Size(col) && !IsZero(value)) {
		cout << "TPZSkylMatrix::PutVal Size" << Size(col);
		cout.flush();
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Index out of range");
	} else if(index >= Size(col)) return 1;
	fElem[col][index] = value;
	//  delete[]newVet;
	this->fDecomposed = 0;
	return( 1 );
}
#endif

template<class TVar>
void TPZSkylMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                                  const TVar alpha,const TVar beta ,const int opt) const {
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	
	if (this->fDecomposed != ENoDecompose) {
        //		DebugStop();
	}
	if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__," <matrixs with incompatible dimensions>" );
	if(z.Rows() != x.Rows() || z.Cols() != x.Cols()) z.Redim(x.Rows(),x.Cols());
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		cout << "x.Cols = " << x.Cols() << " y.Cols()"<< y.Cols() << " z.Cols() " << z.Cols() << " x.Rows() " << x.Rows() << " y.Rows() "<< y.Rows() << " z.Rows() "<< z.Rows() << endl;
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__," incompatible dimensions\n");
	}
	this->PrepareZ(y,z,beta,opt);
	int64_t rows = this->Rows();
	int64_t xcols = x.Cols();
	int64_t ic, r;
	for (ic = 0; ic < xcols; ic++) {
		for( r = 0 ; r < rows ; r++ ) {
			int64_t offset = Size(r);
			TVar val = 0.;
			const TVar *p = &x.g((r-offset+1),ic);
			TVar *diag = fElem[r] + offset-1;
			TVar *diaglast = fElem[r];
			while( diag > diaglast ) {
				val += *diag-- * *p;
				p ++;
			}
			if( diag == diaglast ) val += *diag * *p;
			z(r,ic) += val*alpha;
			TVar *zp = &z((r-offset+1),ic);
			val = x.g(r,ic);
			diag = fElem[r] + offset-1;
			while( diag > diaglast ) {
				*zp += alpha * *diag-- * val;
				zp ++;
			}
		}
	}
}

template<class TVar>
void TPZSkylMatrix<TVar>::SolveSOR(int64_t & numiterations,const TPZFMatrix<TVar> &F,
                                   TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &scratch,const REAL overrelax,
                                   REAL &tol,const int FromCurrent,const int direction)  {
	
	if(residual == &F) {
		cout << "TPZMatrix::SolveSOR called with residual and F equal, no solution\n";
		return;
	}
	REAL res = 2.*tol+1.;
	if(residual) res = Norm(*residual);
	if(!FromCurrent) {
		result.Zero();
	}
    TVar over = overrelax;
	int64_t r = this->Dim();
	int64_t c = F.Cols();
	int64_t i,ifirst = 0, ilast = r, iinc = 1;
	if(direction == -1) {
		ifirst = r-1;
		ilast = 0;
		iinc = -1;
	}
	int64_t it;
	for(it=0; it<numiterations && res > tol; it++) {
		res = 0.;
		scratch = F;
		for(int64_t ic=0; ic<c; ic++) {
			if(direction == 1) {
				//
				// compute the upper triangular part first and put into the scractch vector
				//
				for(i=ifirst; i!=ilast; i+= iinc) {
					//TPZColuna *mydiag = &fDiag[i];
					int64_t offset = Size(i);
					TVar val;
					TVar *diag;
					TVar *diaglast = fElem[i];
					TVar *scratchp = &scratch(i-offset+1,ic);
					val = result(i,ic);
					diag = fElem[i] + offset-1;
					int64_t lastid = diag-diaglast;
					int64_t id;
					for(id=0; id<=lastid; id++) *(scratchp+id) -= *(diag-id) * val;
					/* codeguard fix
					 while( diag >= diaglast ) *scratchp++ -= *diag-- * val;
					 */
				}
				//
				// perform the SOR operation
				//
				for(i=ifirst; i!=ilast; i+= iinc) {
					//TPZColuna *mydiag = &fDiag[i];
					int64_t offset = Size(i);
					TVar val = scratch(i,ic);
					TVar *p = &result(i-offset+1,ic);
					TVar *diag = fElem[i] + offset-1;
					TVar *diaglast = fElem[i];
					while( diag > diaglast ) val -= *diag-- * *p++;
					res += abs(val*val);
					result(i,ic) += val*over/ (*diag);
				}
			} else {
				//
				// the direction is upward
				//
				// put the lower triangular part of the multiplication into the scratch vector
				//
				for(i=ifirst; i!=ilast; i+= iinc) {
					//TPZColuna *mydiag = &fDiag[i];
					int64_t offset = Size(i);
					TVar val = scratch(i,ic);
					TVar *p = &result(i-offset+1,ic);
					TVar *diag = fElem[i] + offset-1;
					TVar *diaglast = fElem[i];
					while( diag > diaglast ) val -= *diag-- * *p++;
					//					res += val*val;
					scratch(i,ic) = val;
				}
				//
				// perform the SOR operation
				//
				for(i=ifirst; i!=ilast; i+= iinc) {
					//TPZColuna *mydiag = &fDiag[i];
					int64_t offset = Size(i);
					//	TVar val = scratch(i,ic);
					TVar *diag;
					TVar *diaglast = fElem[i];
					TVar *scratchp = &scratch(i-offset+1,ic);
					//val= result(i,ic);
					TVar val = scratch(i,ic);
					val -= *diaglast * result(i,ic);
					res += abs(val*val);
					val = over * val / *diaglast;
					result(i,ic) += val;
					val = result(i,ic);
					diag = fElem[i] + offset-1;
					while( diag > diaglast ) *scratchp++ -= *diag-- * val;
				}
			}
		}
		res = sqrt(res);
	}
	if(residual) {
		this->Residual(result,F,*residual);
	}
	numiterations = it;
	tol = res;
}

#ifdef SKYLMATRIX_GETVAL_OPT1
#pragma message ( "warning: Using experimental version of TPZSkylMatrix<TVAr>::GetVal(...)" )
template<class TVar>
const TVar &
TPZSkylMatrix<TVar>::GetVal(const int64_t r,const int64_t c ) const
{
    if (r > c) return GetVal(c,r);
    unsigned dim = this->Dim();
    //EBORIN: Do we really need to do this? May only when running debug version.
    if(r >= dim || c >= dim  || r < 0 || c < 0) {
        cout << "TPZSkylMatrix::GetVal index out of range row = " << r
        << " col = " << c << endl;
        return this->gZero;
    }
    
    // Indice do vetor coluna.
    int64_t index   = c - r;
    if ( index < Size(c) ) {
        return (fElem[c][index]);
    }
    else {
        if(this->gZero != TVar(0.)) {
            cout << "TPZSkylMatrix gZero = " << this->gZero << endl;
            DebugStop();
        }
        return(this->gZero );
    }
}
#else
/**************/
/*** GetVal ***/
template<class TVar>
const TVar &
TPZSkylMatrix<TVar>::GetVal(const int64_t r,const int64_t c ) const
{
	// inicializando row e col para trabalhar com a triangular superior
	int64_t row(r),col(c);
	if ( row > col )
		this->Swap( &row, &col );
	
	if(row >= this->Dim() || col >= this->Dim() || row < 0 || col<0) {
		cout << "TPZSkylMatrix::GetVal index out of range row = " << row
		<< " col = " << col << endl;
		return this->gZero;
	}
	// Indice do vetor coluna.
	int64_t index   = col - row;
	//TPZColuna *pCol = &fDiag[col];
	
	if ( index < Size(col) )
		return( fElem[col][index] );
	else {
		if(this->gZero != TVar(0.)) {
			cout << "TPZSkylMatrix gZero = " << this->gZero << endl;
			DebugStop();
		}
		return(this->gZero );
	}
}
#endif
/******** Operacoes com matrizes SKY LINE  ********/

/******************/
/*** Operator = ***/

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator=(const TPZSkylMatrix<TVar> &A )
{
	Clear();
	Copy( A );
	return( *this );
}

/******************/
/*** Operator + ***/

template<class TVar>
TPZSkylMatrix<TVar>
TPZSkylMatrix<TVar>::operator+(const TPZSkylMatrix<TVar> &A) const
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"<incompatible dimensions>" );
	
	TPZVec<int64_t> skylinesize(this->Dim());
	ComputeMaxSkyline(*this,A,skylinesize);
	TPZSkylMatrix res( this->fRow, skylinesize );
	
	TVar *elemMax;
	TVar *elemMin;
	int64_t  sizeMax;
	int64_t  sizeMin;
	
	for ( int64_t col = 0; col < this->Dim(); col++ )
    {
		// Define o tamanho e os elementos da maior e da menor
		//  coluna.
		if ( Size(col) > A.Size(col) )
		{
			sizeMax = Size(col);
			elemMax = fElem[col];
			sizeMin = A.Size(col);
			elemMin = A.fElem[col];
		}
		else
		{
			sizeMax = A.Size(col);
			elemMax = A.fElem[col];
			sizeMin = Size(col);
			elemMin = fElem[col];
		}
		
		// Inicializa coluna da matriz resultado.
		
		// Efetua a SOMA.
		TVar *dest = res.fElem[col];
		int64_t i;
		for ( i = 0; i < sizeMin; i++ )
			*dest++ = (*elemMax++) + (*elemMin++);
		for ( ; i < sizeMax; i++ )
			*dest++ = *elemMax++;
    }
	
	return( res );
}

/******************/
/*** Operator - ***/

template<class TVar>
TPZSkylMatrix<TVar>
TPZSkylMatrix<TVar>::operator-(const TPZSkylMatrix<TVar> &A ) const
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator-( TPZSkylMatrix ) <incompatible dimensions>" );
	
	TPZVec<int64_t> skylinesize(this->Dim());
	ComputeMaxSkyline(*this,A,skylinesize);
	TPZSkylMatrix<TVar> res( this->fRow, skylinesize );
	
	for ( int64_t col = 0; col < this->fRow; col++ )
    {
		// Define o tamanho e os elementos das colunas das 2 matrizes.
		int64_t  sizeThis  = Size(col);
		TVar *elemThis = fElem[col];
		int64_t  sizeA     = A.Size(col);
		TVar *elemA    = A.fElem[col];
		
		// Inicializa coluna da matriz resultado.
		
		// Efetua a SUBTRACAO.
		TVar *dest = res.fElem[col];
		int64_t i;
		for ( i = 0; (i < sizeThis) && (i < sizeA); i++ ) *dest++ = (*elemThis++) - (*elemA++);
		if ( i == sizeA ) for ( ; i < sizeThis; i++ ) *dest++ = *elemThis++;
		else for ( ; i < sizeA; i++ ) *dest++ = -(*elemA++);
    }
	
	return( res );
}

template<class TVar>
void TPZSkylMatrix<TVar>::AddKel(TPZFMatrix<TVar>&elmat,
                                 TPZVec<int64_t> &source,
                                 TPZVec<int64_t> &destination)
{
    int64_t nelem = source.NElements();
    int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
    for(icoef=0; icoef<nelem; icoef++) {
        ieq = destination[icoef];
        ieqs = source[icoef];
        for(jcoef=icoef; jcoef<nelem; jcoef++) {
            jeq = destination[jcoef];
            jeqs = source[jcoef];
            int64_t row(ieq), col(jeq);
            // invertendo linha-coluna para triangular superior
            if (row > col)
                this->Swap(&row, &col);
#ifdef PZDEBUG
            // checando limites
            if(row >= this->Dim() || col >= this->Dim()) {
                cout << "TPZSkylMatrix::GetVal index out of range row = " <<
                row << " col = " << col << endl;
                DebugStop();
            }
#endif
            // indice do vetor coluna
            int64_t index = col - row;
#ifdef PZDEBUG
            // checando limite da coluna
            if (index >= Size(col)) {
                cerr << "Try TPZSkylMatrix gZero." << endl;
                DebugStop();
            }
#endif
            // executando contribuição
            fElem[col][index] += elmat(ieqs,jeqs);

        }
    }
}

template<>
void TPZSkylMatrix<double>::AddKel(TPZFMatrix<double>&elmat,
                                 TPZVec<int64_t> &source,
                                 TPZVec<int64_t> &destination)
{
    int64_t nelem = source.NElements();
    int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
    for(icoef=0; icoef<nelem; icoef++) {
        ieq = destination[icoef];
        ieqs = source[icoef];
        for(jcoef=icoef; jcoef<nelem; jcoef++) {
            jeq = destination[jcoef];
            jeqs = source[jcoef];
            int64_t row(ieq), col(jeq);
            // invertendo linha-coluna para triangular superior
            if (row > col)
                this->Swap(&row, &col);
#ifdef PZDEBUG
            // checando limites
            if(row >= this->Dim() || col >= this->Dim()) {
                cout << "TPZSkylMatrix::GetVal index out of range row = " <<
                row << " col = " << col << endl;
                DebugStop();
            }
#endif
            // indice do vetor coluna
            int64_t index = col - row;
#ifdef PZDEBUG
            // checando limite da coluna
            if (index >= Size(col)) {
                std::cout << "Skyline wrongly configured " << " row " << row << " col " << col << " Size(col) " << Size(col) << std::endl;
                cerr << "Try TPZSkylMatrix gZero." << endl;
                std::cout << destination << std::endl;
                DebugStop();
            }
#endif
            // executando contribuição
#pragma omp atomic
            fElem[col][index] += elmat(ieqs,jeqs);
            
        }
    }
}


template<>
void TPZSkylMatrix<float>::AddKel(TPZFMatrix<float>&elmat,
                                   TPZVec<int64_t> &source,
                                   TPZVec<int64_t> &destination)
{
    int64_t nelem = source.NElements();
    int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
    for(icoef=0; icoef<nelem; icoef++) {
        ieq = destination[icoef];
        ieqs = source[icoef];
        for(jcoef=icoef; jcoef<nelem; jcoef++) {
            jeq = destination[jcoef];
            jeqs = source[jcoef];
            int64_t row(ieq), col(jeq);
            // invertendo linha-coluna para triangular superior
            if (row > col)
                this->Swap(&row, &col);
#ifdef PZDEBUG
            // checando limites
            if(row >= this->Dim() || col >= this->Dim()) {
                cout << "TPZSkylMatrix::GetVal index out of range row = " <<
                row << " col = " << col << endl;
                DebugStop();
            }
#endif
            // indice do vetor coluna
            int64_t index = col - row;
#ifdef PZDEBUG
            // checando limite da coluna
            if (index >= Size(col)) {
                cerr << "Try TPZSkylMatrix gZero." << endl;
                DebugStop();
            }
#endif
            // executando contribuição
#pragma omp atomic
            fElem[col][index] += elmat(ieqs,jeqs);
            
        }
    }
}

/*******************/
/*** Operator += ***/

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator+=(const TPZSkylMatrix<TVar> &A )
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator+=( TPZSkylMatrix ) <incompatible dimensions>" );
	
	TPZSkylMatrix res((*this)+A);
	*this = res;
	return *this;
}

/*******************/
/*** Operator -= ***/

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator-=(const TPZSkylMatrix<TVar> &A )
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator-=( TPZSkylMatrix ) <incompatible dimensions>" );
	
	TPZSkylMatrix res(*this-A);
	*this = res;
	return *this;
}



/******** Operacoes com valores NUMERICOS ********/
//
// As operacoes com valores numericos sao efetuadas apenas nos
// elementos alocados. Em especial, as operacoes A = 0.0 e A *= 0.0
// desalocam todos os elementos da matriz.
//


/*****************************/
/*** Operator * ( TVar ) ***/

template<class TVar>
TPZSkylMatrix<TVar>
TPZSkylMatrix<TVar>::operator*(const TVar value ) const
{
	TPZSkylMatrix res( *this );
	
	for ( int64_t col = 0; col < this->Dim(); col++ )
    {
		// Aloca nova coluna para o resultado.
		int64_t colsize = Size(col);
		// Efetua a SOMA.
		TVar *elemRes  = res.fElem[col];
		for ( int64_t i = 0; i < colsize; i++ )
			*elemRes++ *= value;
    }
	
	return( res );
}





/******************************/
/*** Operator *= ( TVar ) ***/

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator*=(const TVar value )
{
	if ( IsZero( value ) )
    {
		Zero();
		return( *this );
    }
	
	int64_t col, colmax = this->Dim();
	for (col=0; col<colmax; col++ )
    {
		// Efetua a MULTIPLICACAO.
		TVar *elem = fElem[col];
		TVar *end  = fElem[col+1];
		while ( elem < end ) *elem++ *= value;
    }
	
	this->fDecomposed = 0;
	return( *this );
}



/**************/
/*** Resize ***/
//
// Muda as dimensoes da matriz, mas matem seus valores antigos. Novas
// posicoes sao criadas com ZEROS.
//
template<class TVar>
int TPZSkylMatrix<TVar>::Resize( int64_t newDim ,int64_t ) {
	if ( newDim == this->Dim() )
		return( 1 );
	
	fElem.Resize(newDim+1);
	// Cria nova matrix.
	
	// Copia os elementos para a nova matriz.
	int64_t min = MIN( newDim, this->Dim() );
	int64_t i;
	for ( i = min+1; i <= newDim; i++ )
		fElem[i] = fElem[i-1];
	
	// Zera as posicoes que sobrarem (se sobrarem)
	fStorage.Resize(fElem[newDim]-fElem[0]);
	this->fRow = this->fCol = newDim;
	this->fDecomposed = 0;
	return( 1 );
}



/*************/
/*** Redim ***/
//
// Muda as dimensoes da matriz e ZERA seus elementos.
//
template<class TVar>
int
TPZSkylMatrix<TVar>::Redim( int64_t newDim , int64_t)
{
	if ( newDim == this->Dim() )
    {
		Zero();
		return( 1 );
    }
	
	Clear();
	fElem.Resize(newDim);
	fElem.Fill(0);
	this->fRow = this->fCol = newDim;
	this->fDecomposed = 0;
	return( 1 );
}

template<>
int
TPZSkylMatrix<std::complex<float> >::Decompose_Cholesky(std::list<int64_t> &singular)
{
    DebugStop();
    return -1;
}
template<>
int
TPZSkylMatrix<std::complex<double> >::Decompose_Cholesky(std::list<int64_t> &singular)
{
    DebugStop();
    return -1;
}
template<>
int
TPZSkylMatrix<std::complex<long double> >::Decompose_Cholesky(std::list<int64_t> &singular)
{
    DebugStop();
    return -1;
}
/**************************/
/*** Decompose Cholesky ***/
template<class TVar>
int
TPZSkylMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular)
{
	if(this->fDecomposed == ECholesky) return 1;
	if (  this->fDecomposed )  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
    
#ifdef DUMP_BEFORE_DECOMPOSE
	dump_matrix(this, "TPZSkylMatrix::Decompose_Cholesky(singular)");
#endif
    
	singular.clear();
	TVar pivot;
	TVar Tol;
	ZeroTolerance(Tol);
	int64_t dimension = this->Dim();
	/*  if(Dim() > 100) {
	 cout << "\nTPZSkylMatrix Cholesky decomposition Dim = " << Dim() << endl;
	 cout.flush();
	 }*/
	for ( int64_t k = 0; k < dimension; k++ )
	{
		/*    if(!(k%100) && Dim() > 100) {
		 cout <<  k << ' ';
		 cout.flush();
		 }
		 if(!(k%1000)) cout << endl;*/
		if ( Size(k) == 0 )	return( 0 );
		
		// Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
		//
		TVar sum = 0.0;
		TVar *elem_k = fElem[k]+1;
		TVar *end_k  = fElem[k]+Size(k);
#pragma clang loop vectorize_width(2)
		for ( ; elem_k < end_k; elem_k++ ) sum += (*elem_k) * (*elem_k);
		
		// Faz A(k,k) = sqrt( A(k,k) - sum ).
		//
		pivot = fElem[k][0] - sum;
        //		if ( pivot < ((TVar)1.e-9) ) {
		if(pivot < Tol) {
			singular.push_back(k);
            std::cout << __FUNCTION__ << " Singular equation pivot " << pivot << " k " << k << std::endl;
			pivot = 1.;
		}
		// A matriz nao e' definida positiva.
		
		pivot = fElem[k][0] = sqrt( pivot );
		
		// Loop para i = k+1 ... Dim().
		//
		int64_t i=k+1;
		for ( int64_t j = 2; i<dimension; j++,i++ ) {
			// Se tiverem elementos na linha 'i' cuja coluna e'
			//  menor do que 'K'...
			if ( Size(i) > j ) {
				// Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1,..,k-1.
				sum = 0.0;
				TVar *elem_i = &fElem[i][j];
				TVar *end_i  = fElem[i+1];
				elem_k = &(fElem[k][1]);
				// Vectorizable loop
				unsigned max_l = end_i - elem_i;
				unsigned tmp = end_k - elem_k;
				if (tmp < max_l) max_l = tmp;
#pragma clang loop vectorize_width(2)
				for(unsigned l=0; l<max_l; l++) 
                    sum += (*elem_i++) * (*elem_k++);
				// Faz A(i,k) = (A(i,k) - sum) / A(k,k)
				fElem[i][j-1] = (fElem[i][j-1] -sum) / pivot;
			} else if ( Size(i) == j ) fElem[i][j-1] /= pivot;
			
			// Se nao tiverem estes elementos, sum = 0.0.
			
			// Se nao existir nem o elemento A(i,k), nao faz nada.
		}
	}
	
	if(this->Rows() && fabs(GetVal(this->Rows()-1,this->Rows()-1)) < fabs((TVar)1.e-15))
	{
		singular.push_back(this->Rows()-1);
		PutVal(this->Rows()-1,this->Rows()-1,1.);
	}
	this->fDecomposed  = ECholesky;
	this->fDefPositive = 1;
	return( 1 );
}

template<>
int TPZSkylMatrix<std::complex<float> >::Decompose_Cholesky()
{
    DebugStop();
    return -1;
}
template<>
int TPZSkylMatrix<std::complex<double> >::Decompose_Cholesky()
{
    DebugStop();
    return -1;
}
template<>
int TPZSkylMatrix<std::complex<long double> >::Decompose_Cholesky()
{
    DebugStop();
    return -1;
}

clarg::argBool clk_mig("-skl_chk_pm", "Migrate Skyline pages before Cholesky Decomposition");
clarg::argBool clk_rea("-skl_chk_rea", "Reallocate Skyline data before Cholesky Decomposition");

template<class TVar>
int
TPZSkylMatrix<TVar>::Decompose_Cholesky()
{
    if(this->fDecomposed == ECholesky) return 1;
    if (this->fDecomposed )  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	
#ifdef DUMP_BEFORE_DECOMPOSE
    dump_matrix(this, "TPZSkylMatrix::Decompose_Cholesky()");
#endif
    
    TVar pivot;
    TVar minpivot = 10000.;
    int64_t dimension = this->Dim();
    /*  if(Dim() > 100) {
     cout << "\nTPZSkylMatrix Cholesky decomposition Dim = " << Dim() << endl;
     cout.flush();
     }*/
    
    if (clk_mig.was_set())
        MigratePages();
    if (clk_rea.was_set())
        ReallocForNuma();
    
	//	#define DECOMPOSE_CHOLESKY_OPT2 // EBORIN: Optimization 2 -- See bellow
#ifdef DECOMPOSE_CHOLESKY_OPT2
//#pragma message ( "warning: Using experimental (last_col check) version of TPZSkylMatrix<TVar>::Decompose_Cholesky()" )
    PZError << "warning: Using experimental (last_col check) version of TPZSkylMatrix<TVar>::Decompose_Cholesky()" << std::endl;
	TPZVec<int64_t> last_col(dimension);
	{
        int64_t y = dimension-1;
        for (int64_t k=(dimension-1); k>=0; k--) {
            int64_t min_row = k-Size(k)+1;
            while(y>=min_row) {
                last_col[y--] = k;
            }
        } 
	}
#endif
    
	for ( int64_t k = 0; k < dimension; k++ )
    {
		/*      if(!(k%100) && Dim() > 100) {
		 cout <<  k << ' ';
		 cout.flush();
		 }
		 if(!(k%1000)) cout << endl;*/
		if ( Size(k) == 0 )	return( 0 );
		
		// Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
		//
        
		TVar sum = 0.0;
		TVar *elem_k = fElem[k]+1;
		TVar *end_k  = fElem[k]+Size(k);
#pragma clang loop vectorize_width(2)
		for ( ; elem_k < end_k; elem_k++ ) sum += (*elem_k) * (*elem_k);
		
		// Faz A(k,k) = sqrt( A(k,k) - sum ).
		//
		pivot = fElem[k][0] - sum;
        minpivot = minpivot < pivot ? minpivot : pivot;
		if ( pivot < ((TVar)0.) || IsZero(pivot) ) {
			cout << "TPZSkylMatrix::DecomposeCholesky a matrix nao e positiva definida" << pivot << endl;
			return( 0 );
		}
		// A matriz nao e' definida positiva.
		
		pivot = fElem[k][0] = sqrt( pivot );
		
		// Loop para i = k+1 ... Dim().
		//
		int64_t i=k+1;
#ifdef DECOMPOSE_CHOLESKY_OPT2
		//EBORIN: Este laco computa os elementos da linha k e pode ser computado em paralelo...
		// Cada iteracao computa L[k,i] em funcao de A[k,i], L[0...k-1,k] x L[0...k-1,i]
		// - Apenas a porcao L[a...k-1,k] e L[b...k-1,i] nao zero 
		//   (representada na skyline) precisa ser computada => Numero variavel de 
		//   multiplicacoes no laco interno
		// - L[a...k-1,k] e reaproveitada (i-k) vezes.
		// Idea: Substituir i<dimension por i<last_col(k), onde last_col(k) é a última
		//   coluna da matriz que possi dados na linha k.
		int64_t max_i = last_col[k];
		for ( int64_t j = 2; i <= max_i; j++,i++ ) {
#else
            for ( int64_t j = 2; i<dimension; j++,i++ ) {
#endif
                // Se tiverem elementos na linha 'i' cuja coluna e'
                //  menor do que 'K'...
                if ( Size(i) > j ) {
                    // Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1,..,k-1.
                    sum = 0.0;
                    TVar *elem_i = &fElem[i][j];
                    TVar *end_i  = fElem[i+1];
                    elem_k = &(fElem[k][1]);
                    // Vectorizable loop
                    unsigned max_l = end_i - elem_i;
                    unsigned tmp = end_k - elem_k;
                    if (tmp < max_l) max_l = tmp;
#pragma clang loop vectorize_width(2)
                    for(unsigned l=0; l<max_l; l++)
                        sum += (*elem_i++) * (*elem_k++);
                    // Faz A(i,k) = (A(i,k) - sum) / A(k,k)
                    fElem[i][j-1] = (fElem[i][j-1] -sum) / pivot;
                } else if ( Size(i) == j ) fElem[i][j-1] /= pivot;
                
                // Se nao tiverem estes elementos, sum = 0.0.
                
                // Se nao existir nem o elemento A(i,k), nao faz nada.
            }
            
        }
        
    
    
    this->fDecomposed  = ECholesky;
    this->fDefPositive = 1;
#ifdef PZDEBUG
        std::cout << __PRETTY_FUNCTION__ << " minpivot " << minpivot << std::endl;
#endif
    return( 1 );
}

template<>
int TPZSkylMatrix<std::complex<float> >::Decompose_Cholesky_blk(int64_t blk_sz)
{
    DebugStop();
    return -1;
}
template<>
int TPZSkylMatrix<std::complex<double> >::Decompose_Cholesky_blk(int64_t blk_sz)
{
    DebugStop();
    return -1;
}
template<>
int TPZSkylMatrix<std::complex<long double> >::Decompose_Cholesky_blk(int64_t blk_sz)
{
    DebugStop();
    return -1;
}


template<class TVar>
int
TPZSkylMatrix<TVar>::Decompose_Cholesky_blk(int64_t blk_sz)
{
    if(this->fDecomposed == ECholesky) 
        return 1;
    if (this->fDecomposed )  
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
    
    TVar pivot = (TVar)0.;
    TVar minpivot = 10000.;
    int64_t dimension = this->Dim();
    
    for (int64_t blk_i = 0; blk_i < dimension; blk_i += blk_sz) {
        int64_t first_row = blk_i;
        int64_t first_col = blk_i;
        int64_t last_row  = blk_i + MIN(dimension,blk_i+blk_sz);
        
        // Compute band inner nodes
        for ( int64_t j = first_col; j < dimension; j++) {
            
            if (Size(j) == 0)   
                return( 0 ); // Size of column cannot be zero.
            
            int64_t first_i = MAX(first_row, (j+1)-Size(j));
            int64_t last_i = MIN(last_row, j); // j-1 becase we do not need to compute u(j,j)
            for ( int64_t i = first_i; i < last_i; i++) {
                
                // Compute u(i,j) = (a_ij - SUM_k_1_to_i-1 (u_ki * u_kj) ) / uii 
                
                //	TVar  u_ii = fElem[i][0];
                int64_t I = j-i; // fElem[j][I] = A(i,j) I = j-i
                TVar* u_ij = &fElem[j][I];
                
                TVar sum = 0.0;
                TVar *elem_kj = u_ij+1;
                TVar *end_kj  = fElem[j+1];
                TVar *elem_ki = &fElem[i][1];
                TVar *end_ki  = fElem[i+1];
                
                // Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1,..,k-1.
                sum = 0.0;
                
                unsigned max_l = end_kj - elem_kj;
                unsigned tmp = end_ki - elem_ki;
                if (tmp < max_l) max_l = tmp;
#pragma clang loop vectorize_width(2)
                for(unsigned l=0; l<max_l; l++)
                    sum += (*elem_kj++) * (*elem_ki++);
                
                *u_ij = (*u_ij - sum) / pivot;
            } 
            
            // After computing all the elements of this column, compute the diagonal (ujj)
            if (last_i == j)
            {
                TVar sum = 0.0;
                TVar* u_jj = &fElem[j][0];
                TVar *elem_k = fElem[j]+1;
                TVar *end_k  = fElem[j+1];
#pragma clang loop vectorize_width(2)
                for ( ; elem_k < end_k; elem_k++ ) sum += (*elem_k) * (*elem_k);
                pivot = *u_jj - sum;
                
                minpivot = minpivot < pivot ? minpivot : pivot;
                
                if ( pivot < ((TVar)0.) || IsZero(pivot) ) {
                    cout << "TPZSkylMatrix::DecomposeCholesky a matrix nao e positiva definida" << pivot << endl;
                    return( 0 );
                }
                // Valor do elemento da diagonal k,k
                *u_jj = sqrt( pivot );
            }
        }
    }	
    this->fDecomposed  = ECholesky;
    this->fDefPositive = 1;
    return( 1 );
}

/**********************/
/*** Decompose LDLt ***/
template<class TVar>
int
TPZSkylMatrix<TVar>::Decompose_LDLt(std::list<int64_t> &singular)
{
    if( this->fDecomposed == ELDLt) return 1;
    if ( this->fDecomposed )
    {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed with different decomposition>" );
    }
    
#ifdef DUMP_BEFORE_DECOMPOSE
    dump_matrix(this, "TPZSkylMatrix::Decompose_LDLt(singular)");
#endif
    
    singular.clear();
    
    // Third try
    TVar *elj,*ell;
    int64_t j,l,minj,minl,minrow,dimension = this->Dim();
    TVar sum;
    j = 1;
    while(j < dimension) {
        minj = j-Size(j)+1;
        l = minj;
        while(l <= j) {
            minl = l-Size(l)+1;
            minrow = (minj<minl)? minl:minj;
            int64_t k = minrow;
            //			DiagkPtr = fDiag+minrow;
            elj = fElem[j]+j-minrow;
            ell = fElem[l]+l-minrow;
            sum = 0.;
            //EBORIN: 
            // Is this a hot spot?
            while(k < l) {
                sum += *elj-- * *ell-- * *(fElem[k++]);
            }
            *elj -= sum;
            if(ell != elj) *elj /= *ell;
            else if(IsZero(*elj)) {
                singular.push_back(l);
                *elj = 1.;
            }
            l++;
        }
        j++;
    }
    
    if(this->Rows() && IsZero(GetVal(this->Rows()-1,this->Rows()-1)))
    {
        singular.push_back(this->Rows()-1);
        PutVal(this->Rows()-1,this->Rows()-1,1.);
    }
    this->fDecomposed  = ELDLt;
    this->fDefPositive = 0;
    //if(Dim() > 100) cout << endl;
    return( 1 );
}

template<class TVar>
int
TPZSkylMatrix<TVar>::Decompose_LDLt()
{
    
    if( this->fDecomposed == ELDLt) return 1;
    if (  this->fDecomposed )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed with different decomposition>" );
    
#ifdef DUMP_BEFORE_DECOMPOSE
    dump_matrix(this, "TPZSkylMatrix::Decompose_LDLt()");
#endif
    
    // Third try
    TVar *elj,*ell;
    int64_t j,l,minj,minl,minrow,dimension = this->Dim();
    TPZVec<TVar> diag(dimension);
    for(j=0; j<dimension; j++)
    {
        diag[j] = *fElem[j];
    }
    
    //std::cout << "TPZSkylMatrix<TVar>::Decompose_LDLt: dimension = " << dimension  << std::endl;
    
    TVar sum;
    j = 1;
    while(j < dimension) {
        /*    if(!(j%100) && Dim() > 100) {
         cout <<  j << ' ';
         cout.flush();
         }
         if(!(j%1000)) cout << endl;*/
        minj = j-Size(j)+1;
        l = minj;
        while(l <= j) {
            minl = l-Size(l)+1;
            minrow = (minj<minl)? minl:minj;
            int64_t k = minrow;
            //			DiagkPtr = fDiag+minrow;
            elj = fElem[j]+j-minrow;
            ell = fElem[l]+l-minrow;
            TVar *diagptr = &diag[k];
            sum = 0.;
            while(k < l) {
                //		  sum += *elj-- * *ell-- * *(fElem[k++]);
                //EBORIN: trocar *diagptr++ por *diagptr-- ajuda na vetorização?
                sum += *elj-- * *ell-- * *diagptr++;
                k++;
            }
            *elj -= sum;
            if(ell != elj) *elj /= *ell;
            else if(IsZero(*elj)) {
#ifdef LOG4CXX
                std::stringstream sout;
                sout << "col = " << j << " diagonal " << *elj;
                LOGPZ_DEBUG(logger,sout.str())
#endif
                
                *diagptr = *elj;
                cout << "TPZSkylMatrix pivot = " << *elj << endl;
                cout << "TPZSkylMatrix::DecomposeLDLt zero pivot\n";
                cout << "j = " << j << " l = " << l << endl;
            }
            else
            {
                *diagptr = *elj;
            }
            l++;
        }
        j++;
    }
    this->fDecomposed  = ELDLt;
    this->fDefPositive = 0;
    //if(Dim() > 100) cout << endl;
    return( 1 );
}

/*********************/
/*** Subst Forward ***/
//
//  Faz Ax = b, onde A e' triangular inferior.
//
template<class TVar>
int
TPZSkylMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_Forward not decomposed with cholesky");
    
#ifdef DUMP_BEFORE_SUBST
    dump_matrices(this, B, "TPZSkylMatrix::Subst_Forward(B)");
#endif
    
    //	std::cout << "SubstForward this " << (void *) this << " neq " << Dim() << " normb " << Norm(*B) << std::endl;
    int64_t dimension=this->Dim();
    for ( int64_t j = 0; j < B->Cols(); j++ )
    {
        int64_t k=0;
        while (k<dimension && (*B)(k,j) == TVar(0)) {
            k++;
        }
        //		std::cout << "kstart " << k << std::endl;
        for (; k < dimension; k++ )
        {
            // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
            //
            TVar sum = 0.0;
            TVar *elem_ki = fElem[k]+1;
            TVar *end_ki  = fElem[k+1];
            TVar* BPtr = &(*B)(k,j);   //(k-1,j)
            //	for ( int i = k-1; elem_ki < end_ki; i-- )
            //	  sum += (*elem_ki++) * B->GetVal( i, j );
            
            //EBORIN:
            // Is this a hot-spot?
            // Is it vectorized?
            while(elem_ki < end_ki) sum += (*elem_ki++) * (*--BPtr);//(*BPtr--)
            // Faz B[k,j] = (B[k,j] - sum) / A[k,k].
            //
            //	B->PutVal( k, j, (B->GetVal(k, j) - sum) / row_k->pElem[0] );
            BPtr = &(*B)(k,j);
            *BPtr-= sum;
            *BPtr /= fElem[k][0];
        }
    }
    
    return( 1 );
}

/*********************/
/*** Subst Backward ***/
//
//  Faz Ax = b, onde A e' triangular superior.
//
template<class TVar>
int
TPZSkylMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_Backward not decomposed with cholesky");
    
#ifdef DUMP_BEFORE_SUBST
    dump_matrices(this, B, "TPZSkylMatrix::Subst_Backward(B)");
#endif
    int64_t Dimension = this->Dim();
    if(!Dimension) return 1;	// nothing to do
    int64_t j;
    for ( j = 0; j < B->Cols(); j++ )
    {
        int64_t k = Dimension-1;
        while (k>0 && (*B)(k,j) == TVar(0.)) {
            k--;
        }
        //		std::cout << "kstart " << k << std::endl;
        
        for (;k > 0; k-- )
        {
            // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
            //
            TVar val;
            TVar *elem_ki = fElem[k];
            TVar *end_ki  = fElem[k+1];
            TVar *BPtr = &(*B)(k,j);
            *BPtr /= *elem_ki++;
            val = *BPtr;
            //	BPtr;
            // substract the column of the skyline matrix from the vector.
            //EBORIN:
            // Is this a hot-spot?
            // Is it vectorized?
            while(elem_ki < end_ki) *--BPtr -= (*elem_ki++) * val;
        }
    }
    for( j = 0; j< B->Cols(); j++) (*B)(0,j) /= fElem[0][0];
    return( 1 );
}

/***********************/
/*** Subst L Forward ***/
//
//  Faz a "Forward substitution" assumindo que os elementos
//   da diagonal sao todos iguais a 1.
//
template<class TVar>
int
TPZSkylMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar> *B ) const {
    if ( (B->Rows() != this->Dim()) || (this->fDecomposed != ELDLt && this->fDecomposed != ELU) )
        return( 0 );
    
    int64_t dimension =this->Dim();
    for ( int64_t k = 0; k < dimension; k++ ) {
        for ( int64_t j = 0; j < B->Cols(); j++ ) {
            // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
            //
            TVar sum = 0.0;
            TVar *elem_ki = fElem[k]+1;
            TVar *end_ki  = fElem[k+1];
            TVar *BPtr = &(*B)(k,j);
            while(elem_ki < end_ki) sum += (*elem_ki++) * (*--BPtr);
            
            // Faz B[k,j] = (B[k,j] - sum) / A[k,k].
            //
            //	B->PutVal( k, j, (B->GetVal(k, j) - sum) / row_k->pElem[0] );
            BPtr = &(*B)(k,j);
            *BPtr-= sum;
        }
    }
    return( 1 );
}

/******************/
/*** Subst Diag ***/
//
//  Faz Ax = b, sendo que A e' assumida ser uma matriz diagonal.
//
template<class TVar>
int
TPZSkylMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar> *B ) const {
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ELDLt) return( 0 );
    int64_t dimension = this->Dim();
    if (!dimension) {
        return 1;
    }
    for ( int64_t j = 0; j < B->Cols(); j++ ) {
        TVar *BPtr = &(*B)(0,j);
        int64_t k=0;
        while(k < dimension) *BPtr++ /= *(fElem[k++]);
    }
    return( 1 );
}

/*********************/
/*** Subst Backward ***/
//
//  Faz Ax = b, onde A e' triangular superior.
//
template<class TVar>
int
TPZSkylMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed || this->fDecomposed == ECholesky)
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_LBackward not decomposed properly");
    
    int64_t Dimension = this->Dim();
    for ( int64_t k = Dimension-1; k > 0; k-- ) {
        for ( int64_t j = 0; j < B->Cols(); j++ ) {
            // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
            //
            //		TVar val = 0.0;
            TVar *elem_ki = fElem[k]+1;
            TVar *end_ki  = fElem[k+1];
            TVar *BPtr = &(*B)(k,j);
            // substract the column of the skyline matrix from the vector.
            
            TVar val = *BPtr;
            while(elem_ki < end_ki) *--BPtr -= (*elem_ki++) * val;
        }
    }
    return( 1 );
}


/**************************** PRIVATE ****************************/


/***************/
/*** PutZero ***/

template<class TVar>
int
TPZSkylMatrix<TVar>::Zero()
{
    
    fStorage.Fill(0.);
    this->fDecomposed = 0;
    this->fDefPositive = 0;
    return( 1 );
}

/*************/
/*** CLear ***/

template<class TVar>
int
TPZSkylMatrix<TVar>::Clear()
{
    this->fStorage.Resize(0);
    //	fStorage.Shrink();
    this->fElem.Resize(0);
    this->fRow = this->fCol = 0;
    this->fDecomposed = 0;
    return( 1 );
}

/************/
/*** Copy ***/

template<class TVar>
void
TPZSkylMatrix<TVar>::Copy(const TPZSkylMatrix<TVar> &A )
{
    int64_t dimension = A.Dim();
    this->fRow = this->fCol = dimension;
    fElem.Resize(A.fElem.NElements());
    fStorage = A.fStorage;
    int64_t i;
    TVar *firstp = 0;
    if(fStorage.NElements()) firstp = &fStorage[0];
    for(i=0; i<=dimension; i++)
        fElem[i]=firstp+(A.fElem[i]-A.fElem[0]);
    this->fDecomposed  = A.fDecomposed;
    this->fDefPositive = A.fDefPositive;
    
}

template<class TVar>
TPZSkylMatrix<TVar>
TPZSkylMatrix<TVar>::operator-() const { return operator*(-1.0); }

template <class TVar>
void TPZSkylMatrix<TVar>::Read(TPZStream &buf, void *context )
{
    TPZMatrix<TVar>::Read(buf, context);
    buf.Read( fStorage);
    TPZVec<int64_t> skyl(this->Rows()+1,0);
    buf.Read( skyl);
    TVar *ptr = 0;
    if (this->Rows()) {
        ptr = &fStorage[0];
    }
    fElem.Resize(this->Rows()+1);
    for (int64_t i=0; i<this->Rows()+1; i++) {
        fElem[i] = skyl[i] + ptr;
    }
}

template <class TVar>
void TPZSkylMatrix<TVar>::Write( TPZStream &buf, int withclassid ) const
{
    TPZMatrix<TVar>::Write(buf,withclassid);
    buf.Write( fStorage);
    TPZVec<int64_t> skyl(this->Rows()+1,0);
    TVar *ptr = 0;
    if (this->Rows()) {
        ptr = &fStorage[0];
    }
    for (int64_t i=0; i<this->Rows()+1; i++) {
        skyl[i] = fElem[i] - ptr;
    }
    buf.Write( skyl);
}


template<class TVar>
void TPZSkylMatrix<TVar>::DecomposeColumn(int64_t col, int64_t prevcol){
    TVar *ptrprev;     //Pointer to prev column
    TVar *ptrcol;      //Pointer to col column
    int64_t skprev, skcol; //prev and col Skyline height respectively
    int64_t minline;
    
    skprev = SkyHeight(prevcol);
    skcol = SkyHeight(col);
    
    ptrprev = Diag(prevcol);
    ptrcol = Diag(col);
    
    if((prevcol-skprev) > (col-skcol)){
        minline = prevcol - skprev;
    }else
    {
        minline = col - skcol;
    }
    if(minline > prevcol) {
        cout << "error condition\n";
        cout.flush();
        return;
    }
    TVar *run1 = ptrprev + (prevcol-minline);
    TVar *run2 = ptrcol + (col-minline);
    TVar sum = 0;
    
    while(run1 != ptrprev) sum += (*run1--)*(*run2--);
    *run2-=sum;
    if(run1 != run2){
        *run2 /= *run1;
    }else{
        *run2=sqrt(*run2);
    }
    
}
template<>
void TPZSkylMatrix<std::complex<float> >::DecomposeColumn(int64_t col, int64_t prevcol,std::list<int64_t> &singular)
{
    DebugStop();
}
template<>
void TPZSkylMatrix<std::complex<double> >::DecomposeColumn(int64_t col, int64_t prevcol,std::list<int64_t> &singular)
{
    DebugStop();
}
template<>
void TPZSkylMatrix<std::complex<long double> >::DecomposeColumn(int64_t col, int64_t prevcol,std::list<int64_t> &singular)
{
    DebugStop();
}

template<class TVar>
void TPZSkylMatrix<TVar>::DecomposeColumn(int64_t col, int64_t prevcol,std::list<int64_t> &singular){
    TVar *ptrprev;     //Pointer to prev column
    TVar *ptrcol;      //Pointer to col column
    int64_t skprev, skcol; //prev and col Skyline height respectively
    int64_t minline;
    
    skprev = SkyHeight(prevcol);
    skcol = SkyHeight(col);
    
    ptrprev = Diag(prevcol);
    ptrcol = Diag(col);
    
    if((prevcol-skprev) > (col-skcol)){
        minline = prevcol - skprev;
    }else
    {
        minline = col - skcol;
    }
    if(minline > prevcol) {
        cout << "error condition\n";
        cout.flush();
        return;
    }
    TVar *run1 = ptrprev + (prevcol-minline);
    TVar *run2 = ptrcol + (col-minline);
    TVar sum = 0;
    
    while(run1 != ptrprev) sum += (*run1--)*(*run2--);
    *run2-=sum;
    if(run1 != run2){
        *run2 /= *run1;
    }else{
        TVar pivot = *run2;
        if ( fabs(pivot) < fabs((TVar)1.e-10) ) {
#ifdef LOG4CXX
            std::stringstream sout;
            sout << "equation " << col << " is singular pivot " << pivot;
            LOGPZ_WARN(logger,sout.str())
#endif
            singular.push_back(col);
            pivot = 1.;
        }
        
        *run2=sqrt(pivot);
    }
    
}

template<>
void TPZSkylMatrix<std::complex<float> >::DecomposeColumn2(int64_t col, int64_t prevcol)
{
    DebugStop();
}
template<>
void TPZSkylMatrix<std::complex<double> >::DecomposeColumn2(int64_t col, int64_t prevcol)
{
    DebugStop();
}

template<>
void TPZSkylMatrix<std::complex<long double> >::DecomposeColumn2(int64_t col, int64_t prevcol)
{
    DebugStop();
}


template<class TVar>
void TPZSkylMatrix<TVar>::DecomposeColumn2(int64_t col, int64_t prevcol){
    
    //Cholesky Decomposition
    TVar *ptrprev;     //Pointer to prev column
    TVar *ptrcol;      //Pointer to col column
    int64_t skprev, skcol; //prev and col Skyline height respectively
    int64_t minline;
    
    skprev = SkyHeight(prevcol);
    skcol = SkyHeight(col);
    
    ptrprev = Diag(prevcol);
    ptrcol = Diag(col);
    
    if((prevcol-skprev) > (col-skcol)){
        minline = prevcol - skprev;
    }else
    {
        minline = col - skcol;
    }
    if(minline > prevcol) {
        cout << "error condition\n";
        cout.flush();
        return;
    }
    TVar *run1 = ptrprev + 1;
    TVar *run2 = ptrcol + (col-prevcol)+1;
    TVar *lastptr = ptrprev + prevcol-minline+1;
    TVar sum = 0;
    TVar *modify = ptrcol+(col-prevcol);
#ifndef BLAS
    while(run1 != lastptr) sum += (*run1++)*(*run2++);
    
#else
    int64_t n=lastptr-run1;
    sum = cblas_ddot(n,run1,1,run2,1);
#endif
    *modify-=sum;
    if(col != prevcol){
        *modify /= *ptrprev;
    }else{
        if ( fabs(*modify) < fabs((TVar)1.e-25) ) {
            cout << "TPZSkylMatrix::DecomposeCholesky a matrix nao e positiva definida" << *modify << endl;
            *modify = ((TVar)1.e-10);
        }
        *modify=sqrt(*modify);
    }
}

template <class TVar>
void TPZSkylMatrix<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric) {
    if (nrow != ncol || !symmetric)
    {
        DebugStop();
    }
    TPZMatrix<TVar>::Redim(nrow,ncol);
    TPZVec<int64_t> skyline(nrow);
    fElem.resize(nrow+1);
    fElem.Fill(0);
    for (int64_t i=0; i<nrow; i++) {
        skyline[i]=(i*(rand()))/RAND_MAX;
    }
//    std::cout << "skyline " << skyline << std::endl;
    InitializeElem(skyline,fStorage,fElem);
    
    for (int64_t i=0; i<nrow; i++) {
        TVar sum = 0.;
        for (int64_t j=skyline[i]; j<i; j++) {
            TVar val = ((TVar)rand())/((TVar)RAND_MAX);
            PutVal(j,i,val);
            sum += fabs(val);
        }
        PutVal(i,i,sum+(TVar)1.);
    }
    for (int64_t i=0; i<nrow; i++) {
        TVar sum = (TVar)0.;
        for (int64_t j=0; j<nrow; j++) {
            TVar val;
            val = GetVal(j,i);
            sum += fabs(val);
        }
        PutVal(i,i,sum+(TVar)1.);
    }

}

template class TPZSkylMatrix<float>;
template class TPZSkylMatrix<std::complex<float> >;

template class TPZSkylMatrix<double>;
template class TPZSkylMatrix<std::complex<double> >;

#ifndef BORLAND
template class TPZRestoreClass<TPZSkylMatrix<float>>;
template class TPZRestoreClass<TPZSkylMatrix<double>>;
template class TPZRestoreClass<TPZSkylMatrix<long double>>;

template class TPZRestoreClass<TPZSkylMatrix<std::complex<float>>>;
template class TPZRestoreClass<TPZSkylMatrix<std::complex<double>>>;
template class TPZRestoreClass<TPZSkylMatrix<std::complex<long double>>>;
#endif

template class TPZSkylMatrix<long double>;
template class TPZSkylMatrix<std::complex<long double> >;

template class TPZSkylMatrix<TPZFlopCounter>;

#endif // USING_NEW_SKYLMAT

#if (defined DUMP_BEFORE_DECOMPOSE) || (defined DUMP_BEFORE_SUBST)

pthread_mutex_t dump_matrix_mutex = PTHREAD_MUTEX_INITIALIZER;
unsigned matrix_unique_id = 0;
clarg::argString dm_prefix("-dm_prefix", 
                           "Filename prefix for matrices dumped before decompose/subst", 
                           "matrix_");

template<class TVar>
void dump_matrix(TPZMatrix<TVar>* m, const char* fn_annotation)
{
    if (!dm_prefix.was_set()) return;
    PZ_PTHREAD_MUTEX_LOCK(&dump_matrix_mutex, "dump_matrix");
    std::stringstream fname;
    fname << dm_prefix.get_value() << fn_annotation << "_" << matrix_unique_id++ << ".bin";
    std::cout << "Dump matrix before... (file: " << fname << ")" << std::endl;
    TPZBFileStream fs;
    fs.OpenWrite(fname.str());
    m->Write(fs, 0);
    std::cout << "Dump matrix before... [Done]" << std::endl;
    PZ_PTHREAD_MUTEX_UNLOCK(&dump_matrix_mutex, "dump_matrix");
}

template<class TVar>
void dump_matrices(const TPZMatrix<TVar>* a, const TPZMatrix<TVar>* b, const char* fn_annotation)
{
    if (!dm_prefix.was_set()) return;
    PZ_PTHREAD_MUTEX_LOCK(&dump_matrix_mutex, "dump_matrices");
    std::stringstream fname;
    fname << dm_prefix.get_value() << fn_annotation << "_" << matrix_unique_id++ << ".bin";
    std::cout << "Dump matrix before... (file: " << fname << ")" << std::endl;
    TPZBFileStream fs;
    fs.OpenWrite(fname.str());
    a->Write(fs, 0);
    b->Write(fs, 0);
    std::cout << "Dump matrix before... [Done]" << std::endl;
    PZ_PTHREAD_MUTEX_UNLOCK(&dump_matrix_mutex, "dump_matrices");
}
#endif


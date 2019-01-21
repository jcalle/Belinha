#ifndef TSKYLMATH
#define TSKYLMATH

#ifdef USING_NEW_SKYLMAT

/**
 * @file 
 * @brief Contains a revised version of the TPZSkyline class which
 * implements a skyline storage format.
 * @ingroup matrix
 * @author Edson Borin
 */
#include "pzmatrix.h"
#include "pzvec.h"
//#include "pzmanvector.h" // TODO: verify if this include file is required.

//TODO: why we need the OOPARLIB ifdef? 
#ifdef OOPARLIB
#include "pzsaveable.h"
#include "pzmatdefs.h"
#endif

#include "pzfmatrix.h"

/*
 The fStorate and fElem arrays store the columns data and the pointers to the
 columns in the following way.
 
 Let the following matrix be a 6x6 skyline matrix.
 
 A B D
 C E G
 F H 
 I J
 K 
 L
 
 fStorage is an array that stores the columns data in the following way:	   
 
 A B C D E F G H I J K L
 
 and fElem is an array of pointers to the first element of each column.
 
 fStorage = A B C D E F G H I J K L
 ^ ^   ^     ^     ^   ^  ^
 | |   |     |     |   |  |
 | +-+ +-+   |   +-+ +-+ ++
 |   |   |   |   |   |   |
 fElem    = 0   1   3   6   9  11  12
 
 
 M[r,c] = (fElem[c+1]-1) - (c-r) 
 |                |           
 Diag.            Dist to Diag.  
 
 Ensure (c-r) < Size(c) == fElem[c+1] - fElemn[c-1]
 
 Size(c) == fElem[c+1]-fElem[c]
 
 M[r,c] = (fElem[c+1]-1) - (c-r) 
 =  fElem[c] + Size(c) - 1 - c + r
 =  fElem[c] + (Size(c) - 1 - c + r)
 =  fElem[c] + (r + Size(c) - c - 1)
 =  fElem[c] + (r + (Size(c)-1) - c)
 =  fElem[c] + r + ((Size(c)-1) - c)
 
 Hence:
 element index = r + Size(c) -1 -c 
 
 
 */


/**
 * @brief Implements a skyline storage format. A Skyline matrix is symmetric,
 *        hence square. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSkylMatrix : public TPZMatrix<TVar>
{
public:
    
    TPZSkylMatrix() : TPZMatrix<TVar>(0,0), fElem(0), fStorage(0) 
    { }
    
    virtual ~TPZSkylMatrix() { Clear(); }
    
    /**
     @brief Constructs a dim x dim skyline matrix with zero elements. 
     */ 
    TPZSkylMatrix(const int64_t dim);
    
    //TODO: Verify the descriptions...
    /**
     @brief Constructs a skyline matrix of dimension dim. The skyline array
     indicates the minimum row number which will be accessed by each equation.
     */
    TPZSkylMatrix(const int64_t dim ,const TPZVec<int64_t> &skyline);
    
    TPZSkylMatrix(const TPZSkylMatrix<TVar> &A ) : 
    TPZRegisterClassId(&TPZSkylMatrix::ClassId), TPZMatrix<TVar>(A), fElem(0), fStorage(0)  
    { Copy(A); }
	
    CLONEDEF(TPZSkylMatrix)
    
    virtual int64_t MemoryFootprint() const {
        return (sizeof(TVar*)*fElem.size() + 
                sizeof(TVar)*fStorage.size());
    }	
    
    /**
     @brief Modify the skyline of the matrix, throwing away its values.  The
     skyline array indicates the minimum row number which will be accessed by
     each equation
     */
    void SetSkyline(const TPZVec<int64_t> &skyline);
	
    friend class TPZSkylMatrix<float>;
    friend class TPZSkylMatrix<double>;
    
    /// copy the values from a matrix with a different precision
    template<class TVar2>
    void CopyFrom(TPZSkylMatrix<TVar2> &orig)
    {
        TPZMatrix<TVar>::CopyFrom(orig);
        TPZVec<TVar2> &origstorage = orig.Storage();
        TPZVec<TVar2> &origelem = orig.Elem();
        int64_t nel = origstorage.size();
        fElem.resize(origelem.size());
        fStorage.resize(nel);
        for (int64_t el=0; el<nel; el++) {
            fStorage[el] = origstorage[el];
        }
        int64_t size_el = fElem.size();
        TVar *first = &fStorage[0];
        TVar2 *first_orig = &origstorage[0];
        for (int64_t el=0; el<size_el; el++) {
            fElem[el] = first+(origelem[el]-first_orig);
        }
        
    }
    

    const TPZVec<TVar> &Elem()
    {
        return fElem;
    }
    
    const TPZVec<TVar. &Storage()
    {
        return fStorage;
    }
    /**
     @brief Returns the height of the skyline for a given column.
     */
    int64_t SkyHeight(int64_t col) const { return fElem[col+1]-fElem[col] - 1; }
	
    /**
     @brief Adds a skyline matrix B (with same structure of this).  
     *this += k * B
     */
    void AddSameStruct(TPZSkylMatrix<TVar> &B, double k = 1.);
	
    /** @brief Skyline matrices are always simetric. */
    virtual int IsSimetric() const {return 1;}
	
    /**
     @brief Updates the values of the matrix based on the values of matrix
     mat. It suposes the matrix mat is a TPZSkylMatrix.
     */
    virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat);
	
    /**
     @brief Updates the value of element [row,col]. If element is not inside the
     skyline, it generates an Index out of range error. Otherwise, it returns 1.
     */
    int PutVal(const int64_t row,const int64_t col,const TVar &element );
    
    /**
     @brief Returns the value of element [row,col]
     */
    const TVar& GetVal(const int64_t row,const int64_t col ) const;
	
    TVar& operator()(const int64_t row, const int64_t col);
    
    virtual TVar& s(const int64_t row, const int64_t col) { return operator()(row, col); }
	
    TVar &operator()(const int64_t row) { return operator()(row, row); }
    
    virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                         const TVar alpha,const TVar beta ,const int opt = 0) const ;
    
    TPZSkylMatrix& operator= (const TPZSkylMatrix<TVar> &A );
	
    TPZSkylMatrix operator+  (const TPZSkylMatrix<TVar> &A ) const;
    TPZSkylMatrix operator-  (const TPZSkylMatrix<TVar> &A ) const;
	
    TPZSkylMatrix &operator+=(const TPZSkylMatrix<TVar> &A );
    TPZSkylMatrix &operator-=(const TPZSkylMatrix<TVar> &A );
	
    TPZSkylMatrix operator*  (const TVar v ) const;
    TPZSkylMatrix &operator*=( TVar value );
	
    TPZSkylMatrix operator-() const { return operator*(-1.0); }
    
	
    /**
     @brief Resize the skyline matrix, but keep it's original elements. The
     second parameter is the size of the columns.
     */
    int Resize(const int64_t newDim ,const int64_t );
    
    /**
     @brief Resize the skyline matrix and replace the values by zeros.  The
     second parameter is the size of the columns.
     */
    int Redim(const int64_t newDim ,const int64_t );
    int Redim(const int64_t newDim) {return Redim(newDim,newDim);}
	
    /**
     * @brief Add a contribution of a stiffness matrix
     * @param elmat Element matrix to be contributed
     * @param sourceindex Contains source indexes on current matrix
     * @param destinationindex Contains destine indexes on current matrix
     */
    
    void AddKel(TPZFMatrix<TVar>&elmat, TPZVec<int64_t> &sourceindex,  TPZVec<int64_t> &destinationindex);

    /**
     @brief Replace the values by zeros.
     */
    int Zero();
	
    /**
     @brief Methods to Solve Linear Equations Systems. 
     */
    // @{
    virtual void SolveSOR(int64_t &numiterations,const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
                          TPZFMatrix<TVar> *residual,TPZFMatrix<TVar> &scratch,const REAL overrelax, REAL &tol,
                          const int FromCurrent = 0,const int direction = 1) ;
    
    int Decompose_Cholesky();  // Faz A = GGt.
    int Decompose_Cholesky_blk(int64_t blk_sz);
    
    int Decompose_LDLt    ();  // Faz A = LDLt.
    int Decompose_Cholesky(std::list<int64_t> &singular);  // Faz A = GGt.
    int Decompose_LDLt    (std::list<int64_t> &singular);  // Faz A = LDLt.
    
    int Subst_Forward  ( TPZFMatrix<TVar> *b ) const;
    int Subst_Backward ( TPZFMatrix<TVar> *b ) const;
    int Subst_LForward ( TPZFMatrix<TVar> *b ) const;
    int Subst_LBackward( TPZFMatrix<TVar> *b ) const;
    int Subst_Diag     ( TPZFMatrix<TVar> *b ) const;
    // @}
	
    virtual void AutoFill(int64_t nrow, int64_t ncol, int symmetric);
	
    public:
virtual int ClassId() const;

    
    /**
     * @brief Unpacks the object structure from a stream of bytes
     * @param buf The buffer containing the object in a packed form
     * @param context 
     */
    virtual void Read(TPZStream &buf, void *context );
    
    /**
     * @brief Packs the object structure in a stream of bytes
     * @param buf Buffer which will receive the bytes
     * @param withclassid
     */
    
    virtual void Write( TPZStream &buf, int withclassid );
    virtual void Write( TPZStream &buf, int withclassid ) const;
    
    virtual std::string ClassName() const { return( "TPZSkylMatrix"); }
    
    int64_t Size(const int64_t column) const
    {
#ifdef PZDEBUG
        if (column < 0 || column >= this->Rows()) {
            DebugStop();
        }
#endif
        return fElem[column+1]-fElem[column];
    }
    
protected:
    
    /**
     @brief This method returns a pointer to the diagonal element of the matrix
     of the col column
     */
    TVar* Diag(int64_t col) { return fElem[col+1]-1;}
	
    void DecomposeColumn(int64_t col, int64_t prevcol);
    void DecomposeColumn(int64_t col, int64_t prevcol, std::list<int64_t> &singular);
    void DecomposeColumn2(int64_t col, int64_t prevcol);
    
private:
	
    int  Clear();
    void Copy (const TPZSkylMatrix<TVar> & );
    static int64_t NumElements(const TPZVec<int64_t> &skyline);
    static void InitializeElem(const TPZVec<int64_t> &skyline, TPZVec<TVar> &storage, TPZVec<TVar *> &elem);
    
    /**
     @brief Computes the highest skyline of both objects
     */
    static void ComputeMaxSkyline(const TPZSkylMatrix<TVar> &first, 
                                  const TPZSkylMatrix<TVar> &second, 
                                  TPZVec<int64_t> &res);
    
    void MigratePages() {
        fElem.MigratePages();
        fStorage.MigratePages();
    }
    
    void ReallocForNuma() {
        fElem.ReallocForNuma();
        TVar* old_start = &fStorage[0];
        fStorage.ReallocForNuma();
        TVar* new_start = &fStorage[0];
        for (int i=0; i<fElem.size(); i++) {
            fElem[i] = new_start + (fElem[i]-old_start);
        }
    }	
protected:
    /**
     @brief Storage to keep the first elements to each equation. 
     fElem is of size number of equation+1. 
     fElem[i] is the first element of the skyline of equation i. 
     fElem[Rows()] is one element beyond the last equation.
     */
    TPZVec<TVar *> fElem;
    
private:
    /**
     @brief fStorage is a single array that contains all the data of the skyline
     matrix.
     */
    TPZVec<TVar> fStorage;
};

/**
 @brief Implements iterative sum over N steps.
 */
template<class TVar,int N>
inline TVar TemplateSum(const TVar *p1, const TVar *p2){
    return *p1* *p2 + TemplateSum<N-1>(p1+1,p2+1);
	
}

/**
 @brief Implements product of the values into p1 and p2.
 */
template<>
inline double TemplateSum<double,1>(const double *p1, const double *p2){
    return *p1 * *p2;
}

#else // USING_NEW_SKYLMAT

/**
 * @file
 * @brief Contains TPZSkyline class which implements a skyline storage format.
 */

#include "pzmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

#include "pzfmatrix.h"

/**@note  Esta classe gerencia matrizes do tipo SkyLine. Todas
 matrizes SkyLine sao simetricas (e portanto quadradas).
 */

/**
 * @brief Implements a skyline storage format. A Skyline matrix is symmetric so square. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSkylMatrix : public TPZMatrix<TVar>
{
public:
	TPZSkylMatrix() : TPZRegisterClassId(&TPZSkylMatrix::ClassId),TPZMatrix<TVar>(0,0),fElem(0),fStorage(0) { }
	TPZSkylMatrix(const int64_t dim);
	/**
     @brief Construct a skyline matrix of dimension dim
     skyline indicates the minimum row number which will be accessed by each equation
	 */
	TPZSkylMatrix(const int64_t dim ,const TPZVec<int64_t> &skyline);
	TPZSkylMatrix(const TPZSkylMatrix<TVar> &A ) : TPZRegisterClassId(&TPZSkylMatrix::ClassId),TPZMatrix<TVar>(A), fElem(0), fStorage(0)  { Copy(A); }
	
	CLONEDEF(TPZSkylMatrix)
    
	virtual int64_t MemoryFootprint() const {
        return (sizeof(TVar*)*fElem.size() +
                sizeof(TVar)*fStorage.size());
	}	
    
	/**
     @brief modify the skyline of the matrix, throwing away its values
     skyline indicates the minimum row number which will be accessed by each equation
	 */
	void SetSkyline(const TPZVec<int64_t> &skyline);
	
	/**
     @brief return the height of the skyline for a given column
	 */
	int64_t SkyHeight(int64_t col) { return fElem[col+1]-fElem[col] - 1; }
	
	/** @brief Add a skyline matrix B with same structure of this
	 *  It makes this += k * B
	 */
	void AddSameStruct(TPZSkylMatrix<TVar> &B, double k = 1.);
	
	/** @brief declare the object as simetric matrix*/
	virtual int IsSimetric() const {return 1;}
	
    /** @brief destructor of the skyline matrix */
	virtual ~TPZSkylMatrix() { Clear(); }
    
    /**
	 * @brief Updates the values of the matrix based on the values of the matrix
	 */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat);
    
    friend class TPZSkylMatrix<float>;
    friend class TPZSkylMatrix<double>;
    
    /// copy the values from a matrix with a different precision
    template<class TVar2>
    void CopyFrom(TPZSkylMatrix<TVar2> &orig)
    {
        TPZMatrix<TVar>::CopyFrom(orig);
        int64_t nel = orig.fStorage.size();
        fElem.resize(orig.fElem.size());
        fStorage.resize(nel);
        for (int64_t el=0; el<nel; el++) {
            fStorage[el] = orig.fStorage[el];
        }
        int64_t size_el = fElem.size();
        TVar *first = &fStorage[0];
        TVar2 *first_orig = &orig.fStorage[0];
        for (int64_t el=0; el<size_el; el++) {
            fElem[el] = first+(orig.fElem[el]-first_orig);
        }
        
    }
    
    
	int    PutVal(const int64_t row,const int64_t col,const TVar &element );
	const TVar &GetVal(const int64_t row,const int64_t col ) const;
	
	
	TVar &operator()(const int64_t row, const int64_t col);
	virtual TVar &s(const int64_t row, const int64_t col);
	
	
	TVar &operator()(const int64_t row);
	
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha,const TVar beta ,const int opt = 0) const ;
	// Operadores com matrizes SKY LINE.
	TPZSkylMatrix &operator= (const TPZSkylMatrix<TVar> &A );
	//TPZSkylMatrix &operator= (TTempMat<TPZSkylMatrix> A);
	
	TPZSkylMatrix operator+  (const TPZSkylMatrix<TVar> &A ) const;
	TPZSkylMatrix operator-  (const TPZSkylMatrix<TVar> &A ) const;
	
	TPZSkylMatrix &operator+=(const TPZSkylMatrix<TVar> &A );
	TPZSkylMatrix &operator-=(const TPZSkylMatrix<TVar> &A );
	
	// Operadores com valores NUMERICOS.
	TPZSkylMatrix operator*  (const TVar v ) const;
	TPZSkylMatrix &operator*=( TVar value );
	
	TPZSkylMatrix operator-() const;// { return operator*(-1.0); }
	
	// Redimensiona a matriz, mas mantem seus elementos.
	// o segundo parametro � o tamanho das colunas
	int Resize(const int64_t newDim ,const int64_t );
	
	// Redimensiona a matriz e ZERA seus elementos.
	// o segundo parametro � o tamanho das colunas
	int Redim(const int64_t newDim ,const int64_t );
	int Redim(const int64_t newDim) {return Redim(newDim,newDim);}
	
	// Zera os Elementos da matriz
	int Zero();
    
    /**
     * @brief Add a contribution of a stiffness matrix
     * @param elmat Element matrix to be contributed
     * @param sourceindex Contains source indexes on current matrix
     * @param destinationindex Contains destine indexes on current matrix
     */

    void AddKel(TPZFMatrix<TVar>&elmat, TPZVec<int64_t> &sourceindex,  TPZVec<int64_t> &destinationindex);
	
	
	/*** @brief To Solve Linear Equations ***/
	// @{
	virtual void SolveSOR(int64_t &numiterations,const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
						  TPZFMatrix<TVar> *residual,TPZFMatrix<TVar> &scratch,const REAL overrelax, REAL &tol,
						  const int FromCurrent = 0,const int direction = 1) ;
	
	
	int Decompose_Cholesky();  // Faz A = GGt.
	int Decompose_Cholesky_blk(int64_t blk_sz);
    
	int Decompose_LDLt    ();  // Faz A = LDLt.
	int Decompose_Cholesky(std::list<int64_t> &singular);  // Faz A = GGt.
	int Decompose_LDLt    (std::list<int64_t> &singular);  // Faz A = LDLt.
	
	int Subst_Forward  ( TPZFMatrix<TVar> *b ) const;
	int Subst_Backward ( TPZFMatrix<TVar> *b ) const;
	int Subst_LForward ( TPZFMatrix<TVar> *b ) const;
	int Subst_LBackward( TPZFMatrix<TVar> *b ) const;
	int Subst_Diag     ( TPZFMatrix<TVar> *b ) const;
	// @}
	
	//void TestSpeed(int col, int prevcol);
	virtual void AutoFill(int64_t nrow, int64_t ncol, int symmetric) ;
	
	public:
virtual int ClassId() const;

    /**
	 * @brief Unpacks the object structure from a stream of bytes
	 * @param buf The buffer containing the object in a packed form
	 * @param context 
	 */
	virtual void  Read(TPZStream &buf, void *context );
	/**
	 * @brief Packs the object structure in a stream of bytes
	 * @param buf Buffer which will receive the bytes
	 * @param withclassid
	 */
	virtual void Write( TPZStream &buf, int withclassid ) const;
	
    
	virtual std::string ClassName() const   { return( "TPZSkylMatrix"); }
    
    
	int64_t Size(const int64_t column) const 
	{
#ifdef PZDEBUG
		if (column < 0 || column >= this->Rows()) {
			DebugStop();
		}
#endif
		return fElem[column+1]-fElem[column];
	}
	
    int64_t GetNelemts(){
        return fStorage.size();
    }
    
protected:
	
	/**
     @brief This method returns a pointer to the diagonal element of the matrix of the col column
	 */
	TVar *Diag(int64_t col) { return fElem[col];}
	
	void DecomposeColumn(int64_t col, int64_t prevcol);
	void DecomposeColumn(int64_t col, int64_t prevcol, std::list<int64_t> &singular);
	
	void DecomposeColumn2(int64_t col, int64_t prevcol);
private:
	
	// Aloca uma nova coluna. 'fDiag[col].pElem' deve ser NULL.
	
	//static int  Error(const char *msg1,const char* msg2="" );
	int  Clear();
	void Copy (const TPZSkylMatrix<TVar> & );
    
	static int64_t NumElements(const TPZVec<int64_t> &skyline);
	static void InitializeElem(const TPZVec<int64_t> &skyline, TPZVec<TVar> &storage, TPZVec<TVar *> &elem);
	/**
     @brief Computes the highest skyline of both objects
	 */
	static void ComputeMaxSkyline(const TPZSkylMatrix<TVar> &first, const TPZSkylMatrix<TVar> &second, TPZVec<int64_t> &res);
	
	void MigratePages() {
        fElem.MigratePages();
        fStorage.MigratePages();
	}
    
	void ReallocForNuma() {
        fElem.ReallocForNuma();
        TVar* old_start = &fStorage[0];
        fStorage.ReallocForNuma();
        TVar* new_start = &fStorage[0];
        for (int i=0; i<fElem.size(); i++) {
            fElem[i] = new_start + (fElem[i]-old_start);
        }
	}
    
protected:
	/** @brief Storage to keep the first elements to each equation
	 *
	 * fElem is of size number of equation+1
	 * fElem[i] is the first element of the skyline of equation i
	 * fElem[Rows()] is one element beyond the last equation
	 */
	TPZVec<TVar *> fElem;
private:
	/**
     @brief fStorage is a unique vector which contains all the data of the skyline matrix
	 */
	TPZVec<TVar> fStorage;
};

template<class TVar>
int TPZSkylMatrix<TVar>::ClassId() const{
    return Hash("TPZSkylMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}

/** @brief Implements iterative sum over N steps */
template<class TVar,int N>
inline TVar TemplateSum(const TVar *p1, const TVar *p2){
	return *p1* *p2 + TemplateSum<N-1>(p1+1,p2+1);
	
}
/** @brief Implements product of the values into p1 and p2 */
template<>
inline double TemplateSum<double,1>(const double *p1, const double *p2){
	return *p1 * *p2;
}

#endif // USING_NEW_SKYLMAT

#endif // TSKYLMATH

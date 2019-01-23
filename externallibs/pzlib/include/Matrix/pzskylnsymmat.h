//---------------------------------------------------------------------------

#ifndef pzskylnsymmatH
#define pzskylnsymmatH

//
// Author: Nathan Shauer.
//
// File:   pzskylnsymmat
//
// Class:  TPZSkylMatrix
//
// Obs.: This class manages non symmetric skylyne type matrix
//
//
// Versao: 12 / 2011.
//


#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"


/**@note  Esta classe gerencia matrizes do tipo SkyLine nao simetricas.
*/

/**
 * @brief Implements a skyline storage format
 */
template<class TVar=REAL>
class TPZSkylNSymMatrix : public TPZMatrix<TVar>
{
 public:
  TPZSkylNSymMatrix() : TPZRegisterClassId(&TPZSkylNSymMatrix::ClassId),
    TPZMatrix<TVar>(0,0),fElem(0), fElemb(0), fStorage(0), fStorageb(0) { }
    
  TPZSkylNSymMatrix(const int64_t nrow, const int64_t ncol);
  /**
     Construct a skyline matrix of dimension dim
     skyline indicates the minimum row number which will be accessed by each equation
  */
  TPZSkylNSymMatrix(const int64_t dim ,const TPZVec<int64_t> &skyline);
  
  TPZSkylNSymMatrix(const TPZSkylNSymMatrix &A ) : 
    TPZRegisterClassId(&TPZSkylNSymMatrix::ClassId),TPZMatrix<TVar>(A), fElem(0), fElemb(0), fStorage(0), fStorageb(0)  
    { 
        Copy(A); 
    }

  CLONEDEF(TPZSkylNSymMatrix)
  /**
     modify the skyline of the matrix, throwing away its values
     skyline indicates the minimum row number which will be accessed by each equation
  */
  void SetSkyline(const TPZVec<int64_t> &skyline);

  /**
     return the height of the skyline for a given column
  */
  int SkyHeight(int64_t col) { return fElem[col+1]-fElem[col] - 1; }

  /** Add a skyline matrix B with same structure of this
   *  It makes this += k * B
   */
  //void AddSameStruct(TPZSkylNSymMatrix &B, double k = 1.);

  /**declare the object as non-symmetric matrix*/
  virtual int IsSimetric() const {return 0;}

  virtual ~TPZSkylNSymMatrix() { Clear(); }

  int PutVal(const int64_t row,const int64_t col,const TVar &element );

  const TVar &GetVal(const int64_t row,const int64_t col ) const;


  /// Pega o valor na diagonal ou parte de cima da diagonal
  const TVar &GetValSup(const int64_t row,const int64_t col ) const;

  /// Pega o valor abaixo da diagonal (below)
  const TVar &GetValB(const int64_t row,const int64_t col ) const;


  TVar &operator()(const int64_t row, const int64_t col);
  virtual TVar &s(const int64_t row, const int64_t col);


  TVar &operator()(const int64_t row);

  virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
		       const TVar	alpha,const TVar beta ,const int opt = 0) const ;
    
    
    /** @brief Updates the values of the matrix based on the values of the matrix */
    virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat);

  // Operadores com matrizes SKY LINE.
  //TPZSkylNSymMatrix &operator= (const TPZSkylNSymMatrix &A );
  //TPZSkylMatrix &operator= (TTempMat<TPZSkylMatrix> A);

  //TPZSkylNSymMatrix operator+  (const TPZSkylNSymMatrix &A ) const;
  //TPZSkylNSymMatrix operator-  (const TPZSkylNSymMatrix &A ) const;

  //TPZSkylNSymMatrix &operator+=(const TPZSkylNSymMatrix &A );
  //TPZSkylNSymMatrix &operator-=(const TPZSkylNSymMatrix &A );

  // Operadores com valores NUMERICOS.
  //TPZSkylNSymMatrix operator*  (const REAL v ) const;
  //TPZSkylNSymMatrix &operator*=( REAL v );

  //TPZSkylNSymMatrix operator-() const;// { return operator*(-1.0); }

  // Redimensiona a matriz, mas mantem seus elementos.
  // o segundo parametro � o tamanho das colunas
  //int Resize(const int newDim ,const int );

  // Redimensiona a matriz e ZERA seus elementos.
  // o segundo parametro � o tamanho das colunas
  //int Redim(const int newDim ,const int );
  //int Redim(const int newDim) {return Redim(newDim,newDim);}

  // Zera os Elementos da matriz
  //int Zero();


  /*** Resolucao de sistemas ***/

  int Decompose_LU();  // Faz A = LU.

  //virtual void SolveSOR(int &numiterations,const TPZFMatrix &F, TPZFMatrix &result,
	//		TPZFMatrix *residual,TPZFMatrix &scratch,const REAL overrelax, REAL &tol,
	//		const int FromCurrent = 0,const int direction = 1) ;


  //int Subst_Forward  ( TPZFMatrix *b ) const;
  int Subst_Backward ( TPZFMatrix<TVar> *b ) const;
  int Subst_LForward ( TPZFMatrix<TVar> *b ) const;
  //int Subst_LBackward( TPZFMatrix *b ) const;
  //int Subst_Diag     ( TPZFMatrix *b ) const;

  //void TestSpeed(int col, int prevcol);
	
	/**
	 *@brief Return the id of the matrix defined pzmatrixid.h
	 */
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
	virtual void Write( TPZStream &buf, int withclassid );
    
    /** @brief Fill matrix storage with randomic values */
	/** This method use GetVal and PutVal which are implemented by each type matrices */
	void AutoFill(int64_t nrow, int64_t ncol, int symmetric);


 protected:

  /**
     This method returns a pointer to the diagonal element of the matrix of the col column
  */
  TVar *Diag(int64_t col) { return fElem[col];}

  //void DecomposeColumn(int col, int prevcol);
	//void DecomposeColumn(int col, int prevcol, std::list<int64_t> &singular);

  //void DecomposeColumn2(int col, int prevcol);
 private:

  // Aloca uma nova coluna. 'fDiag[col].pElem' deve ser NULL.

//static int  Error(const char *msg1,const char* msg2="" );
  int  Clear();
  void Copy (const TPZSkylNSymMatrix & );
  int Size(const int64_t column) const {return fElem[column+1]-fElem[column];}
  static int64_t NumElements(const TPZVec<int64_t> &skyline);
  static void InitializeElem(const TPZVec<int64_t> &skyline, TPZManVector<TVar> &storage, TPZVec<TVar *> &elem);
  /**
     Computes the highest skyline of both objects
  */
  static void ComputeMaxSkyline(const TPZSkylNSymMatrix &first, const TPZSkylNSymMatrix &second, TPZVec<int> &res);
	
	/** @brief Zeroes the matrix */
	virtual int Zero(){
		fStorage.Fill(0.);
        fStorageb.Fill(0.);
		return 1;
    }


protected:
  /**
     fElem is of size number of equation+1
     fElem[i] is the first element of the skyline of equation i
     fElem[Rows()] is one element beyond the last equation
  */
  TPZVec<TVar *> fElem;

  /**
     fElemb storages skyline Below diagonal
     fElemb is of size number of equation+1
     fElemb[i] is the first element of the skyline of equation i
     fElemb[Rows()] is one element beyond the last equation
  */
  TPZVec<TVar *> fElemb;

private:
  /**
     fStorage is a unique vector which contains all the data of the skyline matrix
  */
  TPZManVector<TVar> fStorage;

  /**
     fStorage is a unique vector which contains all the data of the skyline matrix below diagonal
  */
  TPZManVector<TVar> fStorageb;
};

/*
template<int N>
inline REAL TemplateSum(const REAL *p1, const REAL *p2){
  return *p1* *p2 + TemplateSum<N-1>(p1+1,p2+1);

}


template<>
inline REAL TemplateSum<1>(const REAL *p1, const REAL *p2){
  return *p1 * *p2;
}

*/

//---------------------------------------------------------------------------

template<class TVar>
int TPZSkylNSymMatrix<TVar>::ClassId() const{
    return Hash("TPZSkylNSymMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}
#endif

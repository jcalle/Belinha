#ifndef ELCT2DISCHH
#define ELCT2DISCHH

#include "pzelct2d.h"

class TCompElT2dTGD : public TPZCompElT2d {

  /** Connect index to new internal connect to discontinuous values */
  int fConnect;
  /** To store the total mass of the variables contained in this element*/
  TPZVec<REAL> fMeanValue;
  /** Coefficient of the diffussive term */
  REAL fCFLDiffussion;

 protected:
  /** Set fConnectDisc[i] 1 to indicate the connect i is discontinuous and
      set the correct positions to point to shape functions*/
  void SetConnectDiscontinuous(int i);
  void SetConnectContinuous(int i,TPZFMatrix &values);

 public:

  static TPZCompEl *CreateElDisc(TPZGeoElT2d *gel,TPZCompMesh &mesh,int &index);

  TCompElT2dTGD(TPZCompMesh &mesh,TPZGeoElT2d *ref,int &index);
  ~TCompElT2dTGD();

  int CanBeDiscontinuous() { return 1; }

  /** To set and return the index of the connect i until to internal connect */
  void SetConnectIndex(int i,int connectindex);
  int ConnectIndex(int i);

  /**returns the number of shapefunctions associated with a connect*/
  int NConnectShapeF(int iconnect);
  int NShapeF();
  
  /**return the number of connects along side iside*/
  int NSideConnects(int iside);
  /**returns the local connect number of c along side*/
  int SideConnectLocId(int c, int side);

  int NConnects() { return 8; }

  /**sets the interpolation order of side to order
  This method only updates the datastructure of the element and
  updates the blocksize of the associated connect object */
  void SetSideOrder(int side, int order);

  void SetCoefDiffussion(REAL coef) { fCFLDiffussion = coef; }
	void IncrementCoefDiffussion();

  /** Return 1 whether not exist neighbours or the connect in neighbours are discontinuous*/
  int NeighbourDiscontinuous(int icon);

  int IsConnectContinuous(const int side);

  /** Compute the values of the shape functions over the point */
  void Shape(TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi);
  /** Compute the values of the shape function on side over the point */
  void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi);

  void Print(ostream &out);

  void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);

  void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
  void CalcRhs(TPZElementMatrix &ef);
  void EvaluateError(void (*fp)(TPZVec<REAL>&loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
        REAL &true_error,REAL &L2_error,TPZBlock* /*flux*/,REAL &estimate);
};


/** The class TCompElT2dWI has the abstraction of the computational discontinuous
triangular element enable to has interface over its common side boundary
with other (only one) discontinuous computational element with interface*/
class TCompElT2dWIT : public TCompElT2dTGD {

 public:
  static TPZCompEl *CreateElDiscWI(TPZGeoElT2d *gel,TPZCompMesh &mesh,int &index);

  /** Constructors and destructor*/
  TCompElT2dWIT(TPZCompMesh &mesh,TPZGeoElT2d *ref,int &index);
  ~TCompElT2dWIT();

  int CanToHaveInterface() { return 1; }

  /** Returns index of the interface over the side*/
  int Interface(int side);
  /** Set index of the interface over side */
  void SetInterface(int side, int index);
  /** Creates one interface element over side of the element if it is possible*/
  int CreateInterface(int side,TPZCompMesh &mesh);

  void Print(ostream &out);
  void DeleteInterfaces();
	
 protected:
  /** Vector of interface indexes */
  int fInterface[7];
};

#endif



#ifndef ELCT2DISCHH
#define ELCT2DISCHH

#include "pzelct2d.h"

class TCompElT2dGD : public TPZCompElT2d {

  /** Vector of masks to discontinuous sides */
  int fConnectDisc[7];
  /** Vector of initial positions of the shape functions associated with each connect.
      The internal connect is associated with the first position of the discontinuous connects*/
  int fPositions[8];
  /** Connect index to new internal connect to discontinuous values */
  int fConnect;
  /** To store the total mass of the variables contained in this element*/
  TPZVec<REAL> fMeanValue;

 protected:
  /** Set fConnectDisc[i] 1 to indicate the connect i is discontinuous and
      set the correct positions to point to shape functions*/
  void SetConnectDiscontinuous(int i);
  void SetConnectContinuous(int i,TPZFMatrix<STATE> &values);
  /**To redimensioning the block of the icon connect that is being actived
     in continuous and copying the adequated values from discontinuous connect*/
  int StablizingConnectContinuous(int icon,TPZStack<TPZCompElSide> &elvec);

 public:

  static TPZCompEl *CreateElDisc(TPZGeoElT2d *gel,TPZCompMesh &mesh,int64_t &index);

  TCompElT2dGD(TPZCompMesh &mesh,TPZGeoElT2d *ref,int64_t &index);
  /**Constructor to discontinuous element but with continuous connect cont[i]*/
  TCompElT2dGD(TPZCompMesh &mesh,TPZGeoElT2d *ref,int64_t &index,TPZVec<int> &cont);
  ~TCompElT2dGD();

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

  /**Ordening the positions into matrix phi of the shape functions to connects*/
  void SetPositions();

  /** To transform one continuous connect in discontinuous or viceverse */
  void MakeDiscontinuous(TPZVec<int> &disc);
  void MakeConnectDiscontinuous(int icon);
  void MakeConnectContinuous(int icon);
  /** Return 1 whether not exist neighbours or the connect in neighbours are discontinuous*/
  int NeighbourDiscontinuous(int icon);

  int IsConnectContinuous(const int side);

  /** Compute the values of the shape functions over the point */
  void Shape(TPZVec<REAL> &point,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi);
  /** Compute the values of the shape function on side over the point */
  void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi);

  /** Ordening the values of the shape functions in correpondence to connects */
  void OrdeningPhi(TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi);

  void Print(std::ostream &out);

  void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);

  void EvaluateError(void (*fp)(TPZVec<REAL>&loc,TPZVec<REAL> &val,TPZFMatrix<STATE> &deriv),
        REAL &true_error,REAL &L2_error,TPZBlock<STATE>* /*flux*/,REAL &estimate);
};


/** The class TCompElT2dWI has the abstraction of the computational discontinuous
triangular element enable to has interface over its common side boundary
with other (only one) discontinuous computational element with interface*/
class TCompElT2dWI : public TCompElT2dGD {

 public:
  static TPZCompEl *CreateElDiscWI(TPZGeoElT2d *gel,TPZCompMesh &mesh,int64_t &index);

  /** Constructors and destructor*/
  TCompElT2dWI(TPZCompMesh &mesh,TPZGeoElT2d *ref,int64_t &index);
  TCompElT2dWI(TPZCompMesh &mesh,TPZGeoElT2d *ref,int64_t &index,TPZVec<int> &continuous);
  ~TCompElT2dWI();

  int CanToHaveInterface() { return 1; }
  void Divide(int64_t index,TPZVec<int64_t> &pv,int interpolatesolution = 0);
  /** Returns index of the interface over the side*/
  int Interface(int side);
  /** Set index of the interface over side */
  void SetInterface(int side, int64_t index);
  /** Creates one interface element over side of the element if it is possible*/
  int CreateInterface(int side,TPZCompMesh &mesh);

  void Print(std::ostream &out);
  void DeleteInterfaces();
	
 protected:
  /** Vector of interface indexes */
  int64_t fInterface[7];
};

#endif



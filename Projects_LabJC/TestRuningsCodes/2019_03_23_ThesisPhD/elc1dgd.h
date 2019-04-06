#ifndef PZELC1DISC_H
#define PZELC1DISC_H

#include "pzelc1d.h"

class TCompEl1dGD : public TPZCompEl1d {

  /** Vector of masks to discontinuous connects */
  int fConnectDisc[3];
  /** Vector of initial positions of the shape functions associated with each connect
      The internal connect has referenced with the first position of the discontinuous connects */
  int fPositions[4];
  /** Connect index to new internal connect to discontinuous values */
  int fConnect;

 protected:
  /** Set fConnectDisc[i] 1 to indicate the connect i is discontinuous and
      set the correct positions to point to shape functions*/
  void SetConnectDiscontinuous(int i);
  void SetConnectContinuous(int i,TPZFMatrix<STATE> &values);
  /**To redimensioning the block of the icon connect that is being actived
     in continuous and copying the adequated values from discontinuous connect*/
  int StablizingConnectContinuous(int icon,TPZStack<TPZCompElSide> &elvec);

 public:

  static TPZCompEl *CreateElDisc(TPZGeoEl1d *gel,TPZCompMesh &mesh, int64_t &index);

  /** Constructors and destructor*/
  TCompEl1dGD(TPZCompMesh &mesh,TPZGeoEl1d *ref, int64_t &index);
  TCompEl1dGD(TPZCompMesh &mesh,TPZGeoEl1d *ref, int64_t &index,TPZVec<int> &continuous);
  ~TCompEl1dGD();

  int CanBeDiscontinuous() { return 1; }

  /**To set and return the index of the connect i until to internal connect */
  void SetConnectIndex(int i, int64_t connectindex);
  int64_t ConnectIndex(int i);

  /**Returns the number of shapefunctions associated with a connect*/
  virtual int NConnectShapeF(int iconnect);
  int NShapeF();

  /**return the number of connects along side iside*/
  virtual int NSideConnects(int iside);
  /**returns the local connect number of c along side*/
  virtual int SideConnectLocId(int c, int side);

  int NConnects() { return 4; }

  /**Ordening the positions into matrix phi of the shape functions to connects*/
  void SetPositions();
  
  /**sets the interpolation order of side to order
  This method only updates the datastructure of the element and
  updates the blocksize of the associated connect object */
  void SetSideOrder(int side, int order);

  /** To transform one continuous connect in discontinuous or viceverse */
  void MakeDiscontinuous(TPZVec<int> &disc);
  /** Set a connect icon continuous*/
  void MakeConnectContinuous(int icon);
  /** Set a connect icon discontinuous*/
  void MakeConnectDiscontinuous(int icon);
  /** Return 1 whether not exist neighbours or the connect in neighbours are discontinuous*/
  int NeighbourDiscontinuous(int icon);

  int IsConnectContinuous(const int i);

  /** Compute the values of the shape functions over the point */
  void Shape(TPZVec<REAL> &point,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi);
  /** Compute the values of the shape function on side over the point */
  void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi);

  /** Ordening the values of the shape functions in correpondence to connects */
  void OrdeningPhi(TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphix);

  virtual void Print(std::ostream &out = std::cout);

  void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension);

/*  void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
  void CalcRhs(TPZElementMatrix &ef);
	*/
  void EvaluateError(void (*fp)(TPZVec<REAL>&loc,TPZVec<REAL> &val,TPZFMatrix<STATE> &deriv),
        REAL &true_error,REAL &L2_error,TPZBlock<STATE>* /*flux*/,REAL &estimate);

};

/** The class TCompEl1dWI has the abstraction of the computational discontinuous
one-dimensional element enable to has interface over its common side boundary
with other (only one) discontinuous computational element with interface*/
class TCompEl1dWI : public TCompEl1dGD {

 public:
  static TPZCompEl *CreateElDiscWI(TPZGeoEl1d *gel,TPZCompMesh &mesh, int64_t &index);

  /** Constructors and destructor*/
  TCompEl1dWI(TPZCompMesh &mesh,TPZGeoEl1d *ref, int64_t &index);
  TCompEl1dWI(TPZCompMesh &mesh,TPZGeoEl1d *ref,
	  int64_t &index,TPZVec<int> &continuous);
  ~TCompEl1dWI();
  void Divide(int64_t index,TPZVec<int64_t> &pv,int interpolatesolution = 0);
  int CanToHaveInterface() { return 1; }

  /** Returns index of the interface over the side*/
  int Interface(int side);
  /** Set index of the interface over side */
  void SetInterface(int side, int64_t index);
  /** Creates one interface element over side of the element if it is possible*/
  int64_t CreateInterface(int side,TPZCompMesh &mesh);
  void DeleteInterfaces();
	
  void Print(std::ostream &out= std::cout);

protected:
  /** Vector of interface indexes */
	int64_t fInterface[3];
};

#endif

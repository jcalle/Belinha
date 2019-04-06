//HEADER FILE FOR CLASS ELEM1D

#ifndef ELCALC1DHPP
#define ELCALC1DHPP

#include <iostream>
#include <string.h>
#include "pzintel.h"
#include "pzquad.h"
#include "pzshapelinear.h"
#include "pzelctemp.h"

class TPZMaterial;
class TPZIntRule;
struct TPZElementMatrix;
//class TPZMatrix;
//class TPZFMatrix;
//class TPZBlock;
class TPZConnect;
class TPZBndCond;
class TPZGeoEl1d;
class TPZGraphMesh;

/**
TPZCompEl1d implements a one-dimensional computational element
*/
class TPZCompEl1d : public TPZIntelGen<pzshape::TPZShapeLinear> {

 protected:
  /** Integration rule for the element*/
  TPZInt1d 	fIntRule;
  /** Vector of connect indices*/
  int64_t fConnectIndexes[3];
  /** Side interpolation order*/
  int fSideOrder;
  /** Side interpolation order liked by the user*/
  int fPreferredSideOrder;

 public:

  TPZCompEl1d(TPZCompMesh &mesh, TPZGeoEl1d *ref, int64_t &index);

  TPZCompEl1d(TPZCompMesh &mesh, TPZGeoEl1d *ref, int64_t &index, int NoConnectCreation);

  virtual ~TPZCompEl1d();

  MElementType Type() { return EOned;}

  void VarRange(int var, double &min, double &max) {
    min = 0.;
    max = 0.;
  }

  /** Sets the connect index of connect i to the element*/
  virtual void SetConnectIndex(int i, int64_t connectindex);

  /**return the index of the ith connectivity of the element*/
  virtual int64_t ConnectIndex(int i);

  /**return the number of connects along side iside*/
  virtual int NSideConnects(int iside);

  /**return the number of corner connects of the element*/
  virtual int NCornerConnects() { return 2; }

  /**return the number of connect objects of the element*/
  virtual int NConnects() { return 3; }

  /**returns the number of shapefunctions associated with a connect*/
  virtual int NConnectShapeF(int iconnect);

  /**returns the local connect number of c along side*/
  virtual int SideConnectLocId(int c, int side);

  /**returns the interpolation dimension of the element*/
  virtual int Dimension() { return 1; }

  /**computes the shape function set at the point x*/
  virtual void Shape(TPZVec<REAL> &x,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi);

  /**returns a reference to an integrationrule suitable for integrating
     the interior of the element*/
  virtual TPZIntPoints &GetIntegrationRule() { return fIntRule;}

  /**Sets the interpolation order for the interior of the element*/
  virtual void SetInterpolationOrder(TPZVec<int> &ord);

  /**Identifies the interpolation order on the interior of the element*/
  virtual void GetInterpolationOrder(TPZVec<int> &ord);

  void SetIntegrationRule(TPZInt1d &int1d) { fIntRule = int1d;}

  void SetIntegrationRule(int order);//Cedric 16/03/99 
  
  /**return the preferred order of the polynomial along side iside*/
  virtual int PreferredSideOrder(int iside);

  /**Sets the preferred interpolation order along a side
  This method only updates the datastructure of the element
  In order to change the interpolation order of an element, use the method PRefine*/
  virtual void SetPreferredSideOrder(int side, int order);

  /**sets the interpolation order of side to order*/
  virtual void SetSideOrder(int side, int order);

  /**returns the actual interpolation order of the polynomial along the side*/
  virtual int SideOrder(int side);

  /**transform a point in the parameter space of the side into a point in the space
     of the master element*/
  virtual void SideParameterToElement(int side, TPZVec<REAL> &par, TPZVec<REAL> &point);

  /**transform a point in the parameter space of the master element into a point in the
     space of the side*/
  virtual void ElementToSideParameter(int side, TPZVec<REAL> &point, TPZVec<REAL> &par);

  /**compute the values of the shape function of the side*/
  virtual void SideShapeFunction(int side, TPZVec<REAL> &point, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &dphi);

  /**allocates dynamically an integration rule adequate for the side
     the caller to the method needs to call the delete method!*/
  virtual TPZIntPoints *CreateSideIntegrationRule(int side);

  void Load() {
	  std::cout << "Precisa implementar.\n";
  }

  void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);

	/** Jorge 09/06/2001
	 * Returns the transformation which transform a point from the side to the interior of the element
	 */
	TPZTransform<STATE> TransformSideToElement(int side);
};

inline TPZTransform<STATE> TPZCompEl1d::TransformSideToElement(int side) {
	return pzshape::TPZShapeLinear::TransformSideToElement(side);
}

#endif


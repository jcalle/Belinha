//HEADER FILE FOR CLASS ELCQ2D

#ifndef PZELCQ2DHPP
#define PZELCQ2DHPP

#include <iostream>
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzshtmat.h"
#include "pzfmatrix.h"
#include "pzshapequad.h"
#include "pzelctemp.h"

struct TPZElementMatrix;
class TPZBndCond;
class TPZConnect;
class TPZMaterial;
class TPZGeoEl;
class TPZGeoElQ2d;
template<class T>
class TPZVec;
class TPZCompMesh;
class TPZIntRule;
class TPZIntQuad;
class TPZShortMatrix;

/** Header file for the computational quadrilateral element class*/
class TPZCompElQ2d : public TPZIntelGen<pzshape::TPZShapeQuad> {

  TPZIntQuad  fIntRule;
  int fConnectIndexes[9];
  int fSideOrder[5];//não nos cantos
  int fPreferredSideOrder[5];

 public:

  TPZCompElQ2d(TPZCompMesh &mesh,TPZGeoElQ2d *ref, int64_t &index);
  TPZCompElQ2d(TPZCompMesh &mesh,TPZGeoElQ2d *ref, int64_t &index,int noconnects);
  ~TPZCompElQ2d() {
  	if(Reference()) {
      if(Reference()->Reference()) {
         RemoveSideRestraintsII(EDelete);
      }
      Reference()->ResetReference();
   }
}

  MElementType Type() { return EQuadrilateral; }

  void VarRange(int var,double &min,double &max);

  int NConnects() { return 9; }

  void SetConnectIndex(int i,int connectindex);

  int ConnectIndex(int i);

  int NConnectShapeF(int connect);

  int Dimension() { return 2; }

  int NCornerConnects() { return 4; }

  /**return the number of dof nodes along side iside*/
  virtual int NSideConnects(int iside);

  virtual int SideConnectLocId(int node, int side);

  TPZIntPoints &GetIntegrationRule() { return fIntRule; }

  /**Sets the interpolation order for the interior of the element*/
  virtual void SetInterpolationOrder(TPZVec<int> &ord);

  /**Identifies the interpolation order on the interior of the element*/
  virtual void GetInterpolationOrder(TPZVec<int> &ord);

  void SetIntegrationRule(TPZIntQuad &intquad) { fIntRule = intquad;}

  void SetIntegrationRule(int ord);

  /**allocates dynamically an integration rule adequate for the side
     the caller to the method needs to call the delete method!*/
  virtual TPZIntPoints *CreateSideIntegrationRule(int side);

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
  virtual void SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);

  /**transform a point in the parameter space of the master element into a point in the
     space of the side*/
  virtual void ElementToSideParameter(int side, TPZVec<REAL> &point, TPZVec<REAL> &par);

  void CornerShape(TPZVec<REAL> &pt, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &dphi);

  /**compute the values of the shape function of the side*/
  virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi);

  void Shape(TPZVec<REAL> &pt, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &dphi);

  void Load();

  void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);

virtual void Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol);

  /**
   * Create a computational mesh which can be used for qualifying purposes
   */
  static TPZCompMesh *CreateMesh();

	/** Jorge 09/06/2001
	 * Returns the transformation which transform a point from the side to the interior of the element
	 */
	TPZTransform<STATE> TransformSideToElement(int side);
};

inline TPZTransform<STATE> TPZCompElQ2d::TransformSideToElement(int side) {
	return pzshape::TPZShapeQuad::TransformSideToElement(side);
}


#endif


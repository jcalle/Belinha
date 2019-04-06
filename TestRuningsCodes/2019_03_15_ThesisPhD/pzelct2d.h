//HEADER FILE FOR CLASS ELEM1D

#ifndef PZCOMPELT2DHPP
#define PZCOMPELT2DHPP

#include <iostream>
#include <string.h>
#include "pzintel.h"
#include "pzquad.h"
#include "pzerror.h"
#include "pzgeoel.h"
#include "pzelctemp.h"
#include "pzshapetriang.h"

struct 	TPZElementMatrix;
class 	TPZConnect;
class    TPZBndCond;
class		TPZMat2dLin;
class		TPZGeoElT2d;
class		TPZCompMesh;
class    TPZGraphGrid;


class TPZCompElT2d : public TPZIntelGen<pzshape::TPZShapeTriang> { 	// header file for the two dimensional triangular
  // element class
 protected:
  TPZIntTriang 	fIntRule;	// integration rule for the element
  int fConnectIndexes[7];
  int	fSideOrder[4];		// interpolation order along the three sides
  int fPreferredSideOrder[4];

 public:

  TPZCompElT2d(TPZCompMesh &mesh,TPZGeoElT2d *ref, int64_t &index);
  TPZCompElT2d(TPZCompMesh &mesh,TPZGeoElT2d *ref, int64_t &index,int NoConnectCreation);
  ~TPZCompElT2d() {
  	if(Reference()) {
      if(Reference()->Reference()) {
         RemoveSideRestraintsII(EDelete);
      }
      Reference()->ResetReference();
   }
  }

  MElementType Type() { return ETriangle; }
  void VarRange(int var,double &min,double &max);

  int NConnects() { return 7; }

  void SetConnectIndex(int i,int connectindex);

  int ConnectIndex(int i);

  int NConnectShapeF(int connect);

  int Dimension() { return 2; }

  int NCornerConnects() { return 3; }

  /**return the number of dof nodes along side iside*/
  virtual int NSideConnects(int iside);

  virtual int SideConnectLocId(int node, int side);

  /**returns a reference to an integrationrule suitable for integrating
     the interior of the element*/
  virtual TPZIntPoints &GetIntegrationRule() { return fIntRule;}

  /**Sets the interpolation order for the interior of the element*/
  virtual void SetInterpolationOrder(TPZVec<int> &ord);

  /**Identifies the interpolation order on the interior of the element*/
  virtual void GetInterpolationOrder(TPZVec<int> &ord);

  void SetIntegrationRule(TPZIntTriang &inttriang) { fIntRule = inttriang;}

  void SetIntegrationRule(int order);//Cedric

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
  virtual void SideShapeFunction(int side, TPZVec<REAL> &point, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &dphi);

  void Shape(TPZVec<REAL> &pt, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &dphi);

  void Load();

  void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);

  /**
   * Create a computational mesh which can be used for qualifying purposes
   */
  static TPZCompMesh *CreateMesh();

	/** Jorge 09/06/2001
	 * Returns the transformation which transform a point from the side to the interior of the element
	 */
	TPZTransform<STATE> TransformSideToElement(int side);
};

inline TPZTransform<STATE> TPZCompElT2d::TransformSideToElement(int side) {
	return pzshape::TPZShapeTriang::TransformSideToElement(side);
}

inline void TPZCompElT2d::VarRange(int var,double &min,double &max) {
	PZError << "TCompElT2d::VarRange is not defined.\n";
   if(var>-1) max = min = 0.;
}

inline void TPZCompElT2d::Load() {
	PZError << "TCompElT2d::Load is called.\n";
}


#endif


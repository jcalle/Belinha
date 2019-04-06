////////////////////////////////////////////////////////////////////////////////
//  Elemento tridimensional cubico
////////////////////////////////////////////////////////////////////////////////
#ifndef ELCC3DHPP
#define ELCC3DHPP

#include <iostream>
#include "pzquad.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzshapecube.h"
#include "pzelctemp.h"

template<class T> class TPZFMatrix;
struct TPZElementMatrix;
class TPZShortMatrix;
class TPZBndCond;
class TPZConnect;
class TPZMaterial;
class TPZGeoEl;
class TPZGeoElC3d;
template<class T> class TPZMatrix;
class TPZMat3dLin;
class TPZCompMesh;
class TPZIntRule;
class TPZIntQuad;
class TPZIntCube3D;
class TPZGraphMesh;

class TPZCompElC3d : public TPZIntelGen<pzshape::TPZShapeCube> {	// header file for the computational element class

  TPZIntCube3D fIntRule;
  int fSideOrder[19];	//ordem das funcoes de interpolacao sobre as 6 caras e 12 lados e interior
  int fPreferredSideOrder[19];
  int fConnectIndexes[27];

  static TPZShortMatrix gEqNumbers;

protected:

  void EqNumber(TPZShortMatrix &mat);	// returns the equation number associated
				//	with the shapefunctions in both directions

  virtual void Shape(TPZVec<REAL> &pt, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &dphi);

  virtual void NormalVector(int side,TPZVec<REAL> &int_point,
			    TPZVec<REAL> &normal,TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &norm_r3);

  void Shape3dNew(TPZVec<REAL> &pt, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &dphi);

  void CreateSideNodes();

 public:

   TPZCompElC3d(TPZCompMesh &mesh,TPZGeoElC3d *ref, int64_t &index);
   TPZCompElC3d(TPZCompMesh &mesh,TPZGeoElC3d *ref, int64_t &index,int noconnects);
  ~TPZCompElC3d() {
      if(Reference()) {
         if(Reference()->Reference()) {
            RemoveSideRestraintsII(EDelete);
         }
         Reference()->ResetReference();
      }
   }

  MElementType Type() { return ECube; }

  virtual void Print(std::ostream & out = std::cout);

  virtual int PreferredSideOrder(int iside);

  virtual int SideOrder(int side);

  void FaceOrder(int face,int &ord1,int &ord2);//Cedric 19/04/98

  virtual void SetSideOrder(int side, int order);

  int NConnects() { return 27; }//estas podem ser
  int NCornerConnects() { return 8; }//substituidas por esta
  int  NSideConnects(int side);

  int Dimension() { return 3;}

  int  SideConnectLocId(int c,int side);

/**creates a cube graphics element and corresponding nodes
   and appends the objects to the maps of the grafgrid*/
  //virtual void CreateGraphEl(TPZGraphMesh &graphmesh);

  // level comparison of elements element and neighbour
  virtual short CompareLevel(TPZCompElC3d &element,TPZCompElC3d &neighbour);

  void FaceParameters(int face,TPZVec<REAL> &pt,REAL &pr1,REAL &pr2,TPZVec<REAL> &coef1,TPZVec<REAL> &coef2);

  void FaceIdsCube(int face,TPZVec<int> &ids,TPZVec<int> &id,int &id0,int &id1);

  /**returns a reference to an integrationrule suitable for integrating
     the interior of the element*/
  TPZIntPoints &GetIntegrationRule() { return fIntRule;}

  void SetIntegrationRule(TPZIntCube3D &intcube) { fIntRule = intcube;}

/**scientific routines*/
/**projects the flux function
   on the finite element space*/
  //virtual void ProjectFlux(TPZElementMatrix &ek, TPZElementMatrix &ef);
void NodeFaceIds(TPZVec<int> &ids,int face);
void VarRange(int var,REAL &min,REAL &max);
void Load();
void SetConnectIndex(int i,int connectindex);
int  ConnectIndex(int i);
int  NConnectShapeF(int side);
void SetInterpolationOrder(TPZVec<int> &ord);
void GetInterpolationOrder(TPZVec<int> &ord);
TPZIntPoints *CreateSideIntegrationRule(int side);
void SetPreferredSideOrder(int side, int order);
void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi);
void SetIntegrationRule(int order);
/**creates a cube graphics element and corresponding nodes
   and appends the objects to the maps of the grafgrid*/
void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension);
static int ShapeFaceId[6][2];
static int FaceConnectLocId[6][9];
static REAL transfor[8][2][2];
static REAL facetransfor[6][2][3];
static int In[18][3];
static int Or[6][2];
static int FaceSons[6][4];
static int FaceNodes[6][4];
static int SideNodes[12][2];
static int MidSideNodes[19][2];
static int FaceSides[6][4];
static REAL MidCoord[19][3];
static REAL MasterCoord[8][3];
static int CornerSons[8][8];
static int InNeigh[8][19][3];

	/** Jorge 09/06/2001
	 * Returns the transformation which transform a point from the side to the interior of the element
	 */
	TPZTransform<REAL> TransformSideToElement(int side);
};

inline TPZTransform<REAL> TPZCompElC3d::TransformSideToElement(int side) {
	return pzshape::TPZShapeCube::TransformSideToElement(side);
}

#endif



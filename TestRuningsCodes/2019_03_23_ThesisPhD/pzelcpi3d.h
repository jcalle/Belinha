////////////////////////////////////////////////////////////////////////////////
//  Elemento tridimensional pir�mide
////////////////////////////////////////////////////////////////////////////////
#ifndef ELCPI3DHPP
#define ELCPI3DHPP

#include <iostream>
#include "pzquad.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzshapepiram.h"
#include "pzelctemp.h"

struct TPZElementMatrix;
class TPZShortMatrix;
class TPZBndCond;
class TPZConnect;
class TPZMaterial;
class TPZGeoEl;
class TPZGeoElT3d;
class TPZGeoElPi3d;
class TPZMat3dLin;
class TPZCompMesh;
class TPZIntRule;
class TPZIntPyram3D;
class TPZGraphMesh;

class TPZCompElPi3d : public TPZIntelGen<pzshape::TPZShapePiram> {	// header file for the computational element class

  TPZIntPyram3D fIntRule;
  int fSideOrder[14];	//ordem das funcoes de interpolacao sobre as 6 faces e 12 lados
  int fPreferredSideOrder[14];//lados + faces + centro (- cantos)
  int fConnectIndexes[19];
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

   TPZCompElPi3d(TPZCompMesh &mesh,TPZGeoElPi3d *ref,int64_t &index);
   TPZCompElPi3d(TPZCompMesh &mesh,TPZGeoElPi3d *ref,int64_t &index,int noconnects);
  ~TPZCompElPi3d() {
      if(Reference()) {
         if(Reference()->Reference()) {
            RemoveSideRestraintsII(EDelete);
         }
         Reference()->ResetReference();
      }
   }

  MElementType Type() { return EPiramide; }

  virtual void Print(std::ostream & out = std::cout);

  virtual int PreferredSideOrder(int iside);

  virtual int SideOrder(int side);

  void FaceOrder(int face,int &ord1,int &ord2);//Cedric 01/06/98

  virtual void SetSideOrder(int side, int order);

  int NConnects() { return 19; }//estas podem ser

  int NCornerConnects() { return 5; }//substituidas por esta

  int NSideConnects(int side);

  int Dimension() { return 3;}

  int  SideConnectLocId(int c,int side);

/**creates a cube graphics element and corresponding nodes
   and appends the objects to the maps of the grafgrid*/
  //virtual void CreateGraphEl(TPZGraphMesh &graphmesh);

  // level comparison of elements element and neighbour
  virtual short CompareLevel(TPZCompElPi3d &element,TPZCompElPi3d &neighbour);

  void FaceParameters(int face,TPZVec<REAL> &pt,REAL &pr1,REAL &pr2,TPZVec<REAL> &coef1,TPZVec<REAL> &coef2);

  void FaceIdsPiram(int face,TPZVec<int> &ids,TPZVec<int> &id,int &id0,int &id1);
  /**returns a reference to an integrationrule suitable for integrating
     the interior of the element*/
  TPZIntPoints &GetIntegrationRule() { return fIntRule;}
  void SetIntegrationRule(TPZIntPyram3D &intpiram) { fIntRule = intpiram;}
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
static int ShapeFaceId[5][4];
static int FaceConnectLocId[5][9];
static int FaceSons[5][4];
static int FaceNodes[5][4];
static int SideNodes[8][2];
static int MidSideNodes[9][2];
static int FaceSides[5][4];
static REAL MidCoord[9][3];
static REAL MasterCoord[5][3];
static int CornerSons[10][5];
static int InNeigh[10][18][3];
static int FaceInRib[5][4];
static int MiddleFace[4];

	/** Jorge 09/06/2001
	 * Returns the transformation which transform a point from the side to the interior of the element
	 */
	TPZTransform<STATE> TransformSideToElement(int side);
};

inline TPZTransform<STATE> TPZCompElPi3d::TransformSideToElement(int side) {
	return pzshape::TPZShapePiram::TransformSideToElement(side);
}

#endif



////////////////////////////////////////////////////////////////////////////////
//  Elemento tridimensional tetraedro
////////////////////////////////////////////////////////////////////////////////
#ifndef ELCT3DHPP
#define ELCT3DHPP

#include <iostream>
#include "pzquad.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzshapetetra.h"
#include "pzelctemp.h"

struct TPZElementMatrix;
class TPZShortMatrix;
class TPZBndCond;
class TPZConnect;
class TPZMaterial;
class TPZGeoEl;
class TPZGeoElT3d;
class TPZMat3dLin;
class TPZCompMesh;
class TPZIntRule;
class TPZIntTetra3D;
class TPZGraphMesh;


class TPZCompElT3d : public TPZIntelGen<pzshape::TPZShapeTetra> {	// header file for the computational element class

  TPZIntTetra3D fIntRule;
  int fSideOrder[11];	//ordem das funcoes de interpolacao sobre as 6 faces e 12 lados
  int fPreferredSideOrder[11];//lados + faces + centro
  int fConnectIndexes[15];
  static TPZShortMatrix gEqNumbers;

protected:

  void EqNumber(TPZShortMatrix &mat);	// returns the equation number associated
				//	with the shapefunctions in both directions

  virtual void Shape(TPZVec<REAL> &pt, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &dphi);

  virtual void NormalVector(int side,TPZVec<REAL> &int_point,
			    TPZVec<REAL> &normal,TPZFMatrix<STATE> &axes, TPZFMatrix<STATE> &norm_r3);

  void Shape3dNew(TPZVec<REAL> &pt, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &dphi);

  void CreateSideNodes();

 public:

   TPZCompElT3d(TPZCompMesh &mesh,TPZGeoElT3d *ref, int64_t &index);
   TPZCompElT3d(TPZCompMesh &mesh,TPZGeoElT3d *ref, int64_t &index,int noconnects);
  ~TPZCompElT3d() {
      if(Reference()) {
         if(Reference()->Reference()) {
            RemoveSideRestraintsII(EDelete);
         }
         Reference()->ResetReference();
      }
   }

  MElementType Type() { return ETetraedro; }

  virtual void Print(std::ostream & out = std::cout);

  virtual int PreferredSideOrder(int iside);

  virtual int SideOrder(int side);

  void FaceOrder(int face,int &ord1,int &ord2);//Cedric 01/06/98

  virtual void SetSideOrder(int side, int order);

  int NConnects() { return 15; }//estas podem ser

  int NCornerConnects() { return 4; }//substituidas por esta

  int NSideConnects(int side);

  int Dimension() { return 3;}

  int  SideConnectLocId(int c,int side);

  /**creates a cube graphics element and corresponding nodes
   and appends the objects to the maps of the grafgrid*/
  //virtual void CreateGraphEl(TPZGraphMesh &graphmesh);

  // Divide the computational element
  virtual short CompareLevel(TPZCompElT3d &element,TPZCompElT3d &neighbour);
  //REAL EstimateError(TFMatrix & TenMed,REAL & NormaSigEF2);// calcula a tensão media nos nós da malha computacional
  //void MedStressEF(DoubleAVec & point,DoubleAVec & sigmaEF);
  void FaceParameters(int face,TPZVec<REAL> &pt,REAL &pr1,REAL &pr2,TPZVec<REAL> &coef1,TPZVec<REAL> &coef2);

  /**returns a reference to an integrationrule suitable for integrating
     the interior of the element*/
  TPZIntPoints &GetIntegrationRule() { return fIntRule;}
  void SetIntegrationRule(TPZIntTetra3D &inttetra) { fIntRule = inttetra;}
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
static int ShapeFaceId[4][3];
static int FaceConnectLocId[4][7];
static int FaceSons[4][4];
static int FaceNodes[4][3];
static int SideNodes[6][2];
static int MidSideNodes[6][2];
static int FaceSides[4][3];
static REAL MidCoord[6][3];
static REAL MasterCoord[4][3];
static int CornerSons[6][5];
static int InNeigh[6][16][3];
static int MiddleFace[4];

	/** Jorge 09/06/2001
	 * Returns the transformation which transform a point from the side to the interior of the element
	 */
	TPZTransform<STATE> TransformSideToElement(int side);
};

inline TPZTransform<STATE> TPZCompElT3d::TransformSideToElement(int side) {
	return pzshape::TPZShapeTetra::TransformSideToElement(side);
}

#endif




#include "pzelcq2d.h"
#include "pzfmatrix.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzelgq2d.h"
#include "pzcmesh.h"
#include "TPZMaterial.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzgraphel.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzgraphelq2dd.h"
#include "pzbndcond.h"
#include "pzgeoelbc.h"

#include "pzmat2dlin.h"
#include "pzelasmat.h"
//#include "../graf/grafgrid.h"
//#include "../graf/grafel.h"
//#include "../graf/grafnode.h"


TPZCompElQ2d::TPZCompElQ2d(TPZCompMesh &mesh,TPZGeoElQ2d *ref, int64_t &index) :
	TPZIntelGen<pzshape::TPZShapeQuad>(mesh,ref,index), fIntRule(2,2) {
  int i;
  for(i=0;i<5;i++) {
    fSideOrder[i] = mesh.GetDefaultOrder();
    fPreferredSideOrder[i] = mesh.GetDefaultOrder();
  }
//  RemoveSideRestraintsII(EInsert);
  ref->SetReference(this);
  for(i=0;i<4;i++) {
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
  }
  for(;i<9;i++) {
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
	 IdentifySideOrder(i);
  }

  TPZVec<int> order(2,2*fSideOrder[8]+2);
  //  TPZVec<int> order(2,1);
  fIntRule.SetOrder(order);
}

TPZCompElQ2d::TPZCompElQ2d(TPZCompMesh &mesh,TPZGeoElQ2d *ref, int64_t &index,int /*noconnects*/) :
	TPZIntelGen<pzshape::TPZShapeQuad>(mesh,ref,index), fIntRule(2,2) {
  int i;
  for(i=0;i<5;i++) {
    fSideOrder[i] = mesh.GetDefaultOrder();
    fPreferredSideOrder[i] = mesh.GetDefaultOrder();
  }
//  RemoveSideRestraintsII(EInsert);
}

void TPZCompElQ2d::VarRange(int var,double &min,double &max) {
	PZError << "TPZCompElQ2d::VarRange is not defined.\n";
   if(var>-1) max = min = 0.;
}

void TPZCompElQ2d::Load() {
	PZError << "TPZCompElQ2d::Load is called.\n";
}

void TPZCompElQ2d::SetConnectIndex(int i,int connectindex) {
  if(i>-1 && i<9) fConnectIndexes[i] = connectindex;
  else {
    PZError << "TPZCompElQ2d::SetConnectIndex. Bad parameter i.\n";
    PZError.flush();
  }
}

int TPZCompElQ2d::ConnectIndex(int i) {
	if(i>8 || i<0) {
   	PZError << "TCompElT2d::ConnectIndex. Bad parameter i.\n";
      return -1;
   }
   return fConnectIndexes[i];
}

int TPZCompElQ2d::NConnectShapeF(int side) {
  switch(side) {
  case 0:
  case 1:
  case 2:
  case 3:
    return 1;
  case 4:
  case 5:
  case 6:
  case 7:
	  return fSideOrder[side-4]-1;
  case 8:
     return (fSideOrder[side-4]-1)*(fSideOrder[side-4]-1);//Cedric Modified 24/02/99
  default:
    PZError << "TPZCompElQ2d::NConnectShapeF, bad parameter iconnect " << side << std::endl;
    return 0;
  }
}

int TPZCompElQ2d::NSideConnects(int side) {
  if(side<0 || side>8) {
    PZError << "TPZCompElT2d::NSideConnects. Bad parameter i.\n";
    return 0;
  }
  if(side<4) return 1;
  if(side<8) return 3;
  return 9;//Cedric
}

/**It do not verify the values of the c*/
int TPZCompElQ2d::SideConnectLocId(int c,int side) {
  switch(side) {
  case 0:
  case 1:
  case 2:
  case 3:
    return side;
  case 4:
  case 5:
  case 6:
  case 7:
    if(!c) return side-4;
    if(c==1) return (side-3)%4;
    if(c==2) return side;
  case 8:
    return c;
  default:
    PZError << "TPZCompElT2d::SideConnectLocId, connect = " << c << std::endl;
    return -1;
  }
}

void TPZCompElQ2d::SetInterpolationOrder(TPZVec<int> &ord) {
  if(ord.NElements()!=5)
    PZError << "TPZCompElT2d::SetInterpolationOrder. ord has bad number of elements.\n";
  for(int i=0;i<5;i++) fPreferredSideOrder[i] = ord[i];
}

void TPZCompElQ2d::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(5);
  for(int i=0;i<5;i++) ord[i] = fSideOrder[i];
}

TPZIntPoints *TPZCompElQ2d::CreateSideIntegrationRule(int side) {  // or TPZInt...(2*ord-1) ???
   if(side<4 || side>8) return new TPZInt1d(0);
   if(side<8) return new TPZInt1d(2*SideOrder(side));
   return new TPZIntQuad(2*fSideOrder[4],2*fSideOrder[4]);
}

int TPZCompElQ2d::PreferredSideOrder(int side) {
  if(side<4) return 0;
  if(side<9) return fPreferredSideOrder[side-4];
  return 0;
}

void TPZCompElQ2d::SetPreferredSideOrder(int side, int order) {
	if(side < 4) return;
   if(side < 9) fPreferredSideOrder[side-4] = order;
   return;
}

void TPZCompElQ2d::SetSideOrder(int side, int order) {
	if(side<0 || side>8) PZError << "TPZCompElQ2d::SetSideOrder. Bad paramenter side.\n";
	if(side>3) {
		fSideOrder[side-4] = order;
		TPZConnect &c = Connect(side);
		int seqnum = c.SequenceNumber();
		int nvar = 1;
		TPZMaterial *mat = Material();
		if(mat) nvar = mat->NStateVariables();
		Mesh()->Block().Set(seqnum,NConnectShapeF(side)*nvar);
	}
}

int TPZCompElQ2d::SideOrder(int side) {
  if(side<4 || side>8) return 0;
  return fSideOrder[side-4];
}

void TPZCompElQ2d::SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {
  switch(side) {
  case 0:
    point[0]=point[1]=-1.;
    return;
  case 1:
    point[0]= 1.;
    point[1]=-1.;
    return;
  case 2:
    point[0]=point[1]=1.;
    return;
  case 3:
    point[0]=-1.;
    point[1]=1.;
    return;
  case 4:
  case 6:
    point[0] = par[0];
	 point[1] = (side==4) ? -1.0 : 1.0;
    return;
  case 5:
  case 7:
    //point[0] = ( side == 1 ) ? 1.0 : -1.0;
    point[0] = ( side == 5 ) ? 1.0 : -1.0;//Cedric 04/05/99
    point[1] = par[0];
    return;
  case 8:
    point[0]=par[0];
    point[1]=par[1];
    return;
  default:
    PZError << "TPZCompElQ2d::SideParameterToElement. Bad paramenter side.\n";
    return;
  }
}

void TPZCompElQ2d::ElementToSideParameter(int side,TPZVec<REAL> &point,TPZVec<REAL> &par) {
  switch(side) {
  case 0:
  case 1:
  case 2:
  case 3:
    return;
  case 4:
    par[0] = point[0];
    return;
  case 6:
    par[0] = -point[0];
    return;
  case 5:
    par[0] = point[1];
    return;
  case 7:
    par[0] = -point[1];
    return;
  case 8:
    par[0]=point[0];
    par[1]=point[1];
    return;
  default:
    PZError << "TCompElQ2d::ElementToSideParameter. Bad paramenter side.\n";
    return;
  }
}

void TPZCompElQ2d::CornerShape(TPZVec<REAL> &pt,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {
  pzshape::TPZShapeQuad::ShapeCorner(pt,phi,dphi);
}

void TPZCompElQ2d::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {  //    id ???
  if(side<0 || side>8) PZError << "TPZCompElQ2d::SideShapeFunction. Bad paramenter side.\n";
  else if(side==8) Shape(point,phi,dphi);
  else if(side<4) {
  		phi(0,0)=1.;
  } else {
    TPZVec<int64_t> id(2);
	TPZManVector<int> faceorder(1);
	faceorder[0] = SideOrder(side);
	id[0] = Reference()->NodeIndex(side-4);
    id[1] = Reference()->NodeIndex((side-3)%4);
    pzshape::TPZShapeLinear::Shape(point,id,faceorder,phi,dphi);
  }
}

void TPZCompElQ2d::Shape(TPZVec<REAL> &x,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {
  TPZVec<int64_t> id(4);//Cedric
  for(int i=0;i<4;i++) id[i] = Reference()->NodeIndex(i);
  TPZManVector<int> ord(5);//Cedric
  for (int ii = 0; ii < 5; ii++) ord[ii] = fSideOrder[ii];
  pzshape::TPZShapeQuad::Shape(x,id,ord,phi,dphi);
}

//Cedric 19/03/99
void TPZCompElQ2d::SetIntegrationRule(int order) {
	TPZIntQuad intquad(order,order);
   SetIntegrationRule(intquad);
}

void TPZCompElQ2d::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
  if(dimension == 2 && Material()->Id() > 0) new TPZGraphElQ2dd(this,&grmesh);
}

void TPZCompElQ2d::Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol) {

  if(var == -1) {
    //    sol[0] = SideOrder(8);
    
    if(fabs(qsi[0]+1.) < 1.e-6) {
      // I am on the left side
      if(fabs(qsi[1]+1.) < 1.e-6) {
	// I am on the bottom side
	sol[0] = (SideOrder(4)+SideOrder(7))/2.;
      } else if(fabs(qsi[1]-1.) < 1.e-6) {
	// I am on the top side
	sol[0] = (SideOrder(6)+SideOrder(7))/2.;
      } else {
	// I am in the middle
	sol[0] = SideOrder(7);
      }
    } else if(fabs(qsi[0]-1.) < 1.e-6) {
      // I am on the right side
      if(fabs(qsi[1]+1.) < 1.e-6) {
	// I am on the bottom side
	sol[0] = (SideOrder(4)+SideOrder(5))/2.;
      } else if(fabs(qsi[1]-1.) < 1.e-6) {
	// I am on the top side
	sol[0] = (SideOrder(5)+SideOrder(6))/2.;
      } else {
	// I am in the middle
	sol[0] = SideOrder(5);
      }
    } else {
      // I am in the middle (horizontally)
      if(fabs(qsi[1]+1.) < 1.e-6) {
	// I am on the bottom side
	sol[0] = SideOrder(4);
      } else if(fabs(qsi[1]-1.) < 1.e-6) {
	// I am on the top side
	sol[0] = SideOrder(6);
      } else {
	// I am in the middle
	sol[0] = SideOrder(8);
      }
    }
  } else {
    TPZInterpolatedElement::Solution(qsi,var,sol);
  }
}

TPZCompMesh *TPZCompElQ2d::CreateMesh() {
   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->SetName("A simple square");
   firstmesh->NodeVec().Resize(4);
   TPZVec<REAL> coord(2);
   coord[0] = 0;
   coord[1] = 0;
   firstmesh->NodeVec()[0].Initialize(coord,*firstmesh);
   coord[0] = 3;
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);
   coord[0] = 3;
   coord[1] = 3;
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);
   coord[0] = 0;
   firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);
   TPZVec<int64_t> nodeindexes(4);
   nodeindexes[0] = 0;
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;
   nodeindexes[3] = 3;
    //elementos geometricos
   TPZGeoEl *elg0 = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);

   firstmesh->BuildConnectivity();

   TPZVec<TPZGeoEl *> sub;
   elg0->Divide(sub);
   sub[0]->Divide(sub);
   TPZGeoElBC(elg0,4,-1);
   TPZGeoElBC(elg0,2,-1);
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   secondmesh->SetName("A simple computational mesh");

//    TPZMat2dLin *mat2d = new TPZMat2dLin(1);
//    TPZFMatrix xk(1,1,1.),xc(1,1,1.),xf(1,1,1.);
//    //   xk(0,1) = xk(1,0) = xc(0,1) = xc(1,0) = 0.;
//    mat->SetMaterial(xk,xc,xf);
   TPZElasticityMaterial *mat = new TPZElasticityMaterial(1,2.,0.3,1.,1.);
   TPZFMatrix<STATE> val1(2,2,0.),val2(2,1,0.);
   TPZBndCond *bc = mat->CreateBC(mat,-1,0,val1,val2);
   secondmesh->InsertMaterialObject(mat);
   secondmesh->InsertMaterialObject(bc);
   secondmesh->AutoBuild();

   secondmesh->AdjustBoundaryElements();

   return secondmesh;
}

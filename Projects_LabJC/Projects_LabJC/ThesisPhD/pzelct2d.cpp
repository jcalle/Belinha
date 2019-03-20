#include "pzelct2d.h"
#include "pzelcq2d.h"
#include "pzerror.h"
#include "pzelgt2d.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzconnect.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzmat2dlin.h"
#include "pzcmesh.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
//#include "pzgraphmesh.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "pzgraphel.h"
#include "pzgeoelbc.h"
//#include "pzgraphnode.h"


TPZCompElT2d::TPZCompElT2d(TPZCompMesh &mesh,TPZGeoElT2d *ref, int64_t &index) :
	TPZIntelGen<pzshape::TPZShapeTriang>(mesh,ref,index), fIntRule(1) {

  int i;
  for(i=0;i<4;i++) {
    fSideOrder[i] = mesh.GetDefaultOrder();
    fPreferredSideOrder[i] = mesh.GetDefaultOrder();
  }
  // Jorge 19/5/99
  /** Compatibilizing the continuity over the neighboard elements */
  for(i=0;i<7;i++) {
    TPZStack<TPZCompElSide> elvec;
    TPZCompElSide thisside(this,i);
    thisside.EqualLevelElementList(elvec,0,0);
    /** Asking if exist neighbour : If this neighboard is continuous all neighboards are continuous */
  }

//  RemoveSideRestraintsII(EInsert);
  ref->SetReference(this);
  for(i=0;i<7;i++) {
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
  }
  for(i=3;i<7;i++)
    IdentifySideOrder(i);
  TPZVec<int> order(2,2*fSideOrder[3]);
  fIntRule.SetOrder(order);
}

TPZCompElT2d::TPZCompElT2d(TPZCompMesh &mesh,TPZGeoElT2d *ref, int64_t &index,int /*noconnects*/) :
	TPZIntelGen<pzshape::TPZShapeTriang>(mesh,ref,index), fIntRule(1) {
  for(int i=0;i<4;i++) {
    fSideOrder[i] = mesh.GetDefaultOrder();
    fPreferredSideOrder[i] = mesh.GetDefaultOrder();
  }
//  RemoveSideRestraintsII(EInsert);
}

void TPZCompElT2d::SetConnectIndex(int i,int connectindex) {
  if(i>-1 && i<7) fConnectIndexes[i] = connectindex;
  else {
    PZError << "TPZCompElT2d::SetConnectIndex. Bad parameter i.\n";
    PZError.flush();
  }
}

int TPZCompElT2d::ConnectIndex(int i) {
	if(i<0 || i>6) {
   	PZError << "TPZCompElT2d::ConnectIndex. Bad parameter i.\n";
      return -1;
   }
   return fConnectIndexes[i];
}

int TPZCompElT2d::NConnectShapeF(int side) {
  switch(side) {
  case 0:
  case 1:
  case 2:
    return 1;
  case 3:
  case 4:
  case 5:
    return fSideOrder[side-3]-1;
  case 6:
    return (fSideOrder[3]-2) < 0 ? 0 : ((fSideOrder[3]-2)*(fSideOrder[3]-1))/2;
  default:
    PZError << "TPZCompElT2d::NConnectShapeF, bad parameter iconnect " << side << std::endl;
    return 0;
  }
}

int TPZCompElT2d::NSideConnects(int side) {
  if(side<0 || side>6) {
    PZError << "TPZCompElT2d::NSideConnects. Bad parameter i.\n";
    return 0;
  }
  if(side<3) return 1;
  if(side<6) return 3;
  return 7;
}

/**It do not verify the values of the c*/
int TPZCompElT2d::SideConnectLocId(int c, int side) {
  switch(side) {
  case 0:
  case 1:
  case 2:
    return side;
  case 3:
  case 4:
  case 5:
    if(!c) return side-3;
    if(c==1) return (side-2)%3;
    if(c==2) return side;
  case 6:
    return c;
  default:
    PZError << "TPZCompElT2d::SideConnectLocId, connect = " << c << std::endl;
    return -1;
  }
}

void TPZCompElT2d::SetInterpolationOrder(TPZVec<int> &ord) {
  if(ord.NElements()!=4)
    PZError << "TPZCompElT2d::SetInterpolationOrder. ord has bad number of elements.\n";
  for(int i=0;i<4;i++) fPreferredSideOrder[i] = ord[i];
}

void TPZCompElT2d::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(4);
  for(int i=0;i<4;i++) ord[i] = fSideOrder[i];
}

TPZIntPoints *TPZCompElT2d::CreateSideIntegrationRule(int side) {  // or TPZInt...(2*ord-1) ???
   if(side<3 || side>6) return new TPZInt1d(0);
   if(side<6) return new TPZInt1d(2*SideOrder(side));
   return new TPZIntTriang(2*fSideOrder[3]);
}

int TPZCompElT2d::PreferredSideOrder(int side) {
  if(side<3) return 0;
  if(side<7) return fPreferredSideOrder[side-3];
  return 0;
}

void TPZCompElT2d::SetPreferredSideOrder(int side, int order) {
	if(side < 3) return;
   if(side < 7) fPreferredSideOrder[side-3] = order;
   return;
}


void TPZCompElT2d::SetSideOrder(int side, int order) {
  if(side<0 || side>6) PZError << "TPZCompElT2d::SetSideOrder. Bad parameter side.\n";
  if(side>2) {
  		fSideOrder[side-3] = order;
		TPZConnect &c = Connect(side);
		int seqnum = c.SequenceNumber();
		int nvar = 1;
		TPZMaterial *mat = Material();
		if(mat) nvar = mat->NStateVariables();
		Mesh()->Block().Set(seqnum,NConnectShapeF(side)*nvar);
  }
}

int TPZCompElT2d::SideOrder(int side) {
  if(side<3 || side>6) return 0;
  return fSideOrder[side-3];
}

void TPZCompElT2d::SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {
  point[0] = point[1] = 0.;
  switch(side) {
  case 0:
    return;
  case 1:
    point[0]=1.;
    return;
  case 2:
    point[1]=1.;
    return;
  case 3:
    point[0] = (1.+par[0])*.5;
	 return;
  case 4:
    point[0] = (1.-par[0])*.5;
	 point[1] = (1.+par[0])*.5;
    return;
  case 5:
    point[1] = (1.-par[0])*.5;
    return;
  case 6:
    point[0]=par[0];
    point[1]=par[1];
    return;
  default:
    PZError << "TPZCompElT2d::SideParameterToElement. Bad paramenter side.\n";
  }
}

void TPZCompElT2d::ElementToSideParameter(int side,TPZVec<REAL> &point,TPZVec<REAL> &par) {
  switch(side) {
  case 0:
  case 1:
  case 2:
    return;
  case 3:
    par[0] = 2.*point[0]-1.;
    return;
  case 4:
    par[0] = point[1]-point[0];
    return;
  case 5:
    par[0] = 1.- 2*point[1];
    return;
  case 6:
    par[0]=point[0];
    par[1]=point[1];
    return;
  default:
    PZError << "TPZCompElT2d::ElementToSideParameter. Bad side = " << side << ".\n";
  }
}

void TPZCompElT2d::CornerShape(TPZVec<REAL> &pt,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {
  pzshape::TPZShapeTriang::ShapeCorner(pt,phi,dphi);
}

void TPZCompElT2d::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {  //    id ???
  if(side<0 || side>6) PZError << "TPZCompElT2d::SideShapeFunction. Bad paramenter side.\n";
  else if(side==6) Shape(point,phi,dphi);
  else if(side<3) {
  	phi(0,0) = 1.;
  } else {
    TPZVec<int64_t> id(2);
	TPZManVector<int> faceorder(1);
	faceorder[0] = SideOrder(side);
	id[0] = Reference()->NodeIndex(side-3);
    id[1] = Reference()->NodeIndex((side-2)%3);
    pzshape::TPZShapeLinear::Shape(point,id,faceorder,phi,dphi);
  }
}

void TPZCompElT2d::Shape(TPZVec<REAL> &x,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {
  TPZVec<int64_t> id(3);
  for(int i=0;i<3;i++) id[i] = Reference()->NodeIndex(i);
  TPZManVector<int> ord(4);
  for (int j = 0; j < 4; j++) ord[j] = fSideOrder[j];
  pzshape::TPZShapeTriang::Shape(x,id,ord,phi,dphi);
}
//Cedric 16/03/99
void TPZCompElT2d::SetIntegrationRule(int order) {
	TPZIntTriang inttriang(order);
   SetIntegrationRule(inttriang);
}

void TPZCompElT2d::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
	  if(dimension == 2) new TPZGraphElTd(this,&grmesh);
}

TPZCompMesh *TPZCompElT2d::CreateMesh() {
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
   TPZVec<int64_t> nodeindexes(3);
   nodeindexes[0] = 0;
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;
    //elementos geometricos
   TPZGeoEl *elg0 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   nodeindexes[0] = 0;
   nodeindexes[1] = 2;
   nodeindexes[2] = 3;
   //elg0 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);

   firstmesh->BuildConnectivity();

   TPZVec<TPZGeoEl *> sub;
   elg0->Divide(sub);
   sub[0]->Divide(sub);
   TPZGeoElBC(elg0,4,-1);
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   secondmesh->SetName("A simple computational mesh");

   TPZMat2dLin *mat = new TPZMat2dLin(1);
   TPZFMatrix<STATE> xk(1,1,1.),xc(1,1,1.),xf(1,1,1.);
   //   xk(0,1) = xk(1,0) = xc(0,1) = xc(1,0) = 0.;
   mat->SetMaterial(xk,xc,xf);
   TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
   TPZBndCond *bc = mat->CreateBC(mat,-1,0,val1,val2);
   secondmesh->InsertMaterialObject(mat);
   secondmesh->InsertMaterialObject((TPZMaterial *)bc);
   secondmesh->AutoBuild();



   return secondmesh;
}

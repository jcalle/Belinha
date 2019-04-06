//METHODS DEFINITION FOR CLASS ELEMQ2D

#include "pzelgq2d.h"
#include "pzelg1d.h"
#include "pzelcq2d.h"
#include "pzelgpoint.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgnode.h"
#include "pzshtmat.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzshapequad.h"

static TPZCompEl *CreateEl(TPZGeoElQ2d *gel,TPZCompMesh &mesh, int64_t &index) {
  return new TPZCompElQ2d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElQ2d::fp)(TPZGeoElQ2d *,TPZCompMesh &, int64_t &) = CreateEl;

TPZGeoElQ2d::TPZGeoElQ2d(int id,TPZVec<int64_t> &nodeindexes,int matid,TPZGeoMesh &mesh):
	TPZGeoElRefLess<pzgeom::TPZGeoQuad>(id,nodeindexes,matid,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=4) {
    PZError << "TPZGeoElQ2d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<4;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<4;i++) fSubEl[i] = 0;
}

TPZGeoElQ2d::TPZGeoElQ2d(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) :
	TPZGeoElRefLess<pzgeom::TPZGeoQuad>(nodeindexes,matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=4) {
    PZError << "TPZGeoElQ2d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<4;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<4;i++) fSubEl[i] = 0;
}

TPZGeoElQ2d::~TPZGeoElQ2d() {}

void TPZGeoElQ2d::Shape(TPZVec<REAL> &param,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {

   REAL x=param[0], y=param[1];

   phi(0,0) = .25*(1.-x)*(1.-y);
   phi(1,0) = .25*(1.+x)*(1.-y);
   phi(2,0) = .25*(1.+x)*(1.+y);
   phi(3,0) = .25*(1.-x)*(1.+y);

   dphi(0,0) = .25*(y-1.);
   dphi(1,0) = .25*(x-1.);

   dphi(0,1) = .25*(1.-y);
   dphi(1,1) =-.25*(1.+x);

   dphi(0,2) = .25*(1.+y);
   dphi(1,2) = .25*(1.+x);

   dphi(0,3) =-.25*(1.+y);
   dphi(1,3) = .25*(1.-x);


/*   REAL spacephi1[3],spacephi2[3],spacedphi1[3],spacedphi2[3];
   TPZFMatrix phi1(2,1,spacephi1,3);
   TPZFMatrix phi2(2,1,spacephi2,3);
   TPZFMatrix dphi1(1,2,spacedphi1,3);
   TPZFMatrix dphi2(1,2,spacedphi2,3);
   REAL x=param[0], y=param[1];

   TPZGeoEl::Shape1d(x,2,phi1,dphi1);
   TPZGeoEl::Shape1d(y,2,phi2,dphi2);

   phi(0,0) = phi1(0,0)*phi2(0,0);
   dphi(0,0) = dphi1(0,0)*phi2(0,0);
   dphi(1,0) = phi1(0,0)*dphi2(0,0);
   phi(1,0) = phi1(1,0)*phi2(0,0);
   dphi(0,1) = dphi1(0,1)*phi2(0,0);
   dphi(1,1) = phi1(1,0)*dphi2(0,0);
   phi(2,0) = phi1(1,0)*phi2(1,0);
   dphi(0,2) = dphi1(0,1)*phi2(1,0);
   dphi(1,2) = phi1(1,0)*dphi2(0,1);
   phi(3,0) = phi1(0,0)*phi2(1,0);
   dphi(0,3) = dphi1(0,0)*phi2(1,0);
   dphi(1,3) = phi1(0,0)*dphi2(0,1);*/
}

TPZGeoElQ2d *TPZGeoElQ2d::CreateGeoEl(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElQ2d(nodeindexes,matind,mesh);
}

int TPZGeoElQ2d::NNodes() {
  return 4;
}
int64_t TPZGeoElQ2d::NodeIndex(int node) {
  if(node<0 || node>3) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElQ2d::NSideNodes(int side) {
  if(side<0 || side>8) {
    PZError << "TPZGeoElQ2d::NSideNodes. Bad parameter side.\n";
    return 0;
  }
  if(side<4) return 1;
  else if (side==8) return 4;
  else return 2;
}

int64_t TPZGeoElQ2d::SideNodeIndex(int side,int node) {
  if(side<0 || side>8) {
    PZError << "TPZGeoElQ2d::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  if(side==8 && node>=0 && node<9) return fNodeIndexes[node];//um dos cantos
  if(side < 4) return fNodeIndexes[side];
  //4,5,6,7
  side-=4;
  if(node==0) return fNodeIndexes[side];      //0,1,2,3
  if(node==1) return fNodeIndexes[(side+1)%4];//1,2,3,0
  return -1;
}

void TPZGeoElQ2d::MidSideNodeIndex(int side, int64_t &index) {
  index = -1;
  if(side<0 || side>8) {
    PZError << "TPZGeoElQ2d::MidSideNodeIndex. Bad parameter side = " << side <<std::endl;
    return;
  }
  //lados 0,1,2
	if(side<4) {//o nó medio do lado 0 é o 0 etc.
		index=fNodeIndexes[side];
		return;
	}
	if(HasSubElement(0) && side < 8) {//caso tenha filhos, se não tiver retorna -1
	    index=((TPZGeoElQ2d *) fSubEl[side-4])->fNodeIndexes[(side-3)%4];
	}
	if(HasSubElement(0) && side == 8) {
	    index=((TPZGeoElQ2d *) fSubEl[0])->fNodeIndexes[2];
	}

}

void TPZGeoElQ2d::NewMidSideNode(int side,int64_t &index) {
  MidSideNodeIndex(side,index);
  if(index < 0) {
    TPZGeoElSide gelside = Neighbour(side);
    while(gelside.Element()) {
      gelside.Element()->MidSideNodeIndex(gelside.Side(),index);
      if(index!=-1) return;
      gelside = gelside.Neighbour();
      if(gelside.Element()==this) break;
    }
    TPZVec<REAL> par(2,0.);//(0,0)
    TPZVec<REAL> coord(3,0.);//(0,0,0)
    //lados do elemento neste caso são : side = 0,1,2,3
    if     (side==4) par[1] = -1.;//(0,-1)
    else if(side==5) par[0] =  1.;//(1,0)
    else if(side==6) par[1] =  1.;//(0,1)
    else if(side==7) par[0] = -1.;//(-1,0)
    X(par,coord);
    index = Mesh()->NodeVec().AllocateNewElement();
    Mesh()->NodeVec()[index].Initialize(coord,*Mesh());
  }
}

/**Determine the coordinate of the center of the element*/
int64_t TPZGeoElQ2d::CenterIndex() {
  TPZGeoElSide nghsd = Neighbour(8);
  while(nghsd.Exists() && !nghsd.HasSubElement() && nghsd.Element() != this) nghsd = nghsd.Neighbour();
// Philippe 15/4/99
  if(nghsd.Exists() && nghsd.HasSubElement()) {
	  int64_t index = -1;
	  nghsd.Element()->MidSideNodeIndex(nghsd.Side(),index);
	  return index;
  }
  TPZVec<REAL> coord(3);
  TPZVec<REAL> param(2,0);
  X(param,coord);//It determine the centroid of the quadrilateral element.
  int indexcenter = Mesh()->NodeVec().AllocateNewElement();
  Mesh()->NodeVec()[indexcenter].Initialize(coord,*Mesh());
  return indexcenter;
  // We can not to take the center type of the master cell because of characteristic
  // point is not preserved for the transformation.
}

int TPZGeoElQ2d::SideDimension(int side) {

	if (side<0 || side>8) {
   	PZError << "TPZGeoElQ2d::SideDimension called with side " << side <<std::endl;
      return 0;
   }
   if(side<4) return 0;
   else if(side==8) return 2;
   return 1;
}

TPZGeoElSide TPZGeoElQ2d::HigherDimensionSides(int side,int targetdimension) {
//targetdimension deve ser 1 ou 2
//se side = 0,1,2,3 targetdimension deve ser 1
//se side = 4,5,6,7 targetdimension deve ser 2
  if( (side<0 || side>7) || (targetdimension != 1 && targetdimension != 2) ) {
    PZError << "TPZGeoElQ2d::HigherDimensionSides called with side = " << side
	    << " targetdimension = " << targetdimension <<std::endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
//side = 0,1,2,3,4,5,6,7
  switch(targetdimension) {//=1,2
	  case 1:
       if       (this == father->SubElement(0)) {
       	 if(side==1) return TPZGeoElSide(this,4);
       	 if(side==3) return TPZGeoElSide(this,7);
       } else if(this == father->SubElement(1)) {
       	 if(side==0) return TPZGeoElSide(this,4);
       	 if(side==2) return TPZGeoElSide(this,5);
       } else if(this == father->SubElement(2)) {
       	 if(side==1) return TPZGeoElSide(this,5);
       	 if(side==3) return TPZGeoElSide(this,6);
       } else if(this == father->SubElement(3)) {
       	 if(side==2) return TPZGeoElSide(this,6);
       	 if(side==0) return TPZGeoElSide(this,7);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
  	  case 2:
       	return TPZGeoElSide(this,8);
  }//switch
  return TPZGeoElSide();
}

void TPZGeoElQ2d::LowerDimensionSides(int side,TPZStack<TPZGeoElSide> &smallsides) {
	if (side < 4) return;
   int i;
   if(side < 8) {//side = 4,5,6,7
   	 smallsides.Push(TPZGeoElSide(this,side-4));
  	    smallsides.Push(TPZGeoElSide(this,(side+1)%4));
   } else if(side==8) {
//   	for (i=0;i<9;i++) smallsides.Push(TPZGeoElSide(this,i));
// Philippe 17/5/99
   	for (i=0;i<8;i++) smallsides.Push(TPZGeoElSide(this,i));
   }
}

void TPZGeoElQ2d::SideMasterCo(int side,TPZFMatrix<STATE> &coord) {
  if(side<0 || side>8) {
    PZError << "TPZGeoElQ2d::SideMasterCo. Bad parameter side.\n";
    return;
  }
  int row = coord.Rows();
  if(side==8) {//side 8 tem dimensão 2
    coord.Redim(row,4);//2x4
    coord(0,1) = coord(1,2) = 1.;// [ {0,1,0,-1},{-1,0,1,0} ]
    coord(0,3) = coord(1,0) =-1.;//ou (0,-1),(1,0),(0,1),(-1,0) nós 0,1,2,3
    return;
  }
  if(side<3) {         //side=0,1,2,3 tem dimensão 0
    coord.Redim(row,1);//2x1 um unico nó
    if     (side==0) {coord(0,0)=-1.;coord(1,0)=-1.;}//(-1,-1) : 0
    else if(side==1) coord(0,0)= 1.;//(1,-1) : 1
    else if(side==2) coord(1,0)= 1.;//(1,1) : 2
    else if(side==3) coord(0,0)=-1.;//(-1,1) : 3
    return;
  }
  coord.Redim(row,2);  //side = 4,5,6,7 tem dimensão 1
  coord.Zero();        //2x2  [{0,0},{0,0}]
  if     (side==4) {coord(0,0)=-1.;coord(0,1)= 1.;coord(1,0)=-1.;coord(1,1)=-1.;}//(-1,-1),(1,-1) : 0,1
  else if(side==5) {coord(0,0)= 1.;coord(0,1)= 1.;coord(1,1)= 1.;}//(1,-1),(1,1) : 1,2
  else if(side==6) {coord(1,0)= 1.;coord(0,1)=-1.;}//(1,1),(-1,1) : 2,3
  else if(side==7) {coord(0,0)=-1.;coord(1,1)=-1.;}//(1,1),(-1,1) : 3,0
}

void TPZGeoElQ2d::Jacobian(TPZVec<REAL> &param,TPZFMatrix<STATE> &jacobian,TPZFMatrix<STATE> &axes,REAL &detjac,TPZFMatrix<STATE> &jacinv){

  int nnodes = NNodes();
#ifdef DEBUG
  if (nnodes != 4) {
    PZError << "TPZGeoElQ2d.jacobian only implemented for"
      " 4, 8 or 9 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 2 || param[0] < -1. || param[0] > 1. ||
     param[1] < -1. || param[1] > 1.) {
    PZError << "TPZGeoElQ2d.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << "\n";
    return;
  }
#endif
  REAL spacephi[4];
  TPZFMatrix<STATE> phi(4,1,spacephi,4);
  REAL spacedphi[8];
  TPZFMatrix<STATE> dphi(2,4,spacedphi,8);
  Shape(param,phi,dphi);
  int i,j;
  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      jacobian(i,j)=0.;

  TPZGeoNode *np;
  TPZVec<REAL> V1(3,0.),V2(3,0.),V2til(3,0.),V3(3,0.);
  REAL V1Norm=0.,V1V2=0.,V2tilNorm=0.;
  for(i=0;i<nnodes;i++) {
    np = NodePtr(i);
    for(j=0;j<3;j++) {
      V1[j] += np->Coord(j)*dphi(0,i);
      V2[j] += np->Coord(j)*dphi(1,i);
    }
  }
  for(i=0;i<3;i++) {
    V1Norm += V1[i]*V1[i];
    V1V2 += V1[i]*V2[i];
  }
  V1Norm = sqrt(V1Norm);
  for(i=0;i<3;i++) {
    V1[i] /= V1Norm;
    V2til[i] = V2[i] - V1V2*V1[i]/V1Norm;
    V2tilNorm += V2til[i]*V2til[i];
  }
  V2tilNorm = sqrt(V2tilNorm);
  jacobian(0,0) = V1Norm;
  jacobian(0,1) = V1V2/V1Norm;
  jacobian(1,1) = V2tilNorm;
  for(i=0;i<3;i++) {
    axes(0,i) = V1[i];
    axes(1,i) = V2til[i]/V2tilNorm;
  }
  detjac = jacobian(0,0)*jacobian(1,1)-jacobian(1,0)*jacobian(0,1);
  jacinv(0,0) = +jacobian(1,1)/detjac;
  jacinv(1,1) = +jacobian(0,0)/detjac;
  jacinv(0,1) = -jacobian(0,1)/detjac;
  jacinv(1,0) = -jacobian(1,0)/detjac;

  axes(2,0) = axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1);
  axes(2,1) = -axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0);
  axes(2,2) = axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0);
}

void TPZGeoElQ2d::X(TPZVec<REAL> & loc,TPZVec<REAL> &result){
  REAL spacephi[4],spacedphi[8];
  int i,j;
  TPZFMatrix<STATE> phi(4,1,spacephi,4);
  TPZFMatrix<STATE> dphi(2,4,spacedphi,8);
  Shape(loc,phi,dphi);
  for(i=0;i<3;i++) {
    result[i] = 0.0;
    for(j=0;j<4;j++)
      result[i] += phi(j,0)*NodePtr(j)->Coord(i);
  }
}
/**It's necessary to define the normal vector to side 4, that is the orthogonal
   vector to the surface*/
void TPZGeoElQ2d::NormalVector(int side,TPZVec<REAL> &param,TPZVec<REAL> &normal,
			     TPZFMatrix<STATE> &axes,TPZFMatrix<STATE> &jac1d) {
  int nnodes = NNodes();
#ifdef DEBUG
  if (nnodes != 4) {
    PZError << "TPZGeoElQ2d.NormalVector, only implemented for"
      " 4 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 2 || param[0] < -1. || param[0] > 1. ||
     param[1] < -1. || param[1] > 1.) {
    PZError << "TPZGeoElQ2d.NormalVector, fl out of range : "
      " point.NElements() = " << param.NElements() <<
      "\npoint[0] = " << param[0] << " point[1] = " << param[1] << "\n";
    return;
  }
  if(normal.NElements() != 3) {
    PZError << "elgq2d.NormalVector normal.capacity() = " << normal.NElements() <<
      "\n";
    return;
  }
  if(side < 0 || side >= 4) {
    PZError << "TPZGeoElQ2d.jacobian invalid side : "
      " side = " << side << "\n";
    return;
  }
#endif

  REAL spacephi[4],spacedphi[8];
//  TPZFMatrix phi(4,1,spacedphi,4);
// Philippe 31;3;99
  TPZFMatrix<STATE> phi(4,1,spacephi,4);
  TPZFMatrix<STATE> dphi(2,4,spacedphi,8);
  Shape(param,phi,dphi);
  TPZGeoNode *np;
  TPZVec<REAL> t(3,0.);
  int i,j,ider = 0;
  if(side==1 || side==3) ider = 1;

  for(i=0;i<nnodes;i++) {
    np = NodePtr(i);
    for(j=0;j<3;j++)
      t[j] += np->Coord(j)*dphi(ider,i);
  }

  //      note that ||t|| != 1 , ||t|| = |J|

  jac1d(0,0) = sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );

  // consistent axes computation Philippe 17/4/97
  TPZVec<REAL> V1(3,0.),V2(3,0.),V2til(3,0.),V3(3,0.),V1til(3,0.);
  REAL V1Norm=0.,V2Norm=0.,V1V2=0.,V2tilNorm=0.,V1tilNorm =0.;
  for(i=0;i<nnodes;i++) {
    np = NodePtr(i);
    for(j=0;j<3;j++) {
      V1[j] += np->Coord(j)*dphi(0,i);
      V2[j] += np->Coord(j)*dphi(1,i);
    }
  }
  for(j=0;j<3;j++) {
    V1Norm += V1[j]*V1[j];
    V2Norm += V2[j]*V2[j];
    V1V2 += V1[j]*V2[j];
  }
  V1Norm = sqrt(V1Norm);
  V2Norm = sqrt(V2Norm);
  for(j=0;j<3;j++) {
    V1[j] /= V1Norm;
    V2[j] /= V2Norm;
    V2til[j] = V2[j] - V1V2*V1[j]/V1Norm/V2Norm;
    V1til[j] = V1[j] - V1V2*V2[j]/V1Norm/V2Norm;
    V2tilNorm += V2til[j]*V2til[j];
    V1tilNorm += V1til[j]*V1til[j];
  }
  V2tilNorm = sqrt(V2tilNorm);
  V1tilNorm = sqrt(V1tilNorm);
  for(j=0;j<3;j++) {
    axes(0,j) = V1[j];
    axes(1,j) = V2til[j]/V2tilNorm;
  }
  switch(side) {
  case 0:
  case 2:
    for(i=0;i<3;i++)
      normal[i] = V2til[i]/V2tilNorm;
    break;
  case 1:
  case 3:
    for(i=0;i<3;i++)
      normal[i] = V1til[i]/V1tilNorm;
    break;
  }
  if(side == 0 || side == 3) for(i=0;i<3;i++) normal[i] *= -1.;

  axes(2,0) = axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1);
  axes(2,1) = -axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0);
  axes(2,2) = axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0);
  return;
}

/** TO SUBDIVISION
********************************************************************************
  Into Divides is necesary to consider the connectivity with the all neighboards*/
void TPZGeoElQ2d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {
   //If exist fSubCell return this sons
   int i;
   if(HasSubElement(0)) {
      SubElVec.Resize(4);
      for(i=0;i<4;i++) SubElVec[i] = fSubEl[i];
      return;
   }

   int j, sub, matid = MaterialId();
	int64_t index;
   int64_t np[4][4];//guarda conectividades dos 4 subelementos

   for(j=0;j<4;j++) np[j][j]=NodeIndex(j);
   for(sub=0;sub<4;sub++) {
      NewMidSideNode(sub+4,index);
      j=(sub+1)%4;
      np[sub][j]=np[j][sub]=index;
   }
   np[0][2] = np[2][0] = np[1][3] = np[3][1] = CenterIndex();

   // creating new subelements
   for(i=0;i<4;i++) {
   	TPZManVector<int64_t> npvec(4);//np[i][0],np[i][1],np[i][2],np[i][3]
	for (j = 0; j < 4; j++) npvec[j] = np[i][j];
      fSubEl[i] = CreateGeoEl(npvec,matid,*Mesh());
   }

   SubElVec.Resize(4);
   for(sub=0;sub<4;sub++) {
      SubElVec[sub] = fSubEl[sub];
      SubElVec[sub]->SetFather(this);
   }
   //nós do lado do atual
   for(i=0;i<4;i++) {       //side do atual                    viz,side
      fSubEl[i]->SetNeighbour((i+1)%4,TPZGeoElSide(fSubEl[(i+1)%4],i));
      fSubEl[(i+1)%4]->SetNeighbour(i,TPZGeoElSide(fSubEl[i],(i+1)%4));
   }
   //conectiv do medio dos lados dos subs
   fSubEl[0]->SetNeighbour(5,TPZGeoElSide(fSubEl[1],7));
   fSubEl[0]->SetNeighbour(6,TPZGeoElSide(fSubEl[3],4));
   fSubEl[1]->SetNeighbour(6,TPZGeoElSide(fSubEl[2],4));
   fSubEl[1]->SetNeighbour(7,TPZGeoElSide(fSubEl[0],5));
   fSubEl[2]->SetNeighbour(4,TPZGeoElSide(fSubEl[1],6));
   fSubEl[2]->SetNeighbour(7,TPZGeoElSide(fSubEl[3],5));
   fSubEl[3]->SetNeighbour(4,TPZGeoElSide(fSubEl[0],6));
   fSubEl[3]->SetNeighbour(5,TPZGeoElSide(fSubEl[2],7));
   //nó central do atual
   for(i=0;i<4;i++) fSubEl[i]->SetNeighbour((i+2)%4,TPZGeoElSide(fSubEl[(i+1)%4],(i+3)%4));
   //procura-se um viz pela face que seja dividido
   TPZGeoElSide dividedneighbour = Neighbour(8);//BuildConnectivity(..) inicializou
   while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
   	if(dividedneighbour.HasSubElement()) break;
      dividedneighbour = dividedneighbour.Neighbour();
   }

   if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
   // we found a neighbour for everybody!!
      TPZManVector<int64_t> nodes(4);
      TPZStack<TPZGeoElSide> neighsub;
      for(i=0;i<4;i++) nodes[i] = NodeIndex(i);//nós globais na ordem e sentido local atual
      dividedneighbour.GetSubElements2(neighsub);//4 subs viz na ordem dos subs do atual
      nodes.Resize(1);
      //side central
      nodes[0] = fSubEl[0]->NodeIndex(2);//certo em ambos sentidos de numeração local
      int locside = neighsub[0].Element()->WhichSide(nodes);//lado do viz do sub
      TPZGeoElSide sub(fSubEl[0],2);//subel
      sub.SetConnectivity(TPZGeoElSide(neighsub[0].Element(),locside));
      //sides pontuais (cantos)
      for(i=0; i<4; i++) {
         nodes[0] = NodeIndex(i);
         int locside = neighsub[i].Element()->WhichSide(nodes);//lado do viz
         TPZGeoElSide sub(fSubEl[i],i);//subel
         sub.SetConnectivity(TPZGeoElSide(neighsub[i].Element(),locside));
         int iplus = (i+1)%4;
         nodes[0] = fSubEl[i]->NodeIndex(iplus);
         locside = neighsub[i].Element()->WhichSide(nodes);//o nó pertence a ambos elementos
         sub = TPZGeoElSide(fSubEl[i],iplus);
         sub.SetConnectivity(TPZGeoElSide(neighsub[i].Element(),locside));
      }
		//side = 4,5,6,7
      nodes.Resize(2);
      for(i=0; i<4; i++) {
      	for(j=4; j<8; j++) {
			int badside = i+6;
			if(i>1) badside = i+2;
			if(j== badside) continue;
         	nodes[0] = fSubEl[i]->SideNodeIndex(j,0);
         	nodes[1] = fSubEl[i]->SideNodeIndex(j,1);//indices 0,1 dos nós do lado j do sub i
            int locside = neighsub[i].Element()->WhichSide(nodes);//lado do viz do subelemento
            TPZGeoElSide sub(fSubEl[i],j);
            sub.SetConnectivity(TPZGeoElSide(neighsub[i].Element(),locside));
         }
      }
      for(i=0; i<4; i++) {//side de face para os 4 subs
      	TPZGeoElSide sub(fSubEl[i],8);
         sub.SetConnectivity(neighsub[i]);
      }
   	return;
   }
	//procura-se um viz pelo lado do atual, que seja dividido
   for(i=4; i<8; i++) {//(side 4 : subs 0,1),(side 5 : subs 1,2),
   	dividedneighbour = Neighbour(i);//(side 6 : subs 2,3),(side 7 : subs 3,0)
      while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
      // we found a divided neighbour for side i
         int iplus = (i-3)%4;
         TPZManVector<int64_t> nodes(2);
         TPZStack<TPZGeoElSide> neighsub;
         nodes[0] = SideNodeIndex(i,0);//fSubEl[i-4]->
         nodes[1] = SideNodeIndex(i,1);//fSubEl[i-4]->
         dividedneighbour.GetSubElements2(neighsub);
         TPZGeoElSide sub(fSubEl[i-4],i);
         sub.SetConnectivity(neighsub[0]);
         sub = TPZGeoElSide(fSubEl[iplus],i);
         sub.SetConnectivity(neighsub[1]);
         nodes.Resize(1);
         nodes[0] = fSubEl[i-4]->SideNodeIndex(iplus,0);
         int locside = neighsub[0].Element()->WhichSide(nodes);
         sub = TPZGeoElSide(fSubEl[i-4],iplus);
         sub.SetConnectivity(TPZGeoElSide(neighsub[0].Element(),locside));
      }
   }
   for(i=0; i<4; i++) {
   	dividedneighbour = Neighbour(i);
      while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
      // we found a divided neighbour for side i
         TPZManVector<int64_t> nodes(1);
         TPZStack<TPZGeoElSide> neighsub;
         nodes[0] = SideNodeIndex(i,0);
         dividedneighbour.GetSubElements2(neighsub);
         TPZGeoElSide sub(fSubEl[i],i);
         sub.SetConnectivity(neighsub[0]);
      }
   }
}


int TPZGeoElQ2d::NSubElements() {
  return 4;
}

int TPZGeoElQ2d::NSideSubElements(int side) {
  if(side < 0 || side > 8) {
    PZError << "TPZGeoElQ2d::NSideSubElements called for side " << side <<std::endl;
    return 0;
  }
  if(side==8) return 4;
  if(side<4) return 1;
  return 2;//side = 4,5,6,7
}

TPZGeoElSide TPZGeoElQ2d::SideSubElement(int side,int position) {
   if (position<0 ||position>3 || side <0 ||side>8) {
   	PZError << "TPZGeoElQ2d::SideSubElement called with position " << position << " side " << side <<std::endl;
      return TPZGeoElSide();
   }
   if(side==8) return TPZGeoElSide(SubElement(position),8);
   if(side<4) {
      if(position!=0) {
         PZError << "TPZGeoElQ2d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
         return TPZGeoElSide(SubElement(side),side);
      }
   }
   //side = 4,5,6,7
   if(position==0) return TPZGeoElSide(SubElement(side-4),side);
   else return TPZGeoElSide(SubElement((side-3)%4),side);

}

void TPZGeoElQ2d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
  if(!fSubEl[0]) {
    sub.Resize(0);
    return;
  }
  if(side < 0 || side > 8) {
    PZError << "TPZGeoElQ2d::SideSubElements called for side " << side <<std::endl;
    return;
  }
  if(side==8) {//8
    for(int i=0;i<4;i++) sub[i] = fSubEl[i];
    return;
  }
  //0,1,2,3
  if(side<4) {
    sub[0]=fSubEl[side];
    return;
  }
  //4,5,6,7
  side-=4;
  sub.Resize(2);
  sub[0] = fSubEl[side];      //0,1,2,3
  sub[1] = fSubEl[(side+1)%4];//1,2,3,0
}

TPZGeoElSide TPZGeoElQ2d::Father(int side) {
	TPZGeoEl *fFather = TPZGeoEl::Father();

   if(!fFather) return TPZGeoElSide();
   int whichsub = -1;
   int i;
   for(i=0; i<4; i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) {//equivale a is = 4  ou is > 3
	   PZError << "TPZGeoElQ2d::Father. fFather isn't father of this element.\n";
   	return TPZGeoElSide();
   }
   //agora o atual elemento é o filho numero whichsub < 4
   if(whichsub == side || side==8) return TPZGeoElSide(fFather,side);//side = 0,1,2,3
   //side = 4,5,6,7
 	if(whichsub == 0 && (side==4 || side==7)) return TPZGeoElSide(fFather,side);//ou
 	if(whichsub == 1 && (side==4 || side==5)) return TPZGeoElSide(fFather,side);//ou side==3 || side==4
 	if(whichsub == 2 && (side==5 || side==6)) return TPZGeoElSide(fFather,side);//ou side==4 || side==5
 	if(whichsub == 3 && (side==6 || side==7)) return TPZGeoElSide(fFather,side);//ou side==4 || side==5
 	//if(wichsub == 3) return TPZGeoElSide();//é feito pelo seguinte caso

//   PZError << "TPZGeoElQ2d::Father. fFather isn't father of this element along the given side.\n";
   return TPZGeoElSide();//inclui os outros casos

}

void TPZGeoElQ2d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {
// REVISAR
   int nsub = NSideSubElements(side);
   sub.Resize(nsub);
   if(!nsub) return;

   int i,j,k;
   if(nsub==1) {//side = 0,1,2,3
   	if(fSubEl[side]->NodeIndex(side)!=refnode[0]) {
      	PZError << "TPZGeoElQ2d::GetSubElement subelement does not contain refnode" <<std::endl;
         return;
      }
	   sub[0]=TPZGeoElSide(fSubEl[side],side);
   	return;
   }
   //nsub = 4
   //if(nsub == 4) nsub = 3; se nsub = 4 => side = 8
   for(i=0;i<nsub;i++) {
   	TPZGeoElSide sidesub = SideSubElement(side,i);
      TPZGeoEl *subel = sidesub.Element();
		for(k=0; k<nsub; k++) {
		   for(j=0;j<4;j++) {
			   if(subel->NodeIndex(j)==refnode[k]) {
            	sub[k] = SideSubElement(side,i);
            }
         }
      }
   }
   //if(side == 8) sub[3] = TPZGeoElSide(fSubEl[3],8);
   return;
}

/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
void TPZGeoElQ2d::BuildTransform(int side,TPZGeoEl *father,TPZTransform<STATE> &t) {
  if(father == this) return;
  TPZGeoEl *fFather = TPZGeoEl::Father();
  if(!fFather) {
	   std::cout << "TPZGeoElQ2d::BuildTransform called for inconsistent parameters\n";
     return;
   }
   int whichsub = -1;
   int i;
   for(i=0; i<4; i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) return;
   int dim = SideDimension(side);
   TPZTransform<STATE> tloc(dim);
   REAL store[6];
   TPZFMatrix<STATE> mult(dim,dim,store,4);
   TPZFMatrix<STATE> sum(dim,1,store+4,2);
   mult.Zero();
   sum.Zero();

   TPZGeoEl *locfather;
   if(father == this) {
	   mult(0,0) = mult(1,1) = 1.;
      locfather = fFather;
   } else {
	   locfather = fFather;//pai do lemento atual
   }
   if(!locfather) {
   	PZError << "TPZGeoElT2d::BuildTransform could not identify the father element\n";
	   return;
   }

   if(side == 8) {//face para face
      mult(0,0) = 0.5;
      mult(1,1) = 0.5;
      sum(0,0) = -0.5;
      sum(1,0) = -0.5;
      switch(whichsub) {//o atual é o filho numero whichsub
         case 0:
            break;
         case 1:
            sum(0,0) *= -1.;
            break;
         case 2:
         	sum(0,0) *= -1.;
            sum(1,0) *= -1.;
            break;
         case 3:
            sum(1,0) *= -1.;
            break;
      }
   } else if(side > 3 && side < 8) {
   	mult(0,0) = 0.5;
      if(whichsub == side-4) sum(0,0) = -0.5;
      else sum(0,0) = 0.5;
   }

   tloc.SetMatrix(mult,sum);
   t = tloc.Multiply(t);
   if(locfather != father) locfather->BuildTransform2(side,father,t);
}

TPZTransform<STATE> TPZGeoElQ2d::SideToSideTransform(int sidefrom,int sideto) {
   if(sideto != 8 && (sidefrom > 3 || sidefrom < 0)) {
      PZError << "TPZGeoElQ2d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   int dimfrom = SideDimension(sidefrom);
   int dimto = SideDimension(sideto);
   if(dimto <= dimfrom){
      PZError << "TPZGeoElQ2d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   TPZTransform<STATE> trans(dimto,dimfrom);
   //daqui dimto > dimfrom
   switch(sidefrom) {
    	//sidefrom de dim = 1
      case 4:
      	trans.Mult()(0,0) =  1.;
         trans.Mult()(1,0) =  0.;
         trans.Sum()(0,0)  =  0.;
         trans.Sum()(1,0)  = -1.;
         break;
      case 5:
         trans.Mult()(0,0) =  0.;
         trans.Mult()(1,0) =  1.;
         trans.Sum()(0,0)  =  1;
         trans.Sum()(1,0)  =  0.;
         break;
      case 6:
         trans.Mult()(0,0) = -1.;
         trans.Mult()(1,0) =  0.;
         trans.Sum()(0,0)  =  0.;
         trans.Sum()(1,0)  =  1.;
         break;
      case 7:
         trans.Mult()(0,0) =  0.;
         trans.Mult()(1,0) = -1.;
         trans.Sum()(0,0)  = -1.;
         trans.Sum()(1,0)  =  0.;
         break;
      //sidefrom de dim = 0
      case 0:
      	//sideto = 4,5 de dim = 1
      	if(sideto == 4) trans.Sum()(0,0) = -1.;
         else if(sideto == 7) trans.Sum()(0,0) = 1.;
         //sideto =6 de dim = 2
         else if(sideto == 8) {
         	trans.Sum()(0,0) = -1.;
            trans.Sum()(1,0) = -1.;
         } else {
            PZError << "TPZGeoElQ2d:SideToSideTransform sidefrom = " << sidefrom <<
            " sideto = " << sideto <<std::endl;
            return trans;
         }
         break;
      case 1:
      	if(sideto == 4) trans.Sum()(0,0) = 1.;
         else if(sideto == 5) trans.Sum()(0,0) = -1.;
         else if(sideto == 8) {
         	trans.Sum()(0,0) = 1.;
            trans.Sum()(1,0) =-1.;
         } else {
            PZError << "TPZGeoElQ2d:SideToSideTransform sidefrom = " << sidefrom <<
            " sideto = " << sideto <<std::endl;
            return trans;
         }
         break;
      case 2:
      	if(sideto == 5) trans.Sum()(0,0) = 1.;
         else if(sideto == 6) trans.Sum()(0,0) = -1.;
         else if(sideto == 8) {
         	trans.Sum()(0,0) = 1.;
            trans.Sum()(1,0) = 1.;
         } else {
            PZError << "TPZGeoElQ2d:SideToSideTransform sidefrom = " << sidefrom <<
            " sideto = " << sideto <<std::endl;
            return trans;
         }
         break;
      case 3:
      	if(sideto == 6) trans.Sum()(0,0) = 1.;
         else if(sideto == 7) trans.Sum()(0,0) = -1.;
         else if(sideto == 8) {
         	trans.Sum()(0,0) =-1.;
            trans.Sum()(1,0) = 1.;
         } else {
            PZError << "TPZGeoElQ2d:SideToSideTransform sidefrom = " << sidefrom <<
            " sideto = " << sideto <<std::endl;
            return trans;
         }
         break;

   }
   return trans;
}

TPZCompEl *TPZGeoElQ2d::CreateBCCompEl(int side,int bc,TPZCompMesh &cmesh) {
   if(side==8) {//8
      TPZManVector<int64_t> nodes(4);
	  for (int j = 0; j < 4; j++)
		  nodes[j] = fNodeIndexes[j];
      TPZGeoElQ2d *gel = CreateGeoEl(nodes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,0));
      TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(this,1));
      TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(this,2));
      TPZGeoElSide(gel,3).SetConnectivity(TPZGeoElSide(this,3));
      TPZGeoElSide(gel,4).SetConnectivity(TPZGeoElSide(this,4));
      TPZGeoElSide(gel,5).SetConnectivity(TPZGeoElSide(this,5));
      TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(this,6));
      TPZGeoElSide(gel,7).SetConnectivity(TPZGeoElSide(this,7));
      TPZGeoElSide(gel,8).SetConnectivity(TPZGeoElSide(this,8));
	  int64_t index;
      return gel->CreateCompEl(cmesh,index);
   }
   else if(side>-1 && side<4) {//side = 0,1,2,3
      TPZBndCond *bcptr = (TPZBndCond *) cmesh.FindMaterial(bc);
      if(!bcptr) {
         PZError << "TPZGeoElQ2d::CreateBCCompEl has no bc.\n";
         return 0;
      }
      TPZGeoElPoint *gel;

      TPZManVector<int64_t> node(1);
	  node[0] = fNodeIndexes[side];
	  gel = new TPZGeoElPoint(node,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,side));
	  int64_t index;
      return gel->CreateCompEl(cmesh,index);

   } else if (side > 3 && side < 8) {//side = 4,5,6,7
      TPZManVector<int64_t> nodes(2);
//Philippe 23/4/99
//      nodes[0] = side-4;//0,1,2,3
//      nodes[1] = (side-3)%4;//1,2,3,0
      nodes[0] = NodeIndex(side-4);//0,1,2,3
      nodes[1] = NodeIndex((side-3)%4);//1,2,3,0
      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,side-4));
      TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(this,(side-3)%4));
      TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(this,side));
      int64_t index;
      return gel->CreateCompEl(cmesh,index);
   } else PZError << "TPZGeoElQ2d::CreateBCCompEl. Side = " << side <<std::endl;
   return 0;
}

TPZGeoElSide TPZGeoElQ2d::Father2(int side){//Augusto:09/01/01
	TPZGeoEl *fFather = TPZGeoEl::Father();

		if(side<0 || side>8 || fFather==0){
		PZError << "TPZGeoEQ2d::Father2 called error" <<std::endl;
	}
	if(fFather->SubElement(0)==this){//fSubEl[0]
		if(side==0) return TPZGeoElSide(fFather,0);
		if(side==1 || side==4)	return TPZGeoElSide(fFather,4);
		if(side==3 || side==7)	return TPZGeoElSide(fFather,7);
		return TPZGeoElSide(fFather,8);
	}
	if(fFather->SubElement(1)==this){//fSubEl[1]
		if(side==1) return TPZGeoElSide(fFather,1);
		if(side==0 || side==4)	return TPZGeoElSide(fFather,4);
		if(side==2 || side==5)	return TPZGeoElSide(fFather,5);
		return TPZGeoElSide(fFather,8);
	}
	if(fFather->SubElement(2)==this){//fSubEl[2]
		if(side==2) return TPZGeoElSide(fFather,2);
		if(side==1 || side==5)	return TPZGeoElSide(fFather,5);
		if(side==3 || side==6)	return TPZGeoElSide(fFather,6);
		return TPZGeoElSide(fFather,8);
	}
	if(fFather->SubElement(3)==this){//fSubEl[3]
		if(side==3) return TPZGeoElSide(fFather,3);
		if(side==2 || side==6)	return TPZGeoElSide(fFather,6);
		if(side==0 || side==7)	return TPZGeoElSide(fFather,7);
		return TPZGeoElSide(fFather,8);
	}
    return TPZGeoElSide();
}

static int subeldata[9][9][2] =
{
	{{0,0}},/*side=0 {isub0{0,1},isub1{0,1},isub2{0,1},...}*/
	{{1,1}},/*side=1*/
	{{2,2}},/*side=2*/
	{{3,3}},/*side=3*/
	{{0,4},{0,1},{1,4}},/*side=4*/
	{{1,5},{1,2},{2,5}},/*side=5*/
	{{2,6},{2,3},{3,6}},/*side=6*/
	{{3,7},{3,0},{0,7}},/*side=7*/
	{{0,2},{0,8},{1,8},{2,8},{3,8},{0,5},{0,6},{2,4},{2,7}}/*side=8*/
};

static int nsubeldata[9] = {1,1,1,1,3,3,3,3,9};


void TPZGeoElQ2d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){//Augusto:09/01/01

  subel.Resize(0);
	if(side<0 || side>8 || !HasSubElement(side)){
		PZError << "TPZGeoElQ2d::GetSubElements2 called error" <<std::endl;
		return;
	}
	int nsub = nsubeldata[side];
	int isub;
	for(isub=0; isub<nsub; isub++) subel.Push(TPZGeoElSide(fSubEl[subeldata[side][isub][0]],
		                                                            subeldata[side][isub][1]));
}

static REAL buildt[4][9][3][2] = {
/*S0*/
     {/*0*/{{0.,0.},{0.,0.},{-1.,-1.}},//por colunas
      /*1*/{{0.,0.},{0.,0.},{0.,0.}},
      /*2*/{{0.,0.},{0.,0.},{0.,0.}},
      /*3*/{{0.,0.},{0.,0.},{0.,0.}},
      /*4*/{{.5,0.},{0.,0.},{-.5,0.}},
      /*5*/{{0.,.5},{0.,0.},{0.,-.5}},
      /*6*/{{-.5,0.},{0.,0.},{-.5,0.}},
      /*7*/{{.5,0.},{0.,0.},{.5,0.}},
      /*8*/{{.5,0.},{0.,.5},{-.5,-.5}}},
/*S1*/
     {/*0*/{{0.,0.},{0.,0.},{0.,0.}},
      /*1*/{{0.,0.},{0.,0.},{1.,-1.}},
      /*2*/{{0.,0.},{0.,0.},{0.,0.}},
      /*3*/{{0.,0.},{0.,0.},{0.,0.}},
      /*4*/{{.5,0.},{0.,0.},{.5,0.}},
      /*5*/{{.5,0.},{0.,0.},{-.5,0.}},
      /*6*/{{-.5,0.},{0.,0.},{.5,0.}},
      /*7*/{{0.,-.5},{0.,0.},{0.,-.5}},
      /*8*/{{.5,0},{0.,.5},{.5,-.5}}},
/*S2*/
     {/*0*/{{0.,0.},{0.,0.},{0.,0.}},
      /*1*/{{0.,0.},{0.,0.},{0.,0.}},
      /*2*/{{0.,0.},{0.,0.},{1.,1.}},
      /*3*/{{0.,0.},{0.,0.},{0.,0.}},
      /*4*/{{.5,0.},{0.,0.},{.5,0.}},
      /*5*/{{.5,0.},{0.,0.},{.5,0.}},
      /*6*/{{.5,0.},{0.,0.},{-.5,0.}},
      /*7*/{{0.,-.5},{0.,0.},{0.,.5}},
      /*8*/{{.5,0.},{0.,.5},{.5,.5}}},
/*S3*/
     {/*0*/{{0.,0.},{0.,0.},{0.,0.}},
      /*1*/{{0.,0.},{0.,0.},{0.,0.}},
      /*2*/{{0.,0.},{0.,0.},{0.,0.}},
      /*3*/{{0.,0.},{0.,0.},{-1.,1.}},
      /*4*/{{.5,0.},{0.,0.},{-.5,0.}},
      /*5*/{{0.,.5},{0.,0.},{0.,.5}},
      /*6*/{{.5,0.},{0.,0.},{.5,0.}},
      /*7*/{{.5,0.},{0.,0.},{-.5,0.}},
      /*8*/{{.5,0.},{0.,.5},{-.5,.5}}}
};

TPZTransform<STATE> TPZGeoElQ2d::BuildTransform2(int side, TPZGeoEl * /*father*/){//Augusto:09/01/01


	if(side<0 || side>9 || !Father(side).Element()){
  	PZError << "TPZGeoElQ2d::BuildTransform2 side out of range or father null\n";
    return TPZTransform<STATE>(0,0);
  }
  TPZTransform<STATE> trans(2,2);
  int son = WhichSubel();
  trans.Mult()(0,0) = buildt[son][side][0][0];
  trans.Mult()(1,0) = buildt[son][side][0][1];
  trans.Mult()(0,1) = buildt[son][side][1][0];
  trans.Mult()(1,1) = buildt[son][side][1][1];
  trans.Sum() (0,0) = buildt[son][side][2][0];
  trans.Sum() (1,0) = buildt[son][side][2][1];

  return trans;
}

static REAL MidSideNode[9][3] = {
/*00*/{-1.,-1.},/*01*/{ 1.,-1.},/*02*/{1.,1.},
/*03*/{-1., 1.},/*04*/{ 0.,-1.},/*05*/{1.,0.},
/*06*/{ 0., 1.},/*07*/{-1., 0.},/*08*/{0.,0.} };

int TPZGeoElQ2d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3), pss(3), pf(3), pfs(3);
  //point son, point side son, point father, point side father : elemento mestre
  pss[1] = 0.; pss[2] = 0.;//1d e 2d
  pfs[1] = 0.; pfs[2] = 0.;
  pf[1] = 0.; pf[2] = 0.;
  for(sn=0;sn<4;sn++){
    TPZGeoEl *son = subs[sn];
    for(sd=0;sd<9;sd++){
      ps[0] = MidSideNode[sd][0];//element
      ps[1] = MidSideNode[sd][1];//master point
      ps[2] = MidSideNode[sd][2];// = 0
      if(son->WhichSide(ps) != sd) std::cout << "Lado nao bate\n";
      TPZTransform<STATE> telsd = pzshape::TPZShapeQuad::TransformElementToSide(sd);//2x2
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform<STATE> t;
	  son->BuildTransform2(sd, gel,t);
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = pzshape::TPZShapeQuad::TransformSideToElement(sdfat);//2x2
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(8).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao furada\n";
        PZError << "son    = " << (son->Id()) <<std::endl;
        PZError << "father = " << ((son->Father2(8).Element())->Id()) <<std::endl;
        PZError << "side   = " << sd <<std::endl <<std::endl;
        int ok;
		std::cin >> ok;
      } else {
		  std::cout << "Transformacao OK!\n";
		  std::cout << "Filho/lado : " << son->Id() << "/" << sd <<std::endl;
		  std::cout << "Pai : " << son->Father2(8).Element()->Id() <<std::endl <<std::endl;
      }
    }
  }
  return 1;
}

REAL TPZGeoElQ2d::Mesure(int dim) {
  if(dim!=2) return 0.;
  REAL fMesure = 0.;
  if(IsZero(Volume())) {
    TPZGeoNode &nod1 = Mesh()->NodeVec()[fNodeIndexes[0]];
    REAL xx, yy, x0 = nod1.Coord(0), y0 = nod1.Coord(1);
    TPZGeoNode &nod2 = Mesh()->NodeVec()[fNodeIndexes[1]];
    xx = x0 - nod2.Coord(0);
    yy = nod2.Coord(1) - y0;
    TPZGeoNode &nod3 = Mesh()->NodeVec()[fNodeIndexes[2]];
    xx *= (nod3.Coord(1) - y0);
    yy *= (nod3.Coord(0) - x0);
    fMesure = fabs(xx + yy);
  }
  return fMesure;
}

void TPZGeoElQ2d::Center(TPZVec<REAL> &center) {
  center[0] = center[1] = 0.;
}

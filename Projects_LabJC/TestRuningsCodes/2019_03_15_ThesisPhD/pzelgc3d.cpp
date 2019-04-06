//METHODS DEFINITION FOR CLASS ELEMGC3D
#include "pzelgc3d.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelc1d.h"
#include "pzelgq2d.h"
#include "pzelcc3d.h"
#include "pzshapecube.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgnode.h"
#include "pzshtmat.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include <stdlib.h>

static TPZCompEl *CreateEl(TPZGeoElC3d *gel,TPZCompMesh &mesh, int64_t &index) {
  return new TPZCompElC3d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElC3d::fp)(TPZGeoElC3d *,TPZCompMesh &, int64_t &) = CreateEl;

TPZGeoElC3d::TPZGeoElC3d(int id,TPZVec<int64_t> &nodeindexes,int matid,TPZGeoMesh &mesh):
	TPZGeoElRefLess<pzgeom::TPZGeoCube>(id,nodeindexes,matid,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=8) {
    PZError << "TPZGeoElC3d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<8;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<8;i++) fSubEl[i] = 0;
}

TPZGeoElC3d::TPZGeoElC3d(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) :
	TPZGeoElRefLess<pzgeom::TPZGeoCube>(nodeindexes,matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=8) {
    PZError << "TPZGeoElC3d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<8;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<8;i++) fSubEl[i] = 0;
}

TPZGeoElC3d::~TPZGeoElC3d() {}

void TPZGeoElC3d::Shape(TPZVec<REAL> &pt,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {

      REAL x[2],dx[2],y[2],dy[2],z[2],dz[2];
      x[0] = (1.-pt[0])/2.;
      x[1] = (1.+pt[0])/2.;
      dx[0] = -0.5;
      dx[1] = 0.5;
      y[0] = (1.-pt[1])/2.;
      y[1] = (1.+pt[1])/2.;
      dy[0] = -0.5;
      dy[1] = 0.5;
      z[0] = (1.-pt[2])/2.;
      z[1] = (1.+pt[2])/2.;
      dz[0] = -0.5;
      dz[1] = 0.5;

      phi(0,0) = x[0]*y[0]*z[0];
      phi(1,0) = x[1]*y[0]*z[0];
      phi(2,0) = x[1]*y[1]*z[0];
      phi(3,0) = x[0]*y[1]*z[0];
      phi(4,0) = x[0]*y[0]*z[1];
      phi(5,0) = x[1]*y[0]*z[1];
      phi(6,0) = x[1]*y[1]*z[1];
      phi(7,0) = x[0]*y[1]*z[1];
      dphi(0,0) = dx[0]*y[0]*z[0];
      dphi(1,0) = x[0]*dy[0]*z[0];
      dphi(2,0) = x[0]*y[0]*dz[0];
      dphi(0,1) = dx[1]*y[0]*z[0];
      dphi(1,1) = x[1]*dy[0]*z[0];
      dphi(2,1) = x[1]*y[0]*dz[0];
      dphi(0,2) = dx[1]*y[1]*z[0];
      dphi(1,2) = x[1]*dy[1]*z[0];
      dphi(2,2) = x[1]*y[1]*dz[0];
      dphi(0,3) = dx[0]*y[1]*z[0];
      dphi(1,3) = x[0]*dy[1]*z[0];
      dphi(2,3) = x[0]*y[1]*dz[0];
      dphi(0,4) = dx[0]*y[0]*z[1];
      dphi(1,4) = x[0]*dy[0]*z[1];
      dphi(2,4) = x[0]*y[0]*dz[1];
      dphi(0,5) = dx[1]*y[0]*z[1];
      dphi(1,5) = x[1]*dy[0]*z[1];
      dphi(2,5) = x[1]*y[0]*dz[1];
      dphi(0,6) = dx[1]*y[1]*z[1];
      dphi(1,6) = x[1]*dy[1]*z[1];
      dphi(2,6) = x[1]*y[1]*dz[1];
      dphi(0,7) = dx[0]*y[1]*z[1];
      dphi(1,7) = x[0]*dy[1]*z[1];
      dphi(2,7) = x[0]*y[1]*dz[1];
}

TPZGeoElC3d *TPZGeoElC3d::CreateGeoEl(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElC3d(nodeindexes,matind,mesh);
}

int TPZGeoElC3d::NNodes() {
  return 8;
}
int64_t TPZGeoElC3d::NodeIndex(int node) {
  if(node<0 || node>7) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElC3d::NSideNodes(int side) {
  if(side<0 || side>26) {
    PZError << "TPZGeoElC3d::NSideNodes. Bad parameter side.\n";
    return 0;
  }
  if(side<8) return 1;//cantos
  if(side>7 && side<20) return 2;//lados
  if(side<26) return 4;//faces
  return 8;//centro
}

int64_t TPZGeoElC3d::SideNodeIndex(int side,int node) {
  if(side<0 || side>26 || node<0) {//27 sides
    PZError << "TPZGeoElC3d::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  //0<=side<=26 e node>=0
  //side central
  if(side==26 && node<8) return fNodeIndexes[node];//um dos 8 cantos
  //side de canto
  if(side < 8) return fNodeIndexes[side];//canto side, deve ser node = 0
  //sides 8,..,26
  if(side>7 && side<20 && node<2) {//lados 8 a 19
      side-=8;
      return fNodeIndexes[TPZCompElC3d::SideNodes[side][node]];
  } else if(side>19 && node<4) {//faces 20 a 25
  		side-=20;
	  	return fNodeIndexes[TPZCompElC3d::FaceNodes[side][node]];
  }
  return -1;
}

void TPZGeoElC3d::MidSideNodeIndex(int side, int64_t &index) {
  index = -1;
  if(side<0 || side>26) {
    PZError << "TPZGeoElC3d::MidSideNodeIndex. Bad parameter side = " << side <<std::endl;
    return;
  }
  //sides 0 a 7
	if(side<8) {//o nó medio do lado 0 é o 0 etc.
		index=fNodeIndexes[side];
		return;
	}
   //o nó medio da face é o centro e o nó medio do centro é o centro
   //como nó de algum filho se este existir
   //caso tenha filhos é o canto de algum filho, se não tiver filhos retorna -1
	if(HasSubElement(0)) {
   	side-=8;
	   index=((TPZGeoElC3d *) fSubEl[TPZCompElC3d::MidSideNodes[side][0]])->fNodeIndexes[TPZCompElC3d::MidSideNodes[side][1]];
	}
}

void TPZGeoElC3d::NewMidSideNode(int side, int64_t &index) {
  MidSideNodeIndex(side,index);
  if(index < 0) {
    TPZGeoElSide gelside = Neighbour(side);
    while(gelside.Element()) {
      gelside.Element()->MidSideNodeIndex(gelside.Side(),index);
      if(index!=-1) return;
      gelside = gelside.Neighbour();
      if(gelside.Element()==this) break;
    }
    TPZVec<REAL> par(3,0.);
    TPZVec<REAL> coord(3,0.);
    if(side < 8) {index = -1; return;}
    //aqui side = 8 a 26
    side-=8;//0,1,..,18
    par[0] = TPZCompElC3d::MidCoord[side][0];
    par[1] = TPZCompElC3d::MidCoord[side][1];
    par[2] = TPZCompElC3d::MidCoord[side][2];
    X(par,coord);
    index = Mesh()->NodeVec().AllocateNewElement();
    Mesh()->NodeVec()[index].Initialize(coord,*Mesh());
  }
}

int TPZGeoElC3d::SideDimension(int side) {

	if (side<0 || side>26) {
   	PZError << "TPZGeoElC3d::SideDimension called with side " << side <<std::endl;
      return 0;
   }
  if(side<8) return 0;//cantos
  if(side>7 && side<20) return 1;//lados
  if(side<26) return 2;//faces
  return 3;//centro

   }

TPZGeoElSide TPZGeoElC3d::HigherDimensionSides(int side,int targetdimension) {
//targetdimension deve ser 1 , 2 ou 3
//se side =  0 a  7 targetdimension deve ser 1
//se side =  8 a 19 targetdimension deve ser 2
//se side = 26      targetdimension deve ser 3
  if( (side<0 || side>25) || (targetdimension < 1 || targetdimension > 3) ) {
     PZError << "TPZGeoElC3d::HigherDimensionSides called with side = " << side
	          << " targetdimension = " << targetdimension <<std::endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
  int bestface;//Cedric 26/05/99
  //side = 0 a 25
  switch(targetdimension) {//=1,2
	  case 1:
     	 if(this == father->SubElement(0)) {
       	 if(side==1) return TPZGeoElSide(this,8);
       	 if(side==3) return TPZGeoElSide(this,11);
          if(side==4) return TPZGeoElSide(this,12);
       } else if(this == father->SubElement(1)) {
       	 if(side==0) return TPZGeoElSide(this,8);
       	 if(side==2) return TPZGeoElSide(this,9);
          if(side==5) return TPZGeoElSide(this,13);
       } else if(this == father->SubElement(2)) {
       	 if(side==1) return TPZGeoElSide(this,9);
       	 if(side==3) return TPZGeoElSide(this,10);
          if(side==6) return TPZGeoElSide(this,14);
       } else if(this == father->SubElement(3)) {
       	 if(side==0) return TPZGeoElSide(this,11);
       	 if(side==2) return TPZGeoElSide(this,10);
          if(side==7) return TPZGeoElSide(this,15);
       } else if(this == father->SubElement(4)) {
       	 if(side==0) return TPZGeoElSide(this,12);
       	 if(side==5) return TPZGeoElSide(this,16);
          if(side==7) return TPZGeoElSide(this,19);
       } else if(this == father->SubElement(5)) {
       	 if(side==1) return TPZGeoElSide(this,13);
       	 if(side==4) return TPZGeoElSide(this,16);
          if(side==6) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(6)) {
       	 if(side==2) return TPZGeoElSide(this,14);
       	 if(side==5) return TPZGeoElSide(this,17);
          if(side==7) return TPZGeoElSide(this,18);
       } else if(this == father->SubElement(7)) {
       	 if(side==3) return TPZGeoElSide(this,15);
       	 if(side==4) return TPZGeoElSide(this,19);
          if(side==6) return TPZGeoElSide(this,18);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
  	  case 2:
     	 if(this == father->SubElement(0)) {

          bestface = BestDimensionSideOfTwoFaces(20,21);
          if((side== 1 || side== 8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(20,24);
          if((side== 3 || side==11) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(21,24);
          if((side== 4 || side==12) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 2) return TPZGeoElSide(this,20);
       	 if(side== 9) return TPZGeoElSide(this,20);
       	 if(side==10) return TPZGeoElSide(this,20);
       	 if(side== 5) return TPZGeoElSide(this,21);
       	 if(side==13) return TPZGeoElSide(this,21);
       	 if(side==16) return TPZGeoElSide(this,21);
          if(side== 7) return TPZGeoElSide(this,24);
          if(side==15) return TPZGeoElSide(this,24);
          if(side==19) return TPZGeoElSide(this,24);
       } else if(this == father->SubElement(1)) {

          bestface = BestDimensionSideOfTwoFaces(20,21);
          if((side== 0 || side== 8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(20,22);
          if((side== 2 || side== 9) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(21,22);
          if((side== 5 || side==13) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 3) return TPZGeoElSide(this,20);
       	 if(side==10) return TPZGeoElSide(this,20);
       	 if(side==11) return TPZGeoElSide(this,20);
       	 if(side== 4) return TPZGeoElSide(this,21);
       	 if(side==12) return TPZGeoElSide(this,21);
       	 if(side==16) return TPZGeoElSide(this,21);
          if(side== 6) return TPZGeoElSide(this,22);
          if(side==14) return TPZGeoElSide(this,22);
          if(side==17) return TPZGeoElSide(this,22);
       } else if(this == father->SubElement(2)) {

          bestface = BestDimensionSideOfTwoFaces(20,22);
          if((side== 1 || side== 9) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(20,23);
          if((side== 3 || side==10) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(22,23);
          if((side== 6 || side==14) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 0) return TPZGeoElSide(this,20);
       	 if(side== 8) return TPZGeoElSide(this,20);
       	 if(side==11) return TPZGeoElSide(this,20);
       	 if(side== 5) return TPZGeoElSide(this,22);
       	 if(side==13) return TPZGeoElSide(this,22);
       	 if(side==17) return TPZGeoElSide(this,22);
          if(side== 7) return TPZGeoElSide(this,23);
          if(side==15) return TPZGeoElSide(this,23);
          if(side==18) return TPZGeoElSide(this,23);
       } else if(this == father->SubElement(3)) {

          bestface = BestDimensionSideOfTwoFaces(20,23);
          if((side== 2 || side==10) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(20,24);
          if((side== 0 || side==11) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(23,24);
          if((side== 7 || side==15) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 1) return TPZGeoElSide(this,20);
       	 if(side== 8) return TPZGeoElSide(this,20);
       	 if(side== 9) return TPZGeoElSide(this,20);
          if(side== 6) return TPZGeoElSide(this,23);
          if(side==14) return TPZGeoElSide(this,23);
          if(side==18) return TPZGeoElSide(this,23);
       	 if(side== 4) return TPZGeoElSide(this,24);
       	 if(side==12) return TPZGeoElSide(this,24);
       	 if(side==19) return TPZGeoElSide(this,24);
       } else if(this == father->SubElement(4)) {

          bestface = BestDimensionSideOfTwoFaces(21,24);
          if((side== 0 || side==12) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(21,25);
          if((side== 5 || side==16) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(24,25);
          if((side== 7 || side==19) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 1) return TPZGeoElSide(this,21);
       	 if(side== 8) return TPZGeoElSide(this,21);
       	 if(side==13) return TPZGeoElSide(this,21);
       	 if(side== 3) return TPZGeoElSide(this,24);
       	 if(side==11) return TPZGeoElSide(this,24);
       	 if(side==15) return TPZGeoElSide(this,24);
          if(side== 6) return TPZGeoElSide(this,25);
          if(side==17) return TPZGeoElSide(this,25);
          if(side==18) return TPZGeoElSide(this,25);
       } else if(this == father->SubElement(5)) {

          bestface = BestDimensionSideOfTwoFaces(21,22);
          if((side== 1 || side==13) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(21,25);
          if((side== 4 || side==16) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(22,25);
          if((side== 6 || side==17) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 0) return TPZGeoElSide(this,21);
       	 if(side== 8) return TPZGeoElSide(this,21);
       	 if(side==12) return TPZGeoElSide(this,21);
       	 if(side== 2) return TPZGeoElSide(this,22);
       	 if(side== 9) return TPZGeoElSide(this,22);
       	 if(side==14) return TPZGeoElSide(this,22);
          if(side== 7) return TPZGeoElSide(this,25);
          if(side==18) return TPZGeoElSide(this,25);
          if(side==19) return TPZGeoElSide(this,25);
       } else if(this == father->SubElement(6)) {

          bestface = BestDimensionSideOfTwoFaces(22,23);
          if((side== 2 || side==14) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(22,25);
          if((side== 5 || side==17) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(23,25);
          if((side== 7 || side==18) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 1) return TPZGeoElSide(this,22);
       	 if(side== 9) return TPZGeoElSide(this,22);
       	 if(side==13) return TPZGeoElSide(this,22);
       	 if(side== 3) return TPZGeoElSide(this,23);
       	 if(side==10) return TPZGeoElSide(this,23);
       	 if(side==15) return TPZGeoElSide(this,23);
          if(side== 4) return TPZGeoElSide(this,25);
          if(side==16) return TPZGeoElSide(this,25);
          if(side==19) return TPZGeoElSide(this,25);
       } else if(this == father->SubElement(7)) {

          bestface = BestDimensionSideOfTwoFaces(23,24);
          if((side== 3 || side==15) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(23,25);
          if((side== 6 || side==18) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(24,25);
          if((side== 4 || side==19) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 0) return TPZGeoElSide(this,24);
       	 if(side==11) return TPZGeoElSide(this,24);
       	 if(side==12) return TPZGeoElSide(this,24);
       	 if(side== 2) return TPZGeoElSide(this,23);
       	 if(side==10) return TPZGeoElSide(this,23);
       	 if(side==14) return TPZGeoElSide(this,23);
          if(side== 5) return TPZGeoElSide(this,25);
          if(side==16) return TPZGeoElSide(this,25);
          if(side==17) return TPZGeoElSide(this,25);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
     case 3:
       return TPZGeoElSide(this,26);//0<=side<=25
  }//switch
  return TPZGeoElSide();
}

int TPZGeoElC3d::BestDimensionSideOfTwoFaces(int face1,int face2) {

    TPZGeoElSide grandfather = Father(face1);
    int levnum1=0,levnum2=0;
    while(grandfather.Element()) {
      levnum1++;
      grandfather = grandfather.Father2();//face1
    }
    grandfather = Father(face2);
    while(grandfather.Element()) {
      levnum2++;
      grandfather = grandfather.Father2();//face2
    }
    if(levnum1 > levnum2) return face1;
    if(levnum2 > levnum1) return face2;
    return 0;
}


void TPZGeoElC3d::LowerDimensionSides(int side,TPZStack<TPZGeoElSide> &smallsides) {
	if (side < 8) return;
   int i;
   if(side < 20) {//side = 8 a 19 : entram os cantos dos lados
   	int s = side-8;
   	smallsides.Push(TPZGeoElSide(this,TPZCompElC3d::SideNodes[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElC3d::SideNodes[s][1]));
   } else if(side < 26) {//side = 20 a 25
   	int s = side-20;   //entram cantos e lados da face
   	smallsides.Push(TPZGeoElSide(this,TPZCompElC3d::FaceNodes[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElC3d::FaceNodes[s][1]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElC3d::FaceNodes[s][2]));
      smallsides.Push(TPZGeoElSide(this,TPZCompElC3d::FaceNodes[s][3]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElC3d::FaceSides[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElC3d::FaceSides[s][1]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElC3d::FaceSides[s][2]));
      smallsides.Push(TPZGeoElSide(this,TPZCompElC3d::FaceSides[s][3]));
      smallsides.Push(TPZGeoElSide(this,side));
   } else if(side==26) {//entram todos os cantos, arestas e faces
   	for (i=0;i<26;i++) smallsides.Push(TPZGeoElSide(this,i));
   }
}

void TPZGeoElC3d::Jacobian(TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv){

  //int nnodes = NNodes();
#ifdef DEBUG
  if (nnodes != 8) {
    PZError << "TPZGeoElC3d.jacobian only implemented for"
      " 8 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 3 || param[0] < -1. || param[0] > 1. ||
     param[1] < -1. || param[1] > 1. || param[2] < -1. || param[2] > 1.) {
    PZError << "TPZGeoElC3d.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
    return;
  }
#endif
  REAL spacephi[12];
  TPZFMatrix<STATE> phi(8,1,spacephi,12);
  REAL spacedphi[30];
  TPZFMatrix<STATE> dphi(3,8,spacedphi,30);
  Shape(param,phi,dphi);
  jacobian.Zero();
  TPZGeoNode *np;
  int i,j;
  for(i=0;i<8;i++) {
    np = NodePtr(i);
    for(j=0;j<3;j++) {
      jacobian(0,j) += np->Coord(j)*dphi(0,i);
      jacobian(1,j) += np->Coord(j)*dphi(1,i);
      jacobian(2,j) += np->Coord(j)*dphi(2,i);
    }
  }

  detjac = -jacobian(0,2)*jacobian(1,1)*jacobian(2,0);//-a02 a11 a20
  detjac += jacobian(0,1)*jacobian(1,2)*jacobian(2,0);//+ a01 a12 a20
  detjac += jacobian(0,2)*jacobian(1,0)*jacobian(2,1);//+ a02 a10 a21
  detjac -= jacobian(0,0)*jacobian(1,2)*jacobian(2,1);//- a00 a12 a21
  detjac -= jacobian(0,1)*jacobian(1,0)*jacobian(2,2);//- a01 a10 a22
  detjac += jacobian(0,0)*jacobian(1,1)*jacobian(2,2);//+ a00 a11 a22

  //jacinv.Zero();
  jacinv(0,0) = (-jacobian(1,2)*jacobian(2,1)+jacobian(1,1)*jacobian(2,2))/detjac;//-a12 a21 + a11 a22
  jacinv(0,1) = ( jacobian(0,2)*jacobian(2,1)-jacobian(0,1)*jacobian(2,2))/detjac;//a02 a21 - a01 a22
  jacinv(0,2) = (-jacobian(0,2)*jacobian(1,1)+jacobian(0,1)*jacobian(1,2))/detjac;//-a02 a11 + a01 a12
  jacinv(1,0) = ( jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2))/detjac;//a12 a20 - a10 a22
  jacinv(1,1) = (-jacobian(0,2)*jacobian(2,0)+jacobian(0,0)*jacobian(2,2))/detjac;//-a02 a20 + a00 a22
  jacinv(1,2) = ( jacobian(0,2)*jacobian(1,0)-jacobian(0,0)*jacobian(1,2))/detjac;//a02 a10 - a00 a12
  jacinv(2,0) = (-jacobian(1,1)*jacobian(2,0)+jacobian(1,0)*jacobian(2,1))/detjac;//-a11 a20 + a10 a21
  jacinv(2,1) = ( jacobian(0,1)*jacobian(2,0)-jacobian(0,0)*jacobian(2,1))/detjac;//a01 a20 - a00 a21
  jacinv(2,2) = (-jacobian(0,1)*jacobian(1,0)+jacobian(0,0)*jacobian(1,1))/detjac;//-a01 a10 + a00 a11

  axes.Zero();
  axes(0,0) = 1.;
  axes(1,1) = 1.;
  axes(2,2) = 1.;
}

void TPZGeoElC3d::X(TPZVec<REAL> & loc,TPZVec<REAL> &result){
  REAL spacephi[12],spacedphi[30];
  int i,j;
  TPZFMatrix<STATE> phi(8,1,spacephi,12);
  TPZFMatrix<STATE> dphi(3,8,spacedphi,30);
  Shape(loc,phi,dphi);
  for(j=0;j<3;j++) {
    result[j] = 0.0;
    for(i=0;i<8;i++) result[j] += NodePtr(i)->Coord(j)*phi(i,0);
  }
}
/**It's necessary to define the normal vector to side 4, that is the orthogonal
   vector to the surface*/

void TPZGeoElC3d::NormalVector(int side,TPZVec<REAL> &param,TPZVec<REAL> &normal,
			     TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &jac1d) {

#ifdef DEBUG
  if (nnodes != 8) {
    PZError << "TPZGeoElC3d.NormalVector, only implemented for"
      " 8 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 3 || param[0] < -1. || param[0] > 1. ||
     param[1] < -1. || param[1] > 1. || param[2] < -1. || param[2] > 1.) {
    PZError << "TPZGeoElC3d.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
    return;
  }
  if(normal.NElements() != 3) {
    PZError << "TPZGeoElC3d::NormalVector normal.capacity() = " << normal.NElements() <<
      "\n";
    return;
  }
  if(side < 0 || side > 5) {//6 faces
    PZError << "TPZGeoElC3d.jacobian invalid side : "
      " side = " << side << "\n";
    return;
  }
#endif

  REAL spacephi[12],spacedphi[30];
  TPZFMatrix<STATE> phi(8,1,spacephi,12);
  TPZFMatrix<STATE> dphi(3,8,spacedphi,30);
  Shape(param,phi,dphi);
  TPZGeoNode *np;
  TPZVec<REAL> n(3,0.);
  int i,j,ider;
  if(side==0 || side==5) ider = 2;
  if(side==1 || side==3) ider = 1;
  if(side==2 || side==4) ider = 0;
  for(i=0;i<4;i++) {
    np = NodePtr(TPZCompElC3d::FaceNodes[side][i]);
    for(j=0;j<3;j++)
      n[j] += np->Coord(j)*dphi(ider,i);
  }

  jac1d(0,0) = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );

  for(i=0;i<3;i++) normal[i] = n[i]/jac1d(0,0);

  switch(side) {
     case 0: for(i=0;i<3;i++) normal[i] *= -1.;
		       break;
     case 1: for(i=0;i<3;i++) normal[i] *= -1.;
             break;
     case 4: for(i=0;i<3;i++) normal[i] *= -1.;
  }
  axes.Zero();
  axes(0,0) = 1.0;
  axes(1,1) = 1.0;
  axes(2,2) = 1.0;
  return;
}

/** TO SUBDIVISION
********************************************************************************
  Into Divides is necesary to consider the connectivity with the all neighboards*/
void TPZGeoElC3d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {

   int i;
   if(HasSubElement(0)) {
      SubElVec.Resize(8);
      for(i=0;i<8;i++) SubElVec[i] = fSubEl[i];
      return;//If exist fSubEl return this sons
   }
   int j, sub, matid = MaterialId();
   int64_t index;
   int np[27];//guarda conectividades dos 8 subelementos

   for(j=0;j<8;j++) np[j]=NodeIndex(j);
   for(sub=8;sub<27;sub++) {
      NewMidSideNode(sub,index);
      np[sub] = index;
   }
   // creating new subelements
   for(i=0;i<8;i++) {
	   TPZManVector<int64_t> cornerindexes(8);
   	for(int j=0;j<8;j++) cornerindexes[j] = np[TPZCompElC3d::CornerSons[i][j]];
      fSubEl[i] = CreateGeoEl(cornerindexes,matid,*Mesh());
   }

   SubElVec.Resize(8);
   for(sub=0;sub<8;sub++) {
      SubElVec[sub] = fSubEl[sub];
      SubElVec[sub]->SetFather(this);
   }
   for(i=0;i<8;i++) {//conectividades entre os filhos : viz interna
   	for(j=0;j<19;j++) {        //lado do subel                                          numero do filho viz.             lado do viz.
      	fSubEl[i]->SetNeighbour(TPZCompElC3d::InNeigh[i][j][0],TPZGeoElSide(fSubEl[TPZCompElC3d::InNeigh[i][j][1]],TPZCompElC3d::InNeigh[i][j][2]));
      }
   }
   //vizinhança externa ao elemento atual
   //procura-se um viz pela face que seja dividido
   TPZGeoElSide dividedneighbour;
   for(int face=20;face<26;face++) {
      dividedneighbour = Neighbour(face);//BuildConnectivity(..) inicializou
      while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         //achou-se um viz para alguma face
         TPZManVector<int64_t> nodes(4);
         TPZStack<TPZGeoElSide> neighsub;
         int f = face-20;
         for(i=0;i<4;i++) nodes[i] = NodeIndex(TPZCompElC3d::FaceNodes[f][i]);//nós globais na ordem e sentido local da face atual
         dividedneighbour.GetSubElements2(neighsub);//4 subs viz na ordem dos subs da face atual
         nodes.Resize(1);
         //side central a face
         int fs0 = TPZCompElC3d::FaceSons[f][0];//1o filho da face
         int fn2 = TPZCompElC3d::FaceNodes[f][2];//nó 2 da face
         nodes[0] = fSubEl[fs0]->NodeIndex(fn2);//nó do centro da face como nó do filho nessa face
         int locside = neighsub[0].Element()->WhichSide(nodes);//lado do viz do sub conectado ao nó centro da face
         TPZGeoElSide sub(fSubEl[fs0],fn2);
         //conectividade do nó do centro da face fechando o ciclo entre os filhos dessa face conectados pelo nó
         sub.SetConnectivity(TPZGeoElSide(neighsub[0].Element(),locside));
         //conectividades dos lados dos subelementos interiores a face
         nodes.Resize(2);
         for(i=0; i<4; i++) {//4 subelementos associados à face
            int fsi = TPZCompElC3d::FaceSons[f][i];
            int inside = TPZCompElC3d::FaceSides[f][(i+1)%4];
            nodes[0] = fSubEl[fsi]->SideNodeIndex(inside,0);
            nodes[1] = fSubEl[fsi]->SideNodeIndex(inside,1);
            int locside = neighsub[i].Element()->WhichSide(nodes);
            TPZGeoElSide sub(fSubEl[fsi],inside);
            TPZGeoElSide subneigh(neighsub[i].Element(),locside);
            sub.SetConnectivity(subneigh);
         }
         //side de face para os 4 subs da face f
         for(i=0; i<4; i++) {//4 subelementos associados à face
            TPZGeoElSide sub(fSubEl[TPZCompElC3d::FaceSons[f][i]],face);
            sub.SetConnectivity(neighsub[i]);
         }
      }
   }//fim faces
	//procura-se um viz pelo lado do atual, que seja dividido
   for(int side=8; side<20; side++) {//(side 8 : subs 0,1),(side 9 : subs 1,2) ...
   	int s = side-8;
   	dividedneighbour = Neighbour(side);
      while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         //achou-se um viz pelo lado i e dividido
         int i = TPZCompElC3d::SideNodes[s][0];//nó local 0 do lado side
         int iplus = TPZCompElC3d::SideNodes[s][1];//nó local 1 do lado
         TPZManVector<int64_t> nodes(2);
         TPZStack<TPZGeoElSide> neighsub;
         nodes[0] = SideNodeIndex(side,0);//nó global 0 do lado side
         nodes[1] = SideNodeIndex(side,1);//nó global 1
         dividedneighbour.GetSubElements2(neighsub);//filhos do viz conectados nesses lados
         TPZGeoElSide sub(fSubEl[i],side);
         //conectividade pelo lado side entre os filhos
         sub.SetConnectivity(neighsub[0]);
         sub = TPZGeoElSide(fSubEl[iplus],side);
         sub.SetConnectivity(neighsub[1]);
         nodes.Resize(1);
         nodes[0] = fSubEl[i]->SideNodeIndex(iplus,0);//nó do medio do lado como nó de um dos filhos
         int locside = neighsub[0].Element()->WhichSide(nodes);
         sub = TPZGeoElSide(fSubEl[i],iplus);
         //conectividade para o nó do medio do lado side do atual como nó de um dos filhos
         sub.SetConnectivity(TPZGeoElSide(neighsub[0].Element(),locside));
      }
   }
   //procura-se um viz pelo canto do atual, que seja dividido
   for(int corner=0; corner<8; corner++) {//cantos
   	dividedneighbour = Neighbour(corner);
      while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
        //achou-se um viz por um dos cantos e que esteja dividido
         TPZManVector<int64_t> nodes(1);
         TPZStack<TPZGeoElSide> neighsub;
         nodes[0] = SideNodeIndex(corner,0);
         dividedneighbour.GetSubElements2(neighsub);
         TPZGeoElSide sub(fSubEl[corner],corner);
         sub.SetConnectivity(neighsub[0]);
      }
   }
}


int TPZGeoElC3d::NSubElements() {
  return 8;
}

int TPZGeoElC3d::NSideSubElements(int side) {
  if(side < 0 || side > 26) {
    PZError << "TPZGeoElC3d::NSideSubElements called for side " << side <<std::endl;
    return 0;
  }
  if(side==26) return 8;//centro
  if(side>19 && side<26) return 4;//faces
  if(side>7) return 2;//lados
  return 1;//cantos
}

TPZGeoElSide TPZGeoElC3d::SideSubElement(int side,int position) {
   if (position<0 || position>8 || side <0 ||side>26) {
   	PZError << "TPZGeoElC3d::SideSubElement called with position " << position << " side " << side <<std::endl;
      return TPZGeoElSide();
   }                              //fSubEl[is]
   if(side==26) return TPZGeoElSide(SubElement(position),26);//centro
   if(side<8) {//cantos
      if(position!=0) {
         PZError << "TPZGeoElC3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
         return TPZGeoElSide(SubElement(side),side);
      }
   }
   if(side>7 && side<20) {//lados
       if(position!=0 && position!=1) {
         PZError << "TPZGeoElC3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
      	int s = side-8;
         return TPZGeoElSide(SubElement(TPZCompElC3d::SideNodes[s][position]),side);
      }
   }
   if(side>19) {//faces
       if(position<0 || position>3) {//position!=0 && position!=1 && position!=2 && position!=3
         PZError << "TPZGeoElC3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
      	int s = side-20;
         return TPZGeoElSide(SubElement(TPZCompElC3d::FaceNodes[s][position]),side);
      }
   }
   return TPZGeoElSide();
}

void TPZGeoElC3d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
   if(!fSubEl[0]) {
      sub.Resize(0);
      return;
   }
   if(side < 0 || side > 26) {
      PZError << "TPZGeoElC3d::SideSubElements called for side " << side <<std::endl;
      return;
   }
   if(side==26) {
      sub.Resize(8);
      for(int i=0;i<8;i++) sub[i] = fSubEl[i];
      return;
   }
   if(side<8) {
      sub.Resize(1);
      sub[0]=fSubEl[side];
      return;
   }
   if(side>7 && side<20) {//lados
      int s = side-8;
      sub.Resize(2);
      sub[0] = fSubEl[TPZCompElC3d::SideNodes[s][0]];
      sub[1] = fSubEl[TPZCompElC3d::SideNodes[s][1]];
      return;
   }
   if(side>19) {//faces
      int s = side-20;
      sub.Resize(4);
      sub[0] = fSubEl[TPZCompElC3d::FaceSons[s][0]];
      sub[1] = fSubEl[TPZCompElC3d::FaceSons[s][1]];
      sub[2] = fSubEl[TPZCompElC3d::FaceSons[s][2]];
      sub[3] = fSubEl[TPZCompElC3d::FaceSons[s][3]];
   }
}

TPZGeoElSide TPZGeoElC3d::Father(int side) {

	TPZGeoEl *fFather = TPZGeoEl::Father();
	if (!fFather) return TPZGeoElSide();
	int whichsub = -1;
   int i;
   for(i=0; i<8; i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) {//equivale a is = 4  ou is > 3
	   PZError << "TPZGeoElC3d::Father. fFather isn't father of this element.\n";
   	return TPZGeoElSide();
   }
   //agora o atual elemento é o filho numero whichsub < 8
   if(side == whichsub || side==26) return TPZGeoElSide(fFather,side);//cantos
   //lados
      if(whichsub == 0 && (side== 8 || side==11 || side==12)) return TPZGeoElSide(fFather,side);
      if(whichsub == 1 && (side== 8 || side== 9 || side==13)) return TPZGeoElSide(fFather,side);
      if(whichsub == 2 && (side== 9 || side==10 || side==14)) return TPZGeoElSide(fFather,side);
      if(whichsub == 3 && (side==10 || side==11 || side==15)) return TPZGeoElSide(fFather,side);
      if(whichsub == 4 && (side==12 || side==16 || side==19)) return TPZGeoElSide(fFather,side);
      if(whichsub == 5 && (side==13 || side==17 || side==16)) return TPZGeoElSide(fFather,side);
      if(whichsub == 6 && (side==14 || side==17 || side==18)) return TPZGeoElSide(fFather,side);
      if(whichsub == 7 && (side==15 || side==18 || side==19)) return TPZGeoElSide(fFather,side);
   //faces
   if(whichsub == 0 && (side==20 || side==21 || side==24)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side==20 || side==21 || side==22)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side==20 || side==22 || side==23)) return TPZGeoElSide(fFather,side);
   if(whichsub == 3 && (side==20 || side==23 || side==24)) return TPZGeoElSide(fFather,side);
   if(whichsub == 4 && (side==21 || side==24 || side==25)) return TPZGeoElSide(fFather,side);
   if(whichsub == 5 && (side==21 || side==22 || side==25)) return TPZGeoElSide(fFather,side);
   if(whichsub == 6 && (side==22 || side==23 || side==25)) return TPZGeoElSide(fFather,side);
   if(whichsub == 7 && (side==23 || side==24 || side==25)) return TPZGeoElSide(fFather,side);
   //outro caso
   return TPZGeoElSide();
}

void TPZGeoElC3d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {

   int nsub = NSideSubElements(side);
   if(!nsub) return;
   sub.Resize(nsub);
   int i,j,k;
   if(nsub==1) {//side = 0 a 7
   	if(fSubEl[side]->NodeIndex(side)!=refnode[0]) {
      	PZError << "TPZGeoElC3d::GetSubElement subelement does not contain refnode" <<std::endl;
         return;
      }
	   sub[0]=TPZGeoElSide(fSubEl[side],side);
   	return;
   }
   //int isub=0;
   for(i=0;i<nsub;i++) {
   	TPZGeoElSide sidesub = SideSubElement(side,i);
      TPZGeoEl *subel = sidesub.Element();
		for(k = 0; k < refnode.NElements(); k++) {
		   for(j=0;j<8;j++) {
			   if(subel->NodeIndex(j)==refnode[k]) {
            	sub[k] = SideSubElement(side,i);
            }
         }
      }
   }
   return;
}

/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
void TPZGeoElC3d::BuildTransform(int side,TPZGeoEl *father,TPZTransform<STATE> &t) {
	TPZGeoEl *fFather = TPZGeoEl::Father();

   if(!Father(side).Element()) return;
   int whichsub = -1;
   int i;
   for(i=0; i<8; i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) return;
   int dim = SideDimension(side);
   TPZTransform<STATE> tloc(dim);
   REAL store[12];
   TPZFMatrix<STATE> mult(dim,dim,store,9);
   TPZFMatrix<STATE> sum(dim,1,store+9,3);
   mult.Zero();
   sum.Zero();

   TPZGeoEl *locfather;
   if(father == this) {
	   mult(0,0) = mult(1,1) = mult(2,2) = 1.;//pai para pai = identidade
      locfather = fFather;
   } else {
	   locfather = fFather;//pai do lemento atual
   }
   if(!locfather) {
	   std::cout << "TPZGeoElC3d::BuildTransform could not identify the father element\n";
	   return;
   }

   if(side == 26) {//pai para filho
      mult(0,0) = 0.5;
      mult(1,1) = 0.5;
      mult(2,2) = 0.5;
      sum(0,0) = -0.5;
      sum(1,0) = -0.5;
      sum(2,0) = -0.5;
      switch(whichsub) {//o atual é o filho numero whichsub
         case 0:
            break;          //x<0 y<0 z<0
         case 1:
            sum(0,0) *= -1.;//x>0 y<0 z<0
            break;
         case 2:
            sum(0,0) *= -1.;//x>0 y>0 z<0
            sum(1,0) *= -1.;
            break;
         case 3:
            sum(1,0) *= -1.;//x<0 y>0 z<0
            break;
         case 4:
            sum(2,0) *= -1.;//x<0 y<0 z>0
            break;
         case 5:
         	sum(0,0) *= -1.;//x>0 y<0 z>0
            sum(2,0) *= -1.;
            break;
         case 6:
            sum(0,0) *= -1.;//x>0 y>0 z>0
            sum(1,0) *= -1.;
            sum(2,0) *= -1.;
            break;
         case 7:
            sum(1,0) *= -1.;//x<0 y>0 z>0
            sum(2,0) *= -1.;
      }
   } else
   if(side>19) {//face do pai para face do filho
      whichsub = -1;
      int s = side-20;
      for(i=0; i<4; i++) if(fFather->SubElement(TPZCompElC3d::FaceSons[s][i]) == this) whichsub = i;
      if(whichsub == -1) return;
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
   } else
   if(side>7) {//e side < 20 : arestas
      whichsub = -1;
      int s = side-8;//0,1,2,3,4,5,6,7,8,9,10,11
      for(i=0; i<2; i++)
      	if(Father(side).Element()->SubElement(TPZCompElC3d::SideNodes[s][i]) == this) whichsub = i;
      if(whichsub == -1) return;
   	mult(0,0) = 0.5;
      if(whichsub==0) sum(0,0) = -0.5;
      if(whichsub==1) sum(0,0) =  0.5;
   }
   tloc.SetMatrix(mult,sum);
   t = tloc.Multiply(t);
   if(locfather != father) locfather->BuildTransform2(side,father,t);
}

TPZTransform<STATE> TPZGeoElC3d::SideToSideTransform(int sidefrom,int sideto) {
   if( (sidefrom > 19 &&  sideto < 26) || (sidefrom > 7 &&  sideto < 20) || sideto < 8) {
      PZError << "TPZGeoElC3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   int dimfrom = SideDimension(sidefrom);
   int dimto = SideDimension(sideto);
   if(dimfrom >= dimto){
      PZError << "TPZGeoElC3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   TPZTransform<STATE> trans(dimto,dimfrom);//retorna zerada
   //agora : sidefrom < sideto ou dimfrom < dimto
 if(sideto == 26) {//interior
   switch(sidefrom) {
      //faces para interior
               //row = dimto = 3 e col = dimfrom = 2
      case  20://sidefrom face de dim = 2 : t(3,2)
      	trans.Mult()(0,0) =  1.;//         fMult(row,col) : fMult(3,2)
         trans.Mult()(1,1) =  1.;//         fSum(row,1)    : fSum(3,1)
         trans.Sum()(2,0)  = -1.;//sides : 0123 -> 0123
         break;
      case  21:
      	trans.Mult()(0,0) =  1.;//sides : 0123 -> 0154
         trans.Mult()(2,1) =  1.;//         R2  -> R3
         trans.Sum()(1,0)  = -1.;//         2D -> 3D
         break;
      case  22:
      	trans.Mult()(1,0) =  1.;//sides : 0123 -> 1265
         trans.Mult()(2,1) =  1.;
         trans.Sum()(0,0)  =  1.;
         break;
      case  23:
      	trans.Mult()(0,0) =  1.;//sides : 0123 -> 3267
         trans.Mult()(2,1) =  1.;
         trans.Sum()(1,0)  =  1.;
         break;
      case  24:
      	trans.Mult()(1,0) =  1.;//sides : 0123 -> 0374
         trans.Mult()(2,1) =  1.;
         trans.Sum()(0,0)  = -1.;
         break;
      case  25:
      	trans.Mult()(0,0) =  1.;//sides : 0123 -> 4567
         trans.Mult()(1,1) =  1.;
         trans.Sum()(2,0)  =  1.;
         break;
      //lados para interior
      		  //row = dimto = 3 e col = dimfrom = 1
      case  8://sidefrom lado de dim = 1 : t(3,1)
      	trans.Mult()(0,0) =  1.;//          fMult(row,col) : fMult(3,1)
         trans.Sum()(1,0)  = -1.;//          fSum(row,1)    : fSum(3,1)
         trans.Sum()(2,0)  = -1.;//sides : 01 -> 01
         break;
      case  9:
      	trans.Mult()(1,0) =  1.;//sides : 01 -> 12
         trans.Sum()(0,0)  =  1.;//        R1 -> R3
         trans.Sum()(2,0)  = -1.;//        1D -> 3D
         break;
      case 10:
      	trans.Mult()(0,0) = -1.;//sides : 01 -> 23
         trans.Sum()(1,0)  =  1.;
         trans.Sum()(2,0)  = -1.;
         break;
      case 11:
      	trans.Mult()(1,0) = -1.;//sides : 01 -> 30
         trans.Sum()(0,0)  = -1.;
         trans.Sum()(2,0)  = -1.;
         break;
      case 12:
      	trans.Mult()(2,0) =  1.;//sides : 01 -> 04
         trans.Sum()(0,0)  = -1.;
         trans.Sum()(1,0)  = -1.;
         break;
      case 13:
      	trans.Mult()(2,0) =  1.;//sides : 01 -> 15
         trans.Sum()(0,0)  =  1.;
         trans.Sum()(1,0)  = -1.;
         break;
      case 14:
      	trans.Mult()(2,0) =  1.;//sides : 01 -> 26
         trans.Sum()(0,0)  =  1.;
         trans.Sum()(1,0)  =  1.;
         break;
      case 15:
      	trans.Mult()(2,0) =  1.;//sides : 01 -> 37
         trans.Sum()(0,0)  = -1.;
         trans.Sum()(1,0)  =  1.;
         break;
      case 16:
      	trans.Mult()(0,0) =  1.;//sides : 01 -> 45
         trans.Sum()(1,0)  = -1.;
         trans.Sum()(2,0)  =  1.;
         break;
      case 17:
      	trans.Mult()(1,0) =  1.;//sides : 01 -> 56
         trans.Sum()(0,0)  =  1.;
         trans.Sum()(2,0)  =  1.;
         break;
      case 18:
      	trans.Mult()(0,0) = -1.;//sides : 01 -> 67
         trans.Sum()(1,0)  =  1.;
         trans.Sum()(2,0)  =  1.;
         break;
      case 19:
      	trans.Mult()(1,0) = -1.;//sides : 01 -> 70
         trans.Sum()(0,0)  = -1.;
         trans.Sum()(2,0)  =  1.;
         break;
   	//canto para interior
			//row = dimto = 3 e col = dimfrom = 0
      case 0://sidefrom lado de dim = 0 : t(3,0)
      	trans.Sum()(0,0)  = -1.;//          fMult(row,col) : fMult(3,0)
         trans.Sum()(1,0)  = -1.;//          fSum(row,1)    : fSum(3,1)
         trans.Sum()(2,0)  = -1.;//sides : 0 -> 0
         break;
      case 1:
      	trans.Sum()(0,0)  =  1.;//sides : 0 -> 1
         trans.Sum()(1,0)  = -1.;
         trans.Sum()(2,0)  = -1.;
         break;
      case 2:
      	trans.Sum()(0,0)  =  1.;//sides : 0 -> 2
         trans.Sum()(1,0)  =  1.;
         trans.Sum()(2,0)  = -1.;
         break;
      case 3:
      	trans.Sum()(0,0)  = -1.;//sides : 0 -> 3
         trans.Sum()(1,0)  =  1.;
         trans.Sum()(2,0)  = -1.;
         break;
      case 4:
      	trans.Sum()(0,0)  = -1.;//sides : 0 -> 4
         trans.Sum()(1,0)  = -1.;
         trans.Sum()(2,0)  =  1.;
         break;
      case 5:
      	trans.Sum()(0,0)  =  1.;//sides : 0 -> 5
         trans.Sum()(1,0)  = -1.;
         trans.Sum()(2,0)  =  1.;
         break;
      case 6:
      	trans.Sum()(0,0)  =  1.;//sides : 0 -> 6
         trans.Sum()(1,0)  =  1.;
         trans.Sum()(2,0)  =  1.;
         break;
      case 7:
      	trans.Sum()(0,0)  = -1.;//sides : 0 -> 7
         trans.Sum()(1,0)  =  1.;
         trans.Sum()(2,0)  =  1.;
   }
 }//fim if sideto = 26
 if(sideto > 19 && sideto < 26) {
 	switch(sidefrom) {
   	//lados para faces
   		   //row = dimto = 2 e col = dimfrom = 1
      case  8://sidefrom lado de dim =10 : t(2,1)
      	trans.Mult()(0,0) =  1.;//          fMult(row,col) : fMult(2,1)
         trans.Sum()(1,0)  = -1.;//          fSum(row,1)    : fSum(2,1)
         break;                  //sides :  8 -> 20
      case  9:
      	if(sideto==20) {
            trans.Mult()(1,0) =  1.;
            trans.Sum()(0,0)  =  1.;
            break;
         } else
         if(sideto==22) {
            trans.Mult()(0,0) =  1.;
            trans.Sum()(1,0)  = -1.;
            break;
         }
      case 10:
      	if(sideto==20) {//8-9-10-11
            trans.Mult()(0,0) = -1.;//sides : 10 -> 20
            trans.Sum()(1,0)  =  1.;
            break;
         } else
         if(sideto==23) {//10-14-18-15
            trans.Mult()(0,0) = -1.;//sides : 10 -> 20
            trans.Sum()(1,0)  = -1.;
            break;
         }
      case 11:
      	if(sideto==20) {
            trans.Mult()(1,0) = -1.;
            trans.Sum()(0,0)  = -1.;
            break;
         } else
         if(sideto==24) {//11-15-19-12
            trans.Mult()(0,0) = -1.;
            trans.Sum()(1,0)  = -1.;
            break;
         }
      case 12:
         trans.Mult()(1,0) =  1.;
         trans.Sum()(0,0)  = -1.;
         break;
      case 13:
      	if(sideto==21) {//8-13-16-12
            trans.Mult()(1,0) =  1.;
            trans.Sum()(0,0)  =  1.;
            break;
         } else
         if(sideto==22) {//9-14-17-13
            trans.Mult()(1,0) =  1.;
            trans.Sum()(0,0)  = -1.;
            break;
         }
      case 14:
         trans.Mult()(1,0) =  1.;
         trans.Sum()(0,0)  =  1.;
         break;
      case 15:
      	if(sideto==23) {
            trans.Mult()(1,0) =  1.;
            trans.Sum()(0,0)  = -1.;
            break;
         } else
         if(sideto==24) {
            trans.Mult()(1,0) =  1.;
            trans.Sum()(0,0)  =  1.;
            break;
         }
      case 16:
      	if(sideto==21) {
            trans.Mult()(0,0) =  1.;
            trans.Sum()(1,0)  =  1.;
            break;
         } else
         if(sideto==25) {//16-17-18-19
            trans.Mult()(0,0) =  1.;
            trans.Sum()(1,0)  = -1.;
            break;
         }
      case 17:
      	if(sideto==22) {
            trans.Mult()(0,0) =  1.;
            trans.Sum()(1,0)  =  1.;
            break;
         } else
         if(sideto==25) {
            trans.Mult()(1,0) =  1.;
            trans.Sum()(0,0)  =  1.;
            break;
         }
      case 18:
         trans.Mult()(0,0) = -1.;
         trans.Sum()(1,0)  =  1.;
         break;
      case 19:
      	if(sideto==24) {
            trans.Mult()(0,0) = -1.;
            trans.Sum()(1,0)  =  1.;
            break;
         } else
         if(sideto==25) {
            trans.Mult()(1,0) = -1.;
            trans.Sum()(0,0)  = -1.;
            break;
         }
   	//cantos para faces
   		   //row = dimto = 2 e col = dimfrom = 0
      case  0://sidefrom lado de dim =0 : t(2,0)
      	trans.Sum()(0,0) = -1.;//          fMult(row,col) : fMult(2,0)
         trans.Sum()(1,0) = -1.;//          fSum(row,1)    : fSum(2,1)
         break;                 //sides :  0 -> 20,21,24
      case  1:
      	if(sideto==20 || sideto==21) {//sides :  1 -> 20,21
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) = -1.;
            break;
         } else
         if(sideto==22) {//sides :  1 -> 22
            trans.Sum()(0,0) = -1.;
            trans.Sum()(1,0) = -1.;
            break;
         }
      case  2:
      	if(sideto==20) {//sides :  2 -> 20
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) =  1.;
            break;
         } else
         if(sideto==22 || sideto==23) {//sides :  2 -> 22,23
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) = -1.;
            break;
         }
      case  3:
       	if(sideto==20) {//sides :  3 -> 20
            trans.Sum()(0,0) = -1.;
            trans.Sum()(1,0) =  1.;
            break;
         } else
         if(sideto==23) {//sides :  3 -> 23
            trans.Sum()(0,0) = -1.;
            trans.Sum()(1,0) = -1.;
            break;
         } else
         if(sideto==24) {//sides :  3 -> 24
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) = -1.;
            break;
         }
      case  4:
      	if(sideto==21 || sideto==24) {
            trans.Sum()(0,0) = -1.;
            trans.Sum()(1,0) =  1.;
            break;
         } else
         if(sideto==25) {
            trans.Sum()(0,0) = -1.;
            trans.Sum()(1,0) = -1.;
            break;
         }
      case  5:
       	if(sideto==21) {
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) =  1.;
            break;
         } else
         if(sideto==22) {
            trans.Sum()(0,0) = -1.;
            trans.Sum()(1,0) =  1.;
            break;
         } else
         if(sideto==25) {
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) = -1.;
            break;
         }
      case  6:
         trans.Sum()(0,0) =  1.;
         trans.Sum()(1,0) =  1.;
         break;
      case  7:
      	if(sideto==23 || sideto==25) {
            trans.Sum()(0,0) = -1.;
            trans.Sum()(1,0) =  1.;
            break;
         } else
         if(sideto==24) {
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) =  1.;
            break;
         }
      PZError << "TPZGeoElC3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
   }//switch
 }//fim lados
 //cantos para faces
 if(sideto > 7 && sideto < 20) {//8 a 19
 	switch(sidefrom) {
   	//canto para lado
   		   //row = dimto = 1 e col = dimfrom = 0
            //sidefrom canto de dim =0 : t(1,0)
      case  0://fMult(row,col) : fMult(1,0) , fSum(row,1) : fSum(1,1)
         if(sideto==11)               trans.Sum()(0,0)  =  1.;
         if(sideto== 8 || sideto==12) trans.Sum()(0,0)  = -1.;
         //sides : 8 ->  0 , 11 -> 0  , 12 -> 0  ou
         break;//  0 ->  8 , 0  -> 11 , 0  -> 12
      case  1:
      	if(sideto== 8)               trans.Sum()(0,0)  =  1.;
         if(sideto== 9 || sideto==13) trans.Sum()(0,0)  = -1.;
      	break;
      case  2:
      	if(sideto== 9)               trans.Sum()(0,0)  =  1.;
         if(sideto==10 || sideto==14) trans.Sum()(0,0)  = -1.;
      	break;
      case  3:
      	if(sideto==10)               trans.Sum()(0,0)  =  1.;
         if(sideto==11 || sideto==15) trans.Sum()(0,0)  = -1.;
      	break;
      case  4:
         if(sideto==16)               trans.Sum()(0,0)  = -1.;
         if(sideto==12 || sideto==19) trans.Sum()(0,0)  =  1.;
      	break;
      case  5:
      	if(sideto==13 || sideto==16) trans.Sum()(0,0)  =  1.;
         if(sideto==17)               trans.Sum()(0,0)  = -1.;
      	break;
      case  6:
         if(sideto==18)               trans.Sum()(0,0)  = -1.;
      	if(sideto==14 || sideto==17) trans.Sum()(0,0)  =  1.;
      	break;
      case  7:
      	if(sideto==15 || sideto==18) trans.Sum()(0,0)  =  1.;
         if(sideto==19)               trans.Sum()(0,0)  = -1.;
      	break;
  	}//switch
 }//if cantos
   return trans;
}

TPZCompEl *TPZGeoElC3d::CreateBCCompEl(int side,int bc,TPZCompMesh &cmesh) {
	if(side<0 || side>26) return 0;

   if(side==26) {
     std::cout << "TPZGeoElC3d::CreateBCCompEl with side = 26 not implemented\n";
      return 0;
   }
   int64_t index;
	if(side<8) {
      TPZMaterial *bcptr = cmesh.FindMaterial(bc);
      if(!bcptr) {
      PZError << "TPZGeoElC3d::CreateBCCompEl has no bc.\n";
      return 0;
      }
      TPZCompEl *cel = Reference();
      if(!cel) {
      PZError << "TPZGeoElC3d::CreateBCCompEl has no computational element\n";
      return 0;
      }
	  TPZManVector<int64_t> nodeindexes(1);
      TPZGeoElPoint *gel;
      nodeindexes[0] = fNodeIndexes[side];
      gel = new TPZGeoElPoint(nodeindexes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,side));
      return gel->CreateCompEl(cmesh,index);
   } else if (side > 7 && side < 20) {//side = 8 a 19 : lados
      TPZManVector<int64_t> nodes(2);
      int s = side-8;
      nodes[0] = NodeIndex(TPZCompElC3d::SideNodes[s][0]);
      nodes[1] = NodeIndex(TPZCompElC3d::SideNodes[s][1]);
      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,TPZCompElC3d::SideNodes[s][0]));
      TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(this,TPZCompElC3d::SideNodes[s][1]));
      TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(this,side));
      return gel->CreateCompEl(cmesh,index);
   } else if (side > 19) {//side = 20 a 25 : faces
      TPZManVector<int64_t> nodes(4);
      int s = side-20;
      nodes[0] = NodeIndex(TPZCompElC3d::FaceNodes[s][0]);
      nodes[1] = NodeIndex(TPZCompElC3d::FaceNodes[s][1]);
      nodes[2] = NodeIndex(TPZCompElC3d::FaceNodes[s][2]);
      nodes[3] = NodeIndex(TPZCompElC3d::FaceNodes[s][3]);
      TPZGeoElQ2d *gel = new TPZGeoElQ2d(nodes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,TPZCompElC3d::FaceNodes[s][0]));
      TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(this,TPZCompElC3d::FaceNodes[s][1]));
      TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(this,TPZCompElC3d::FaceNodes[s][2]));
      TPZGeoElSide(gel,3).SetConnectivity(TPZGeoElSide(this,TPZCompElC3d::FaceNodes[s][3]));
      TPZGeoElSide(gel,4).SetConnectivity(TPZGeoElSide(this,TPZCompElC3d::FaceSides[s][0]));
      TPZGeoElSide(gel,5).SetConnectivity(TPZGeoElSide(this,TPZCompElC3d::FaceSides[s][1]));
      TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(this,TPZCompElC3d::FaceSides[s][2]));
      TPZGeoElSide(gel,7).SetConnectivity(TPZGeoElSide(this,TPZCompElC3d::FaceSides[s][3]));
      TPZGeoElSide(gel,8).SetConnectivity(TPZGeoElSide(this,side));
      return gel->CreateCompEl(cmesh,index);
   } else PZError << "TPZGeoElC3d::CreateBCCompEl. Side = " << side <<std::endl;
   return 0;
}

void TPZGeoElC3d::NodeFaceIds(TPZVec<int> &ids,int face) {

	ids.Resize(4,-1);
   switch(face) {
   case 20:
         ids[0] = NodeIndex(0);
         ids[1] = NodeIndex(1);
         ids[2] = NodeIndex(2);
         ids[3] = NodeIndex(3);
         return;
   case 21:
         ids[0] = NodeIndex(0);
         ids[1] = NodeIndex(1);
         ids[2] = NodeIndex(5);
         ids[3] = NodeIndex(4);
         return;
   case 22:
         ids[0] = NodeIndex(1);
         ids[1] = NodeIndex(2);
         ids[2] = NodeIndex(6);
         ids[3] = NodeIndex(5);
         return;
   case 23:
         ids[0] = NodeIndex(3);
         ids[1] = NodeIndex(2);
         ids[2] = NodeIndex(6);
         ids[3] = NodeIndex(7);
         return;
   case 24:
         ids[0] = NodeIndex(0);
         ids[1] = NodeIndex(3);
         ids[2] = NodeIndex(7);
         ids[3] = NodeIndex(4);
         return;
   case 25:
         ids[0] = NodeIndex(4);
         ids[1] = NodeIndex(5);
         ids[2] = NodeIndex(6);
         ids[3] = NodeIndex(7);
         return;
   default :
	   std::cout << "TPZCompElC3d::NodeFaceIds bad side , side = " << face <<std::endl;
   }
}

static int fatherside[8][27] = {
/*00*/{0,8,20,11,12,21,26,24,8,20,20,11,12,21,26,24,21,26,26,24,20,21,26,26,24,26,26},
/*01*/{8,1,9,20,21,13,22,26,8,9,20,20,21,13,22,26,21,22,26,26,20,21,22,26,26,26,26},
/*02*/{20,9,2,10,26,22,14,23,20,9,10,20,26,22,14,23,26,22,23,26,20,26,22,23,26,26,26},
/*03*/{11,20,10,3,24,26,23,15,20,20,10,11,24,26,23,15,26,26,23,24,20,26,26,23,24,26,26},
/*04*/{12,21,26,24,4,16,25,19,21,26,26,24,12,21,26,24,16,25,25,19,26,21,26,26,24,25,26},
/*05*/{21,13,22,26,16,5,17,25,21,22,26,26,21,13,22,26,16,17,25,25,26,21,22,26,26,25,26},
/*06*/{26,22,14,23,25,17,6,18,26,22,23,26,26,22,14,23,25,17,18,25,26,26,22,23,26,25,26},
/*07*/{24,26,23,15,19,25,18,7,26,26,23,24,24,26,23,15,25,25,18,19,26,26,26,23,24,25,26},
};

TPZGeoElSide TPZGeoElC3d::Father2(int side){//Augusto:09/01/01

	if(side<0 || side>26){
		PZError << "TPZGeoElC3d::Father2 called error" <<std::endl;
        return TPZGeoElSide();
	}
	int subelindex = WhichSubel();
	if(fatherside[subelindex][side]<0){
		PZError << "TPZGeoElC3d::Father2 called with index error\n";
		return TPZGeoElSide();
	}
	return Father(side);
}
  
static int subeldata[27][26][2] = {//Augusto:09/01/01
/*00*/{{0,0}},
/*01*/{{1,1}},
/*02*/{{2,2}},
/*03*/{{3,3}},
/*04*/{{4,4}},
/*05*/{{5,5}},
/*06*/{{6,6}},
/*07*/{{7,7}},
/*08*/{{0,8},{0,1},{1,8}},
/*09*/{{1,9},{1,2},{2,9}},
/*10*/{{2,10},{2,3},{3,10}},
/*11*/{{0,11},{0,3},{3,11}},
/*12*/{{0,12},{0,4},{4,12}},
/*13*/{{1,13},{1,5},{5,13}},
/*14*/{{2,14},{2,6},{6,14}},
/*15*/{{3,15},{3,7},{7,15}},
/*16*/{{4,16},{4,5},{5,16}},
/*17*/{{5,17},{5,6},{6,17}},
/*18*/{{6,18},{6,7},{7,18}},
/*19*/{{4,19},{4,7},{7,19}},
/*20*/{{0,2},{1,20},{2,20},{3,20},{0,9},{0,10},{2,8},{2,11}},
/*21*/{{0,21},{1,21},{4,21},{5,21},{0,13},{0,16},{5,8},{5,12}},
/*22*/{{1,22},{2,22},{5,22},{6,22},{1,14},{1,17},{6,9},{6,13}},
/*23*/{{2,23},{3,23},{6,23},{7,23},{2,15},{2,18},{7,10},{7,14}},
/*24*/{{0,24},{3,24},{4,24},{7,24},{0,15},{0,19},{7,11},{7,12}},
/*25*/{{4,25},{5,25},{6,25},{7,25},{4,17},{4,18},{6,16},{6,19}},
/*26*/{{0,26},{1,26},{2,26},{3,26},{4,26},{5,26},{6,26},{7,26},{0,22},
       {0,23},{0,25},{1,25},{2,21},{2,24},{2,25},{3,25},{4,22},{4,23},
       {6,21},{6,24},{0,14},{0,17},{0,18},{6,8},{6,11},{6,12}}
};

static int nsubeldata[27] = {1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,3,8,8,8,8,8,8,26};


void TPZGeoElC3d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){

	  subel.Resize(0);
    if(side<0 || side>26 || !HasSubElement(side)){
       PZError << "TPZGeoElC3d::GetSubelements2 called with error arguments\n";
       return;
    }
    int nsub = nsubeldata[side];
    for(int i=0;i<nsub;i++)
    		subel.Push(TPZGeoElSide(fSubEl[subeldata[side][i][0]],
                                       subeldata[side][i][1]));
}

static REAL buildt[8][27][4][3] = {//por colunas
/*S0*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{-1.,-1.,-1.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*07*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*09*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*10*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*13*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*14*/{{0.,0.,0.5},{0.,0.,0.},{0.,0.,0.},{0.,0.,-0.5}},
      /*15*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*16*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*17*/{{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-0.5,0.}},
      /*18*/{{-0.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-0.5,0.,0.}},
      /*19*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*21*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*22*/{{0.,0.5,0.},{0.,0.,0.5},{0.,0.,0.},{0.,-0.5,-0.5}},
      /*23*/{{0.5,0.,0.},{0.,0.,0.5},{0.,0.,0.},{-0.5,0.,-0.5}},
      /*24*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*25*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{-0.5,-0.5,0.}},
      /*26*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{-.5,-.5,-.5}}},
/*S1*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1,-1,-1}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*07*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*10*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*11*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
  		/*12*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*13*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*14*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*15*/{{0.,0.,0.5},{0.,0.,0.},{0.,0.,0.},{0.,0.,-0.5}},
      /*16*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
	  	/*17*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*18*/{{-0.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*19*/{{0.,-0.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-0.5,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*21*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*22*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*23*/{{0.5,0.,0.},{0.,0.,0.5},{0.,0.,0.},{0.5,0.,-0.5}},
      /*24*/{{0.,0.5,0.},{0.,0.,0.5},{0.,0.,0.},{0.,-0.5,-0.5}},
      /*25*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,-0.5,0.}},
      /*26*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,-.5,-.5}}},
/*S2*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1,1,-1}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*07*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*11*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*12*/{{0.,0.,0.5},{0.,0.,0.},{0.,0.,0.},{0.,0.,-0.5}},
      /*13*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*14*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*15*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*16*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*17*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*18*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*19*/{{0.,-0.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*21*/{{0.5,0.,0.},{0.,0.,0.5},{0.,0.,0.},{0.5,0.,-0.5}},
      /*22*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*23*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*24*/{{0.,0.5,0.},{0.,0.,0.5},{0.,0.,0.},{0.,0.5,-0.5}},
      /*25*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.5,0.}},
      /*26*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,.5,-.5}}},
/*S3*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{-1,1,-1}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*07*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*09*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*12*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*13*/{{0.,0.,.5},{0.,0.,0.},{0.,0.,0.},{0.,0.,-.5}},
      /*14*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*15*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*16*/{{0.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-0.5,0.,0.}},
      /*17*/{{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*18*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*19*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*21*/{{0.5,0.,0.},{0.,0.,0.5},{0.,0.,0.},{-0.5,0.,-0.5}},
      /*22*/{{0.,0.5,0.},{0.,0.,0.5},{0.,0.,0.},{0.,0.5,-0.5}},
      /*23*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*24*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*25*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{-0.5,0.5,0.}},
      /*26*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{-.5,.5,-.5}}},
/*S4*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{-1,-1,1}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*07*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*09*/{{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-0.5,0.}},
      /*10*/{{-0.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-0.5,0.,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*13*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*14*/{{0.,0.,0.5},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.5}},
      /*15*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*16*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*17*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*18*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*19*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*20*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{-0.5,-0.5,0.}},
      /*21*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*22*/{{0.,0.5,0.},{0.,0.,0.5},{0.,0.,0.},{0.,-0.5,0.5}},
      /*23*/{{0.5,0.,0.},{0.,0.,0.5},{0.,0.,0.},{-0.5,0.,0.5}},
      /*24*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*25*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*26*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{-.5,-.5,.5}}},
/*S5*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1,-1,1}},
      /*06*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*07*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*10*/{{-0.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*11*/{{0.,-0.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-0.5,0.}},
      /*12*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*13*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*14*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*15*/{{0.,0.,0.5},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.5}},
      /*16*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*17*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*18*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*19*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*20*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,-0.5,0.}},
      /*21*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*22*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*23*/{{0.5,0.,0.},{0.,0.,0.5},{0.,0.,0.},{0.5,0.,0.5}},
      /*24*/{{0.,0.5,0.},{0.,0.,0.5},{0.,0.,0.},{0.,-0.5,0.5}},
      /*25*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*26*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,-.5,.5}}},
/*S6*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1,1,1}},
      /*07*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*08*/{{0.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*11*/{{0.,-0.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*12*/{{0.,0.,0.5},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.5}},
      /*13*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*14*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*15*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*16*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*17*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*18*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*19*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*20*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.5,0.}},
      /*21*/{{0.5,0.,0.},{0.,0.,0.5},{0.,0.,0.},{0.5,0.,0.5}},
      /*22*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*23*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*24*/{{0.,0.5,0.},{0.,0.,0.5},{0.,0.,0.},{0.,0.5,0.5}},
      /*25*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*26*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,.5,.5}}},
/*S7*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*07*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{-1,1,1}},
      /*08*/{{0.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-0.5,0.,0.}},
      /*09*/{{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*12*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*13*/{{0.,0.,0.5},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.5}},
      /*14*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*15*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*16*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*17*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*18*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*19*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*20*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{-0.5,0.5,0.}},
      /*21*/{{0.5,-0.5,0.5},{0.,0.,0.5},{0.,0.,0.},{-0.5,0.,0.5}},
      /*22*/{{0.,0.5,0.},{0.,0.,0.5},{0.,0.,0.},{0.,0.5,0.5}},
      /*23*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*24*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*25*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*26*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{-.5,.5,.5}}}
};

TPZTransform<STATE> TPZGeoElC3d::BuildTransform2(int side, TPZGeoEl * /*father*/){//Augusto:09/01/01

	if(side<0 || side>26 || !Father(side).Element()){
  	PZError << "TPZGeoElC3d::BuildTransform2 side out of range or father null\n";
    return TPZTransform<STATE>(0,0);
  }
  TPZTransform<STATE> trans(3,3);
  int son = WhichSubel();
  trans.Mult()(0,0) = buildt[son][side][0][0];
  trans.Mult()(1,0) = buildt[son][side][0][1];
  trans.Mult()(2,0) = buildt[son][side][0][2];
  trans.Mult()(0,1) = buildt[son][side][1][0];
  trans.Mult()(1,1) = buildt[son][side][1][1];
  trans.Mult()(2,1) = buildt[son][side][1][2];
  trans.Mult()(0,2) = buildt[son][side][2][0];
  trans.Mult()(1,2) = buildt[son][side][2][1];
  trans.Mult()(2,2) = buildt[son][side][2][2];
  trans.Sum() (0,0) = buildt[son][side][3][0];
  trans.Sum() (1,0) = buildt[son][side][3][1];

  trans.Sum() (2,0) = buildt[son][side][3][2];

  return trans;
}

static REAL MidSideNode[27][3] = {
/*00*/{-1.,-1.,-1.},/*01*/{1.,-1.,-1.},/*02*/{1.,1.,-1.},/*03*/{-1.,1.,-1.},
/*04*/{-1.,-1., 1.},/*05*/{1.,-1., 1.},/*06*/{1.,1., 1.},/*07*/{-1.,1., 1.},
/*08*/{ 0.,-1.,-1.},/*09*/{1., 0.,-1.},/*10*/{0.,1.,-1.},/*11*/{-1.,0.,-1.},
/*12*/{-1.,-1., 0.},/*13*/{1.,-1., 0.},/*14*/{1.,1., 0.},/*15*/{-1.,1., 0.},
/*16*/{ 0.,-1., 1.},/*17*/{1., 0., 1.},/*18*/{0.,1., 1.},/*19*/{-1.,0., 1.},
/*20*/{ 0., 0.,-1.},/*21*/{0.,-1., 0.},/*22*/{1.,0., 0.},/*23*/{ 0.,1., 0.},
/*24*/{-1., 0., 0.},/*25*/{0., 0., 1.},/*26*/{0.,0., 0.} };

int TPZGeoElC3d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3),pss(3),pf(3),pfs(3);
                    //point son, point side son, point father, point side father : elemento mestre
  for(sn=0;sn<8;sn++){
    TPZGeoEl *son = subs[sn];
    for(sd=0;sd<27;sd++){
      ps[0] = MidSideNode[sd][0];//element
      ps[1] = MidSideNode[sd][1];//master point
      ps[2] = MidSideNode[sd][2];
      if(son->WhichSide(ps) != sd)std::cout << "Lado nao bate\n";
      TPZTransform<STATE> telsd = pzshape::TPZShapeCube::TransformElementToSide(sd);//2x2
      telsd.Apply(ps,pss);//son element -> son side 
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform<STATE> t;
	  son->BuildTransform2(sd, gel,t);
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = pzshape::TPZShapeCube::TransformSideToElement(sdfat);//2x2
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(26).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao errada\n";
        PZError << "son    = " << (son->Id()) <<std::endl;
        PZError << "father = " << ((son->Father2(26).Element())->Id()) <<std::endl;
        PZError << "side   = " << sd <<std::endl <<std::endl;
        int ok;
		std::cin >> ok;
      } else {
		  std::cout << "Transformacao OK!\n";
		  std::cout << "Filho/lado : " << son->Id() << "/" << sd <<std::endl;
       std::cout << "Pai : " << son->Father2(26).Element()->Id() <<std::endl <<std::endl;
      }
    }
  }
  return 1;
}


//METHODS DEFINITION FOR CLASS ELEMGTD
#include "pzelgt3d.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelc1d.h"
#include "pzelgt2d.h"
#include "pzelgpi3d.h"
#include "pzelct3d.h"
#include "pzelcpi3d.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgnode.h"
#include "pzshtmat.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include <stdlib.h>

static TPZCompEl *CreateEl(TPZGeoElT3d *gel,TPZCompMesh &mesh, int64_t &index) {
  return new TPZCompElT3d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElT3d::fp)(TPZGeoElT3d *,TPZCompMesh &, int64_t &) = CreateEl;

TPZGeoElT3d::TPZGeoElT3d(int id,TPZVec<int64_t> &nodeindexes,int matid,TPZGeoMesh &mesh):
	TPZGeoElRefLess<pzgeom::TPZGeoTetrahedra>(id,nodeindexes,matid,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=4) {
    PZError << "TPZGeoElT3d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<4;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<6;i++) fSubEl[i] = 0;
}

TPZGeoElT3d::TPZGeoElT3d(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) :
	TPZGeoElRefLess<pzgeom::TPZGeoTetrahedra>(nodeindexes,matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=4) {
    PZError << "TPZGeoElT3d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<4;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<6;i++) fSubEl[i] = 0;
}

TPZGeoElT3d::~TPZGeoElT3d() {}

void TPZGeoElT3d::Shape(TPZVec<REAL> &pt,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {
  phi(0,0)  = 1-pt[0]-pt[1]-pt[2];
  phi(1,0)  = pt[0];
  phi(2,0)  = pt[1];
  phi(3,0)  = pt[2];

  dphi(0,0) = -1.0;
  dphi(1,0) = -1.0;
  dphi(2,0) = -1.0;
  dphi(0,1) =  1.0;
  dphi(1,1) =  0.0;
  dphi(2,1) =  0.0;
  dphi(0,2) =  0.0;
  dphi(1,2) =  1.0;
  dphi(2,2) =  0.0;
  dphi(0,3) =  0.0;
  dphi(1,3) =  0.0;
  dphi(2,3) =  1.0;
}

TPZGeoElT3d *TPZGeoElT3d::CreateGeoEl(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElT3d(nodeindexes,matind,mesh);
}

int TPZGeoElT3d::NNodes() {
  return 4;
}
int64_t TPZGeoElT3d::NodeIndex(int node) {
  if(node<0 || node>3) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElT3d::NSideNodes(int side) {
  if(side<0 || side>14) {
    PZError << "TPZGeoElT3d::NSideNodes. Bad parameter side.\n";
    return 0;
  }
  if(side<4) return 1; //cantos
  if(side>3 && side<10) return 2;//lados
  if(side<14) return 3;//faces
  return 4;//centro
}

int64_t TPZGeoElT3d::SideNodeIndex(int side,int node) {
  if(side<0 || side>14 || node<0) {//15 sides
    PZError << "TPZGeoElT3d::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  //0<=side<=14 e node>=0
  //side central
  if(side==14 && node<4) return fNodeIndexes[node];//um dos 4 cantos
  //side de canto
  if(side < 4) return fNodeIndexes[side];//canto side, deve ser node = 0
  //sides 4 a 9
  if(side>3 && side<10 && node<2) {//lados 4 a 9
      side-=4;
      return fNodeIndexes[TPZCompElT3d::SideNodes[side][node]];
  } else if(side>9 && node<3) {//faces 10 a 13
  		side-=10;
	  	return fNodeIndexes[TPZCompElT3d::FaceNodes[side][node]];
  }
  return -1;
}

void TPZGeoElT3d::MidSideNodeIndex(int side, int64_t &index) {
  index = -1;
  if(side<0 || side>14) {
    PZError << "TPZGeoElT3d::MidSideNodeIndex. Bad parameter side = " << side <<std::endl;
    return;
  }
  //sides 0 a 3
	if(side<4) {//o nó medio do lado 0 é o 0 etc.
		index=fNodeIndexes[side];
		return;
	}
   //o nó medio da face é o centro da face e o nó medio do centro é o centro
   //como nó de algum filho se este existir
   //caso tenha filhos é o canto de algum filho, se não tiver filhos retorna -1
	if(HasSubElement(0)) {
   	side-=4;
	   index=((TPZGeoElT3d *) fSubEl[TPZCompElT3d::MidSideNodes[side][0]])->fNodeIndexes[TPZCompElT3d::MidSideNodes[side][1]];
	}
}

void TPZGeoElT3d::NewMidSideNode(int side, int64_t &index) {
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
    if(side < 4) {index = -1; return;}
    //aqui side = 4 a 14
    side-=4;//0 a 10
    par[0] = TPZCompElT3d::MidCoord[side][0];
    par[1] = TPZCompElT3d::MidCoord[side][1];
    par[2] = TPZCompElT3d::MidCoord[side][2];
    X(par,coord);
    index = Mesh()->NodeVec().AllocateNewElement();
    Mesh()->NodeVec()[index].Initialize(coord,*Mesh());
  }
}

int TPZGeoElT3d::SideDimension(int side) {

	if (side<0 || side>14) {
   	PZError << "TPZGeoElT3d::SideDimension called with side " << side <<std::endl;
      return 0;
   }
  if(side<4) return 0;//cantos
  if(side>3 && side<10) return 1;//lados
  if(side<14) return 2;//faces
  return 3;//centro

   }

TPZGeoElSide TPZGeoElT3d::HigherDimensionSides(int side,int targetdimension) {
//targetdimension deve ser 1 , 2 ou 3
//se side =  0 a 3  targetdimension deve ser 1
//se side =  4 a 9  targetdimension deve ser 2
//se side = 14      targetdimension deve ser 3
  if( (side<0 || side>14) || (targetdimension < 1 || targetdimension > 3) ) {
     PZError << "TPZGeoElT3d::HigherDimensionSides called with side = " << side
	          << " targetdimension = " << targetdimension <<std::endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
  int bestface;
  //side = 0 a 14
  switch(targetdimension) {//=1,2
	  case 1:
       if(father->NSides() == 19) {//o pai é uma pirâmide
          if(this == father->SubElement(6)) {
          	 if(side==0) return TPZGeoElSide(father->SubElement(4), 9);
             if(side==1) return TPZGeoElSide(father->SubElement(0), 5);
             if(side==2) return TPZGeoElSide(this,5);//canto do mesmo id
          } else if(this == father->SubElement(7)) {
             if(side==1) return TPZGeoElSide(father->SubElement(4),10);
             if(side==2) return TPZGeoElSide(father->SubElement(4),11);
             if(side==3) return TPZGeoElSide(this,7);

          } else if(this == father->SubElement(8)) {
             if(side==1) return TPZGeoElSide(this,5);
             if(side==2) return TPZGeoElSide(father->SubElement(3), 7);

          } else if(this == father->SubElement(9)) {
          	 if(side==0) return TPZGeoElSide(this,7);
             if(side==3) return TPZGeoElSide(father->SubElement(0),8);
          }
       	return TPZGeoElSide();
       }
     	 if(this == father->SubElement(0)) {
       	 if(side==1) return TPZGeoElSide(this,4);
       	 if(side==2) return TPZGeoElSide(this,6);
          if(side==3) return TPZGeoElSide(this,7);
       } else if(this == father->SubElement(1)) {
       	 if(side==0) return TPZGeoElSide(this,4);
       	 if(side==2) return TPZGeoElSide(this,5);
          if(side==3) return TPZGeoElSide(this,8);
       } else if(this == father->SubElement(2)) {
       	 if(side==0) return TPZGeoElSide(this,6);
       	 if(side==1) return TPZGeoElSide(this,5);
          if(side==3) return TPZGeoElSide(this,9);
       } else if(this == father->SubElement(3)) {
       	 if(side==0) return TPZGeoElSide(this,7);
       	 if(side==1) return TPZGeoElSide(this,8);
          if(side==2) return TPZGeoElSide(this,9);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
  	  case 2:
     	 if(this == father->SubElement(0)) {

          bestface = BestDimensionSideOfTwoFaces(10,11);
          if((side==1 || side==4) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(10,13);
          if((side==2 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(11,13);
          if((side==3 || side==7) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 5) return TPZGeoElSide(this,10);
       	 if(side== 8) return TPZGeoElSide(this,11);
       	 if(side== 9) return TPZGeoElSide(this,13);
       } else if(this == father->SubElement(1)) {

          bestface = BestDimensionSideOfTwoFaces(10,11);
          if((side==0 || side==4) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(10,12);
          if((side==2 || side==5) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(11,12);
          if((side==3 || side==8) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 6) return TPZGeoElSide(this,10);
       	 if(side== 7) return TPZGeoElSide(this,11);
       	 if(side== 9) return TPZGeoElSide(this,12);
       } else if(this == father->SubElement(2)) {

          bestface = BestDimensionSideOfTwoFaces(10,12);
          if((side==1 || side==5) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(10,13);
          if((side==0 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(12,13);
          if((side==2 || side==9) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 4) return TPZGeoElSide(this,10);
       	 if(side== 7) return TPZGeoElSide(this,13);
       	 if(side== 8) return TPZGeoElSide(this,12);
       } else if(this == father->SubElement(3)) {

          bestface = BestDimensionSideOfTwoFaces(11,13);
          if((side==0 || side==7) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(11,12);
          if((side==1 || side==8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(12,13);
          if((side==2 || side==9) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side== 4) return TPZGeoElSide(this,11);
       	 if(side== 5) return TPZGeoElSide(this,12);
       	 if(side== 6) return TPZGeoElSide(this,13);
       }
       //o pai é uma pirâmide
       if(father->NSides() == 19) {
          if(this == father->SubElement(6)) {

             TPZGeoElPi3d *sub = (TPZGeoElPi3d *) father->SubElement(4);//pois primeior subira pelo canto para aresta e logo de aresta para face
             bestface = sub->BestDimensionSideOfTwoFaces(14,17);
             if(side==0 && bestface) return TPZGeoElSide(sub,bestface);
             sub = (TPZGeoElPi3d *) father->SubElement(0);
             bestface = sub->BestDimensionSideOfTwoFaces(13,14);
             if(side==1 && bestface) return TPZGeoElSide(sub,bestface);//para side=5 a aresta do sub 0 nao é 5, entao nao pode, outro irmao o considerara
             if(side==2) return TPZGeoElSide(sub,13);//para side 5 as faces do atual que contem este side
             if(side==4 || side==7 || side==8) return TPZGeoElSide(this,11);

          } else if(this == father->SubElement(7)) {
             TPZGeoElPi3d *sub = (TPZGeoElPi3d *) father->SubElement(4);//1
             bestface = sub->BestDimensionSideOfTwoFaces(14,15);
             if(side==1 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(15,16);
             if(side==2 && bestface) return TPZGeoElSide(sub,bestface);
             sub = (TPZGeoElPi3d *) father->SubElement(1);
             if(side==3 || side==7) return TPZGeoElSide(sub,13);
             if(side==4 || side==5 || side==6) return TPZGeoElSide(this,10);

          } else if(this == father->SubElement(8)) {

             TPZGeoElPi3d *sub = (TPZGeoElPi3d *) father->SubElement(3);//2
             if(side==1) return TPZGeoElSide(sub,13);
             bestface = sub->BestDimensionSideOfTwoFaces(13,16);
             if(side==2 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==6 || side==7 || side==9) return TPZGeoElSide(this,13);
          } else if(this == father->SubElement(9)) {

             TPZGeoElPi3d *sub = (TPZGeoElPi3d *) father->SubElement(0);
             bestface = sub->BestDimensionSideOfTwoFaces(13,17);
             if(side==3 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==7) return TPZGeoElSide(sub,13);
             if(side==5 || side==8 || side==9) return TPZGeoElSide(this,12);
          }
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
     case 3:
       return TPZGeoElSide(this,14);//0<=side<=14
  }//switch
  return TPZGeoElSide();
}

int TPZGeoElT3d::BestDimensionSideOfTwoFaces(int face1,int face2) {

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


void TPZGeoElT3d::LowerDimensionSides(int side,TPZStack<TPZGeoElSide> &smallsides) {
	if (side < 4) return;
   int i;
   if(side < 10) {//side = 4 a 9 : entram cantos dos lados
   	int s = side-4;
   	smallsides.Push(TPZGeoElSide(this,TPZCompElT3d::SideNodes[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElT3d::SideNodes[s][1]));
   } else if(side < 14) {//entram cantos e lados da face
   	int s = side-10;
   	smallsides.Push(TPZGeoElSide(this,TPZCompElT3d::FaceNodes[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElT3d::FaceNodes[s][1]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElT3d::FaceNodes[s][2]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElT3d::FaceSides[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElT3d::FaceSides[s][1]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElT3d::FaceSides[s][2]));
      smallsides.Push(TPZGeoElSide(this,side));
   } else if(side==14) {
   	for (i=0;i<14;i++) smallsides.Push(TPZGeoElSide(this,i));
   }
}

void TPZGeoElT3d::Jacobian(TPZVec<REAL> &param,TPZFMatrix<STATE> &jacobian,TPZFMatrix<STATE> &axes,REAL &detjac,TPZFMatrix<STATE> &jacinv){

  //int nnodes = NNodes();
#ifdef DEBUG
  if (nnodes != 4) {
    PZError << "TPZGeoElT3d.jacobian only implemented for"
      " 4 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 3 || param[0] < 0. || param[0] > 1. ||
     param[1] < 0. || param[1] > 1. || param[2] < 0. || param[2] > 1.) {
    PZError << "TPZGeoElT3d.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
    return;
  }
#endif
  REAL spacephi[10];
  TPZFMatrix<STATE> phi(4,1,spacephi,10);
  REAL spacedphi[20];
  TPZFMatrix<STATE> dphi(3,4,spacedphi,20);
  Shape(param,phi,dphi);
  jacobian.Zero();
  TPZGeoNode *np;
  int i,j;
  for(i=0;i<4;i++) {
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

void TPZGeoElT3d::X(TPZVec<REAL> & loc,TPZVec<REAL> &result){
  REAL spacephi[10],spacedphi[20];
  int i,j;
  TPZFMatrix<STATE> phi(4,1,spacephi,10);
  TPZFMatrix<STATE> dphi(3,4,spacedphi,20);
  Shape(loc,phi,dphi);
  for(j=0;j<3;j++) {
    result[j] = 0.0;
    for(i=0;i<4;i++) result[j] += NodePtr(i)->Coord(j)*phi(i,0);
  }
}
/**It's necessary to define the normal vector to side 4, that is the orthogonal
   vector to the surface*/

void TPZGeoElT3d::NormalVector(int /*side*/,TPZVec<REAL> &/*param*/,TPZVec<REAL> &/*normal*/,
			     TPZFMatrix<STATE> &/*axes*/,TPZFMatrix<STATE> &/*jac1d*/) {
/*
#ifdef DEBUG
  if (nnodes != 8) {
    PZError << "TPZGeoElT3d.NormalVector, only implemented for"
      " 8 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 3 || param[0] < -1. || param[0] > 1. ||
     param[1] < -1. || param[1] > 1. || param[2] < -1. || param[2] > 1.) {
    PZError << "TPZGeoElT3d.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
    return;
  }
  if(normal.NElements() != 3) {
    PZError << "TPZGeoElT3d::NormalVector normal.capacity() = " << normal.NElements() <<
      "\n";
    return;
  }
  if(side < 0 || side > 5) {//6 faces
    PZError << "TPZGeoElT3d.jacobian invalid side : "
      " side = " << side << "\n";
    return;
  }
#endif

  REAL spacephi[12],spacedphi[30];
  TPZFMatrix phi(8,1,spacephi,12);
  TPZFMatrix dphi(3,8,spacedphi,30);
  Shape(param,phi,dphi);
  TPZGeoNode *np;
  TPZVec<REAL> n(3,0.);
  int i,j,ider;
  if(side==0 || side==5) ider = 2;
  if(side==1 || side==3) ider = 1;
  if(side==2 || side==4) ider = 0;
  for(i=0;i<4;i++) {
    np = NodePtr(TPZCompElT3d::FaceNodes[side][i]);
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
  return;  */
}

/** TO SUBDIVISION
********************************************************************************
  Into Divides is necesary to consider the connectivity with the all neighboards*/
void TPZGeoElT3d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {

   int i;
   if(HasSubElement(0)) {
      SubElVec.Resize(6);
      for(i=0;i<6;i++) SubElVec[i] = fSubEl[i];
      return;//If exist fSubEl return this sons
   }
   int j, sub, matid = MaterialId();
   int64_t index;
   int64_t np[10];//guarda conectividades dos 6 subelementos

   for(j=0;j<4;j++) np[j]=NodeIndex(j);
   for(sub=4;sub<10;sub++) {
      NewMidSideNode(sub,index);
      np[sub] = index;
   }
   // creating new subelements
   for(i=0;i<4;i++) {
	   TPZManVector<int64_t> cornerindexes(4);
   	for(int j=0;j<4;j++) cornerindexes[j] = np[TPZCompElT3d::CornerSons[i][j]];
      fSubEl[i] = CreateGeoEl(cornerindexes,matid,*Mesh());
   }
   for(;i<6;i++) {
      TPZManVector<int64_t> cornerindexes(5);
      for(int j=0;j<5;j++) cornerindexes[j] = np[TPZCompElT3d::CornerSons[i][j]];
      TPZGeoElPi3d *subpi=0;
      fSubEl[i] = subpi->CreateGeoEl(cornerindexes,matid,*Mesh());
   }
   if(SubElVec.NElements()!=6) SubElVec.Resize(6);
   for(sub=0;sub<6;sub++) {
      SubElVec[sub] = fSubEl[sub];
      SubElVec[sub]->SetFather(this);
   }
   for(i=0;i<6;i++) {//conectividades entre os filhos : viz interna
   	for(j=0;j<16;j++) {
      	int elside = TPZCompElT3d::InNeigh[i][j][0];//lado do subel
         if(elside == -1) break;
         int k = TPZCompElT3d::InNeigh[i][j][1];//número do filho viz
         int l = TPZCompElT3d::InNeigh[i][j][2];//lado do viz.
         TPZGeoElSide neighside(fSubEl[k],l);
      	fSubEl[i]->SetNeighbour(elside,neighside);
      }
   }
   //vizinhança externa ao elemento atual
   //procura-se um viz pela face que seja dividido
   TPZGeoElSide dividedneighbour;
   for(int face=10;face<14;face++) {
      dividedneighbour = Neighbour(face);//BuildConnectivity(..) inicializou
      while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         //achou-se um viz para alguma face
         TPZManVector<int64_t> nodes(3);
         TPZStack<TPZGeoElSide> neighsub;
         int f = face-10;
         for(i=0;i<3;i++) nodes[i] = NodeIndex(TPZCompElT3d::FaceNodes[f][i]);//nós globais na ordem e sentido local da face atual
         dividedneighbour.GetSubElements2(neighsub);//4 subs viz na ordem dos subs da face atual
         nodes.Resize(1);
         //conectividades dos lados dos subelementos interiores a face
         nodes.Resize(2);
         for(i=0; i<3; i++) {//somente 3 arestas interiores à face
            int fsi = TPZCompElT3d::FaceSons[f][i];
            int inside = TPZCompElT3d::FaceSides[f][(i+1)%3];//
            nodes[0] = fSubEl[fsi]->SideNodeIndex(inside,0);
            nodes[1] = fSubEl[fsi]->SideNodeIndex(inside,1);
            int locside = neighsub[i].Element()->WhichSide(nodes);
            TPZGeoElSide sub(fSubEl[fsi],inside);
            TPZGeoElSide subneigh(neighsub[i].Element(),locside);
            sub.SetConnectivity(subneigh);
         }
         //side de face para os 4 subs da face f
         for(i=0; i<4; i++) {
            if(i==3) {
            	TPZGeoElSide sub = SideSubElement(face,i);
               sub.SetConnectivity(neighsub[i]);
               continue;
            }
            TPZGeoElSide sub(fSubEl[TPZCompElT3d::FaceSons[f][i]],face);
            sub.SetConnectivity(neighsub[i]);
         }
      }
   }//fim faces
	//procura-se um viz pelo lado do atual, que seja dividido
   for(int side=4; side<10; side++) {//(side 8 : subs 0,1),(side 9 : subs 1,2) ...
   	int s = side-4;
   	dividedneighbour = Neighbour(side);
      while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         //achou-se um viz pelo lado i e dividido
         int i = TPZCompElT3d::SideNodes[s][0];//nó local 0 do lado side
         int iplus = TPZCompElT3d::SideNodes[s][1];//nó local 1 do lado
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
   for(int corner=0; corner<4; corner++) {//cantos
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


int TPZGeoElT3d::NSubElements() {
  return 6;
}

int TPZGeoElT3d::NSideSubElements(int side) {
  if(side < 0 || side > 14) {
    PZError << "TPZGeoElT3d::NSideSubElements called for side " << side <<std::endl;
    return 0;
  }
  if(side==14) return 6;//centro
  if(side>9 && side<14) return 4;//faces
  if(side>3) return 2;//lados
  return 1;//cantos
}

TPZGeoElSide TPZGeoElT3d::SideSubElement(int side,int position) {
   if (position<0 || position>6 || side <0 ||side>14) {
   	PZError << "TPZGeoElT3d::SideSubElement called with position " << position << " side " << side <<std::endl;
      return TPZGeoElSide();
   }                              //fSubEl[is]
   if(side==14) {
      if(position > 3) return TPZGeoElSide(SubElement(position),18);//centro : pirâmides
      return TPZGeoElSide(SubElement(position),14);//centro : tetraedros
	}
   if(side<4) {//cantos
      if(position!=0) {
         PZError << "TPZGeoElT3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
         return TPZGeoElSide(SubElement(side),side);
      }
   }
   if(side>3 && side<10) {//lados
       if(position!=0 && position!=1) {
         PZError << "TPZGeoElT3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
      	int s = side-4;
         return TPZGeoElSide(SubElement(TPZCompElT3d::SideNodes[s][position]),side);
      }
   }
   if(side>9) {//faces
       if(position<0 || position>3) {//position!=0 && position!=1 && position!=2 && position!=3
         PZError << "TPZGeoElT3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
      	int s = side-10;//s = 0,1,2,3
         if(s==0 && position == 3) return TPZGeoElSide(SubElement(TPZCompElT3d::FaceSons[s][3]),15);//lado do subelemento
         if(s==1 && position == 3) return TPZGeoElSide(SubElement(TPZCompElT3d::FaceSons[s][3]),14);//interior à face
         if(s==2 && position == 3) return TPZGeoElSide(SubElement(TPZCompElT3d::FaceSons[s][3]),17);//nestes casos uma
         if(s==3 && position == 3) return TPZGeoElSide(SubElement(TPZCompElT3d::FaceSons[s][3]),16);//piramide
         return TPZGeoElSide(SubElement(TPZCompElT3d::FaceSons[s][position]),side);
      }
   }
   return TPZGeoElSide();
}

void TPZGeoElT3d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
   if(!fSubEl[0]) {
      sub.Resize(0);
      return;
   }
   if(side < 0 || side > 14) {
      PZError << "TPZGeoElT3d::SideSubElements called for side " << side <<std::endl;
      return;
   }
   if(side==14) {
      sub.Resize(6);
      for(int i=0;i<6;i++) sub[i] = fSubEl[i];
      return;
   }
   if(side<4) {
      sub.Resize(1);
      sub[0]=fSubEl[side];
      return;
   }
   if(side>3 && side<10) {//lados
      int s = side-4;
      sub.Resize(2);
      sub[0] = fSubEl[TPZCompElT3d::SideNodes[s][0]];
      sub[1] = fSubEl[TPZCompElT3d::SideNodes[s][1]];
      return;
   }
   if(side>9) {//faces
      int s = side-10;
      sub.Resize(4);
      sub[0] = fSubEl[TPZCompElT3d::FaceSons[s][0]];
      sub[1] = fSubEl[TPZCompElT3d::FaceSons[s][1]];
      sub[2] = fSubEl[TPZCompElT3d::FaceSons[s][2]];
      sub[3] = fSubEl[TPZCompElT3d::FaceSons[s][3]];
   }
}

TPZGeoElSide TPZGeoElT3d::Father(int side) {
	TPZGeoEl *fFather = TPZGeoEl::Father();

   if(!fFather) return TPZGeoElSide();
   int whichsub = -1;
   int i,nsubel = 6;
   if(fFather->NSides() == 19) nsubel = 10;//pai é pirâmide
   for(i=0; i<nsubel; i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) {
	   PZError << "TPZGeoElT3d::Father. fFather isn't father of this element.\n";
   	return TPZGeoElSide();
   }
   if(fFather->NSides() == 19) {//o pai é uma pirâmide
      if(whichsub == 6 &&  side==11) return TPZGeoElSide(fFather,14);
      if(whichsub == 7 &&  side==10) return TPZGeoElSide(fFather,15);
      if(whichsub == 8 &&  side==13) return TPZGeoElSide(fFather,16);
      if(whichsub == 9 &&  side==12) return TPZGeoElSide(fFather,17);
      if(side==14) return TPZGeoElSide(fFather,18);//pai do tetraedro pelo interior
      return TPZGeoElSide();
   }
   //agora o atual elemento é o filho numero whichsub < 6
   //os filhos interiores não tém pai associados a seus cantos
   if((side<4 && side == whichsub) || side==14) return TPZGeoElSide(fFather,side);//cantos
   //lados
   if(whichsub == 0 && (side== 4 || side== 6 || side== 7)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side== 4 || side== 5 || side== 8)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side== 5 || side== 6 || side== 9)) return TPZGeoElSide(fFather,side);
   if(whichsub == 3 && (side== 7 || side== 8 || side== 9)) return TPZGeoElSide(fFather,side);
   //faces
   if(whichsub == 0 && (side==10 || side==11 || side==13)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side==10 || side==11 || side==12)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side==10 || side==12 || side==13)) return TPZGeoElSide(fFather,side);
   if(whichsub == 3 && (side==11 || side==12 || side==13)) return TPZGeoElSide(fFather,side);
   if(whichsub == 4 && (side==14 || side==16))             return TPZGeoElSide(fFather,side-3);
   if(whichsub == 5 && (side==15 || side==17))             return TPZGeoElSide(fFather,side-5);
   //outro caso
   return TPZGeoElSide();
}

void TPZGeoElT3d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {

   int nsub = NSideSubElements(side);
   if(!nsub) return;
   sub.Resize(nsub);
   int i,j,k;
   if(nsub==1) {//side = 0 a 3
   	if(fSubEl[side]->NodeIndex(side)!=refnode[0]) {
      	PZError << "TPZGeoElT3d::GetSubElement subelement does not contain refnode" <<std::endl;
         return;
      }
	   sub[0]=TPZGeoElSide(fSubEl[side],side);
   	return;
   }
   //int isub=0;
   for(i=0;i<nsub;i++) {
   	TPZGeoElSide sidesub = SideSubElement(side,i);
      TPZGeoEl *subel = sidesub.Element();
		for(k = 0; k < refnode.NElements(); k++) {//k<nsub?
		   for(j=0;j<4;j++) {//4 vantos do pai
			   if(subel->NodeIndex(j)==refnode[k]) {
            	sub[k] = SideSubElement(side,i);//sub[k]?
            }
         }
      }
   }
   if(side > 9 && side < 14) sub[3] = SideSubElement(side,3);
   if(side == 14) {
   	sub[4] = SideSubElement(side,4);
      sub[5] = SideSubElement(side,5);
   }
   return;
}

/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
//transforma tetraedro para tetraedro ou tetraedro para pirâmide
//subelementos 0,1,2,3 do tetraedro -> tetraedro
//subelementos 6,7,8,9 da tetraedro -> pirâmide
//transformação entre lado filho e lado pai, o side é o side do filho
//o lado do filho esta contido no lado do pai
void TPZGeoElT3d::BuildTransform(int side,TPZGeoEl *father,TPZTransform<STATE> &t) {
	TPZGeoEl *fFather = TPZGeoEl::Father();
	if(!fFather || side > 14) return;
   int whichsub = -1;
   int i,nsides = fFather->NSides();
   if(nsides == 15)//o pai é um tetraedro
   	for(i=0;i<4;i++)  if(fFather->SubElement(i) == this) whichsub = i;
   if(nsides == 19)//o pai é uma pirâmide
   	for(i=6;i<10;i++) if(fFather->SubElement(i) == this) whichsub = i;
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
	   std::cout << "TPZGeoElT3d::BuildTransform could not identify the father element\n";
	   return;
   }

   if(side==14) {//pai para filho ou filho mestre para pai mestre
      mult(0,0) = 0.5;
      mult(1,1) = 0.5;
      mult(2,2) = 0.5;
      switch(whichsub) {//o atual é o tetraedro filho numero whichsub
         case 0:
            break;
         case 1:
            sum(0,0) = .5;
            break;
         case 2:
            sum(1,0) = .5;
            break;
         case 3:
         	sum(2,0) = .5;
            break;
         case 6:
         	mult.Zero();
            mult(0,0) =  0.5;
            mult(0,1) =  1.;
            mult(0,2) =  0.5;
            mult(1,2) =  0.5;
            mult(2,0) = -0.5;
            mult(2,2) = -0.5;
         	sum(0,0)  = -.5;
            sum(1,0)  = -.5;
            sum(2,0)  =  .5;
            break;
         case 7:
         	mult.Zero();
            mult(0,0) = -0.5;
            mult(0,1) = -0.5;
            mult(0,2) = -1.;
            mult(1,0) = -0.5;
            mult(1,1) =  0.5;
            mult(2,0) =  0.5;
            mult(2,1) =  0.5;
         	sum(0,0)  = 1.;
            break;
         case 8:
         	mult.Zero();
            mult(0,0) =  0.5;
            mult(0,1) =  0.5;
            mult(0,2) =  1.;
            mult(1,0) = -0.5;
            mult(1,1) =  0.5;
            mult(2,0) = -0.5;
            mult(2,1) = -0.5;
         	sum(0,0)  = -.5;
            sum(1,0)  =  .5;
            sum(2,0)  =  .5;
            break;
         case 9:
         	mult.Zero();
            mult(0,0) = -0.5;
            mult(0,1) = -0.5;
            mult(0,2) = -1.;
            mult(1,0) = -0.5;
            mult(1,1) =  0.5;
            mult(2,0) =  0.5;
            mult(2,1) =  0.5;
      }
   } else if(side>9) {//face do pai para face do filho
         whichsub = -1;
         int s = side-10;;
         if(nsides == 15) {
            for(i=0; i<3; i++) if(fFather->SubElement(TPZCompElT3d::FaceSons[s][i]) == this) whichsub = i;
         }
         if(nsides == 19) {
               i = TPZCompElT3d::MiddleFace[s]-13;//face da pirâmide que contem a face do tetraedro atual
            	if(fFather->SubElement(TPZCompElPi3d::FaceSons[i][3]) == this) whichsub = 3;
         }
         if(whichsub == -1) return;
         mult(0,0) = 0.5;
         mult(1,1) = 0.5;
         switch(whichsub) {//o atual é o filho numero whichsub
            case 0:        //ambos são tetraedros para os casos 0,1,2
               break;
            case 1:
               sum(0,0) = 0.5;
               break;
            case 2:
               sum(1,0) = 0.5;
               break;
            case 3://filho tetraedro e pai pirâmide
               if(side==11 || side==13) {//basta com o side do filho já
                  mult(0,1) = 0.5;        //que o pai é pirâmide
                  mult(1,0) =-0.5;
                  mult(1,1) = 0.;
                  sum(1,0)  = 0.5;
               } else
               if(side==10) {
                  mult(0,0) =-0.5;
                  mult(1,0) = 0.5;
                  sum(0,0)  = 0.5;
               } else
               if(side==12) {
                  mult(0,1) = 0.5;
                  mult(1,1) =-0.5;
                  sum(1,0)  = 0.5;
               }
         }
   } else if(side>3) {//4 a 8
      whichsub = -1;
      int s = side-4;
      for(i=0; i<2; i++) {
         if(nsides==15) if(fFather->SubElement(TPZCompElT3d::SideNodes[s][i]) == this) whichsub = i;
         if(nsides==19) if(fFather->SubElement(TPZCompElPi3d::SideNodes[s+1][i]) == this) whichsub = i;
      }
      if(whichsub == -1) return;
   	mult(0,0) = 0.5;
      if(whichsub==0) sum(0,0) = -0.5;
      if(whichsub==1) sum(0,0) =  0.5;
   }

   tloc.SetMatrix(mult,sum);
   t = tloc.Multiply(t);
   if(locfather != father) locfather->BuildTransform2(side,father,t);
}

TPZTransform<STATE> TPZGeoElT3d::SideToSideTransform(int sidefrom,int sideto) {
   if( (sidefrom > 9 &&  sideto < 14) || (sidefrom > 3 &&  sideto < 10) || sideto < 4) {
      PZError << "TPZGeoElT3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   int dimfrom = SideDimension(sidefrom);
   int dimto = SideDimension(sideto);
   if(dimfrom >= dimto){
      PZError << "TPZGeoElT3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   TPZTransform<STATE> trans(dimto,dimfrom);//retorna zerada
   //agora : sidefrom < sideto ou dimfrom < dimto
 if(sideto == 14) {//interior
   switch(sidefrom) {
      //faces para interior
               //row = dimto = 3 e col = dimfrom = 2
      case  10://sidefrom face de dim = 2 : t(3,2)
      	trans.Mult()(0,0) =  1.;//         fMult(row,col) : fMult(3,2)
         trans.Mult()(1,1) =  1.;//         fSum(row,1)    : fSum(3,1)
         break;
      case  11:
      	trans.Mult()(0,0) =  1.;//sides : 0123 -> 0154
         trans.Mult()(2,1) =  1.;//         R2  -> R3
         break;
      case  12:
      	trans.Mult()(0,0) = -1.;
         trans.Mult()(0,1) = -1.;
      	trans.Mult()(1,0) =  1.;//sides : 0123 -> 1265
         trans.Mult()(2,1) =  1.;
         trans.Sum()(0,0)  =  1.;
         break;
      case  13:
      	trans.Mult()(1,0) =  1.;//sides : 0123 -> 3267
         trans.Mult()(2,1) =  1.;
         break;
      //lados para interior
      		  //row = dimto = 3 e col = dimfrom = 1
      case  4://sidefrom lado de dim = 1 : t(3,1)
      	trans.Mult()(0,0) =  .5;//          fMult(row,col) : fMult(3,1)
         trans.Sum()(0,0)  =  .5;//          fSum(row,1)    : fSum(3,1)
         break;
      case  5:
      	trans.Mult()(0,0) = -.5;//sides :
         trans.Mult()(1,0) =  .5;
         trans.Sum()(0,0)  =  .5;//
         trans.Sum()(1,0)  =  .5;//
         break;
      case 6:
      	trans.Mult()(1,0) = -.5;//sides :
         trans.Sum()(1,0)  =  .5;
         break;
      case 7:
      	trans.Mult()(2,0) =  .5;//sides :
         trans.Sum()(2,0)  =  .5;
         break;
      case 8:
      	trans.Mult()(0,0) = -.5;//sides :
         trans.Mult()(2,0) =  .5;
         trans.Sum()(0,0)  =  .5;
         trans.Sum()(2,0)  =  .5;
         break;
      case 9:
      	trans.Mult()(1,0) = -.5;
         trans.Mult()(2,0) =  .5;
         trans.Sum()(1,0)  =  .5;//
         trans.Sum()(2,0)  =  .5;//
         break;
   	//canto para interior
			//row = dimto = 3 e col = dimfrom = 0
      case 0://sidefrom lado de dim = 0 : t(3,0)
      	trans.Sum()(0,0)  = 0.;//          fMult(row,col) : fMult(3,0)
         trans.Sum()(1,0)  = 0.;//          fSum(row,1)    : fSum(3,1)
         trans.Sum()(2,0)  = 0.;
         break;
      case 1:
      	trans.Sum()(0,0)  = 1.;
         trans.Sum()(1,0)  = 0.;
         trans.Sum()(2,0)  = 0.;
         break;
      case 2:
      	trans.Sum()(0,0)  = 0.;
         trans.Sum()(1,0)  = 1.;
         trans.Sum()(2,0)  = 0.;
         break;
      case 3:
      	trans.Sum()(0,0)  = 0.;
         trans.Sum()(1,0)  = 0.;
         trans.Sum()(2,0)  = 1.;
         break;
   }
 }//fim if sideto = 14
 if(sideto > 9 && sideto < 14) {
 	switch(sidefrom) {
   	//lados para faces
   		   //row = dimto = 2 e col = dimfrom = 1
      case  4://sidefrom lado de dim =10 : t(2,1)
      	if(sideto==10 || sideto==11) {
            trans.Mult()(0,0) = .5;//          fMult(row,col) : fMult(2,1)
            trans.Sum()(0,0)  = .5;//          fSum(row,1)    : fSum(2,1)
         }
         break;
      case  5:
      	if(sideto==10) {
         	trans.Mult()(0,0) = -.5;
            trans.Mult()(1,0) =  .5;
            trans.Sum()(0,0)  =  .5;
            trans.Sum()(1,0)  =  .5;
         } else
         if(sideto==12) {
            trans.Mult()(0,0) = .5;
            trans.Sum()(0,0)  = .5;
         }
         break;
      case 6:
      	if(sideto==10) {
            trans.Mult()(1,0) = -.5;
            trans.Sum()(1,0)  =  .5;
         } else
         if(sideto==13) {
            trans.Mult()(0,0) = -.5;
            trans.Sum()(0,0)  =  .5;
         }
         break;
      case 7:
      	if(sideto==11 || sideto==13) {
            trans.Mult()(1,0) = .5;
            trans.Sum()(1,0)  = .5;
         }
         break;
      case 8:
      	if(sideto==11) {
            trans.Mult()(0,0) = -.5;
            trans.Mult()(1,0) =  .5;
            trans.Sum()(0,0)  =  .5;
            trans.Sum()(1,0)  =  .5;
         } else
         if(sideto==12) {
            trans.Mult()(1,0) =  .5;
            trans.Sum()(1,0)  =  .5;
         }
         break;
      case 9:
      	if(sideto==12 || sideto==13) {
            trans.Mult()(0,0) = -.5;
            trans.Mult()(1,0) =  .5;
            trans.Sum()(0,0)  =  .5;
            trans.Sum()(1,0)  =  .5;
         }
         break;
   	//cantos para faces
   		   //row = dimto = 2 e col = dimfrom = 0
      case  0://sidefrom lado de dim =0 : t(2,0)
      	if(sideto==10 || sideto==11|| sideto==13) {
            trans.Sum()(0,0) = 0.;//          fMult(row,col) : fMult(2,0)
            trans.Sum()(1,0) = 0.;//          fSum(row,1)    : fSum(2,1)
         }
         break;
      case  1:
      	if(sideto==10 || sideto==11) {
            trans.Sum()(0,0) = 1.;
            trans.Sum()(1,0) = 0.;
         } else
         if(sideto==12) {
            trans.Sum()(0,0) = 0.;
            trans.Sum()(1,0) = 0.;
         }
         break;
      case  2:
      	if(sideto==10) {
            trans.Sum()(0,0) =  0.;
            trans.Sum()(1,0) =  1.;
         } else
         if(sideto==12 || sideto==13) {
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) =  0.;
         }
         break;
      case  3:
       	if(sideto==11 || sideto==12 || sideto==13) {
            trans.Sum()(0,0) = 0.;
            trans.Sum()(1,0) = 1.;
         }
         break;
      default:
	      PZError << "TPZGeoElT3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
   }//switch
 }//fim lados
 //cantos para faces
 if(sideto > 3 && sideto < 10) {//4 a 9
 	switch(sidefrom) {
   	//canto para lado
   		   //row = dimto = 1 e col = dimfrom = 0
            //sidefrom canto de dim =0 : t(1,0)
      case  0://fMult(row,col) : fMult(1,0) , fSum(row,1) : fSum(1,1)
         if(sideto== 6)              trans.Sum()(0,0)  =  1.;
         if(sideto== 4 || sideto==7) trans.Sum()(0,0)  = -1.;
         break;
      case  1:
      	if(sideto== 4)              trans.Sum()(0,0)  =  1.;
         if(sideto== 5 || sideto==8) trans.Sum()(0,0)  = -1.;
      	break;
      case  2:
      	if(sideto== 5)               trans.Sum()(0,0)  =  1.;
         if(sideto== 6 || sideto== 9) trans.Sum()(0,0)  = -1.;
      	break;
      case  3:
         if(sideto==7 || sideto==8 || sideto==9) trans.Sum()(0,0)  =  1.;
  	}//switch
 }//if cantos
   return trans;
}

TPZCompEl *TPZGeoElT3d::CreateBCCompEl(int side,int bc,TPZCompMesh &cmesh) {
	if(side<0 || side>14) return 0;

   if(side==14) {
	   std::cout << "TPZGeoElT3d::CreateBCCompEl with side = 14 not implemented\n";
      return 0;
   }
   int64_t index;
	if(side<4) {
      TPZMaterial *bcptr = cmesh.FindMaterial(bc);
      if(!bcptr) {
      PZError << "TPZGeoElT3d::CreateBCCompEl has no bc.\n";
      return 0;
      }
      TPZCompEl *cel = Reference();
      if(!cel) {
      PZError << "TPZGeoElT3d::CreateBCCompEl has no computational element\n";
      return 0;
      }
	  TPZManVector<int64_t> nodeindexes(1);
      TPZGeoElPoint *gel;
      nodeindexes[0] = fNodeIndexes[side];
      gel = new TPZGeoElPoint(nodeindexes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,side));
      return gel->CreateCompEl(cmesh,index);
   } else if (side > 3 && side < 10) {//side =4 a 9 : lados
      TPZManVector<int64_t> nodes(2);
      int s = side-4;
      nodes[0] = NodeIndex(TPZCompElT3d::SideNodes[s][0]);
      nodes[1] = NodeIndex(TPZCompElT3d::SideNodes[s][1]);
      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,TPZCompElT3d::SideNodes[s][0]));
      TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(this,TPZCompElT3d::SideNodes[s][1]));
      TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(this,side));
      return gel->CreateCompEl(cmesh,index);
   } else if (side > 9) {//side = 10 a 13 : faces
      TPZManVector<int64_t> nodes(3);
      int s = side-10;
      nodes[0] = NodeIndex(TPZCompElT3d::FaceNodes[s][0]);
      nodes[1] = NodeIndex(TPZCompElT3d::FaceNodes[s][1]);
      nodes[2] = NodeIndex(TPZCompElT3d::FaceNodes[s][2]);
      TPZGeoElT2d *gel = new TPZGeoElT2d(nodes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,TPZCompElT3d::FaceNodes[s][0]));
      TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(this,TPZCompElT3d::FaceNodes[s][1]));
      TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(this,TPZCompElT3d::FaceNodes[s][2]));
      TPZGeoElSide(gel,3).SetConnectivity(TPZGeoElSide(this,TPZCompElT3d::FaceSides[s][0]));
      TPZGeoElSide(gel,4).SetConnectivity(TPZGeoElSide(this,TPZCompElT3d::FaceSides[s][1]));
      TPZGeoElSide(gel,5).SetConnectivity(TPZGeoElSide(this,TPZCompElT3d::FaceSides[s][2]));
      TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(this,side));
      return gel->CreateCompEl(cmesh,index);
   } else PZError << "TPZGeoElT3d::CreateBCCompEl. Side = " << side <<std::endl;
   return 0;
}

void TPZGeoElT3d::NodeFaceIds(TPZVec<int> &ids,int face) {

	ids.Resize(3,-1);
   if((face>-1 && face<4) || (face>9 && face<14)) {
   	if(face>9) face = face-10;
      ids[0] = NodeIndex(TPZCompElT3d::FaceNodes[face][0]);
      ids[1] = NodeIndex(TPZCompElT3d::FaceNodes[face][1]);
      ids[2] = NodeIndex(TPZCompElT3d::FaceNodes[face][2]);
      return;
   }
   std::cout << "TPZCompElT3d::NodeFaceIds bad side , side = " << face <<std::endl;
}
//cada lado do filho esta contido em que lado do pai?
static int fatherside[6][19] = {
/*00*/{0,4,6,7,4,10,6,7,11,13,10,11,14,13,14,-1,-1,-1,-1},
/*01*/{4,1,5,8,4,5,10,11,8,12,10,11,12,14,14,-1,-1,-1,-1},
/*02*/{6,5,2,9,10,5,6,13,12,9,10,14,12,13,14,-1,-1,-1,-1},
/*03*/{7,8,9,3,11,12,13,7,8,9,14,11,12,13,14,-1,-1,-1,-1},
/*04*/{4,8,9,6,7,11,12,13,10,11,11,13,13,14,11,14,13,14,14},
/*05*/{8,4,6,9,5,11,10,13,12,12,10,10,12,14,14,10,14,12,14},
};
static int fatherside2[4][19] = {//tetraedro com pai pirâmide
/*06*/{9,5,13,10,14,13,18,14,14,18,18,14,18,18,18,-1,-1,-1,-1},
/*07*/{6,10,11,13,15,15,15,13,18,18,15,18,18,18,18,-1,-1,-1,-1},
/*08*/{12,13,7,11,18,13,16,16,18,16,18,18,18,16,18,-1,-1,-1,-1},
/*09*/{13,9,12,8,18,17,18,13,17,17,18,18,17,18,18,-1,-1,-1,-1} };

TPZGeoElSide TPZGeoElT3d::Father2(int side){//Augusto:09/01/01
	TPZGeoEl *fFather = TPZGeoEl::Father();

	if(side<0 || side>14 || !fFather){
		PZError << "TPZGeoElT3d::Father2 called error" <<std::endl;
        return TPZGeoElSide();
	}
	int subelindex = WhichSubel();
	if(fatherside[subelindex][side]<0){
		PZError << "TPZGeoElT3d::Father2 called with index error\n";
		return TPZGeoElSide();
	}
	if(fFather->NSides()==19){//pai pirâmide
  	return TPZGeoElSide(fFather,fatherside2[subelindex-6][side]);
  } else return TPZGeoElSide(fFather,fatherside[subelindex][side]);
}

static int subeldata[15][11][2] = {//TAMANHO DISTINTO
//lados do pai:{{conectividades dos filhos}}
/*00*/{{0,0}},//os lados dos filhos formam
/*01*/{{1,1}},//uma particao do lado do pai
/*02*/{{2,2}},
/*03*/{{3,3}},
/*04*/{{0,4},{0,1},{1,4}},
/*05*/{{1,5},{1,2},{2,5}},
/*06*/{{2,6},{2,0},{0,6}},
/*07*/{{0,7},{0,3},{3,7}},
/*08*/{{1,8},{1,3},{3,8}},
/*09*/{{2,9},{2,3},{3,9}},
/*10*/{{0,10},{1,10},{2,10},{5,15},{0,5},{1,6},{2,4}},
/*11*/{{0,11},{1,11},{3,11},{4,14},{0,8},{1,7},{3,4}},
/*12*/{{1,12},{2,12},{3,12},{5,17},{1,9},{2,8},{3,5}},
/*13*/{{0,13},{2,13},{3,13},{4,16},{0,9},{2,7},{3,6}},
/*14*/{{0,14},{1,14},{2,14},{3,14},{4,18},{5,18},{4,13},{0,12},{1,13},{2,11},{3,10}}
};

static int nsubeldata[15] = {1,1,1,1,3,3,3,3,3,3,7,7,7,7,11};


void TPZGeoElT3d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){

   subel.Resize(0);
   if(side<0 || side>14 || !HasSubElement(side)){
      PZError << "TPZGeoElT3d::GetSubelements2 called with error arguments\n";
      return;
   }
   int nsub = nsubeldata[side];
   for(int i=0;i<nsub;i++)
       subel.Push(TPZGeoElSide(fSubEl[subeldata[side][i][0]],
                                       subeldata[side][i][1]));
}

static REAL buildt[6][19][4][3] = {//por colunas
/*S0*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*05*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*08*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*09*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*10*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*11*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*12*/{{-0.5,0.5,0},{-0.5,0,0.5},{0.,0.,0.},{0.5,0,0}},
      /*13*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0,0,0.}},
      /*14*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,0.,0.}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
/*S1*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*06*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*07*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*09*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*10*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*11*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*12*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0,0,0.}},
      /*13*/{{0,0.5,0},{0,0,0.5},{0.,0.,0.},{0.5,0,0}},
      /*14*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,0.,0.}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
/*S2*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,1,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*07*/{{0,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*08*/{{0,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*10*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*11*/{{0.5,0.,0.},{0.,0.,0.5},{0.,0.,0.},{0.,0.5,0.}},
      /*12*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0.5,0,0.}},
      /*13*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0.5,0,0.}},
      /*14*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,.5,0.}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
/*S3*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,1}},
      /*04*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*05*/{{0.25,0,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*06*/{{0.25,0,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*10*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.5}},
      /*11*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*12*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0,0.5,0.}},
      /*13*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0,0.5,0.}},
      /*14*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,0.,.5}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
/*S4*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*08*/{{0.25,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*09*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*10*/{{-0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*11*/{{-0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*12*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*13*/{{-0.25,0.25,0.},{0.,0.,0.25},{0.,0.,0.},{0.25,0.25,0.25}},
      /*14*/{{0.,0.5,0.},{-0.5,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*15*/{{-0.5,0.5,0.},{-0.5,0.,0.},{0.,0.,0.},{0.5,0.,0.5}},
      /*16*/{{0.,0.5,0.},{-0.5,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*17*/{{-0.5,0.5,0.},{-0.5,0.,0.5},{0.,0.,0.},{0.5,0.,0.}},
      /*18*/{{0.,0.,.25},{-.25,.25,0.},{-.25,-.25,.25},{.25,.25,.25}}},
/*S5*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*06*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*07*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*08*/{{-0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*09*/{{0.25,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*10*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*11*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*12*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*13*/{{0.,0.,-0.25},{-0.25,0.25,0.},{0.,0.,0.},{0.25,0.25,0.25}},
      /*14*/{{0.,0.,-0.5},{0.,0.5,-0.5},{0.,0.,0.},{0.5,0.,0.5}},
      /*15*/{{-0.5,0.5,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*16*/{{0.,0.,-0.5},{0.5,0.,-0.5},{0.,0.,0.},{0.,0.5,0.5}},
      /*17*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*18*/{{0.,0.,-.25},{-.25,.25,0.},{.25,.25,-.25},{.25,.25,.25}}}
};

TPZTransform<STATE> TPZGeoElT3d::BuildTransform2(int side, TPZGeoEl *father){//Augusto:09/01/01

	if(side<0 || side>18 || !father){
  	PZError << "TPZGeoElT3d::BuildTransform2 side out of range or father null\n";
    return TPZTransform<STATE>(0,0);
  }
  TPZTransform<STATE> trans(3,3);
  int son = father->WhichSubel();
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

REAL TPZGeoElT3d::MidSideNode[15][3] = {
/*00*/{.0,.0},/*01*/{1.,.0},/*02*/{0.,1.,.0},/*03*/{.0,0.,1.0},/*04*/{.5,.0,.0},
/*05*/{.5,.5},/*06*/{0.,.5},/*07*/{0.,0.,.5},/*08*/{.5,0.,0.5},/*09*/{.0,.5,.5},
/*10*/{1./3.,1./3., 0.  }  ,/*11*/{1./3., .0  ,1./3.},
/*12*/{1./3.,1./3.,1./3.}  ,/*13*/{ 0.  ,1./3.,1./3.},/*14*/{1./4.,1./4.,1./4.} };

int TPZGeoElT3d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3), pss(3), pf(3), pfs(3);
  //point son, point side son, point father, point side father : elemento mestre
  for(sn=0;sn<6;sn++){
    TPZGeoEl *son = subs[sn];
    int nsides = son->NSides();
    for(sd=0;sd<nsides;sd++){
      TPZTransform<STATE> telsd(0,0);
      if(nsides==15){
        ps[0] = MidSideNode[sd][0];//tetraedro
        ps[1] = MidSideNode[sd][1];
        ps[2] = MidSideNode[sd][2];
        telsd = pzshape::TPZShapeTetra::TransformElementToSide(sd);
      } else if(nsides==19){
        ps[0] = TPZGeoElPi3d::MidSideNode[sd][0];//pirâmide
        ps[1] = TPZGeoElPi3d::MidSideNode[sd][1];
        ps[2] = TPZGeoElPi3d::MidSideNode[sd][2];
        telsd = pzshape::TPZShapePiram::TransformElementToSide(sd);
      }
      if(son->WhichSide(ps) != sd) std::cout << "Lado nao bate\n";
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform<STATE> t;
	  gel->BuildTransform2(sd, son,t);//para não furar pirâmide com pai tetraedro
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      /*if(nsides==15)*/
      telsd = pzshape::TPZShapeTetra::TransformSideToElement(sdfat); //else
      //if(nsides==19) telsd = TPZShapePiram::TransformSideToElement(sdfat);
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(nsides-1).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao errada\n";
        PZError << "son    = " << (son->Id()) <<std::endl;
        PZError << "father = " << ((son->Father2(nsides-1).Element())->Id()) <<std::endl;
        PZError << "side   = " << sd <<std::endl <<std::endl;
        int ok;
		std::cin >> ok;
      } else {
		  std::cout << "Transformacao OK!\n";
		  std::cout << "Filho/lado : " << son->Id() << "/" << sd <<std::endl;
		  std::cout << "Pai : " << son->Father2(nsides-1).Element()->Id() <<std::endl <<std::endl;
      }
    }
  }
  return 1;
}


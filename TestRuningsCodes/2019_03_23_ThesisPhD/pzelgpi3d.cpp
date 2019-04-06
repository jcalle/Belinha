//METHODS DEFINITION FOR CLASS ELEMPI3D
#include "pzelgpi3d.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelc1d.h"
#include "pzelgt2d.h"
#include "pzelgq2d.h"
#include "pzelgt3d.h"
#include "pzshapepiram.h"
#include "pzshapetetra.h"
#include "pzelcpi3d.h"
#include "pzelct3d.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgnode.h"
#include "pzshtmat.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include <stdlib.h>

static TPZCompEl *CreateEl(TPZGeoElPi3d *gel,TPZCompMesh &mesh, int64_t &index) {
  return new TPZCompElPi3d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElPi3d::fp)(TPZGeoElPi3d *,TPZCompMesh &, int64_t &) = CreateEl;

TPZGeoElPi3d::TPZGeoElPi3d(int id,TPZVec<int64_t> &nodeindexes,int matid,TPZGeoMesh &mesh):
	TPZGeoElRefLess<pzgeom::TPZGeoPyramid>(id,nodeindexes,matid,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=5) {
    PZError << "TPZGeoElPi3d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<5;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<10;i++) fSubEl[i] = 0;
}

TPZGeoElPi3d::TPZGeoElPi3d(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) :
	TPZGeoElRefLess<pzgeom::TPZGeoPyramid>(nodeindexes,matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=5) {
    PZError << "TPZGeoElPi3d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<5;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<10;i++) fSubEl[i] = 0;
}

TPZGeoElPi3d::~TPZGeoElPi3d() {}

void TPZGeoElPi3d::Shape(TPZVec<REAL> &pt,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {

  /*if(abs(pt[0])<1.e-10 && abs(pt[1])<1.e-10 && pt[2]==1.) {
  	//para testes com transformações geometricas
   //(0,0,1) nunca é um ponto de integração
     phi(0,0)  = 0.;
     phi(1,0)  = 0.;
     phi(2,0)  = 0.;
     phi(3,0)  = 0.;
     phi(4,0)  = 1.;
     for(int i=0;i<5;i++) {
        dphi(0,i) = 0.;
        dphi(1,i) = 0.;
        dphi(2,i) = 0.;
     }
     return;
  }*/
  REAL T0xz = .5*(1.-pt[2]-pt[0]) / (1.-pt[2]);
  REAL T0yz = .5*(1.-pt[2]-pt[1]) / (1.-pt[2]);
  REAL T1xz = .5*(1.-pt[2]+pt[0]) / (1.-pt[2]);
  REAL T1yz = .5*(1.-pt[2]+pt[1]) / (1.-pt[2]);
  REAL lmez = (1.-pt[2]);
  phi(0,0)  = T0xz*T0yz*lmez;
  phi(1,0)  = T1xz*T0yz*lmez;
  phi(2,0)  = T1xz*T1yz*lmez;
  phi(3,0)  = T0xz*T1yz*lmez;
  phi(4,0)  = pt[2];
  REAL lmexmez = 1.-pt[0]-pt[2];
  REAL lmeymez = 1.-pt[1]-pt[2];
  REAL lmaxmez = 1.+pt[0]-pt[2];
  REAL lmaymez = 1.+pt[1]-pt[2];
  dphi(0,0) = -.25*lmeymez / lmez;
  dphi(1,0) = -.25*lmexmez / lmez;
  dphi(2,0) = -.25*(lmeymez+lmexmez-lmexmez*lmeymez/lmez) / lmez;

  dphi(0,1) =  .25*lmeymez / lmez;
  dphi(1,1) = -.25*lmaxmez / lmez;
  dphi(2,1) = -.25*(lmeymez+lmaxmez-lmaxmez*lmeymez/lmez) / lmez;

  dphi(0,2) =  .25*lmaymez / lmez;
  dphi(1,2) =  .25*lmaxmez / lmez;
  dphi(2,2) = -.25*(lmaymez+lmaxmez-lmaxmez*lmaymez/lmez) / lmez;

  dphi(0,3) = -.25*lmaymez / lmez;
  dphi(1,3) =  .25*lmexmez / lmez;
  dphi(2,3) = -.25*(lmaymez+lmexmez-lmexmez*lmaymez/lmez) / lmez;

  dphi(0,4) =  0.0;
  dphi(1,4) =  0.0;
  dphi(2,4) =  1.0;
}

TPZGeoElPi3d *TPZGeoElPi3d::CreateGeoEl(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElPi3d(nodeindexes,matind,mesh);
}

int TPZGeoElPi3d::NNodes() {
  return 5;
}
int64_t TPZGeoElPi3d::NodeIndex(int node) {
  if(node<0 || node>4) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElPi3d::NSideNodes(int side) {
  if(side<0 || side>18) {
    PZError << "TPZGeoElPi3d::NSideNodes. Bad parameter side.\n";
    return 0;
  }
  if(side<5) return 1;//cantos
  if(side>4 && side<13) return 2;//lados
  if(side == 13) return 4;//face quadrilateral
  if(side<18) return 3;//faces triangulares
  return 5;//centro
}

int64_t TPZGeoElPi3d::SideNodeIndex(int side,int node) {
  if(side<0 || side>18 || node<0) {//19 sides
    PZError << "TPZGeoElPi3d::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  //0<=side<=14 e node>=0
  //side central
  if(side==18 && node<5) return fNodeIndexes[node];//um dos 5 cantos
  //side de canto
  if(side < 5) return fNodeIndexes[side];//canto side, deve ser node = 0
  //sides 5 a 12
  if(side>4 && side<13 && node<2) {//lados 5 a 12
      side-=5;
      return fNodeIndexes[TPZCompElPi3d::SideNodes[side][node]];
  } else if(side == 13 && node<4) {//face 13
  		side = 0;
	  	return fNodeIndexes[TPZCompElPi3d::FaceNodes[side][node]];
  } else if(side>13 && node<3) {//faces 14 a 18
  		side-=13;
	  	return fNodeIndexes[TPZCompElPi3d::FaceNodes[side][node]];
  }
  return -1;
}

void TPZGeoElPi3d::MidSideNodeIndex(int side, int64_t &index) {
  index = -1;
  if(side<0 || side>18) {
    PZError << "TPZGeoElPi3d::MidSideNodeIndex. Bad parameter side = " << side <<std::endl;
    return;
  }
  //sides 0 a 3
	if(side<5) {//o nó medio do lado 0 é o 0 etc.
		index=fNodeIndexes[side];
		return;
	}
   //o nó medio da face é o centro da face e o nó medio do centro é o centro
   //como nó de algum filho se este existir
   //caso tenha filhos é o canto de algum filho, se não tiver filhos retorna -1
	if(HasSubElement(0)) {
   	side-=5;
	   index=((TPZGeoElPi3d *) fSubEl[TPZCompElPi3d::MidSideNodes[side][0]])->fNodeIndexes[TPZCompElPi3d::MidSideNodes[side][1]];
	}
}

void TPZGeoElPi3d::NewMidSideNode(int side, int64_t &index) {
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
    if(side < 5) {index = -1; return;}
    //aqui side = 5 a 13
    side-=5;//0 a 8
    par[0] = TPZCompElPi3d::MidCoord[side][0];
    par[1] = TPZCompElPi3d::MidCoord[side][1];
    par[2] = TPZCompElPi3d::MidCoord[side][2];
    X(par,coord);
    index = Mesh()->NodeVec().AllocateNewElement();
    Mesh()->NodeVec()[index].Initialize(coord,*Mesh());
  }
}

int TPZGeoElPi3d::SideDimension(int side) {

	if (side<0 || side>18) {
   	PZError << "TPZGeoElPi3d::SideDimension called with side " << side <<std::endl;
      return 0;
   }
  if(side<5) return 0;//cantos
  if(side>4 && side<13) return 1;//lados
  if(side<18) return 2;//faces
  return 3;//centro

   }

TPZGeoElSide TPZGeoElPi3d::HigherDimensionSides(int side,int targetdimension) {
//targetdimension deve ser 1 , 2 ou 3
//se side =  0 a 4  targetdimension deve ser 1
//se side =  5 a 12  targetdimension deve ser 2
//se side = 18      targetdimension deve ser 3
  if( (side<0 || side>18) || (targetdimension < 1 || targetdimension > 3) ) {
     PZError << "TPZGeoElPi3d::HigherDimensionSides called with side = " << side
	          << " targetdimension = " << targetdimension <<std::endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
  int bestface;
  //side = 0 a 18
  switch(targetdimension) {//=1,2
	  case 1:
       if(father->NSides() == 15) return TPZGeoElSide();//o pai é um tetraedro 
     	 if(this == father->SubElement(0)) {
       	 if(side==1) return TPZGeoElSide(this,5);
       	 if(side==3) return TPZGeoElSide(this,8);
          if(side==4) return TPZGeoElSide(this,9);
       } else if(this == father->SubElement(1)) {
       	 if(side==0) return TPZGeoElSide(this,5);
       	 if(side==2) return TPZGeoElSide(this,6);
          if(side==4) return TPZGeoElSide(this,10);
       } else if(this == father->SubElement(2)) {
       	 if(side==1) return TPZGeoElSide(this,6);
       	 if(side==3) return TPZGeoElSide(this,7);
          if(side==4) return TPZGeoElSide(this,11);
       } else if(this == father->SubElement(3)) {
       	 if(side==0) return TPZGeoElSide(this,8);
       	 if(side==2) return TPZGeoElSide(this,7);
          if(side==4) return TPZGeoElSide(this,12);
       } else if(this == father->SubElement(4)) {
       	 if(side==0) return TPZGeoElSide(this,9);
       	 if(side==1) return TPZGeoElSide(this,10);
          if(side==2) return TPZGeoElSide(this,11);
          if(side==3) return TPZGeoElSide(this,12);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
  	  case 2:
     	 if(this == father->SubElement(0)) {

          bestface = BestDimensionSideOfTwoFaces(13,14);
          if((side==1 || side==5) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(13,17);
          if((side==3 || side==8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(14,17);
          if((side==4 || side==9) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==2 || side==6 || side==7) return TPZGeoElSide(this,13);
       	 if(side==10) return TPZGeoElSide(this,14);
       	 if(side==12) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(1)) {

          bestface = BestDimensionSideOfTwoFaces(13,14);
          if((side==0 || side==5) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(13,15);
          if((side==2 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(14,15);
          if((side==4 || side==10) && bestface) return TPZGeoElSide(this,bestface);
          if(side==3 || side==7 || side==8) return TPZGeoElSide(this,13);
       	 if(side==9) return TPZGeoElSide(this,14);
       	 if(side==11) return TPZGeoElSide(this,15);
       } else if(this == father->SubElement(2)) {

          bestface = BestDimensionSideOfTwoFaces(13,15);
          if((side==1 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(13,16);
          if((side==3 || side==7) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(15,16);
          if((side==4 || side==11) && bestface) return TPZGeoElSide(this,bestface);
          if(side==0 || side==5 || side==8) return TPZGeoElSide(this,13);
       	 if(side==10) return TPZGeoElSide(this,15);
       	 if(side==12) return TPZGeoElSide(this,16);
       } else if(this == father->SubElement(3)) {

          bestface = BestDimensionSideOfTwoFaces(13,17);
          if((side==0 || side==8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(13,16);
          if((side==2 || side==7) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,17);
          if((side==4 || side==12) && bestface) return TPZGeoElSide(this,bestface);
          if(side==1 || side==5 || side==6) return TPZGeoElSide(this,13);
       	 if(side==9) return TPZGeoElSide(this,17);
       	 if(side==11) return TPZGeoElSide(this,16);
       } else if(this == father->SubElement(4)) {

          if(father->NSides() == 15) {//o pai é um tetraedro

          	 TPZGeoElT3d *sub = (TPZGeoElT3d *) father->SubElement(0);
             bestface = sub->BestDimensionSideOfTwoFaces(10,11);
             if(side==0 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(10,13);
             if(side==3 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(11,13);
             if(side==4 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==8)             return TPZGeoElSide(sub,10);
             sub = (TPZGeoElT3d *) father->SubElement(3);
             bestface = sub->BestDimensionSideOfTwoFaces(11,12);
             if(side==1 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(12,13);
             if(side==2 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==6)             return TPZGeoElSide(sub,12);
             if(side== 5 || side== 9 || side==10) return TPZGeoElSide(this,14);
             if(side== 7 || side==11 || side==12) return TPZGeoElSide(this,16);
             return TPZGeoElSide();
          }
          bestface = BestDimensionSideOfTwoFaces(14,17);
          if((side==0 || side==9) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(14,15);
          if((side==1 || side==10) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(15,16);
          if((side==2 || side==11) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,17);
          if((side==3 || side==12) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==5) return TPZGeoElSide(this,14);
       	 if(side==6) return TPZGeoElSide(this,15);
       	 if(side==7) return TPZGeoElSide(this,16);
       	 if(side==8) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(5)) {

          if(father->NSides() == 15) {//o pai é um tetraedro

          	 TPZGeoElT3d *sub = (TPZGeoElT3d *) father->SubElement(1);
             bestface = sub->BestDimensionSideOfTwoFaces(10,11);
             if(side==1 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(11,12);
             if(side==0 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(10,12);
             if(side==4 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==5)             return TPZGeoElSide(sub,11);
             sub = (TPZGeoElT3d *) father->SubElement(2);
             bestface = sub->BestDimensionSideOfTwoFaces(10,13);
             if(side==2 && bestface) return TPZGeoElSide(sub,bestface);
             bestface = sub->BestDimensionSideOfTwoFaces(12,13);
             if(side==3 && bestface) return TPZGeoElSide(sub,bestface);
             if(side==7)             return TPZGeoElSide(sub,13);
             if(side== 8 || side== 9 || side==12) return TPZGeoElSide(this,17);
             if(side== 6 || side==10 || side==11) return TPZGeoElSide(this,15);
             return TPZGeoElSide();
          }
          TPZGeoElPi3d *sub = (TPZGeoElPi3d *) father->SubElement(0);
          bestface = sub->BestDimensionSideOfTwoFaces(14,17);
          if(side==1 && bestface==17) return TPZGeoElSide(father->SubElement(9),12);
          sub = (TPZGeoElPi3d *) father->SubElement(3);
          bestface = sub->BestDimensionSideOfTwoFaces(16,17);
          if(side==2 && bestface==17) return TPZGeoElSide(father->SubElement(9),12);
          sub = (TPZGeoElPi3d *) father->SubElement(2);
          bestface = sub->BestDimensionSideOfTwoFaces(15,16);
          if(side==3 && bestface==16) return TPZGeoElSide(father->SubElement(8),13);
       	 if(side==5) return TPZGeoElSide(father->SubElement(4),14);//(this,14)
       	 if(side==7) return TPZGeoElSide(father->SubElement(4),16);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
     case 3:
       return TPZGeoElSide(this,18);//0<=side<=18
  }//switch
  return TPZGeoElSide();
}

int TPZGeoElPi3d::BestDimensionSideOfTwoFaces(int face1,int face2) {

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


void TPZGeoElPi3d::LowerDimensionSides(int side,TPZStack<TPZGeoElSide> &smallsides) {
	if (side < 5) return;
   int i;
   if(side < 13) {//side = 5 a 12 : entram os cantos dos lados
   	int s = side-5;
   	smallsides.Push(TPZGeoElSide(this,TPZCompElPi3d::SideNodes[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElPi3d::SideNodes[s][1]));
   } else if(side < 18) {//entram cantos e lados da face
   	int s = side-13;
   	smallsides.Push(TPZGeoElSide(this,TPZCompElPi3d::FaceNodes[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElPi3d::FaceNodes[s][1]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElPi3d::FaceNodes[s][2]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElPi3d::FaceSides[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElPi3d::FaceSides[s][1]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElPi3d::FaceSides[s][2]));
      if(s==0) {//face quadrilateral
         smallsides.Push(TPZGeoElSide(this,TPZCompElPi3d::FaceNodes[s][3]));
         smallsides.Push(TPZGeoElSide(this,TPZCompElPi3d::FaceSides[s][3]));
      }
      smallsides.Push(TPZGeoElSide(this,side));
   } else if(side==18) {//entram todos os cantos, arestas e faces
   	for (i=0;i<18;i++) smallsides.Push(TPZGeoElSide(this,i));
   }
}

void TPZGeoElPi3d::Jacobian(TPZVec<REAL> &param,TPZFMatrix<STATE> &jacobian,TPZFMatrix<STATE> &axes,REAL &detjac,TPZFMatrix<STATE> &jacinv){

#ifdef DEBUG
  int nnodes = NNodes();
  if (nnodes != 5) {
    PZError << "TPZGeoElPi3d.jacobian only implemented for"
      " 5 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 3 || param[0] < 0. || param[0] > 1. ||
     param[1] < 0. || param[1] > 1. || param[2] < 0. || param[2] > 1.) {
    PZError << "TPZGeoElPi3d.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
    return;
  }
#endif
  REAL spacephi[5];
  TPZFMatrix<STATE> phi(5,1,spacephi,5);
  REAL spacedphi[15];
  TPZFMatrix<STATE> dphi(3,5,spacedphi,15);
  Shape(param,phi,dphi);
  jacobian.Zero();
  TPZGeoNode *np;
  int i,j;
  for(i=0;i<5;i++) {
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

void TPZGeoElPi3d::X(TPZVec<REAL> & loc,TPZVec<REAL> &result){
  REAL spacephi[10],spacedphi[20];
  int i,j;
  TPZFMatrix<STATE> phi(5,1,spacephi,10);
  TPZFMatrix<STATE> dphi(3,5,spacedphi,20);
  Shape(loc,phi,dphi);
  for(j=0;j<3;j++) {
    result[j] = 0.0;
    for(i=0;i<5;i++) result[j] += NodePtr(i)->Coord(j)*phi(i,0);
  }
}
/**It's necessary to define the normal vector to side 4, that is the orthogonal
   vector to the surface*/

void TPZGeoElPi3d::NormalVector(int /*side*/,TPZVec<REAL> &/*param*/,TPZVec<REAL> &/*normal*/,
			     TPZFMatrix<STATE> &/*axes*/,TPZFMatrix<STATE> &/*jac1d*/) {
/*
#ifdef DEBUG
  if (nnodes != 8) {
    PZError << "TPZGeoElPi3d.NormalVector, only implemented for"
      " 8 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 3 || param[0] < -1. || param[0] > 1. ||
     param[1] < -1. || param[1] > 1. || param[2] < -1. || param[2] > 1.) {
    PZError << "TPZGeoElPi3d.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
    return;
  }
  if(normal.NElements() != 3) {
    PZError << "TPZGeoElPi3d::NormalVector normal.capacity() = " << normal.NElements() <<
      "\n";
    return;
  }
  if(side < 0 || side > 5) {//6 faces
    PZError << "TPZGeoElPi3d.jacobian invalid side : "
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
    np = NodePtr(TPZCompElPi3d::FaceNodes[side][i]);
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
void TPZGeoElPi3d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {

   int i;
   if(HasSubElement(0)) {
      SubElVec.Resize(10);
      for(i=0;i<10;i++) SubElVec[i] = fSubEl[i];
      return;//If exist fSubEl return this sons
   }
   int j, sub, matid = MaterialId();
   int64_t index;
   int np[14];//guarda conectividades dos 6 subelementos

   for(j=0;j<5;j++) np[j]=NodeIndex(j);
   for(sub=5;sub<14;sub++) {
      NewMidSideNode(sub,index);
      np[sub] = index;
   }
   // creating new subelements
   for(i=0;i<6;i++) {
	   TPZManVector<int64_t> cornerindexes(5);
   	for(int j=0;j<5;j++) cornerindexes[j] = np[TPZCompElPi3d::CornerSons[i][j]];
      fSubEl[i] = CreateGeoEl(cornerindexes,matid,*Mesh());
   }
   for(;i<10;i++) {
	   TPZManVector<int64_t> cornerindexes(4);
      for(int j=0;j<4;j++) cornerindexes[j] = np[TPZCompElPi3d::CornerSons[i][j]];
      TPZGeoElT3d *subt=0;
      fSubEl[i] = subt->CreateGeoEl(cornerindexes,matid,*Mesh());
   }
   if(SubElVec.NElements()!=10) SubElVec.Resize(10);
   for(sub=0;sub<10;sub++) {
      SubElVec[sub] = fSubEl[sub];
      SubElVec[sub]->SetFather(this);
   }
   for(i=0;i<10;i++) {//conectividades entre os filhos : viz interna
   	for(j=0;j<18;j++) {
      	int elside = TPZCompElPi3d::InNeigh[i][j][0];
         if(elside == -1) break;
      	fSubEl[i]->SetNeighbour(elside,TPZGeoElSide(fSubEl[TPZCompElPi3d::InNeigh[i][j][1]],TPZCompElPi3d::InNeigh[i][j][2]));
      }                              //lado do subel                                          numero do filho viz.             lado do viz.
   }
   //vizinhança externa ao elemento atual
   //procura-se um viz pela face que seja dividido
   TPZGeoElSide dividedneighbour;
   for(int face=13;face<18;face++) {
      dividedneighbour = Neighbour(face);//BuildConnectivity(..) inicializou
      while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         //achou-se um viz para alguma face
         TPZManVector<int64_t> nodes(4);
         TPZStack<TPZGeoElSide> neighsub;
         int f = face-13;
         for(i=0;i<4;i++) nodes[i] = NodeIndex(TPZCompElPi3d::FaceNodes[f][i]);//nós globais na ordem e sentido local da face atual
         if(f > 0) nodes.Resize(3);
         dividedneighbour.GetSubElements2(neighsub);//3 ou 4 subs viz na ordem dos subs da face atual
         if(f == 0) {//só para a face 13
            nodes.Resize(1);
            //geoside central a face
            int fs0 = TPZCompElPi3d::FaceSons[f][0];//1o filho da face
            int fn2 = TPZCompElPi3d::FaceNodes[f][2];//nó 2 da face
            nodes[0] = fSubEl[fs0]->NodeIndex(fn2);//nó do centro da face como nó do filho nessa face
            int locside = neighsub[0].Element()->WhichSide(nodes);//lado do viz do sub conectado ao nó centro da face
            TPZGeoElSide sub(fSubEl[fs0],fn2);
            //conectividade do nó do centro da face fechando o ciclo entre os filhos dessa face conectados pelo nó
            sub.SetConnectivity(TPZGeoElSide(neighsub[0].Element(),locside));
         }
         //conectividades dos lados dos subelementos interiores a face
         nodes.Resize(2);
         for(i=0; i<4; i++) {//4 subelementos associados à face
            int fsi = TPZCompElPi3d::FaceSons[f][i];
            int inside = TPZCompElPi3d::FaceInRib[f][i];
            if(inside == -1) continue;
            nodes[0] = fSubEl[fsi]->SideNodeIndex(inside,0);
            nodes[1] = fSubEl[fsi]->SideNodeIndex(inside,1);
            int locside = neighsub[i].Element()->WhichSide(nodes);
            TPZGeoElSide sub(fSubEl[fsi],inside);
            TPZGeoElSide subneigh(neighsub[i].Element(),locside);
            sub.SetConnectivity(subneigh);
         }
         //geoside de face para os 4 subs da face f
         for(i=0; i<4; i++) {//4 subelementos associados à face
            if(i==3) {
            	TPZGeoElSide sub = SideSubElement(face,i);
               sub.SetConnectivity(neighsub[i]);
               continue;
            }
            TPZGeoElSide sub(fSubEl[TPZCompElPi3d::FaceSons[f][i]],face);
            sub.SetConnectivity(neighsub[i]);
         }
      }
   }//fim faces
	//procura-se um viz pela aresta do atual, que seja dividido
   for(int side=5; side<13; side++) {//(side 5 : subs 0,1),(side 6 : subs 1,2) ...
   	int s = side-5;
   	dividedneighbour = Neighbour(side);
      while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         //achou-se um viz pelo lado i e dividido
		  TPZManVector<int64_t> nodestore(2);
         int i = TPZCompElPi3d::SideNodes[s][0];//nó local 0 do lado side
         int iplus = TPZCompElPi3d::SideNodes[s][1];//nó local 1 do lado
		 TPZStack<TPZGeoElSide> neighsubstore;
		 nodestore[0] = SideNodeIndex(side,0);//nó global 0 do lado side
		 nodestore[1] = SideNodeIndex(side,1);//nó global 1
         dividedneighbour.GetSubElements2(neighsubstore);//filhos do viz conectados nesses lados
         TPZGeoElSide sub(fSubEl[i],side);
         //conectividade pelo lado side entre os filhos
         sub.SetConnectivity(neighsubstore[0]);
         sub = TPZGeoElSide(fSubEl[iplus],side);
         sub.SetConnectivity(neighsubstore[1]);
		 nodestore.Resize(1);
		 nodestore[0] = fSubEl[i]->SideNodeIndex(iplus,0);//nó do medio do lado como nó de um dos filhos
         int locside = neighsubstore[0].Element()->WhichSide(nodestore);
         sub = TPZGeoElSide(fSubEl[i],iplus);
         //conectividade para o nó do medio do lado side do atual como nó de um dos filhos
         sub.SetConnectivity(TPZGeoElSide(neighsubstore[0].Element(),locside));
      }
   }
   //procura-se um viz pelo canto do atual, que seja dividido
   for(int corner=0; corner<5; corner++) {//cantos
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


int TPZGeoElPi3d::NSubElements() {
  return 10;
}

int TPZGeoElPi3d::NSideSubElements(int side) {
  if(side < 0 || side > 18) {
    PZError << "TPZGeoElPi3d::NSideSubElements called for side " << side <<std::endl;
    return 0;
  }
  if(side==18) return 10;//centro
  if(side>12 && side<18) return 4;//faces
  if(side>4) return 2;//lados
  return 1;//cantos
}

TPZGeoElSide TPZGeoElPi3d::SideSubElement(int side,int position) {
   if (position<0 || position>10 || side <0 ||side>18) {
   	PZError << "TPZGeoElPi3d::SideSubElement called with position " << position << " side " << side <<std::endl;
      return TPZGeoElSide();
   }                              //fSubEl[is]
   if(side==18) {
      if(position > 5) return TPZGeoElSide(SubElement(position),14);//centro : tetraedro
      return TPZGeoElSide(SubElement(position),18);//centro : pirâmides
   }
   if(side<5) {//cantos
      if(position!=0) {
         PZError << "TPZGeoElPi3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
         return TPZGeoElSide(SubElement(side),side);
      }
   }
   if(side>4 && side<13) {//lados
       if(position!=0 && position!=1) {
         PZError << "TPZGeoElPi3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
      	int s = side-5;
         return TPZGeoElSide(SubElement(TPZCompElPi3d::SideNodes[s][position]),side);
      }
   }
   if(side>12) {//faces
       if(position<0 || position>4) {//position!=0 && position!=1 && position!=2 && position!=3
         PZError << "TPZGeoElPi3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
      	int s = side-13;//s = 0,1,2,3,4
         if(s==1 && position == 3) return TPZGeoElSide(SubElement(TPZCompElPi3d::FaceSons[s][3]),11);
         if(s==2 && position == 3) return TPZGeoElSide(SubElement(TPZCompElPi3d::FaceSons[s][3]),10);
         if(s==3 && position == 3) return TPZGeoElSide(SubElement(TPZCompElPi3d::FaceSons[s][3]),13);
         if(s==4 && position == 3) return TPZGeoElSide(SubElement(TPZCompElPi3d::FaceSons[s][3]),12);
         return TPZGeoElSide(SubElement(TPZCompElPi3d::FaceSons[s][position]),side);
      }
   }
   return TPZGeoElSide();
}

void TPZGeoElPi3d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
   if(!fSubEl[0]) {
      sub.Resize(0);
      return;
   }
   if(side < 0 || side > 18) {
      PZError << "TPZGeoElPi3d::SideSubElements called for side " << side <<std::endl;
      return;
   }
   if(side==18) {
      sub.Resize(10);
      for(int i=0;i<10;i++) sub[i] = fSubEl[i];
      return;
   }
   if(side<5) {
      sub.Resize(1);
      sub[0]=fSubEl[side];
      return;
   }
   if(side>4 && side<13) {//lados
      int s = side-5;
      sub.Resize(2);
      sub[0] = fSubEl[TPZCompElPi3d::SideNodes[s][0]];
      sub[1] = fSubEl[TPZCompElPi3d::SideNodes[s][1]];
      return;
   }
   if(side>12) {//faces
      int s = side-13;
      sub.Resize(4);
      sub[0] = fSubEl[TPZCompElPi3d::FaceSons[s][0]];
      sub[1] = fSubEl[TPZCompElPi3d::FaceSons[s][1]];
      sub[2] = fSubEl[TPZCompElPi3d::FaceSons[s][2]];
      sub[3] = fSubEl[TPZCompElPi3d::FaceSons[s][3]];
   }
}

TPZGeoElSide TPZGeoElPi3d::Father(int side) {
	TPZGeoEl *fFather = TPZGeoEl::Father();
	if (!fFather) return TPZGeoElSide();

   int whichsub = -1;
   int i,nsides = fFather->NSides();
   if(nsides == 15)//o pai é um tetraedro
   	for(i=4;i<6;i++)  if(fFather->SubElement(i) == this) whichsub = i;
   if(nsides == 19)//o pai é uma pirâmide
   	for(i=0;i<6;i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) {
	   PZError << "TPZGeoElPi3d::Father. fFather isn't father of this element.\n";
   	return TPZGeoElSide();
   }
   if(nsides == 15) {//face da pirâmide cujo pai é um tetraedro
      //if(side >13) {
         if(whichsub == 4 && (side==14 || side==16)) return TPZGeoElSide(fFather,side-3); else
         if(whichsub == 5 && (side==15 || side==17)) return TPZGeoElSide(fFather,side-5);
				 if(side==18) return TPZGeoElSide(fFather,14);//pai tetraedro para o interior do filho pirâmide
      //}
      return TPZGeoElSide();
   }
   //agora o atual elemento é o filho numero whichsub < 10
   //os filhos interiores não tém pai associados a seus cantos
   if((side<5 && side == whichsub) || side==18) return TPZGeoElSide(fFather,side);//cantos
   //lados
   if(whichsub == 0 && (side== 5 || side== 8 || side== 9)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side== 5 || side== 6 || side==10)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side== 6 || side== 7 || side==11)) return TPZGeoElSide(fFather,side);
   if(whichsub == 3 && (side== 7 || side== 8 || side==12)) return TPZGeoElSide(fFather,side);
   if(whichsub == 4 && (side== 9 || side==10 || side==11 || side==12)) return TPZGeoElSide(fFather,side);
   //faces
   if(whichsub == 0 && (side==13 || side==14 || side==17)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side==13 || side==14 || side==15)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side==13 || side==15 || side==16)) return TPZGeoElSide(fFather,side);
   if(whichsub == 3 && (side==13 || side==16 || side==17)) return TPZGeoElSide(fFather,side);
   if(whichsub == 4 && (side==14 || side==15 || side==16 || side==17)) return TPZGeoElSide(fFather,side);
   if(whichsub == 6 &&  side==11)                          return TPZGeoElSide(fFather,11);
   if(whichsub == 7 &&  side==10)                          return TPZGeoElSide(fFather,10);
   if(whichsub == 8 &&  side==13)                          return TPZGeoElSide(fFather,13);
   if(whichsub == 9 &&  side==12)                          return TPZGeoElSide(fFather,12);
   //outro caso
   return TPZGeoElSide();
}

void TPZGeoElPi3d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {

   int nsub = NSideSubElements(side);
   if(!nsub) return;
   sub.Resize(nsub);
   int i,j,k;
   if(nsub==1) {//side = 0 a 4
   	if(fSubEl[side]->NodeIndex(side)!=refnode[0]) {
      	PZError << "TPZGeoElPi3d::GetSubElement subelement does not contain refnode" <<std::endl;
         return;
      }
	   sub[0]=TPZGeoElSide(fSubEl[side],side);
   	return;
   }
   //int isub=0;
   for(i=0;i<nsub;i++) {
   	TPZGeoElSide sidesub = SideSubElement(side,i);
      TPZGeoEl *subel = sidesub.Element();
      int nnod = 4,son;
      for(son=6;son<10;son++) if(fSubEl[son] == subel) break;
      if(son == 10) nnod = 5;
		for(k = 0; k < refnode.NElements(); k++) {//se o subelemento k do thisside tiver o canto refnode[k]
		   for(j=0;j<nnod;j++) {  //este é um subelemento procurado
			   if(subel->NodeIndex(j)==refnode[k]) {
            	sub[k] = SideSubElement(side,i);
            }
         }
      }
   }
   if(side > 13 && side < 18) sub[3] = SideSubElement(side,3);
   if(side == 18) {
   	sub[5] = SideSubElement(side,5);
      sub[6] = SideSubElement(side,6);
   	sub[7] = SideSubElement(side,7);
      sub[8] = SideSubElement(side,8);
      sub[9] = SideSubElement(side,9);
   }
   return;
}

/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
//transforma piramide para piramide ou pirâmide para tetraedro
void TPZGeoElPi3d::BuildTransform(int side,TPZGeoEl *father,TPZTransform<STATE> &t) {
	TPZGeoEl *fFather = TPZGeoEl::Father();
	if(!fFather || side > 18) return;
   int whichsub = -1;
   int i,nsides = fFather->NSides();
   if(nsides == 15)//o pai é um tetraedro
   	for(i=4;i<6;i++)  if(fFather->SubElement(i) == this) whichsub = i;
   if(nsides == 19)//o pai é uma pirâmide
   	for(i=0;i<6;i++) if(fFather->SubElement(i) == this) whichsub = i;
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
	   std::cout << "TPZGeoElPi3d::BuildTransform could not identify the father element\n";
	   return;
   }

   if(side == 18) {//pai mestre para filho
      mult(0,0) = 0.5;//ou transformação entre elemento mestre do filho
      mult(1,1) = 0.5;//para o elemento mestre do pai
      mult(2,2) = 0.5;
      switch(whichsub) {//o atual é o filho numero whichsub
         case 0:
         	sum(0,0) = -.5;
            sum(1,0) = -.5;
            break;
         case 1:
         	sum(0,0) =  .5;
            sum(1,0) = -.5;
            break;
         case 2:
         	sum(0,0) = .5;
            sum(1,0) = .5;
            break;
         case 3:
         	sum(0,0) = -.5;
            sum(1,0) =  .5;
            break;
         case 4:
            if(nsides==15) {//pai tetraedro
               mult.Zero();
               mult(0,1) = -0.25;
               mult(0,2) = -0.25;
               mult(1,1) =  0.25;
               mult(1,2) = -0.25;
               mult(2,0) =  0.25;
               mult(2,2) =  0.25;
               sum(0,0)  = .25;
               sum(1,0)  = .25;
               sum(2,0)  = .25;
            } else {//pai pirâmide
               sum(2,0) = .5;
            }
         case 5:
         	if(nsides==19) {//pai pirâmide
               mult(0,0) *= -1.;
               mult(2,2) *= -1.;
               sum(2,0)   = .5;
            } else {//pai tetraedro
               mult.Zero();
               mult(0,1) = -0.25;
               mult(0,2) =  0.25;
               mult(1,1) =  0.25;
               mult(1,2) =  0.25;
               mult(2,0) = -0.25;
               mult(2,2) = -0.25;
               sum(0,0)  = .25;
               sum(1,0)  = .25;
               sum(2,0)  = .25;
            }
      }
   } else if(side>12) {//face do pai para face do filho
      whichsub = -1;
      int s = side-13;//lado da pirâmide
      if(nsides == 15) {//pai tetraedro
      	i = TPZCompElPi3d::MiddleFace[s-1]-10;//face do tetraedro que contém a face da pirâmide atual
         if(fFather->SubElement(TPZCompElT3d::FaceSons[i][3]) == this) whichsub = 3;
      }
      if(nsides == 19) {//pai pirâmide
         int n = 3;
         if(side==13) n = 4;
         for(i=0; i<n; i++) if(fFather->SubElement(TPZCompElPi3d::FaceSons[s][i]) == this) whichsub = i;
      }
      if(whichsub == -1) return;
      mult(0,0) = 0.5;
      mult(1,1) = 0.5;
   	if(side == 13) {//face quadrilateral
			sum(0,0) = -0.5;
         sum(1,0) = -0.5;
         switch(whichsub) {//o atual é a face numero whichsub do filho dentro da face do pai
            case 0:
               break;
            case 1:
               sum(0,0) = .5;
               break;
            case 2:
               sum(0,0) = .5;
               sum(1,0) = .5;
               break;
            case 3:
               sum(1,0) = .5;
         }
      } else {//face triangular
         switch(whichsub) {//o atual é o filho numero whichsub
            case 0://filho pirâmide e pai pirâmide casos 0,1,2
               break;
            case 1:
               sum(0,0) = 0.5;
               break;
            case 2:
               sum(1,0) = 0.5;
               break;
            case 3://filho piramide e pai tetraedro caso 3
               if(side==14 || side==16) {//basta com o side ja que o side do tetraedro acaba em 13 para as faces
                  mult(0,0) = 0.;
                  mult(0,1) =-0.5;
                  mult(1,0) = 0.5;
                  sum(0,0) = 0.5;
                  sum(1,0) = 0.;
               } else
               if(side==15) {
                  mult(0,0) =-0.5;
                  mult(1,0) = 0.5;
                  sum(0,0) = 0.5;
                  sum(1,0) = 0.;
               } else
               if(side==17) {
                  mult(0,1) = 0.5;
                  mult(1,1) =-0.5;
                  sum(0,0) = 0.;
                  sum(1,0) = 0.5;
               }
         }
      }
   } else if(side>4) {//e side < 13
      whichsub = -1;
      int s = side-5;//0 a 7
      for(i=0; i<2; i++) {
         if(nsides==15) if(fFather->SubElement(TPZCompElT3d::SideNodes[s-1][i]) == this) whichsub = i;
         if(nsides==19) if(fFather->SubElement(TPZCompElPi3d::SideNodes[s][i]) == this) whichsub = i;
      }
      if(whichsub == -1) return;
   	mult(0,0) = 0.5;
      if(whichsub==0) sum(0,0) = -0.5;//subelemento 0 do lado
      if(whichsub==1) sum(0,0) =  0.5;//subelemento 1 do lado
   }
   tloc.SetMatrix(mult,sum);
   t = tloc.Multiply(t);
   if(locfather != father) locfather->BuildTransform2(side,father,t);
}

TPZTransform<STATE> TPZGeoElPi3d::SideToSideTransform(int sidefrom,int sideto) {
   if( (sidefrom > 12 &&  sideto < 18) || (sidefrom > 4 &&  sideto < 13) || sideto < 5) {
      PZError << "TPZGeoElPi3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   int dimfrom = SideDimension(sidefrom);
   int dimto = SideDimension(sideto);
   if(dimfrom >= dimto){
      PZError << "TPZGeoElPi3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   TPZTransform<STATE> trans(dimto,dimfrom);//retorna zerada
   //agora : sidefrom < sideto ou dimfrom < dimto
 if(sideto == 18) {//interior
   switch(sidefrom) {
      //faces para interior!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               //row = dimto = 3 e col = dimfrom = 2
      case  13://sidefrom face de dim = 2 : t(3,2)
      	trans.Mult()(0,0) =  1.;//         fMult(row,col) : fMult(3,2)
         trans.Mult()(1,1) =  1.;//         fSum(row,1)    : fSum(3,1)
         break;
      case  14:
      	trans.Mult()(0,0) =  2.;
         trans.Mult()(0,1) =  1.;
         trans.Mult()(1,1) =  1.;
         trans.Mult()(2,1) =  1.;
         trans.Sum()(0,0)  = -1.;
         trans.Sum()(1,0)  = -1.;
         break;
      case  15:
         trans.Mult()(0,1) = -1.;
         trans.Mult()(1,0) =  2.;
         trans.Mult()(1,1) =  1.;
         trans.Mult()(2,1) =  1.;
         trans.Sum()(0,0)  =  1.;
         trans.Sum()(1,0)  = -1.;
         break;
      case  16:
         trans.Mult()(0,0) =  2.;
         trans.Mult()(0,1) =  1.;
         trans.Mult()(1,1) = -1.;
         trans.Mult()(2,1) =  1.;
         trans.Sum()(0,0)  = -1.;
         trans.Sum()(1,0)  =  1.;
         break;
      case  17:
         trans.Mult()(0,1) =  1.;
         trans.Mult()(1,0) =  2.;
         trans.Mult()(1,1) =  1.;
         trans.Mult()(2,1) =  1.;
         trans.Sum()(0,0)  = -1.;
         trans.Sum()(1,0)  = -1.;
         break;
      //lados para interior!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      		  //row = dimto = 3 e col = dimfrom = 1
      case  5://sidefrom lado de dim = 1 : t(3,1)
      	trans.Mult()(0,0) =  1.;//          fMult(row,col) : fMult(3,1)
         trans.Sum()(1,0)  = -1.;//          fSum(row,1)    : fSum(3,1)
         break;
      case  6:
         trans.Mult()(1,0) =  1.;
         trans.Sum()(0,0)  =  1.;
         break;
      case 7:
      	trans.Mult()(0,0) = -1.;
         trans.Sum()(1,0)  =  1.;
         break;
      case 8:
      	trans.Mult()(1,0) = -1.;
         trans.Sum()(0,0)  = -1.;
         break;
      case 9:
      	trans.Mult()(0,0) =  .5;
         trans.Mult()(1,0) =  .5;
         trans.Mult()(2,0) =  .5;
         trans.Sum()(0,0)  = -.5;
         trans.Sum()(1,0)  = -.5;
         trans.Sum()(2,0)  =  .5;
         break;
      case 10:
      	trans.Mult()(0,0) = -.5;
         trans.Mult()(1,0) =  .5;
         trans.Mult()(2,0) =  .5;
         trans.Sum()(0,0)  =  .5;
         trans.Sum()(1,0)  = -.5;
         trans.Sum()(2,0)  =  .5;
         break;
      case 11:
      	trans.Mult()(0,0) = -.5;
         trans.Mult()(1,0) = -.5;
         trans.Mult()(2,0) =  .5;
         trans.Sum()(0,0)  =  .5;
         trans.Sum()(1,0)  =  .5;
         trans.Sum()(2,0)  =  .5;
         break;
      case 12:
      	trans.Mult()(0,0) =  .5;
         trans.Mult()(1,0) = -.5;
         trans.Mult()(2,0) =  .5;
         trans.Sum()(0,0)  = -.5;
         trans.Sum()(1,0)  =  .5;
         trans.Sum()(2,0)  =  .5;
         break;
   	//canto para interior !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			    //row = dimto = 3 e col = dimfrom = 0
      case 0://sidefrom lado de dim = 0 : t(3,0)
      	trans.Sum()(0,0)  =-1.;//          fMult(row,col) : fMult(3,0)
         trans.Sum()(1,0)  =-1.;//          fSum(row,1)    : fSum(3,1)
         trans.Sum()(2,0)  = 0.;
         break;
      case 1:
      	trans.Sum()(0,0)  = 1.;
         trans.Sum()(1,0)  =-1.;
         trans.Sum()(2,0)  = 0.;
         break;
      case 2:
      	trans.Sum()(0,0)  = 1.;
         trans.Sum()(1,0)  = 1.;
         trans.Sum()(2,0)  = 0.;
         break;
      case 3:
      	trans.Sum()(0,0)  =-1.;
         trans.Sum()(1,0)  = 1.;
         trans.Sum()(2,0)  = 0.;
         break;
      case 4:
      	trans.Sum()(0,0)  = 0.;
         trans.Sum()(1,0)  = 0.;
         trans.Sum()(2,0)  = 1.;
         break;
   }
 }//fim if sideto = 18
 if(sideto > 12 && sideto < 18) {
 	switch(sidefrom) {
   	//lados para faces!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   		   //row = dimto = 2 e col = dimfrom = 1
      case  5://sidefrom lado de dim =10 : t(2,1)
      	if(sideto==13) {
            trans.Mult()(0,0) = 1.;//          fMult(row,col) : fMult(2,1)
            trans.Sum()(1,0)  =-1.;//          fSum(row,1)    : fSum(2,1)
         } else
         if(sideto==14) {
            trans.Mult()(0,0) = .5;
            trans.Sum()(0,0)  = .5;
         }
         break;
      case  6:
      	if(sideto==13) {
            trans.Mult()(1,0) = 1.;
            trans.Sum()(0,0)  = 1.;
         } else
         if(sideto==15) {
            trans.Mult()(0,0) = .5;
            trans.Sum()(0,0)  = .5;
         }
         break;
      case 7:
      	if(sideto==13) {
            trans.Mult()(0,0) = -1.;
            trans.Sum()(1,0)  =  1.;
         } else
         if(sideto==16) {
            trans.Mult()(0,0) = -.5;
            trans.Sum()(0,0)  =  .5;
         }
         break;
      case 8:
      	if(sideto==13) {
            trans.Mult()(1,0) = -1.;
            trans.Sum()(0,0)  = -1.;
         } else
         if(sideto==17) {
            trans.Mult()(0,0) = -.5;
            trans.Sum()(0,0)  =  .5;
         }
         break;
      case 9:
      	if(sideto==14 || sideto==17) {
            trans.Mult()(1,0) = .5;
            trans.Sum()(1,0)  = .5;
         }
         break;
      case 10:
      	if(sideto==14) {
            trans.Mult()(0,0) = -.5;
            trans.Mult()(1,0) =  .5;
            trans.Sum()(0,0)  =  .5;
            trans.Sum()(1,0)  =  .5;
         } else
         if(sideto==15) {
            trans.Mult()(1,0) =  .5;
            trans.Sum()(1,0)  =  .5;
         }
         break;
      case 11:
      	if(sideto==15 || sideto==16) {
            trans.Mult()(0,0) = -.5;
            trans.Mult()(1,0) =  .5;
            trans.Sum()(0,0)  =  .5;
            trans.Sum()(1,0)  =  .5;
         }
         break;
      case 12:
      	if(sideto==17) {
            trans.Mult()(0,0) = -.5;
            trans.Mult()(1,0) =  .5;
            trans.Sum()(0,0)  =  .5;
            trans.Sum()(1,0)  =  .5;
         } else
         if(sideto==16) {
            trans.Mult()(1,0) =  .5;
            trans.Sum()(1,0)  =  .5;
         }
         break;
   	//cantos para faces!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   		   //row = dimto = 2 e col = dimfrom = 0
      case  0://sidefrom lado de dim =0 : t(2,0)
      	if(sideto==14 || sideto==17) {
            trans.Sum()(0,0) = 0.;//          fMult(row,col) : fMult(2,0)
            trans.Sum()(1,0) = 0.;//          fSum(row,1)    : fSum(2,1)
         } else
         if(sideto==13) {
            trans.Sum()(0,0) =-1.;
            trans.Sum()(1,0) =-1.;
         }
         break;
      case  1:
      	if(sideto==13) {
            trans.Sum()(0,0) = 1.;
            trans.Sum()(1,0) =-1.;
         } else
         if(sideto==14) {
            trans.Sum()(0,0) = 1.;
            trans.Sum()(1,0) = 0.;
         } else
         if(sideto==15) {
            trans.Sum()(0,0) = 0.;
            trans.Sum()(1,0) = 0.;
         }
         break;
      case  2:
      	if(sideto==13) {
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) =  1.;
         } else
         if(sideto==15 || sideto==16) {
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) =  0.;
         }
         break;
      case  3:
      	if(sideto==13) {
            trans.Sum()(0,0) =-1.;
            trans.Sum()(1,0) = 1.;
         } else
         if(sideto==16) {
            trans.Sum()(0,0) = 0.;
            trans.Sum()(1,0) = 0.;
         } else
         if(sideto==17) {
            trans.Sum()(0,0) = 1.;
            trans.Sum()(1,0) = 0.;
         }
         break;
      case 4:
         if(sideto==14 || sideto==15 || sideto==16 || sideto==17) {
            trans.Sum()(0,0) = 0.;
            trans.Sum()(1,0) = 1.;
         }
         break;
      default:
	      PZError << "TPZGeoElPi3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
   }//switch
 }//fim lados cantos para faces
 if(sideto > 4 && sideto < 13) {//5 a 12
 	switch(sidefrom) {
   	//canto para lado !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   		   //row = dimto = 1 e col = dimfrom = 0
            //sidefrom canto de dim =0 : t(1,0)
      case  0://fMult(row,col) : fMult(1,0) , fSum(row,1) : fSum(1,1)
         if(sideto== 8)               trans.Sum()(0,0)  =  1.;
         if(sideto== 5 || sideto== 9) trans.Sum()(0,0)  = -1.;
         break;
      case  1:
      	if(sideto== 5)               trans.Sum()(0,0)  =  1.;
         if(sideto== 6 || sideto==10) trans.Sum()(0,0)  = -1.;
      	break;
      case  2:
      	if(sideto== 6)               trans.Sum()(0,0)  =  1.;
         if(sideto== 7 || sideto==11) trans.Sum()(0,0)  = -1.;
      	break;
      case  3:
      	if(sideto== 7)               trans.Sum()(0,0)  =  1.;
         if(sideto== 8 || sideto==12) trans.Sum()(0,0)  = -1.;
      	break;
      case  4:
         if(sideto==9 || sideto==10 || sideto==11 || sideto==12) trans.Sum()(0,0) = 1.;
  	}//switch
 }//if cantos
   return trans;
}

TPZCompEl *TPZGeoElPi3d::CreateBCCompEl(int side,int bc,TPZCompMesh &cmesh) {
	if(side<0 || side>18) return 0;

   if(side==18) {
     std::cout << "TPZGeoElPi3d::CreateBCCompEl with side = 18 not implemented\n";
      return 0;
   }
   int64_t index;
	if(side<5) {
      TPZMaterial *bcptr = cmesh.FindMaterial(bc);
      if(!bcptr) {
      PZError << "TPZGeoElPi3d::CreateBCCompEl has no bc.\n";
      return 0;
      }
      TPZCompEl *cel = Reference();
      if(!cel) {
      PZError << "TPZGeoElPi3d::CreateBCCompEl has no computational element\n";
      return 0;
      }
	  TPZManVector<int64_t> nodeindexes(1);
      TPZGeoElPoint *gel;
      nodeindexes[0] = fNodeIndexes[side];
      gel = new TPZGeoElPoint(nodeindexes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,side));
      return gel->CreateCompEl(cmesh,index);
   } else if (side > 4 && side < 13) {//side =5 a 12 : lados
      TPZManVector<int64_t> nodes(2);
      int s = side-5;
      nodes[0] = NodeIndex(TPZCompElPi3d::SideNodes[s][0]);
      nodes[1] = NodeIndex(TPZCompElPi3d::SideNodes[s][1]);
      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::SideNodes[s][0]));
      TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::SideNodes[s][1]));
      TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(this,side));
      return gel->CreateCompEl(cmesh,index);
   } else if (side > 12) {//side = 13 a 17 : faces
      TPZManVector<int64_t> nodes(4);//4o = -1 para face triangular
      int s = side-13;
      nodes[0] = NodeIndex(TPZCompElPi3d::FaceNodes[s][0]);
      nodes[1] = NodeIndex(TPZCompElPi3d::FaceNodes[s][1]);
      nodes[2] = NodeIndex(TPZCompElPi3d::FaceNodes[s][2]);
      nodes[3] = NodeIndex(TPZCompElPi3d::FaceNodes[s][3]);
      TPZGeoElT2d *gelt;
      TPZGeoElQ2d *gelq;
      if(side==13) {
      	gelq = new TPZGeoElQ2d(nodes,bc,*Mesh());
         TPZGeoElSide(gelq,0).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceNodes[s][0]));
         TPZGeoElSide(gelq,1).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceNodes[s][1]));
         TPZGeoElSide(gelq,2).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceNodes[s][2]));
         TPZGeoElSide(gelq,3).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceNodes[s][3]));
         TPZGeoElSide(gelq,4).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceSides[s][0]));
         TPZGeoElSide(gelq,5).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceSides[s][1]));
         TPZGeoElSide(gelq,6).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceSides[s][2]));
         TPZGeoElSide(gelq,7).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceSides[s][3]));
         TPZGeoElSide(gelq,8).SetConnectivity(TPZGeoElSide(this,side));
         return gelq->CreateCompEl(cmesh,index);
      } else {
         nodes.Resize(3);
	      gelt = new TPZGeoElT2d(nodes,bc,*Mesh());
         TPZGeoElSide(gelt,0).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceNodes[s][0]));
         TPZGeoElSide(gelt,1).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceNodes[s][1]));
         TPZGeoElSide(gelt,2).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceNodes[s][2]));
         TPZGeoElSide(gelt,3).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceSides[s][0]));
         TPZGeoElSide(gelt,4).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceSides[s][1]));
         TPZGeoElSide(gelt,5).SetConnectivity(TPZGeoElSide(this,TPZCompElPi3d::FaceSides[s][2]));
         TPZGeoElSide(gelt,6).SetConnectivity(TPZGeoElSide(this,side));
         return gelt->CreateCompEl(cmesh,index);
      }
   } else PZError << "TPZGeoElPi3d::CreateBCCompEl. Side = " << side <<std::endl;
   return 0;
}

void TPZGeoElPi3d::NodeFaceIds(TPZVec<int> &ids,int face) {

	ids.Resize(4,-1);
   if((face>-1 && face<5) || (face>12 && face<18)) {
   	if(face>12) face = face-13;
      ids[0] = NodeIndex(TPZCompElPi3d::FaceNodes[face][0]);
      ids[1] = NodeIndex(TPZCompElPi3d::FaceNodes[face][1]);
      ids[2] = NodeIndex(TPZCompElPi3d::FaceNodes[face][2]);
		ids[3] = NodeIndex(TPZCompElPi3d::FaceNodes[face][3]);
      if(face > 0) ids.Resize(3);
      return;
   }
   std::cout << "TPZCompElPi3d::NodeFaceIds bad side , side = " << face <<std::endl;
}


//cada lado do filho {s0,s1,s2,...,s19} está contido num lado do pai: sn=fatside

static int fatherside[10][19] = {
/*00*/{0,5,13,8,9,5,13,13,8,9,14,18,17,13,14,18,18,17,18},
/*01*/{5,1,6,13,10,5,6,13,13,14,10,15,18,13,14,15,18,18,18},
/*02*/{13,6,2,7,11,13,6,7,13,18,15,11,16,13,18,15,16,18,18},
/*03*/{8,13,7,3,12,13,13,7,8,17,18,16,12,13,18,18,16,17,18},
/*04*/{9,10,11,12,4,14,15,16,17,9,10,11,12,18,14,15,16,17,18},
/*05*/{10,9,12,11,13,14,17,16,15,18,18,18,18,18,18,18,18,18,18},
/*06*/{9,5,13,10,14,13,18,14,14,18,18,14,18,18,18,-1,-1,-1,-1},
/*07*/{6,10,11,13,15,15,15,13,18,18,15,18,18,18,18,-1,-1,-1,-1},
/*08*/{12,13,7,11,18,13,16,16,18,16,18,18,18,16,18,-1,-1,-1,-1},
/*09*/{13,9,12,8,18,17,18,13,17,17,18,18,17,18,18,-1,-1,-1,-1},
};
static int fatherside2[2][19] = {//pirâmides filhos do tetraedro
/*04*/{4,8,9,6,7,11,12,13,10,11,11,13,13,14,11,14,13,14,14},
/*05*/{8,4,6,9,5,11,10,13,12,12,10,10,12,14,14,10,14,12,14} };

TPZGeoElSide TPZGeoElPi3d::Father2(int side){//Augusto:09/01/01

	if(side<0 || side>18 || Father(side).Element()==0){
		PZError << "TPZGeoElPi3d::Father2 called error" <<std::endl;
        return TPZGeoElSide();
	}
	int subelindex = WhichSubel();
	if(fatherside[subelindex][side]<0){
		PZError << "TPZGeoElPi3d::Father2 called with index error\n";
		return TPZGeoElSide();
	}
	if(Father(side).Element()->NSides()==15){//pai tetraedro
  	return TPZGeoElSide(Father(side).Element(),fatherside2[subelindex-4][side]);
  } else return TPZGeoElSide(Father(side).Element(),fatherside[subelindex][side]);
}

static int subeldata[19][27][2] = {//CASO DIFERENTE TAMANHO
/*00*/{{0,0}},
/*01*/{{1,1}},
/*02*/{{2,2}},
/*03*/{{3,3}},
/*04*/{{4,4}},
/*05*/{{0,5},{0,1},{1,5}},
/*06*/{{1,6},{1,2},{2,6}},
/*07*/{{2,7},{2,3},{3,7}},
/*08*/{{0,8},{0,3},{3,8}},
/*09*/{{0,9},{0,4},{4,9}},
/*10*/{{1,10},{1,4},{4,10}},
/*11*/{{2,11},{2,4},{4,11}},
/*12*/{{3,12},{3,4},{4,12}},
/*13*/{{0,13},{1,13},{2,13},{3,13},{0,2},{0,6},{0,7},{1,7},{2,8}},
/*14*/{{0,14},{1,14},{4,14},{6,11},{0,10},{1,9},{4,5}},
/*15*/{{1,15},{2,15},{4,15},{7,10},{1,11},{2,10},{4,6}},
/*16*/{{2,16},{3,16},{4,16},{8,13},{2,12},{3,11},{4,7}},
/*17*/{{0,17},{3,17},{4,17},{9,12},{0,12},{3,9},{4,8}},
/*18*/{{0,18},{1,18},{2,18},{3,18},{4,18},{5,18},{6,14},{7,14},{8,14},
       {9,14},{0,11},{1,12},{2,9},{3,10},{0,15},{0,16},{1,16},{1,17},
       {2,14},{2,17},{3,14},{3,15},{4,13},{6,12},{7,13},{8,11},{9,10}}
};

static int nsubeldata[19] = {1,1,1,1,1,3,3,3,3,3,3,3,3,9,7,7,7,7,27};


void TPZGeoElPi3d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){

	 subel.Resize(0);
   if(side<0 || side>18 || !HasSubElement(side)){
      PZError << "TPZGeoElPi3d::GetSubElements2 called with error arguments\n";
      return;
   }
   int nsub = nsubeldata[side];
   for(int i=0;i<nsub;i++)
       subel.Push(TPZGeoElSide(fSubEl[subeldata[side][i][0]],
                                      subeldata[side][i][1]));
}

static REAL buildt[10][19][4][3] = {//por colunas
/*S0*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{-1,-1,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*06*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*07*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*10*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*11*/{{-0.25,-0.25,0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,-0.25,0.25}},
      /*12*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*13*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*14*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*15*/{{0.,1,0.},{-0.5,0.5,0.5},{0.,0.,0.},{0.,-1,0.}},
      /*16*/{{1,0.,0.},{0.5,-0.5,0.5},{0.,0.,0.},{-1,0.,0.}},
      /*17*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{-.5,-.5,0.}}},
/*S1*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1,-1,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*07*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*08*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*09*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*11*/{{-.25,.25,0.},{0.,0.,0.},{0.,0.,0.},{.25,.25,0.}},
      /*12*/{{0.25,-0.25,0.25},{0.,0.,0.},{0.,0.,0.},{0.25,-0.25,0.25}},
      /*13*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*14*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{1,0.,0.},{0.5,-0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,1,0.},{0.5,0.5,0.5},{0.,0.,0.},{0.,-1,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,-.5,0.}}},
/*S2*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1,1,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*08*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*09*/{{0.25,0.25,0.25},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.25}},
      /*10*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*12*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*13*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*14*/{{1,0.,0.},{0.5,0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*16*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*17*/{{0.,1,0.},{0.5,0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,.5,0.}}},
/*S3*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{-1,1,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*06*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*09*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*10*/{{-0.25,0.25,0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,0.25,0.25}},
      /*11*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*13*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*14*/{{1,0.,0.},{0.5,0.5,0.5},{0.,0.,0.},{-1,0.,0.}},
      /*15*/{{0.,1,0.},{-0.5,0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{-.5,.5,0.}}},
/*S4*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,1}},
      /*05*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*08*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*13*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.5}},
      /*14*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*16*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*17*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,0.,.5}}},
/*S5*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{-0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*08*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*09*/{{-0.25,0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{0.25,-0.25,0.25}},
      /*10*/{{0.25,0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,-0.25,0.25}},
      /*11*/{{0.25,-0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,0.25,0.25}},
      /*12*/{{-0.25,-0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.25}},
      /*13*/{{-0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.5}},
      /*14*/{{-1,0.,0.},{-0.5,0.5,-0.5},{0.,0.,0.},{0.5,-0.5,0.5}},
      /*15*/{{0.,1,0.},{0.5,0.5,-0.5},{0.,0.,0.},{-0.5,-0.5,0.5}},
      /*16*/{{-1,0.,0.},{-0.5,-0.5,-0.5},{0.,0.,0.},{0.5,0.5,0.5}},
      /*17*/{{0.,1,0.},{-0.5,0.5,-0.5},{0.,0.,0.},{0.5,-0.5,0.5}},
      /*18*/{{-.5,0.,0.},{0.,.5,0.},{0.,0.,-.5},{0.,0.,.5}}},
/*S6*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.25,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*05*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*06*/{{0.25,0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,-0.25,0.25}},
      /*07*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*08*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*09*/{{0.25,-0.25,0.25},{0.,0.,0.},{0.,0.,0.},{0.25,-0.25,0.25}},
      /*10*/{{0.5,-0.5,-0.5},{0.5,0.5,-0.5},{0.,0.,0.},{-0.5,-0.5,0.5}},
      /*11*/{{0.5,-0.5,0.},{0.5,0.,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*12*/{{0.,1,0.},{0.5,0.5,0.5},{0.,0.,0.},{0.,-1,0.}},
      /*13*/{{0.5,0.5,-0.5},{1,0.,0.},{0.,0.,0.},{-0.5,-0.5,0.5}},
      /*14*/{{.5,-.5,-.5},{.5,.5,-.5},{1.,.0,.0},{-.5,-.5,.5}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
/*S7*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*05*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*06*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*08*/{{-0.25,0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{0.25,-0.25,0.25}},
      /*09*/{{-0.25,-0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.25}},
      /*10*/{{-0.5,0.5,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*11*/{{-0.5,-0.5,0.5},{-1,0.,0.},{0.,0.,0.},{1,0.,0.}},
      /*12*/{{0.,1,0.},{-0.5,0.5,-0.5},{0.,0.,0.},{0.5,-0.5,0.5}},
      /*13*/{{-0.5,0.5,0.5},{-1,0.,0.},{0.,0.,0.},{1,0.,0.}},
      /*14*/{{-.5,-.5,.5},{-.5,.5,.5},{-1.,0.,0.},{1.,0.,0.}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
/*S8*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.25,-0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,0.25,0.25}},
      /*05*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*06*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*07*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*08*/{{0.25,0.25,0.25},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.25}},
      /*09*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*10*/{{0.5,-0.5,-0.5},{0.5,0.5,-0.5},{0.,0.,0.},{-0.5,0.5,0.5}},
      /*11*/{{0.5,-0.5,-0.5},{1,0.,0.},{0.,0.,0.},{-0.5,0.5,0.5}},
      /*12*/{{0.,1,0.},{0.5,0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*13*/{{0.5,-0.5,0.},{0.5,0.,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*14*/{{.5,-.5,-.5},{.5,.5,-.5},{1.,0.,0.},{-.5,.5,.5}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
/*S9*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{-0.25,-0.25,0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,-0.25,0.25}},
      /*05*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*06*/{{0.25,-0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,0.25,0.25}},
      /*07*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*08*/{{0.25,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*09*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*10*/{{-0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*11*/{{-0.5,-0.5,0.5},{-1,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*12*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*13*/{{-0.5,0.5,0.5},{-1,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*14*/{{-.5,-.5,.5},{-.5,.5,.5},{-1.,0.,0.},{0.,0.,0.}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}}
};

TPZTransform<STATE> TPZGeoElPi3d::BuildTransform2(int side, TPZGeoEl *father){//Augusto:09/01/01

	if(side<0 || side>18 || !father){
  	PZError << "TPZGeoElPi3d::BuildTransform2 side out of range or father null\n";
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


REAL TPZGeoElPi3d::MidSideNode[19][3] = {
/*00*/{-1.,-1.},   /*01*/{1.,-1.},   /*02*/{1.,1.},/*03*/{-1.,1.},/*04*/{0.,0.,1.},
/*05*/{ 0.,-1.},   /*06*/{1., 0.},   /*07*/{0.,1.},/*08*/{-1.,0.},
/*09*/{-.5,-.5,.5},/*10*/{.5,-.5,.5},/*11*/{.5,.5,.5},/*12*/{-.5,.5,.5},
/*13*/{0.,  0. ,  0. },/*14*/{  0.  ,-2./3.,1./3.},/*15*/{2./3.,0.,1./3.},
/*16*/{0.,2./3.,1./3.},/*17*/{-2./3.,  0.  ,1./3.},/*18*/{  0. ,0.,1./5.} };


int TPZGeoElPi3d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3),pss(3),pf(3),pfs(3);
                    //point son, point side son, point father, point side father : elemento mestre
  for(sn=0;sn<10;sn++){
    TPZGeoEl *son = subs[sn];
    int nsides = son->NSides();
    for(sd=0;sd<nsides;sd++){
      TPZTransform<STATE> telsd(0,0);
      if(nsides==15){
        ps[0] = TPZGeoElT3d::MidSideNode[sd][0];//tetraedro
        ps[1] = TPZGeoElT3d::MidSideNode[sd][1];
        ps[2] = TPZGeoElT3d::MidSideNode[sd][2];
        telsd = pzshape::TPZShapeTetra::TransformElementToSide(sd);
      } else if(nsides==19){
        ps[0] = MidSideNode[sd][0];//pirâmide
        ps[1] = MidSideNode[sd][1];
        ps[2] = MidSideNode[sd][2];
        telsd = pzshape::TPZShapePiram::TransformElementToSide(sd);
      }
      if(son->WhichSide(ps) != sd)std::cout << "Lado nao bate\n";
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform<STATE> t;
	  gel->BuildTransform2(sd, son,t);//para não furar pirâmide com pai tetraedro
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = pzshape::TPZShapePiram::TransformSideToElement(sdfat); //else
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



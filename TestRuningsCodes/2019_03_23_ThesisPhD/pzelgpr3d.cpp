//METHODS DEFINITION FOR CLASS ELEMPR3D
#include "pzelgpr3d.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelc1d.h"
#include "pzelgt2d.h"
#include "pzelgq2d.h"
#include "pzelcpr3d.h"
#include "pzshapeprism.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgnode.h"
#include "pzshtmat.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include <stdlib.h>

static TPZCompEl *CreateEl(TPZGeoElPr3d *gel,TPZCompMesh &mesh, int64_t &index) {
  return new TPZCompElPr3d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElPr3d::fp)(TPZGeoElPr3d *,TPZCompMesh &, int64_t &) = CreateEl;

TPZGeoElPr3d::TPZGeoElPr3d(int id,TPZVec<int64_t> &nodeindexes,int matid,TPZGeoMesh &mesh):
	TPZGeoElRefLess<pzgeom::TPZGeoPrism>(id,nodeindexes,matid,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=6) {
    PZError << "TPZGeoElPr3d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<6;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<8;i++) fSubEl[i] = 0;
}

TPZGeoElPr3d::TPZGeoElPr3d(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) :
  TPZGeoElRefLess<pzgeom::TPZGeoPrism>(nodeindexes,matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod!=6) {
    PZError << "TPZGeoElPr3d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<6;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<8;i++) fSubEl[i] = 0;
}

TPZGeoElPr3d::~TPZGeoElPr3d() {}

void TPZGeoElPr3d::Shape(TPZVec<REAL> &pt,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {

  phi(0,0)  = .5*(1.-pt[0]-pt[1])*(1.-pt[2]);
  phi(1,0)  = .5*pt[0]*(1.-pt[2]);
  phi(2,0)  = .5*pt[1]*(1.-pt[2]);
  phi(3,0)  = .5*(1.-pt[0]-pt[1])*(1.+pt[2]);
  phi(4,0)  = .5*pt[0]*(1.+pt[2]);
  phi(5,0)  = .5*pt[1]*(1.+pt[2]);

  dphi(0,0) = -.5*(1.-pt[2]);
  dphi(1,0) = -.5*(1.-pt[2]);
  dphi(2,0) = -.5*(1.-pt[0]-pt[1]);

  dphi(0,1) =  .5*(1.-pt[2]);
  dphi(1,1) =  .0;
  dphi(2,1) = -.5*pt[0];

  dphi(0,2) =  .0;
  dphi(1,2) =  .5*(1.-pt[2]);
  dphi(2,2) = -.5*pt[1];

  dphi(0,3) = -.5*(1.+pt[2]);
  dphi(1,3) = -.5*(1.+pt[2]);
  dphi(2,3) =  .5*(1.-pt[0]-pt[1]);

  dphi(0,4) =  .5*(1.+pt[2]);
  dphi(1,4) =  .0;
  dphi(2,4) =  .5*pt[0];

  dphi(0,5) =  .0;
  dphi(1,5) =  .5*(1.+pt[2]);
  dphi(2,5) =  .5*pt[1];
}

TPZGeoElPr3d *TPZGeoElPr3d::CreateGeoEl(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElPr3d(nodeindexes,matind,mesh);
}

int TPZGeoElPr3d::NNodes() {
  return 6;
}
int64_t TPZGeoElPr3d::NodeIndex(int node) {
  if(node<0 || node>6) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElPr3d::NSideNodes(int side) {
  if(side<0 || side>20) {
    PZError << "TPZGeoElPr3d::NSideNodes. Bad parameter side.\n";
    return 0;
  }
  if(side<6) return 1;//cantos
  if(side>5 && side<15) return 2;//lados
  if(side==15 || side==19) return 3;//faceS triangulares
  if(side>15 && side<19) return 4;//faces  quadrilaterais
  return 6;//centro
}

int64_t TPZGeoElPr3d::SideNodeIndex(int side,int node) {
  if(side<0 || side>20 || node<0) {//21 sides
    PZError << "TPZGeoElPr3d::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  //0<=side<=20 e node>=0
  //side central
  if(side==20 && node<6) return fNodeIndexes[node];//um dos 5 cantos
  //side de canto
  if(side < 6) return fNodeIndexes[side];//canto side, deve ser node = 0
  //sides 6 a 14
  if(side>5 && side<15 && node<2) {//arestas 6 a 14
      side-=6;
      return fNodeIndexes[TPZCompElPr3d::SideNodes[side][node]];
  } else if((side == 15 || side == 19) && node<3) {//faces 15 e 19
  		side -= 15;
	  	return fNodeIndexes[TPZCompElPr3d::FaceNodes[side][node]];
  } else if((side>15 && side<19) && node<4) {//faces 16 a 18
  		side-=15;
	  	return fNodeIndexes[TPZCompElPr3d::FaceNodes[side][node]];
  }
  return -1;
}

void TPZGeoElPr3d::MidSideNodeIndex(int side, int64_t &index) {
  index = -1;
  if(side<0 || side>20) {
    PZError << "TPZGeoElPr3d::MidSideNodeIndex. Bad parameter side = " << side <<std::endl;
    return;
  }
  //sides 0 a 5
	if(side<6) {//o nó medio do lado 0 é o 0 etc.
		index=fNodeIndexes[side];
		return;
	}
   //o nó medio da face é o centro da face e o nó medio do centro é o centro
   //como nó de algum filho se este existir
   //caso tenha filhos é o canto de algum filho, se não tiver filhos retorna -1
	if(HasSubElement(0)) {
   	int s = side - 6;//arestas e faces
      if(side==15) return;
      if(side>15) s = s-1;
      int k = TPZCompElPr3d::MidSideNodes[s][0];
      int l = TPZCompElPr3d::MidSideNodes[s][1];
      index=((TPZGeoElPr3d *) fSubEl[k])->fNodeIndexes[l];
	   //index=((TPZGeoElPr3d *) fSubEl[])->fNodeIndexes[TPZCompElPr3d::MidSideNodes[side][1]];
	}
}

void TPZGeoElPr3d::NewMidSideNode(int side, int64_t &index) {
  MidSideNodeIndex(side,index);
  if(index < 0) {
    int k = 0;
    if(side>14) k = 1;
    TPZGeoElSide gelside = Neighbour(side+k);
    while(gelside.Element()) {
      gelside.Element()->MidSideNodeIndex(gelside.Side(),index);
      if(index!=-1) return;
      gelside = gelside.Neighbour();
      if(gelside.Element()==this) break;
    }
    TPZVec<REAL> par(3,0.);
    TPZVec<REAL> coord(3,0.);
    if(side < 6) {index = -1; return;}
    //aqui side = 6 a 19
    side-=6;//0 a 14 : medios das arestas e faces
    par[0] = TPZCompElPr3d::MidCoord[side][0];
    par[1] = TPZCompElPr3d::MidCoord[side][1];
    par[2] = TPZCompElPr3d::MidCoord[side][2];
    X(par,coord);
    index = Mesh()->NodeVec().AllocateNewElement();
    Mesh()->NodeVec()[index].Initialize(coord,*Mesh());
  }
}

int TPZGeoElPr3d::SideDimension(int side) {

	if (side<0 || side>20) {
   	PZError << "TPZGeoElPr3d::SideDimension called with side " << side <<std::endl;
      return 0;
   }
  if(side<6)  return 0;//cantos
  if(side<15) return 1;//arestas 6 a 14
  if(side<20) return 2;//faces 15 a 19
  return 3;//centro

   }

TPZGeoElSide TPZGeoElPr3d::HigherDimensionSides(int side,int targetdimension) {
//targetdimension deve ser 1 , 2 ou 3
//se side =  0 a 4  targetdimension deve ser 1
//se side =  5 a 12  targetdimension deve ser 2
//se side = 18      targetdimension deve ser 3
  if( (side<0 || side>20) || (targetdimension < 1 || targetdimension > 3) ) {
     PZError << "TPZGeoElPr3d::HigherDimensionSides called with side = " << side
	          << " targetdimension = " << targetdimension <<std::endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
  int bestface;
  //side = 0 a 18
  switch(targetdimension) {//=1,2
	  case 1:
     	 if(this == father->SubElement(0)) {
       	 if(side==1) return TPZGeoElSide(this,6);
       	 if(side==2) return TPZGeoElSide(this,8);
          if(side==3) return TPZGeoElSide(this,9);
       } else if(this == father->SubElement(1)) {
       	 if(side==0) return TPZGeoElSide(this,6);
       	 if(side==2) return TPZGeoElSide(this,7);
          if(side==4) return TPZGeoElSide(this,10);
       } else if(this == father->SubElement(2)) {
       	 if(side==0) return TPZGeoElSide(this,8);
       	 if(side==1) return TPZGeoElSide(this,7);
          if(side==5) return TPZGeoElSide(this,11);
       } else if(this == father->SubElement(3)) {
       	 if(side==3) return TPZGeoElSide(father->SubElement(0),8);
       	 if(side==5) return TPZGeoElSide(father->SubElement(0),6);
       	 if(side==4) return TPZGeoElSide(father->SubElement(1),7);
       } else if(this == father->SubElement(4)) {
       	 if(side==0) return TPZGeoElSide(this,9);
       	 if(side==4) return TPZGeoElSide(this,12);
          if(side==5) return TPZGeoElSide(this,14);
       } else if(this == father->SubElement(5)) {
       	 if(side==1) return TPZGeoElSide(this,10);
          if(side==3) return TPZGeoElSide(this,12);
       	 if(side==5) return TPZGeoElSide(this,13);
       } else if(this == father->SubElement(6)) {
       	 if(side==2) return TPZGeoElSide(this,11);
       	 if(side==3) return TPZGeoElSide(this,14);
          if(side==4) return TPZGeoElSide(this,13);
       } else if(this == father->SubElement(7)) {
       	 if(side==0) return TPZGeoElSide(father->SubElement(4),14);
       	 if(side==1) return TPZGeoElSide(father->SubElement(5),13);
       	 if(side==2) return TPZGeoElSide(father->SubElement(4),12);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
  	  case 2:
     	 if(this == father->SubElement(0)) {

          bestface = BestDimensionSideOfTwoFaces(15,16);
          if((side==1 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(15,18);
          if((side==2 || side==8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,18);
          if((side==3 || side==9) && bestface) return TPZGeoElSide(this,bestface);
          if(side==7) return TPZGeoElSide(this,15);
       	 if(side==4 || side==10 || side==12) return TPZGeoElSide(this,16);
          if(side==5 || side==11 || side==14) return TPZGeoElSide(this,18);
       } else if(this == father->SubElement(1)) {

          bestface = BestDimensionSideOfTwoFaces(15,16);
          if((side==0 || side==6) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(15,17);
          if((side==2 || side==7) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,17);
          if((side==4 || side==10) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==8) return TPZGeoElSide(this,15);
          if(side==3 || side== 9 || side==12) return TPZGeoElSide(this,16);
          if(side==5 || side==11 || side==13) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(2)) {

          bestface = BestDimensionSideOfTwoFaces(15,18);
          if((side==0 || side==8) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(15,17);
          if((side==1 || side==7) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(17,18);
          if((side==5 || side==11) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==6) return TPZGeoElSide(this,15);
          if(side==3 || side== 9 || side==14) return TPZGeoElSide(this,18);
          if(side==4 || side==10 || side==13) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(3)) {

          bestface = BestDimensionSideOfTwoFaces(15,18);
          if(side==3 && bestface) return TPZGeoElSide(father->SubElement(0),bestface);
          bestface = BestDimensionSideOfTwoFaces(15,17);
          if(side==4 && bestface) return TPZGeoElSide(father->SubElement(1),bestface);
          bestface = BestDimensionSideOfTwoFaces(15,16);
          if(side==5 && bestface) return TPZGeoElSide(father->SubElement(0),bestface);
          if(side==12 || side==13 || side==14) return TPZGeoElSide(this,19);
          if(side== 0 || side== 9) return TPZGeoElSide(father->SubElement(6),18);//sub 6 dado que em ambos o canto é o zero
          if(side== 1 || side==10) return TPZGeoElSide(father->SubElement(5),17);//o canto é o mesmo
          //if(side== 2 || side==11) return TPZGeoElSide(father->SubElement(1),16);
          //como nao existe um irmao com o mesmo canto => nao existe HigherDimensioSide para este caso
       } else if(this == father->SubElement(4)) {

          bestface = BestDimensionSideOfTwoFaces(16,18);
          if((side==0 || side==9) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,19);
          if((side==4 || side==12) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(18,19);
          if((side==5 || side==14) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==13) return TPZGeoElSide(this,19);
          if(side==1 || side== 6 || side==10) return TPZGeoElSide(this,16);
          if(side==2 || side== 8 || side==11) return TPZGeoElSide(this,18);
       } else if(this == father->SubElement(5)) {

          bestface = BestDimensionSideOfTwoFaces(16,17);
          if((side==1 || side==10) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(16,19);
          if((side==3 || side==12) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(17,19);
          if((side==5 || side==13) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==14) return TPZGeoElSide(this,19);
          if(side==0 || side== 6 || side== 9) return TPZGeoElSide(this,16);
          if(side==2 || side== 7 || side==11) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(6)) {

          bestface = BestDimensionSideOfTwoFaces(17,18);
          if((side==2 || side==11) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(18,19);
          if((side==3 || side==14) && bestface) return TPZGeoElSide(this,bestface);
          bestface = BestDimensionSideOfTwoFaces(17,19);
          if((side==4 || side==13) && bestface) return TPZGeoElSide(this,bestface);
       	 if(side==12) return TPZGeoElSide(this,19);
          if(side==0 || side== 8 || side== 9) return TPZGeoElSide(this,18);
          if(side==1 || side== 7 || side==10) return TPZGeoElSide(this,17);
       } else if(this == father->SubElement(7)) {

          bestface = BestDimensionSideOfTwoFaces(18,19);
          if(side==0 && bestface) return TPZGeoElSide(father->SubElement(4),bestface);
          bestface = BestDimensionSideOfTwoFaces(17,19);
          if(side==1 && bestface) return TPZGeoElSide(father->SubElement(5),bestface);
          bestface = BestDimensionSideOfTwoFaces(16,19);
          if(side==2 && bestface) return TPZGeoElSide(father->SubElement(4),bestface);
          if(side==6 || side==7 || side==8) return TPZGeoElSide(this,15);
          if(side== 3 || side== 9) return TPZGeoElSide(father->SubElement(2),18);//o canto é o mesmo
          //if(side== 5 || side==11) return TPZGeoElSide(father->SubElement(4),16);//para este nao existe um irmao com o mesmo canto logo a transformacao deve ser consertada para este caso
          if(side== 4 || side==10) return TPZGeoElSide(father->SubElement(1),17);//o canto é o mesmo
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
     case 3:
       return TPZGeoElSide(this,20);//0<=side<=20
  }//switch
  return TPZGeoElSide();
}

int TPZGeoElPr3d::BestDimensionSideOfTwoFaces(int face1,int face2) {

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

//feito
void TPZGeoElPr3d::LowerDimensionSides(int side,TPZStack<TPZGeoElSide> &smallsides) {
	if (side < 6) return;
   int i;
   if(side < 15) {//side = 6 a 14 : entram os cantos das arestas
   	int s = side-6;
   	smallsides.Push(TPZGeoElSide(this,TPZCompElPr3d::SideNodes[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElPr3d::SideNodes[s][1]));
   } else if(side < 20) {//entram cantos e arestas da face
   	int s = side-15;
   	smallsides.Push(TPZGeoElSide(this,TPZCompElPr3d::FaceNodes[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElPr3d::FaceNodes[s][1]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElPr3d::FaceNodes[s][2]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElPr3d::FaceSides[s][0]));
  	   smallsides.Push(TPZGeoElSide(this,TPZCompElPr3d::FaceSides[s][1]));
   	smallsides.Push(TPZGeoElSide(this,TPZCompElPr3d::FaceSides[s][2]));
      if(s>0 && s<4) {//face quadrilateral
         smallsides.Push(TPZGeoElSide(this,TPZCompElPr3d::FaceNodes[s][3]));
         smallsides.Push(TPZGeoElSide(this,TPZCompElPr3d::FaceSides[s][3]));
      }
      smallsides.Push(TPZGeoElSide(this,side));
   } else if(side==20) {//entram todos os cantos, arestas e faces
   	for (i=0;i<20;i++) smallsides.Push(TPZGeoElSide(this,i));
   }
}

void TPZGeoElPr3d::SideMasterCo(int /*side*/,TPZFMatrix<STATE> &/*coord*/) {
}

//feito
void TPZGeoElPr3d::Jacobian(TPZVec<REAL> &param,TPZFMatrix<STATE> &jacobian,TPZFMatrix<STATE> &axes,REAL &detjac,TPZFMatrix<STATE> &jacinv){

#ifdef DEBUG
  int nnodes = NNodes();
  if (nnodes != 6) {
    PZError << "TPZGeoElPr3d.jacobian only implemented for"
      " 6 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 3 || param[0] < 0. || param[0] > 1. ||
     param[1] < 0. || param[1] > 1. || param[2] < 0. || param[2] > 1.) {
    PZError << "TPZGeoElPr3d.jacobian. param out of range : "
      " param.NElements() = " << param.NElements() <<
      "\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
    return;
  }
#endif
  REAL spacephi[6];
  TPZFMatrix<STATE> phi(6,1,spacephi,6);
  REAL spacedphi[18];
  TPZFMatrix<STATE> dphi(3,6,spacedphi,18);
  Shape(param,phi,dphi);
  jacobian.Zero();
  TPZGeoNode *np;
  int i,j;
  for(i=0;i<6;i++) {
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

void TPZGeoElPr3d::X(TPZVec<REAL> & loc,TPZVec<REAL> &result){
  REAL spacephi[6],spacedphi[18];
  int i,j;
  TPZFMatrix<STATE> phi(6,1,spacephi,5);
  TPZFMatrix<STATE> dphi(3,6,spacedphi,18);
  Shape(loc,phi,dphi);
  for(j=0;j<3;j++) {
    result[j] = 0.0;
    for(i=0;i<6;i++) result[j] += NodePtr(i)->Coord(j)*phi(i,0);
  }
}
/**It's necessary to define the normal vector to side 4, that is the orthogonal
   vector to the surface*/

void TPZGeoElPr3d::NormalVector(int /*side*/,TPZVec<REAL> &/*param*/,TPZVec<REAL> &/*normal*/,
			     TPZFMatrix<STATE> &/*axes*/,TPZFMatrix <STATE> &/*jac1d*/) {
	std::cout << "TPZGeoElPr3d::NormalVector Not implemented yet.\n";
}

/** TO SUBDIVISION
********************************************************************************
  Into Divides is necesary to consider the connectivity with the all neighboards*/
void TPZGeoElPr3d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {

   int i;
   if(HasSubElement(0)) {
      SubElVec.Resize(8);
      for(i=0;i<8;i++) SubElVec[i] = fSubEl[i];
      return;//If exist fSubEl return this sons
   }
   int j, sub, matid = MaterialId();
   int64_t index;
   int64_t np[18];//guarda conectividades dos 6 subelementos
   //numero de nós medios dos lados : arestas e faces
   for(j=0;j<6;j++) np[j]=NodeIndex(j);
   for(sub=6;sub<18;sub++) {
      NewMidSideNode(sub,index);
      np[sub] = index;
   }
   // creating new subelements
   for(i=0;i<8;i++) {
	   TPZManVector<int64_t> cornerindexes(6);
   	for(int j=0;j<6;j++) cornerindexes[j] = np[TPZCompElPr3d::CornerSons[i][j]];
      fSubEl[i] = CreateGeoEl(cornerindexes,matid,*Mesh());
   }
   SubElVec.Resize(8);
   for(sub=0;sub<8;sub++) {
      SubElVec[sub] = fSubEl[sub];
      SubElVec[sub]->SetFather(this);
   }
   for(i=0;i<8;i++) {//conectividades entre os filhos : viz interna
   	for(j=0;j<19;j++) {
      	int elside = TPZCompElPr3d::InNeigh[i][j][0];
         if(elside == -1) break;
      	fSubEl[i]->SetNeighbour(elside,TPZGeoElSide(fSubEl[TPZCompElPr3d::InNeigh[i][j][1]],TPZCompElPr3d::InNeigh[i][j][2]));
      }                              //lado do subel                                          numero do filho viz.             lado do viz.
   }
   //vizinhança externa ao elemento atual
   //procura-se um viz pela face que seja dividido
   TPZGeoElSide dividedneighbour;
   for(int face=15;face<20;face++) {
      dividedneighbour = Neighbour(face);//BuildConnectivity(..) inicializou
      while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         //achou-se um viz para alguma face
         TPZManVector<int64_t> nodes(4);
         TPZStack<TPZGeoElSide> neighsub;
         int f = face-15;
         for(i=0;i<4;i++) nodes[i] = NodeIndex(TPZCompElPr3d::FaceNodes[f][i]);//nós globais na ordem e sentido local da face atual
         if(f==0 || f==4) nodes.Resize(3);
         dividedneighbour.GetSubElements2(neighsub);//3 ou 4 subs viz na ordem dos subs da face atual
         if(f>0 && f<4) {//faces 16,17 e 18
            nodes.Resize(1);
            //geoside central a face
            int fs0 = TPZCompElPr3d::FaceSons[f][0];//1o filho da face
            int fn2 = TPZCompElPr3d::FaceNodes[f][2];//nó 2 da face
            nodes[0] = fSubEl[fs0]->NodeIndex(fn2);//nó do centro da face como nó do filho nessa face
            int locside = neighsub[0].Element()->WhichSide(nodes);//lado do viz do sub conectado ao nó centro da face
            TPZGeoElSide sub(fSubEl[fs0],fn2);
            //conectividade do nó do centro da face fechando o ciclo entre os filhos dessa face conectados pelo nó
            sub.SetConnectivity(TPZGeoElSide(neighsub[0].Element(),locside));
         }
         //conectividades dos lados dos subelementos interiores a face
         nodes.Resize(2);
         for(i=0; i<4; i++) {//4 subelementos associados à face
            int fsi = TPZCompElPr3d::FaceSons[f][i];
            int inside = TPZCompElPr3d::FaceInRib[f][i];
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
            TPZGeoElSide sub(fSubEl[TPZCompElPr3d::FaceSons[f][i]],face);
            sub.SetConnectivity(neighsub[i]);
         }
      }
   }//fim faces
	//procura-se um viz pela aresta do atual, que seja dividido
   for(int side=6; side<15; side++) {//(side 6 : subs 0,1),(side 6 : subs 1,2) ...
   	int s = side-6;
   	dividedneighbour = Neighbour(side);
      while(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
         //int i = TPZCompElPr3d::SideNodes[s][0];//nó local 0 do lado side
         int iplus = TPZCompElPr3d::SideNodes[s][1];//nó local 1 do lado
         TPZManVector<int64_t> nodes(2);
         TPZStack<TPZGeoElSide> neighsub;
         nodes[0] = SideNodeIndex(side,0);//nó global 0 do lado side
         nodes[1] = SideNodeIndex(side,1);//nó global 1
         dividedneighbour.GetSubElements2(neighsub);//filhos do viz conectados nesses lados
         //TPZGeoElSide sub(fSubEl[i],side);
         TPZGeoElSide sub(fSubEl[TPZCompElPr3d::RibSons[s][0]],side);
         //conectividade pelo lado side entre os filhos
         sub.SetConnectivity(neighsub[0]);
         //sub = TPZGeoElSide(fSubEl[iplus],side);
         sub = TPZGeoElSide(fSubEl[TPZCompElPr3d::RibSons[s][1]],side);
         sub.SetConnectivity(neighsub[1]);
         nodes.Resize(1);
         //nodes[0] = fSubEl[i]->SideNodeIndex(iplus,0);//nó do medio do lado como nó de um dos filhos
         nodes[0] = fSubEl[TPZCompElPr3d::MidSideNodes[s][0]]->SideNodeIndex(TPZCompElPr3d::MidSideNodes[s][1],0);
         int locside = neighsub[0].Element()->WhichSide(nodes);
         //sub = TPZGeoElSide(fSubEl[i],iplus);
         sub = TPZGeoElSide(fSubEl[TPZCompElPr3d::RibSons[s][0]],iplus);
         //conectividade para o nó do medio do lado side do atual como nó de um dos filhos
         sub.SetConnectivity(TPZGeoElSide(neighsub[0].Element(),locside));
      }
   }
   //procura-se um viz pelo canto do atual, que seja dividido
   for(int corner=0; corner<6; corner++) {//cantos
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
         int k = 0;
         if(corner>2) k = 1;
         TPZGeoElSide sub(fSubEl[corner+k],corner);
         sub.SetConnectivity(neighsub[0]);
      }
   }
}


int TPZGeoElPr3d::NSubElements() {
  return 8;
}

int TPZGeoElPr3d::NSideSubElements(int side) {
  if(side < 0 || side > 20) {
    PZError << "TPZGeoElPr3d::NSideSubElements called for side " << side <<std::endl;
    return 0;
  }
  if(side==20) return 8;//centro
  if(side>14 && side<20) return 4;//faces
  if(side>5) return 2;//lados
  return 1;//cantos
}

TPZGeoElSide TPZGeoElPr3d::SideSubElement(int side,int position) {
   if (position<0 || position>8 || side <0 ||side>20) {
   	PZError << "TPZGeoElPr3d::SideSubElement called with position " << position << " side " << side <<std::endl;
      return TPZGeoElSide();
   }
   if(side==20) {
      return TPZGeoElSide(SubElement(position),20);
   }
   if(side<6) {//cantos
      if(position!=0) {
         PZError << "TPZGeoElPr3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
         return TPZGeoElSide(SubElement(side),side);
      }
   }
   if(side>5 && side<15) {//lados
       if(position!=0 && position!=1) {
         PZError << "TPZGeoElPr3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
      	int s = side-6;
         return TPZGeoElSide(SubElement(TPZCompElPr3d::RibSons[s][position]),side);
      }
   }
   if(side>14) {//faces
       if(position<0 || position>4) {//position!=0 && position!=1 && position!=2 && position!=3
         PZError << "TPZGeoElPr3d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
      	int s = side-15;//s = 0,1,2,3,4
         int k = 0;
         if(side==15 && position==3) k =  4;
         if(side==19 && position==3) k = -4;
         return TPZGeoElSide(SubElement(TPZCompElPr3d::FaceSons[s][position]),side+k);
      }
   }
   return TPZGeoElSide();
}

void TPZGeoElPr3d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
   if(!fSubEl[0]) {
      sub.Resize(0);
      return;
   }
   if(side < 0 || side > 20) {
      PZError << "TPZGeoElPr3d::SideSubElements called for side " << side <<std::endl;
      return;
   }
   if(side==20) {
      sub.Resize(10);
      for(int i=0;i<8;i++) sub[i] = fSubEl[i];
      return;
   }
   if(side<6) {
      sub.Resize(1);
      sub[0]=fSubEl[side];
      return;
   }
   if(side>5 && side<15) {//lados
      int s = side-6;
      sub.Resize(2);
      sub[0] = fSubEl[TPZCompElPr3d::SideNodes[s][0]];
      sub[1] = fSubEl[TPZCompElPr3d::SideNodes[s][1]];
      return;
   }
   if(side>14) {//faces
      int s = side-15;
      sub.Resize(4);
      sub[0] = fSubEl[TPZCompElPr3d::FaceSons[s][0]];
      sub[1] = fSubEl[TPZCompElPr3d::FaceSons[s][1]];
      sub[2] = fSubEl[TPZCompElPr3d::FaceSons[s][2]];
      sub[3] = fSubEl[TPZCompElPr3d::FaceSons[s][3]];
   }
}

TPZGeoElSide TPZGeoElPr3d::Father(int side) {

	TPZGeoEl *fFather = TPZGeoEl::Father();
	if (!fFather) return TPZGeoElSide();
   int whichsub = -1;
   int i;
 	for(i=0;i<8;i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) {
	   PZError << "TPZGeoElPr3d::Father. fFather isn't father of this element.\n";
   	return TPZGeoElSide();
   }
   //agora o atual elemento é o filho numero whichsub < 6
   if(side==20) return TPZGeoElSide(fFather,side);
   //os filhos interiores não tém pai associados a seus cantos
   if(side == whichsub && side<3) return TPZGeoElSide(fFather,side);//cantos
   if((side+1) == whichsub && side>2 && side<6) return TPZGeoElSide(fFather,side);//cantos
   //lados
   if(whichsub == 0 && (side== 6 || side== 8 || side== 9)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side== 6 || side== 7 || side==10)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side== 7 || side== 8 || side==11)) return TPZGeoElSide(fFather,side);
   if(whichsub == 4 && (side== 9 || side==12 || side==14)) return TPZGeoElSide(fFather,side);
   if(whichsub == 5 && (side==10 || side==12 || side==13)) return TPZGeoElSide(fFather,side);
   if(whichsub == 6 && (side==11 || side==13 || side==14)) return TPZGeoElSide(fFather,side);
   //para os filhos 3 e 7 nenhuma aresta esta contida nalguma aresta do pai
   //faces
   if(whichsub == 0 && (side==15 || side==16 || side==18)) return TPZGeoElSide(fFather,side);
   if(whichsub == 1 && (side==15 || side==16 || side==17)) return TPZGeoElSide(fFather,side);
   if(whichsub == 2 && (side==15 || side==17 || side==18)) return TPZGeoElSide(fFather,side);
   if(whichsub == 3 &&  side==19)                          return TPZGeoElSide(fFather,15);
   if(whichsub == 4 && (side==16 || side==18 || side==19)) return TPZGeoElSide(fFather,side);
   if(whichsub == 5 && (side==16 || side==17 || side==19)) return TPZGeoElSide(fFather,side);
   if(whichsub == 6 && (side==17 || side==18 || side==19)) return TPZGeoElSide(fFather,side);
   if(whichsub == 7 &&  side==15)                          return TPZGeoElSide(fFather,19);
   //outro caso
   return TPZGeoElSide();
}

void TPZGeoElPr3d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {

   int nsub = NSideSubElements(side);
   if(!nsub) return;
   sub.Resize(nsub);
   int i,j,k;
   if(nsub==1) {//side = 0 a 5
      int s = side;
      if(side>2) s +=1;
   	if(fSubEl[s]->NodeIndex(side)!=refnode[0]) {
      	PZError << "TPZGeoElPr3d::GetSubElement subelement does not contain refnode" <<std::endl;
         return;
      }
	   sub[0]=TPZGeoElSide(fSubEl[s],side);
   	return;
   }
   //int isub=0;
   for(i=0;i<nsub;i++) {
   	TPZGeoElSide sidesub = SideSubElement(side,i);
      TPZGeoEl *subel = sidesub.Element();
		for(k = 0; k < refnode.NElements(); k++) {//se o subelemento k do thisside tiver o canto refnode[k]
		   for(j=0;j<6;j++) {  //este é um subelemento procurado
			   if(subel->NodeIndex(j)==refnode[k]) {
            	sub[k] = SideSubElement(side,i);
            }
         }
      }
   }
   if(side ==15 || side == 19) sub[3] = SideSubElement(side,3);
   if(side == 20) {
      sub[6] = sub[5];
      sub[5] = sub[4];
      sub[4] = sub[3];
   	sub[3] = SideSubElement(side,3);
      sub[7] = SideSubElement(side,7);
   }
   return;

}

/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
//transforma piramide para piramide ou pirâmide para tetraedro
void TPZGeoElPr3d::BuildTransform(int side,TPZGeoEl *father,TPZTransform<STATE> &t) {
	TPZGeoEl *fFather = TPZGeoEl::Father();
   if(!fFather || side > 20) return;
   int whichsub = -1,i;
 	for(i=0;i<8;i++)  if(fFather->SubElement(i) == this) whichsub = i;
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
	   std::cout << "TPZGeoElPr3d::BuildTransform could not identify the father element\n";
	   return;
   }

   if(side == 20) {
      mult(0,0) = 0.5;//transformação entre elemento mestre do filho
      mult(1,1) = 0.5;//para o elemento mestre do pai
      mult(2,2) = 0.5;
      switch(whichsub) {//o atual é o filho numero whichsub
         case 0:
            sum(2,0) = -.5;
            break;
         case 1:
         	sum(0,0) =  .5;
            sum(2,0) = -.5;
            break;
         case 2:
         	sum(1,0) = .5;
            sum(2,0) =-.5;
            break;
         case 3:
            mult(0,1) =  0.5;
            mult(1,1)*= -1.;
            mult(2,2)*= -1.;
            sum(1,0)  =  .5;
            sum(2,0)  = -.5;
            break;
         case 4:
            sum(2,0) =  .5;
            break;
         case 5:
         	sum(0,0) =  .5;
            sum(2,0) =  .5;
            break;
         case 6:
         	sum(1,0) =  .5;
            sum(2,0) =  .5;
            break;
         case 7:
            mult(0,1) =  0.5;
            mult(1,1)*= -1.;
            mult(2,2)*= -1.;
            sum(1,0)  =  .5;
            sum(2,0)  =  .5;
      }
   //face do filho para a face do pai
   } else if(side>14) {//15 a 19
      int s = side-15;//face do prisma
      if(whichsub==3 || whichsub==7) {
         if(side==15) s = 4;
         if(side==19) s = 0;
      }
      whichsub = -1;
      for(i=0; i<4; i++) if(fFather->SubElement(TPZCompElPr3d::FaceSons[s][i]) == this) whichsub = i;
      if(whichsub == -1) return;
      mult(0,0) = 0.5;
      mult(1,1) = 0.5;
      if(side==15 || side==19) {
         switch(whichsub) {//face triangular
            case 0:
               break;
            case 1:
               sum(0,0) = 0.5;
               break;
            case 2:
               sum(1,0) = 0.5;
               break;
            case 3:
               mult(0,0) =  0.5;
               mult(0,1) =  0.5;
               mult(1,0) =  0.;
               mult(1,1) = -0.5;
               sum(0,0)  =  0.;
               sum(1,0)  =  0.5;
               break;
         }
	   } else {//face quadrilateral
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
      }
   } else if(side>5) {//e side < 15 : arestas
      whichsub = -1;
      int s = side-6;//0 a 7
      for(i=0; i<2; i++)
	      if(fFather->SubElement(TPZCompElPr3d::RibSons[s][i]) == this) whichsub = i;
      if(whichsub == -1) return;
   	mult(0,0) = 0.5;
      if(whichsub==0) sum(0,0) = -0.5;//subelemento 0 do lado
      if(whichsub==1) sum(0,0) =  0.5;//subelemento 1 do lado
   }
   tloc.SetMatrix(mult,sum);
   t = tloc.Multiply(t);
   if(locfather != father) locfather->BuildTransform2(side,father,t);
}

TPZTransform<STATE> TPZGeoElPr3d::SideToSideTransform(int sidefrom,int sideto) {
   if( (sidefrom > 14 &&  sideto < 20) || (sidefrom > 5 &&  sideto < 15) || sideto < 6) {
      PZError << "TPZGeoElPr3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   int dimfrom = SideDimension(sidefrom);
   int dimto = SideDimension(sideto);
   if(dimfrom >= dimto){
      PZError << "TPZGeoElPr3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   TPZTransform<STATE> trans(dimto,dimfrom);//retorna zerada
   //agora : sidefrom < sideto ou dimfrom < dimto
 if(sideto == 20) {//interior
   switch(sidefrom) {//face mestre para elemento mestre
      //faces para interior_____________________________________________________
      case  15:////sidefrom face de dim = 2 : t(3,2)
      	trans.Mult()(0,0) =  1.;//         row = dimto = 3 e col = dimfrom = 2
         trans.Mult()(1,1) =  1.;//         fMult(row,col) : fMult(3,2)
         trans.Sum()(2,0) = -1.;
         break;                  //         fSum(row,1)    : fSum(3,1)
      case  16:
      	trans.Mult()(0,0) =  .5;
         trans.Mult()(2,1) =  1.;
         trans.Sum()(0,0)  = .5;
         break;
      case  17:
         trans.Mult()(0,0) = -.5;
         trans.Mult()(1,0) =  .5;
         trans.Mult()(2,1) =  1.;
         trans.Sum()(0,0)  =  .5;
         trans.Sum()(1,0)  =  .5;
         break;
      case  18:
         trans.Mult()(1,0) =  .5;
         trans.Mult()(2,1) =   1.;
         trans.Sum()(1,0)  =  .5;
         break;
      case  19:
         trans.Mult()(0,0) =  1.;
         trans.Mult()(1,1) =  1.;
         trans.Sum()(2,0)  =  1.;
         break;
      //lados para interior_____________________________________________________
      case  6://row = dimto = 3 e col = dimfrom = 1
      	trans.Mult()(0,0) =  .5;//          sidefrom lado de dim = 1 : t(3,1)
         trans.Sum()(0,0)  =  .5;//          fMult(row,col) : fMult(3,1)
         trans.Sum()(2,0)  = -1.;//          fSum(row,1)    : fSum(3,1)
         break;
      case 7:
      	trans.Mult()(0,0) = -.5;
         trans.Mult()(1,0) =  .5;
         trans.Sum()(0,0)  =  .5;
         trans.Sum()(1,0)  =  .5;
         trans.Sum()(2,0)  = -1.;
         break;
      case 8:
         trans.Mult()(1,0) = -.5;
         trans.Sum()(1,0)  =  .5;
         trans.Sum()(2,0)  = -1.;
         break;
      case 9:
      	trans.Mult()(2,0) =  1.;
         break;
      case 10:
      	trans.Mult()(2,0) =  1.;
         trans.Sum()(0,0)  =  1.;
         break;
      case 11:
         trans.Mult()(2,0) = 1.;
         trans.Sum()(1,0)  = 1.;
         break;
      case 12:
      	trans.Mult()(0,0) =  .5;
         trans.Sum()(0,0)  =  .5;
         trans.Sum()(2,0)  =  1.;
         break;
      case 13:
      	trans.Mult()(0,0) = -.5;
         trans.Mult()(1,0) =  .5;
         trans.Sum()(0,0)  =  .5;
         trans.Sum()(1,0)  =  .5;
         trans.Sum()(2,0)  =  1.;
         break;
      case 14:
         trans.Mult()(1,0) = -.5;
         trans.Sum()(1,0)  =  .5;
         trans.Sum()(2,0)  =  1.;
         break;
   	//canto para interior_____________________________________________________
      case 0://row = dimto = 3 e col = dimfrom = 0
      	trans.Sum()(0,0)  =  0.;//          sidefrom lado de dim = 0 : t(3,0)
         trans.Sum()(1,0)  =  0.;//          fMult(row,col) : fMult(3,0)
         trans.Sum()(2,0)  = -1.;//          fSum(row,1)    : fSum(3,1)
         break;
      case 1:
      	trans.Sum()(0,0)  =  1.;
         trans.Sum()(2,0)  = -1.;
         break;
      case 2:
         trans.Sum()(1,0)  =  1.;
         trans.Sum()(2,0)  = -1.;
         break;
      case 3:
      	trans.Sum()(2,0)  =  1.;
         break;
      case 4:
      	trans.Sum()(0,0)  = 1.;
         trans.Sum()(2,0)  = 1.;
         break;
      case 5:
         trans.Sum()(1,0)  = 1.;
         trans.Sum()(2,0)  = 1.;
         break;
   }
 }//fim if sideto = 20
 if(sideto > 14 && sideto < 20) {
 	switch(sidefrom) {
   	//lados para faces________________________________________________________
      case  6://row = dimto = 2 e col = dimfrom = 1
         if(sideto==15) {
            trans.Mult()(0,0) = .5;
            trans.Sum()(0,0)  = .5;
         } else if(sideto==16) {
            trans.Mult()(0,0) =  1.;
            trans.Sum()(1,0)  = -1.;
         }
         break;
      case 7:
      	if(sideto==15) {
            trans.Mult()(0,0) = -.5;
            trans.Mult()(1,0) =  .5;
            trans.Sum()(0,0)  =  .5;
            trans.Sum()(1,0)  =  .5;
         } else if(sideto==17) {
            trans.Mult()(0,0) =  1.;
            trans.Sum()(1,0)  = -1.;
         }
         break;
      case 8:
      	if(sideto==15) {
            trans.Mult()(1,0) = -.5;
            trans.Sum()(1,0)  =  .5;
         } else if(sideto==18) {
            trans.Mult()(0,0) = -1.;
            trans.Sum()(1,0)  = -1.;
         }
         break;
      case 9:
         if(sideto==16 || sideto==18) {
            trans.Mult()(1,0) =  1.;
            trans.Sum()(0,0)  = -1.;
         }
         break;
      case 10:
      	if(sideto==16) {
            trans.Mult()(1,0) =  1.;
            trans.Sum()(0,0)  =  1.;
         } else if(sideto==17) {
            trans.Mult()(1,0) =  1.;
            trans.Sum()(0,0)  = -1.;
         }
         break;
      case 11:
      	if(sideto==17 || sideto==18) {
            trans.Mult()(1,0) =  1.;
            trans.Sum()(0,0)  =  1.;
         }
         break;
      case  12:
      	if(sideto==16) {
            trans.Mult()(0,0) = 1.;
            trans.Sum()(1,0)  = 1.;
         } else if(sideto==19) {
            trans.Mult()(0,0) = .5;
            trans.Sum()(0,0)  = .5;
         }
         break;
      case  13:
      	if(sideto==19) {
            trans.Mult()(0,0) = -.5;
            trans.Mult()(1,0) =  .5;
            trans.Sum()(0,0)  =  .5;
            trans.Sum()(1,0)  =  .5;
         } else if(sideto==17) {
            trans.Mult()(0,0) =  1.;
            trans.Sum()(1,0)  =  1.;
         }
         break;
      case  14:
      	if(sideto==19) {
            trans.Mult()(1,0) = -.5;
            trans.Sum()(1,0)  =  .5;
         } else if(sideto==18) {
            trans.Mult()(0,0) = -1.;
            trans.Sum()(1,0)  =  1.;
         }
         break;
   	//cantos para faces_______________________________________________________
      case  0://row = dimto = 2 e col = dimfrom = 0
      	if(sideto==15) {//sidefrom lado de dim =0 : t(2,0)
            trans.Sum()(0,0) = 0.;
            trans.Sum()(1,0) = 0.;
         } else
         if(sideto==18 || sideto==16) {
            trans.Sum()(0,0) =-1.;//          fMult(row,col) : fMult(2,0)
            trans.Sum()(1,0) =-1.;//          fSum(row,1)    : fSum(2,1)
         }
         //if(sideto==18 && fFather->SubElement(3) == this) {//dado que nao corresponde
         //   trans.Sum()(0,0) =  1.;                      //nenhum irmão com o
         //   trans.Sum()(1,0) =  1.;                      //mesmo canto
         //}
         break;
      case  1:
      	if(sideto==15) {
            trans.Sum()(0,0) = 1.;
         } else
         if(sideto==16) {
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) = -1.;
         } else
         if(sideto==17) {
            trans.Sum()(0,0) = -1.;
            trans.Sum()(1,0) = -1.;
         }
         break;
      case  2:
      	if(sideto==15) {
            trans.Sum()(1,0) =  1.;
         } else
         if(sideto==17 || sideto==18) {
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) = -1.;
         }
         break;
      case  3:
      	if(sideto==19) {
            trans.Sum()(0,0) = 0.;
            trans.Sum()(1,0) = 0.;
         } else
         if(sideto==16 || sideto==18) {
            trans.Sum()(0,0) = -1.;
            trans.Sum()(1,0) =  1.;
         }
         break;
      case 4:
      	if(sideto==19) {
            trans.Sum()(0,0) = 1.;
         } else
         if(sideto==16) {
            trans.Sum()(0,0) =  1.;
            trans.Sum()(1,0) =  1.;
         } else
         if(sideto==17) {
            trans.Sum()(0,0) = -1.;
            trans.Sum()(1,0) =  1.;
         }
         break;
      case  5:
      	if(sideto==19) {
            trans.Sum()(1,0) = 1.;
         } else
         if(sideto==17 || sideto==18) {
            trans.Sum()(0,0) = 1.;
            trans.Sum()(1,0) = 1.;
         }
         break;
      default:
	      PZError << "TPZGeoElPr3d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
   }//switch
 }//fim lados cantos para faces
 if(sideto > 5 && sideto < 15) {//6 a 14
 	switch(sidefrom) {
   	//canto para lado_________________________________________________________
      case  0:
         if(sideto== 8)               trans.Sum()(0,0)  =  1.;//sidefrom canto de dim =0 : t(1,0)
         if(sideto== 6 || sideto== 9) trans.Sum()(0,0)  = -1.;//row = dimto = 1 e col = dimfrom = 0
         break;                                               //fMult(row,col) : fMult(1,0) , fSum(row,1) : fSum(1,1)
      case  1:
      	if(sideto== 6)               trans.Sum()(0,0)  =  1.;
         if(sideto== 7 || sideto==10) trans.Sum()(0,0)  = -1.;
      	break;
      case  2:
         if(sideto== 7)               trans.Sum()(0,0)  =  1.;
         if(sideto== 8 || sideto==11) trans.Sum()(0,0)  = -1.;
      	break;
      case  3:
      	if(sideto== 9 || sideto==14) trans.Sum()(0,0)  =  1.;
         if(sideto==12)               trans.Sum()(0,0)  = -1.;
      	break;
      case  4:
         if(sideto==13)               trans.Sum()(0,0) = -1.;
         if(sideto==10 || sideto==12) trans.Sum()(0,0) =  1.;
         break;
      case  5:
      	if(sideto==14)               trans.Sum()(0,0)  = -1.;
         if(sideto==11 || sideto==13) trans.Sum()(0,0)  =  1.;
  	}//switch
 }//if cantos
   return trans;
}

TPZCompEl *TPZGeoElPr3d::CreateBCCompEl(int side,int bc,TPZCompMesh &cmesh) {
	if(side<0 || side>20) return 0;

   if(side==20) {
	   std::cout << "TPZGeoElPr3d::CreateBCCompEl with side = 18 not implemented\n";
      return 0;
   }
   int64_t index;
	if(side<6) {
      TPZMaterial *bcptr = cmesh.FindMaterial(bc);
      if(!bcptr) {
      PZError << "TPZGeoElPr3d::CreateBCCompEl has no bc.\n";
      return 0;
      }
      TPZCompEl *cel = Reference();
      if(!cel) {
      PZError << "TPZGeoElPr3d::CreateBCCompEl has no computational element\n";
      return 0;
      }
	  TPZManVector<int64_t> nodeindexes(1);
      TPZGeoElPoint *gel;
      nodeindexes[0] = fNodeIndexes[side];
      gel = new TPZGeoElPoint(nodeindexes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,side));
      return gel->CreateCompEl(cmesh,index);
   } else if (side > 5 && side < 15) {//side = 6 a 14 : arestas
      TPZManVector<int64_t> nodes(2);
      int s = side-6;
      nodes[0] = NodeIndex(TPZCompElPr3d::SideNodes[s][0]);
      nodes[1] = NodeIndex(TPZCompElPr3d::SideNodes[s][1]);
      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::SideNodes[s][0]));
      TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::SideNodes[s][1]));
      TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(this,side));
      return gel->CreateCompEl(cmesh,index);
   } else if (side > 14) {//side = 15 a 19 : faces
      TPZManVector<int64_t> nodes(4);//4o = -1 para face triangular
      int s = side-15;
      nodes[0] = NodeIndex(TPZCompElPr3d::FaceNodes[s][0]);
      nodes[1] = NodeIndex(TPZCompElPr3d::FaceNodes[s][1]);
      nodes[2] = NodeIndex(TPZCompElPr3d::FaceNodes[s][2]);
      nodes[3] = NodeIndex(TPZCompElPr3d::FaceNodes[s][3]);
      TPZGeoElT2d *gelt;
      TPZGeoElQ2d *gelq;
      if(side>15 && side<19) {
      	gelq = new TPZGeoElQ2d(nodes,bc,*Mesh());
         TPZGeoElSide(gelq,0).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceNodes[s][0]));
         TPZGeoElSide(gelq,1).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceNodes[s][1]));
         TPZGeoElSide(gelq,2).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceNodes[s][2]));
         TPZGeoElSide(gelq,3).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceNodes[s][3]));
         TPZGeoElSide(gelq,4).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceSides[s][0]));
         TPZGeoElSide(gelq,5).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceSides[s][1]));
         TPZGeoElSide(gelq,6).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceSides[s][2]));
         TPZGeoElSide(gelq,7).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceSides[s][3]));
         TPZGeoElSide(gelq,8).SetConnectivity(TPZGeoElSide(this,side));
         return gelq->CreateCompEl(cmesh,index);
      } else {
         nodes.Resize(3);
	      gelt = new TPZGeoElT2d(nodes,bc,*Mesh());
         TPZGeoElSide(gelt,0).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceNodes[s][0]));
         TPZGeoElSide(gelt,1).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceNodes[s][1]));
         TPZGeoElSide(gelt,2).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceNodes[s][2]));
         TPZGeoElSide(gelt,3).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceSides[s][0]));
         TPZGeoElSide(gelt,4).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceSides[s][1]));
         TPZGeoElSide(gelt,5).SetConnectivity(TPZGeoElSide(this,TPZCompElPr3d::FaceSides[s][2]));
         TPZGeoElSide(gelt,6).SetConnectivity(TPZGeoElSide(this,side));
         return gelt->CreateCompEl(cmesh,index);
      }
   } else PZError << "TPZGeoElPr3d::CreateBCCompEl. Side = " << side <<std::endl;
   return 0;
}

void TPZGeoElPr3d::NodeFaceIds(TPZVec<int> &ids,int face) {

	ids.Resize(4,-1);
   if((face>-1 && face<6) || (face>14 && face<20)) {//cada condicao refer-se à mesma face
   	if(face>14) face = face-15;
      ids[0] = NodeIndex(TPZCompElPr3d::FaceNodes[face][0]);
      ids[1] = NodeIndex(TPZCompElPr3d::FaceNodes[face][1]);
      ids[2] = NodeIndex(TPZCompElPr3d::FaceNodes[face][2]);
		ids[3] = NodeIndex(TPZCompElPr3d::FaceNodes[face][3]);
      if(face==0 && face==4) ids.Resize(3);//face triangular
      return;
   }
   std::cout << "TPZCompElPr3d::NodeFaceIds bad side , side = " << face <<std::endl;
}

static int fatherside[8][21] = {
/*00*/{0,6,8,9,16,18,6,15,8,9,16,18,16,20,18,15,16,20,18,20,20},
/*01*/{6,1,7,16,10,17,6,7,15,16,10,17,16,17,20,15,16,17,20,20,20},
/*02*/{8,7,2,18,17,11,15,7,8,18,17,11,20,17,18,15,20,17,18,20,20},
/*03*/{18,17,16,8,7,6,20,20,20,18,17,16,15,15,15,20,20,20,20,15,20},
/*04*/{9,16,18,3,12,14,16,20,18,9,16,18,12,19,14,20,16,20,18,19,20},
/*05*/{16,10,17,12,4,13,16,17,20,16,10,17,12,13,19,20,16,17,20,19,20},
/*06*/{18,17,11,14,13,5,20,17,18,18,17,11,19,13,14,20,20,17,18,19,20},
/*07*/{14,13,12,18,17,16,19,19,19,18,17,16,20,20,20,19,20,20,20,20,20},
};

TPZGeoElSide TPZGeoElPr3d::Father2(int side){//Augusto:09/01/01
	TPZGeoEl *fFather = TPZGeoEl::Father();

	if(side<0 || side>20 || fFather==0){
		PZError << "TPZGeoElPr3d::Father2 called error" <<std::endl;
		return TPZGeoElSide();
	}
	int subelindex = WhichSubel();
	if(fatherside[subelindex][side]<0){
		PZError << "TPZGeoElPr3d::Father2 called with index error\n";
		return TPZGeoElSide();
	}
	return TPZGeoElSide(fFather,fatherside[subelindex][side]);
}

static int subeldata[21][21][2] = {
/*00*/{{0,0}},
/*01*/{{1,1}},
/*02*/{{2,2}},
/*03*/{{4,3}},
/*04*/{{5,4}},
/*05*/{{6,5}},
/*06*/{{0,6},{0,1},{1,6}},
/*07*/{{1,7},{1,2},{2,7}},
/*08*/{{0,8},{0,2},{2,8}},
/*09*/{{0,9},{0,3},{4,9}},
/*10*/{{1,10},{1,4},{5,10}},
/*11*/{{2,11},{2,5},{6,11}},
/*12*/{{4,12},{4,4},{5,12}},
/*13*/{{5,13},{5,5},{6,13}},
/*14*/{{4,14},{4,5},{6,14}},
/*15*/{{0,15},{1,15},{2,15},{3,19},{0,7},{1,8},{2,6}},
/*16*/{{0,16},{1,16},{4,16},{5,16},{0,10},{0,12},{5,6},{5,9},{0,4}},
/*17*/{{1,17},{2,17},{5,17},{6,17},{1,11},{1,13},{6,7},{6,10},{1,5}},
/*18*/{{0,18},{2,18},{4,18},{6,18},{0,11},{0,14},{6,8},{6,9},{2,3}},
/*19*/{{4,19},{5,19},{6,19},{7,15},{4,13},{5,14},{6,12}},
/*20*/{{0,20},{1,20},{2,20},{3,20},{4,20},{5,20},{6,20},{7,20},{0,17},
       {0,19},{4,17},{1,18},{1,19},{5,18},{2,16},{2,19},{6,16},{3,15},
       {0,13},{1,14},{2,12}}
};

static int nsubeldata[21] = {1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,7,9,9,9,7,21};


void TPZGeoElPr3d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){

    subel.Resize(0);
    if(side<0 || side>20 || !HasSubElement(side)){
       PZError << "TPZGeoElPr3d::GetSubElements2 called with error arguments\n";
       return;
    }
    int nsub = nsubeldata[side];
    for(int i=0;i<nsub;i++)
    		subel.Push(TPZGeoElSide(fSubEl[subeldata[side][i][0]],
                                       subeldata[side][i][1]));
}

static REAL buildt[8][21][4][3] = {//por colunas
/*S0*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,-1}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*07*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*10*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*11*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*12*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*13*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*14*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*17*/{{-0.25,0.25,0.},{0.,0.,0.5},{0.,0.,0.},{0.25,0.25,-0.5}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,0.,-.5}}},
/*S1*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1.,0.,-1.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*08*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*09*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*11*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*12*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*13*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*14*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*16*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*17*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*18*/{{0.,0.25,0.},{0.,0.,0.5},{0.,0.,0.},{0.5,0.25,-0.5}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,0.,-.5}}},
/*S2*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,1.,-1.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*09*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*10*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*12*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*13*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*14*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*16*/{{0.25,0.,0.},{0.,0.,0.5},{0.,0.,0.},{0.25,0.5,-0.5}},
      /*17*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,.5,-.5}}},
/*S3*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*08*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*09*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*10*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*11*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*12*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*13*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*14*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*15*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*16*/{{0.25,0.,0.},{0.,0.,-0.5},{0.,0.,0.},{0.25,0.5,-0.5}},
      /*17*/{{0.,-0.25,0.},{0.,0.,-0.5},{0.,0.,0.},{0.5,0.25,-0.5}},
		  /*18*/{{0.25,-0.25,0},{0,0,0.5},{0.,0.,0.},{0.25,0.25,-0.5}},      /*19*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*20*/{{.5,0.,0.},{.5,-.5,0.},{0.,0.,-.5},{0.,.5,-.5}}},
/*S4*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,1.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*07*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*10*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*11*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*13*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*14*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*17*/{{-0.25,0.25,0.},{0.,0.,0.5},{0.,0.,0.},{0.25,0.25,0.5}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,0.,.5}}},
/*S5*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1.,0.,1.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*08*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*09*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*11*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*13*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*14*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*16*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*17*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*18*/{{0.,0.25,0.},{0.,0.,0.5},{0.,0.,0.},{0.5,0.25,0.5}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,0.,.5}}},
/*S6*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,1.,1.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*09*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*10*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*12*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*13*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*14*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*16*/{{0.25,0.,0.},{0.,0.,0.5},{0.,0.,0.},{0.25,0.5,0.5}},
      /*17*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,.5,.5}}},
/*S7*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*08*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*09*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*10*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*11*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*12*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*13*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*14*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*15*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*16*/{{0.25,0.,0.},{0.,0.,-0.5},{0.,0.,0.},{0.25,0.5,0.5}},
      /*17*/{{0.,-0.25,0.},{0.,0.,-0.5},{0.,0.,0.},{0.5,0.25,0.5}},
      /*18*/{{0.25,-0.25,0.},{0.,0.,-0.5},{0.,0.,0.},{0.25,0.25,0.5}},
      /*19*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*20*/{{.5,0.,0.},{.5,-.5,0.},{0.,0.,-.5},{0.,.5,.5}}}
};

TPZTransform<STATE> TPZGeoElPr3d::BuildTransform2(int side, TPZGeoEl * /*father*/){//Augusto:09/01/01
	TPZGeoEl *fFather = TPZGeoEl::Father();

	if(side<0 || side>20 || !fFather){
  	PZError << "TPZGeoElPr3d::BuildTransform2 side out of range or father null\n";
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

static REAL MidSideNode[21][3] = {
/*00*/{0.,.0,-1.},/*01*/{1.,0.,-1.},/*02*/{.0,1.,-1.},/*03*/{.0,.0, 1.},
/*04*/{1.,.0, 1.},/*05*/{0.,1., 1.},/*06*/{.5,.0,-1.},/*07*/{.5,.5,-1.},
/*08*/{.0,.5,-1.},/*09*/{0.,.0, 0.},/*10*/{1.,.0, 0.},/*11*/{.0,1., 0.},
/*12*/{.5,.0, 1.},/*13*/{.5,.5, 1.},/*14*/{.0,.5, 1.},/*15*/{1./3.,1./3.,-1.},
/*16*/{.5,.0, 0.},/*17*/{.5,.5, 0.},/*18*/{0.,.5, 0.},/*19*/{1./3.,1./3., 1.},
/*20*/{1./3.,1./3.,0.} };

int TPZGeoElPr3d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3),pss(3),pf(3),pfs(3);
                    //point son, point side son, point father, point side father : elemento mestre
  for(sn=0;sn<8;sn++){
    TPZGeoEl *son = subs[sn];
    for(sd=0;sd<21;sd++){
      ps[0] = MidSideNode[sd][0];//element
      ps[1] = MidSideNode[sd][1];//master point
      ps[2] = MidSideNode[sd][2];
      if(son->WhichSide(ps) != sd) std::cout << "Lado nao bate\n";
      TPZTransform<STATE> telsd = pzshape::TPZShapePrism::TransformElementToSide(sd);//2x2
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform<STATE> t;
	  son->BuildTransform2(sd, gel,t);
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = pzshape::TPZShapePrism::TransformSideToElement(sdfat);//2x2
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(20).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao errada\n";
        PZError << "son    = " << (son->Id()) <<std::endl;
        PZError << "father = " << ((son->Father2(20).Element())->Id()) <<std::endl;
        PZError << "side   = " << sd <<std::endl <<std::endl;
        int ok;
		std::cin >> ok;
      } else {
		  std::cout << "Transformacao OK!\n";
		  std::cout << "Filho/lado : " << son->Id() << "/" << sd <<std::endl;
		  std::cout << "Pai : " << son->Father2(20).Element()->Id() <<std::endl <<std::endl;
      }
    }
  }
  return 1;
}


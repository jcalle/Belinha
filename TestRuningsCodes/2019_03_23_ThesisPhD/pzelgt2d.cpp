//METHODS DEFINITION FOR CLASS ELEMT2D

#include "pzelgt2d.h"
#include "pzelct2d.h"
#include "pzelgq2d.h"
#include "pzelg1d.h"
#include "pzshapetriang.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzgnode.h"
#include "pzfmatrix.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzelgpoint.h"

static TPZCompEl *CreateEl(TPZGeoElT2d *gel,TPZCompMesh &mesh, int64_t &index) {
  return new TPZCompElT2d(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElT2d::fp)(TPZGeoElT2d *,TPZCompMesh &, int64_t &) = CreateEl;

/**The vector with nodes has first the three corners nodes. The last three nodes are for middle nodes
   on the sides*/
TPZGeoElT2d::TPZGeoElT2d(int id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh)
  : TPZGeoElRefLess<pzgeom::TPZGeoTriangle>(id,nodeindexes,matind,mesh) {

  int i,nnod = nodeindexes.NElements();
  if(nnod != 3) {
    PZError << "TPZGeoElT2d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<3;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<4;i++) fSubEl[i] = 0;
}

TPZGeoElT2d::TPZGeoElT2d() {
  int i;
  for(i=0;i<3;i++) fNodeIndexes[i] = -1;
  for(i=0;i<4;i++) fSubEl[i] = 0;
}

TPZGeoElT2d::TPZGeoElT2d(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) :
	TPZGeoElRefLess<pzgeom::TPZGeoTriangle>(nodeindexes,matind,mesh) {


  int i,nnod = nodeindexes.NElements();
  if(nnod != 3) {
    PZError << "TPZGeoElT2d::Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }

  for(i=0;i<3;i++) fNodeIndexes[i] = nodeindexes[i];
  for(i=0;i<4;i++) fSubEl[i] = 0;
}

TPZGeoElT2d::~TPZGeoElT2d() { }

void TPZGeoElT2d::Shape(TPZVec<REAL> &param,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {
  REAL x = param[0], y = param[1];
  phi(0,0) = 1.-x-y;
  phi(1,0) = x;
  phi(2,0) = y;
  dphi(0,0) = dphi(1,0) = -1.;
  dphi(0,1) = dphi(1,2) = 1.;
  dphi(1,1) = dphi(0,2) = 0.;
}

TPZGeoElT2d *TPZGeoElT2d::CreateGeoEl(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  return new TPZGeoElT2d(nodeindexes,matind,mesh);
}

int64_t TPZGeoElT2d::NodeIndex(int node) {
  if(node<0 || node>2) return -1;
  return fNodeIndexes[node];
}

int TPZGeoElT2d::NSideNodes(int side) {
  if(side<0 || side>6) {
    PZError << "TPZGeoElT2d::NSideNodes. Bad parameter side.\n";
    return 0;
  }
  if(side<3) return 1;
  if(side<6) return 2;
  return 3;
}

int64_t TPZGeoElT2d::SideNodeIndex(int side,int node) {
  if(side<0 || side>6) {
    PZError << "TPZGeoElT2d::SideNodeIndex. Bad parameter side.\n";
    return -1;
  }
  if(side==6) return fNodeIndexes[node];
  if(side < 3) return fNodeIndexes[side];
  //3,4,5
  side-=3;//cedric
  if(node==0) return fNodeIndexes[side];
  if(node==1) return fNodeIndexes[(side+1)%3];
  return -1;
}

void TPZGeoElT2d::MidSideNodeIndex(int side,int64_t &index) {
  index = -1;
  if(side<0 || side>5) {
    PZError << "TPZGeoElT2d::MidSideNodeIndex. Bad parameter side = " << side <<std::endl;
    return;
  }
  //lados 0,1,2
  if(side<3) {//o nó medio do lado 0 é o 0 etc.
    index=fNodeIndexes[side];
    return;
  }
  // CEDRIC VERIFICAR : feito!!!
  //lados 3,4,5
	if(HasSubElement(0)) {//caso tenha filhos, se não tiver retorna -1
	  switch(side) {
		  case 3:
		    index=((TPZGeoElT2d *) fSubEl[0])->fNodeIndexes[1];
          break;
		  case 4:
		    index=((TPZGeoElT2d *) fSubEl[1])->fNodeIndexes[2];
          break;
		  case 5:
		    index=((TPZGeoElT2d *) fSubEl[2])->fNodeIndexes[0];
          break;
	  }
   }
}

void TPZGeoElT2d::NewMidSideNode(int side,int64_t &index) {
  MidSideNodeIndex(side,index);
  if(index < 0) {
    TPZGeoElSide gelside = Neighbour(side);
    while(gelside.Exists()) {
      gelside.Element()->MidSideNodeIndex(gelside.Side(),index);
      if(index!=-1) return;
      gelside = gelside.Neighbour();
      if(gelside.Element()==this) break;
    }
    TPZVec<REAL> par(2,.5);
    TPZVec<REAL> coord(3,0.);
    //lados do elemento neste caso são : side = 0,1,2
    if(side==3) par[1] = 0.;
    else if(side==5) par[0] = 0.;
    X(par,coord);
    index = Mesh()->NodeVec().AllocateNewElement();
    Mesh()->NodeVec()[index].Initialize(coord,*Mesh());
  }
}

/**Determine the coordinate of the center of the element*/
int64_t TPZGeoElT2d::CenterIndex() {
  TPZVec<REAL> coord(3);
  TPZVec<REAL> param(2,1./3.);
  X(param,coord);  //It determine the centroid of the triangular element.
  int64_t indexcenter = Mesh()->NodeVec().AllocateNewElement();
  Mesh()->NodeVec()[indexcenter].Initialize(coord,*Mesh());
  return indexcenter;
  // We can not to take the center type of the master cell because of characteristic
  // point is not preserved for the transformation.
}

int TPZGeoElT2d::SideDimension(int side) {
  if(side<3 || side>6) return 0;
  if(side==6) return 2;
  return 1;
}

TPZGeoElSide TPZGeoElT2d::HigherDimensionSides(int side,int targetdimension) {
//targetdimension deve ser 1 ou 2
//o lado 6 tem dim. 2 e neste caso targetdimension = 3 logo : 0 <= side <= 5
  if( (side<0 || side>5) || (targetdimension != 1 && targetdimension != 2) ) {
    PZError << "TPZGeoElT2d::HigherDimensionSides called with side = " << side
	    << " targetdimension = " << targetdimension <<std::endl;
    return TPZGeoElSide();//retorna objeto nulo {0,-1}
  }
  TPZGeoEl *father = TPZGeoEl::Father();
  if (!father || Father(side).Exists()) return TPZGeoElSide();
//side = 0,1,2,3,4,5
  switch(targetdimension) {
	  case 1:
       if(this == father->SubElement(0)) {
       	 if(side==1) return TPZGeoElSide(this,3);
       	 if(side==2) return TPZGeoElSide(this,5);
       } else if(this == father->SubElement(1)) {
       	 if(side==0) return TPZGeoElSide(this,3);
       	 if(side==2) return TPZGeoElSide(this,4);
       } else if(this == father->SubElement(2)) {
       	 if(side==0) return TPZGeoElSide(this,5);
       	 if(side==1) return TPZGeoElSide(this,4);
       } else if(this == father->SubElement(3)) {
       	 if(side==0) return TPZGeoElSide(father->SubElement(1),4);
       	 if(side==1) return TPZGeoElSide(father->SubElement(0),5);
       	 if(side==2) return TPZGeoElSide(father->SubElement(0),3);
       }
       return TPZGeoElSide();//retorna objeto nulo {0,-1}
  	  case 2:
       	return TPZGeoElSide(this,6);
  }//switch
  return TPZGeoElSide();
}

void TPZGeoElT2d::LowerDimensionSides(int side,TPZStack<TPZGeoElSide> &smallsides) {
// CEDRIC VERIFICAR : feito!!!
	if (side < 3) return;
  switch(side) {
	  case 3:
   	 smallsides.Push(TPZGeoElSide(this,0));
  	    smallsides.Push(TPZGeoElSide(this,1));
	    break;
	  case 4:
   	 smallsides.Push(TPZGeoElSide(this,1));
  	    smallsides.Push(TPZGeoElSide(this,2));
	    break;
	  case 5:
   	 smallsides.Push(TPZGeoElSide(this,2));
  	    smallsides.Push(TPZGeoElSide(this,0));
	    break;
	  case 6:
   	 smallsides.Push(TPZGeoElSide(this,0));
  	    smallsides.Push(TPZGeoElSide(this,1));
   	 smallsides.Push(TPZGeoElSide(this,2));
  	    smallsides.Push(TPZGeoElSide(this,3));
   	 smallsides.Push(TPZGeoElSide(this,4));
  	    smallsides.Push(TPZGeoElSide(this,5));
//   	 smallsides.Push(TPZGeoElSide(this,6));
//Philippe 17/5/99
  }
}

void TPZGeoElT2d::SideMasterCo(int side,TPZFMatrix<STATE> &coord) {
// CEDRIC VERIFICAR , comparar com 1D : verificado!!!
  if(side<0 || side>6) {
    PZError << "TPZGeoElT2d::SideMasterCo. Bad parameter side.\n";
    return;
  }
  int row = coord.Rows();//2x3
  if(side==6) {
    coord.Redim(row,3);          //side 6 tem dimensão 2
    coord(0,1) = coord(1,2) = 1.;// [ {0,1,0},{0,0,1} ]
    return;                      //ou (0,0),(1,0),(0,1) nós 0,1,2
  }
  if(side<3) {         //side<3 tem dimensão 0
    coord.Redim(row,1);//2x1
    coord.Zero();      //side = 0 => (0,0)
    if(side==1) coord(0,0)=1.;//(1,0)
    else if(side==2) coord(1,0)=1.;//(0,1)
    return;
  }
  coord.Redim(row,2);  //side = 3,4,5 tem dimensão 1
  coord.Zero();        //2x2  [{0,0},{0,0}]
  if(side==3) coord(0,1) = 1.;//(0,0),(1,0)
  else if(side==5) coord(1,0) = 1.;//(0,1),(0,0)
  else if (side == 4) coord(0,0) = coord(1,1) = 1.;//(1,0),(0,1)
}

void TPZGeoElT2d::Jacobian(TPZVec<REAL> &param,TPZFMatrix<STATE> &jacobian,TPZFMatrix<STATE> &axes,REAL &detjac,TPZFMatrix<STATE> &jacinv) {

  REAL spacephi[3],spacedphi[6];
  int i,j;
  TPZFMatrix<STATE> phi(3,1,spacephi,3);
  TPZFMatrix<STATE> dphi(2,3,spacedphi,6);
  jacobian.Zero();
  Shape(param,phi,dphi);

  TPZVec<REAL> V1(3,0.),V2(3,0.),V2til(3,0.),V3(3,0.);
  double V1Norm=0.,V1V2=0.,V2tilNorm=0.;
  TPZGeoNode *np;

  for(i=0;i<3;i++) {
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

void TPZGeoElT2d::X(TPZVec<REAL> &par,TPZVec<REAL> &result){

  REAL spacephi[3],spacedphi[6];
  TPZFMatrix<STATE> phi(3,1,spacephi,3);
  TPZFMatrix<STATE> dphi(2,3,spacedphi,6);
  Shape(par,phi,dphi);

  int i,j;
  for(i=0;i<3;i++) {
    result[i] = 0.0;
    for(j=0;j<3;j++)
      result[i] += phi(j,0)*NodePtr(j)->Coord(i);
  }
}

void TPZGeoElT2d::NormalVector(int side,TPZVec<REAL> &loc,
			     TPZVec<REAL> &normal,TPZFMatrix<STATE> &axes,TPZFMatrix<STATE> &jacside) {
  if(side < 0 || side > 6) {
    PZError << "TPZGeoElT2d.NormalVector invalid side : "
      " side = " << side << "\n";
    return;
  }

  TPZVec<REAL> t(3,0.);
  int i,id,ic,j=(side+3)%3;
  REAL detjac;
  TPZFMatrix<STATE> jacinv(0);
  if(side==6 || side<3) {
    REAL norm = 0.;
    TPZVec<REAL> t1(3,0.);
    for(i=0;i<3;i++) {
      t[i]= NodePtr(j)->Coord(i);
      t1[i]= t[i] - NodePtr((j+1)%3)->Coord(i);
      t[i]-=NodePtr((j+2)%3)->Coord(i);
    }
    if(side==6) {
    	normal[0]=t[1]*t1[2]-t1[1]*t[2];
      normal[1]=t[2]*t1[0]-t1[2]*t[0];
      normal[2]=t[0]*t1[1]-t[1]*t1[0];
      jacside.Redim(2,2);
      Jacobian(loc,jacside,axes,detjac,jacinv);
    }
    else {
      for(i=0;i<3;i++) normal[i]=t[i]+t1[i];
      jacside.Redim(0,0);
      axes.Zero();
      for(i=0;i<3;i++) axes(i,i)=1.;
    }
    for(i=0;i<3;i++) norm += normal[i]*normal[i];
    norm = sqrt(norm);
    for(i=0;i<3;i++) normal[i]/=norm;
    return;
  }
  REAL spacephi[3],spacedphi[6];
  TPZFMatrix<STATE> phi(3,1,spacephi,3);
  TPZFMatrix<STATE> dphi(2,3,spacedphi,6);
  Shape(loc,phi,dphi);

  TPZGeoNode* np;
  REAL ider[2] = {0.,0.}, sq2 = sqrt(2.)/2.;
  switch(j) {
  case 0: ider[0] = 1.;break;
  case 1: ider[0] = -sq2; ider[1] = sq2; break;
  case 2: ider[1] = -1.; break;
  }
  for(i=0;i<3;i++) {
    np = NodePtr(i);
    for(ic=0;ic<3;ic++)
      for(id=0;id<2;id++)
	     t[ic] += ider[id] * (np->Coord(ic)) * dphi(id,i);
  }

  REAL jac1dvar,tnorm;
  tnorm = sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
  if(j != 1) jac1dvar = 0.5*tnorm;
  else jac1dvar = sq2*tnorm;
  jacside(0,0) = jac1dvar;

  TPZVec<REAL> V1(3,0.),V2(3,0.),V2til(3,0.),V3(3,0.),V1til(3,0.),V1Ttil(3,0.);
  REAL V1Norm=0.,V2Norm=0.,V1V2=0.,V2tilNorm=0.,V1tilNorm =0.,tV1 = 0.,V1TtilNorm = 0.;
  for(i=0;i<3;i++) {
    np = NodePtr(i);
    for(ic=0;ic<3;ic++) {
      V1[ic] += np->Coord(ic)*dphi(0,i);
      V2[ic] += np->Coord(ic)*dphi(1,i);
    }
  }
  for(ic=0;ic<3;ic++) {
    V1Norm += V1[ic]*V1[ic];
    V2Norm += V2[ic]*V2[ic];
    V1V2 += V1[ic]*V2[ic];
    tV1 += t[ic]*V1[ic];
  }
  V1Norm = sqrt(V1Norm);
  V2Norm = sqrt(V2Norm);
  for(ic=0;ic<3;ic++) {
    V1[ic] /= V1Norm;
    V2[ic] /= V2Norm;
    //Jorge 14/10/99
    t[ic] /= tnorm;
    V2til[ic] = V2[ic] - V1V2*V1[ic]/V1Norm/V2Norm;
    V1til[ic] = V1[ic] - V1V2*V2[ic]/V1Norm/V2Norm;
    //Jorge 14/10/99
    V1Ttil[ic] = V1[ic] - tV1*t[ic]/V1Norm/tnorm;
    V2tilNorm += V2til[ic]*V2til[ic];
    V1tilNorm += V1til[ic]*V1til[ic];
    V1TtilNorm += V1Ttil[ic]*V1Ttil[ic];
  }
  V2tilNorm = sqrt(V2tilNorm);
  V1tilNorm = sqrt(V1tilNorm);
  V1TtilNorm = sqrt(V1TtilNorm);
  for(ic=0;ic<3;ic++) {
    axes(0,ic) = V1[ic];
    axes(1,ic) = V2til[ic]/V2tilNorm;
  }
  switch(j) {
  case 0:
    normal[0] = -V2til[0]/V2tilNorm;
    normal[1] = -V2til[1]/V2tilNorm;
    normal[2] = -V2til[2]/V2tilNorm;
    break;
  case 1:
    normal[0] = V1Ttil[0]/V1TtilNorm;
    normal[1] = V1Ttil[1]/V1TtilNorm;
    normal[2] = V1Ttil[2]/V1TtilNorm;
    break;
  case 2:
    normal[0] = -V1til[0]/V1tilNorm;
    normal[1] = -V1til[1]/V1tilNorm;
    normal[2] = -V1til[2]/V1tilNorm;
    break;
  }
  axes(2,0) = axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1);
  axes(2,1) = -axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0);
  axes(2,2) = axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0);
}

/** TO SUBDIVISION
********************************************************************************
  Into Divides is necesary to consider the connectivity with the all neighboards*/
void TPZGeoElT2d::Divide(TPZVec<TPZGeoEl *> &SubElVec) {
// VERIFICAR CEDRIC
   //If exist fSubCell return this sons
   int i;
   if(HasSubElement(0)) {
      SubElVec.Resize(4);
      for(i=0;i<4;i++) SubElVec[i] = fSubEl[i];
      return;
   }

   int j, sub, matid = MaterialId();
   int64_t index;
   int np[4][3];
   for(j=0;j<3;j++) np[j][j]=NodeIndex(j);
   for(sub=0;sub<3;sub++) {
      NewMidSideNode(sub+3,index);
      j=(sub+1)%3;
      np[sub][j]=np[j][sub]=index;
   }
   np[3][0] = np[1][2]; np[3][1] = np[0][2]; np[3][2] = np[0][1];
   // creating new subelements
   for(i=0;i<4;i++) {
   	TPZManVector<int64_t> npvec(3);
	for (j = 0; j < 3; j++) npvec[j] = np[i][j];
      fSubEl[i] = CreateGeoEl(npvec,matid,*Mesh());
   }

   SubElVec.Resize(4);
   for(sub=0;sub<4;sub++) {
      SubElVec[sub] = fSubEl[sub];
      SubElVec[sub]->SetFather(this);
   }
/*   for(i=0;i<3;i++) {
     fSubEl[i]->SetNeighbour((i+1)%3,TPZGeoElSide(fSubEl[3],i));
     fSubEl[3]->SetNeighbour(i,TPZGeoElSide(fSubEl[(i+1)%3],i));
     fSubEl[(i+1)%3]->SetNeighbour(i,TPZGeoElSide(fSubEl[i],(i+1)%3));
   }  */
   fSubEl[0]->SetNeighbour(1,TPZGeoElSide(fSubEl[1],0));
   fSubEl[1]->SetNeighbour(0,TPZGeoElSide(fSubEl[3],2));
   fSubEl[3]->SetNeighbour(2,TPZGeoElSide(fSubEl[0],1));//3/2->0/1

   fSubEl[1]->SetNeighbour(2,TPZGeoElSide(fSubEl[2],1));
   fSubEl[2]->SetNeighbour(1,TPZGeoElSide(fSubEl[3],0));
   fSubEl[3]->SetNeighbour(0,TPZGeoElSide(fSubEl[1],2));//3/0->1/2

   fSubEl[2]->SetNeighbour(0,TPZGeoElSide(fSubEl[0],2));
   fSubEl[0]->SetNeighbour(2,TPZGeoElSide(fSubEl[3],1));
   fSubEl[3]->SetNeighbour(1,TPZGeoElSide(fSubEl[2],0));//3/1->2/0
   //conectiv do medio dos lados dos subs
   fSubEl[0]->SetNeighbour(4,TPZGeoElSide(fSubEl[3],4));
   fSubEl[3]->SetNeighbour(4,TPZGeoElSide(fSubEl[0],4));//3/4->0/4
   fSubEl[1]->SetNeighbour(5,TPZGeoElSide(fSubEl[3],5));
   fSubEl[3]->SetNeighbour(5,TPZGeoElSide(fSubEl[1],5));//3/5->1/5
   fSubEl[2]->SetNeighbour(3,TPZGeoElSide(fSubEl[3],3));
   fSubEl[3]->SetNeighbour(3,TPZGeoElSide(fSubEl[2],3));//3/3->2/3
   //procura-se um viz pela face que seja dividido
   TPZGeoElSide dividedneighbour = Neighbour(6);
   if(dividedneighbour.Exists()) {
     while(dividedneighbour.Element() != this) {
       if(dividedneighbour.HasSubElement()) break;
       dividedneighbour = dividedneighbour.Neighbour();
     }
     if(dividedneighbour.Exists() && dividedneighbour.Element() != this) {
   // we found a neighbour for everybody!!
       TPZManVector<int64_t> nodes(3);
       TPZStack<TPZGeoElSide> neighsub;
       nodes[0] = NodeIndex(0);
       nodes[1] = NodeIndex(1);
       nodes[2] = NodeIndex(2);
       dividedneighbour.GetSubElements2(neighsub);

       nodes.Resize(1);
       for(i=0; i<3; i++) {
         nodes[0] = NodeIndex(i);
         TPZGeoEl *gel = neighsub[i].Element();
         int locside = gel->WhichSide(nodes);
         TPZGeoElSide sub(fSubEl[i],i);
         sub.SetConnectivity(TPZGeoElSide(neighsub[i].Element(),locside));
         int iplus = (i+1)%3;
         nodes[0] = fSubEl[i]->NodeIndex(iplus);
         locside = neighsub[i].Element()->WhichSide(nodes);
         sub = TPZGeoElSide(fSubEl[i],iplus);
         sub.SetConnectivity(TPZGeoElSide(neighsub[i].Element(),locside));
       }
       nodes.Resize(2);
       for(i=0; i<3; i++) {
      	for(j=3; j<6; j++) {
           nodes[0] = fSubEl[i]->SideNodeIndex(j,0);
           nodes[1] = fSubEl[i]->SideNodeIndex(j,1);
           int locside = neighsub[i].Element()->WhichSide(nodes);
           TPZGeoElSide sub(fSubEl[i],j);
           sub.SetConnectivity(TPZGeoElSide(neighsub[i].Element(),locside));
         }
       }
       for(i=0; i<4; i++) {
      	TPZGeoElSide sub(fSubEl[i],6);
         sub.SetConnectivity(neighsub[i]);
       }
   	 return;
     }
   }
   for(i=3; i<6; i++) {
   	dividedneighbour = Neighbour(i);
      if(!dividedneighbour.Exists()) continue;
      while(dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Element() != this) {
      // we found a divided neighbour for side i
         int iplus = (i+1)%3;
         TPZManVector<int64_t> nodes(2);
         TPZStack<TPZGeoElSide> neighsub;
         nodes[0] = SideNodeIndex(i,0);
         nodes[1] = SideNodeIndex(i,1);
         dividedneighbour.GetSubElements2(neighsub);
         TPZGeoElSide sub(fSubEl[i-3],i);
         sub.SetConnectivity(neighsub[0]);
         sub = TPZGeoElSide(fSubEl[iplus],i);
         sub.SetConnectivity(neighsub[1]);
         nodes.Resize(1);
         nodes[0] = fSubEl[i-3]->SideNodeIndex(iplus,0);
         int locside = neighsub[0].Element()->WhichSide(nodes);
         sub = TPZGeoElSide(fSubEl[i-3],iplus);
         sub.SetConnectivity(TPZGeoElSide(neighsub[0].Element(),locside));
      }
   }
   for(i=0; i<3; i++) {
   	dividedneighbour = Neighbour(i);
      if(!dividedneighbour.Exists()) continue;
      while(dividedneighbour.Element() != this) {
         if(dividedneighbour.HasSubElement()) break;
         dividedneighbour = dividedneighbour.Neighbour();
      }
      if(dividedneighbour.Element() != this) {
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

int TPZGeoElT2d::NSubElements() {
  return 4;
}

int TPZGeoElT2d::NSideSubElements(int side) {
  if(side < 0 || side > 6) {
    PZError << "TPZGeoElT2d::NSideSubElements called for side " << side <<std::endl;
    return 0;
  }
  if(side==6) return 4;
  if(side<3) return 1;
  return 2;
}

void TPZGeoElT2d::SideSubElements(int side,TPZVec<TPZGeoEl *> &sub) {
  if(!fSubEl[0]) {
    sub.Resize(0);
    return;
  }
  if(side < 0 || side > 6) {
    PZError << "TPZGeoElT2d::SideSubElements called for side " << side <<std::endl;
    return;
  }
  if(side==6) {//6
    for(int i=0;i<4;i++) sub[i] = fSubEl[i];
    return;
  }
  //0,1,2
  if(side<3) {
    sub[0]=fSubEl[side];
    return;
  }
  //3,4,5
  side-=3;
  sub.Resize(2);
  sub[0] = fSubEl[side];
  sub[1] = fSubEl[(side+1)%3];
}

TPZGeoElSide TPZGeoElT2d::SideSubElement(int side,int position) {
   if (position<0 ||position>3 || side <0 ||side>6) {
   	PZError << "TPZGeoElT2d::SideSubElement called with position " << position << " side " << side <<std::endl;
      return TPZGeoElSide();
   }
   if(side==6) return TPZGeoElSide(SubElement(position),6);
   if(side<3) {
      if(position!=0) {
         PZError << "TPZGeoElT2d::SideSubElement called with position " << position << " side " << side <<std::endl;
         return TPZGeoElSide();
      } else {
         return TPZGeoElSide(SubElement(side),side);
      }
   }
   if(position==0) return TPZGeoElSide(SubElement(side-3),side);
   else if(position == 1) return TPZGeoElSide(SubElement((side-2)%3),side);
   PZError << "TPZGeoElT2d::SideSubElement called with position " << position << " side " << side <<std::endl;
   return TPZGeoElSide();
}

void TPZGeoElT2d::GetSubElement(int side,TPZVec<int> &refnode,TPZVec<TPZGeoElSide> &sub) {
  int nsub = NSideSubElements(side);
  sub.Resize(nsub);
  if(!nsub) return;
  if(side<3 && refnode[0]==NodeIndex(side)) {
    sub[0]=TPZGeoElSide(fSubEl[side],side);
    return;
  }
//  if(nsub==1) {  // then side=0 or =1 or =2
//    sub[0]=TPZGeoElSide(fSubEl[side],side);
//    return;
//  }

  int i,j,k,nref=refnode.NElements();
  for(i=0;i<nsub;i++) {
    TPZGeoElSide sidesub = SideSubElement(side,i);
    TPZGeoEl *subel = sidesub.Element();
    for(k=0; k<nref; k++) {
      for(j=0;j<3;j++) {
         if(subel->NodeIndex(j)==refnode[k]) {
            sub[k] = SideSubElement(side,i);
         }
      }
    }
  }
  if(side==6) sub[3] = SideSubElement(side,3);
  return;
}

TPZGeoElSide TPZGeoElT2d::Father(int side) {
	TPZGeoEl *fFather = TPZGeoEl::Father();
   if(!fFather) return TPZGeoElSide();

   int whichsub = -1;
   int i;
   for(i=0; i<4; i++) if(fFather->SubElement(i) == this) whichsub = i;
   if(whichsub == -1) {//equivale a is = 4  ou is > 3
	   PZError << "TPZGeoElT2d::Father. fFather isn't father of this element.\n";
   	return TPZGeoElSide();
   }
   //agora o atual elemento é o filho numero whichsub < 4
   if(whichsub == 3 && side != 6) return TPZGeoElSide();
   if(whichsub == side || side==6) return TPZGeoElSide(fFather,side);//side = 0,1,2
	if(whichsub != side && side<3) return TPZGeoElSide();
   //side = 3,4,5
 	if(whichsub == 0 && side!=4) return TPZGeoElSide(fFather,side);//ou side==3 || side==5
 	if(whichsub == 1 && side!=5) return TPZGeoElSide(fFather,side);//ou side==3 || side==4
 	if(whichsub == 2 && side!=3) return TPZGeoElSide(fFather,side);//ou side==4 || side==5
 	//if(wichsub == 3) return TPZGeoElSide();//é feito pelo seguinte caso

//   PZError << "TPZGeoElT2d::Father. fFather isn't father of this element along the given side.\n";
   return TPZGeoElSide();//inclui os outros casos
}

/**accumulates the transformation of the jacobian which maps the current
   master element space into the space of the master element of the father*/
void TPZGeoElT2d::BuildTransform(int side,TPZGeoEl *father,TPZTransform<STATE> &t) {
  if(this == father) return;
  TPZGeoEl *fFather = TPZGeoEl::Father();

   if(!fFather) {
	   std::cout << "TPZGeoElT2d::BuildTransform called for inconsistent parameters\n";
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
	   locfather = fFather;
   }
   if(!locfather) {
	   std::cout << "TPZGeoElT2d::BuildTransform could not identify the father element\n";
	   return;
   }

   if(side == 6) {
      mult(0,0) = 0.5;
      mult(1,1) = 0.5;
      //  if(typediv==1) {
      switch(whichsub) {
         case 0:
            break;
         case 1:
            sum(0,0) = 0.5;
            break;
         case 2:
            sum(1,0) = 0.5;
            break;
         case 3:
            mult(0,0) = -0.5;
            mult(0,1) = 0.;
            mult(1,0) = 0.;
            mult(1,1) = -0.5;
            sum(0,0) = 0.5;
            sum(1,0) = 0.5;
            break;
      }
   } else if(side <3) {
   	return;
   } else {
   	mult(0,0) = 0.5;
      if(whichsub == side-3) sum(0,0) = -0.5;
      else sum(0,0) = 0.5;
   }

   tloc.SetMatrix(mult,sum);
   t = tloc.Multiply(t);
   if(locfather != father) locfather->BuildTransform2(side,father,t);
}

TPZTransform<STATE> TPZGeoElT2d::SideToSideTransform(int sidefrom,int sideto) {
   if(sideto != 6 && (sidefrom > 2 || sidefrom < 0)) {
      PZError << "TPZGeoElT2d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   int dimfrom = SideDimension(sidefrom);
   int dimto = SideDimension(sideto);
   if(dimto <= dimfrom){
      PZError << "TPZGeoElT2d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
      return TPZTransform<STATE>(0,0);
   }
   TPZTransform<STATE> trans(dimto,dimfrom);
   switch(sidefrom) {
      case 3:
         trans.Mult()(1,0) = 0.;
         trans.Mult()(0,0) = .5;
         trans.Sum()(0,0) = .5;
         trans.Sum()(1,0) = 0.;
         break;
      case 4:
         trans.Mult()(0,0) = -.5;
         trans.Mult()(1,0) = .5;
         trans.Sum()(0,0) = .5;
         trans.Sum()(1,0) = .5;
         break;
      case 5:
         trans.Mult()(0,0) = 0.;
         trans.Mult()(1,0) = -.5;
         trans.Sum()(0,0) = 0.;
         trans.Sum()(1,0) = .5;
         break;
      case 0:
      	if(sideto == 3) trans.Sum()(0,0) = -1.;
         else if(sideto == 5) trans.Sum()(0,0) = 1.;
         else if(sideto == 6) {
         	trans.Sum()(0,0) = 0.;
            trans.Sum()(1,0) = 0.;
         } else {
            PZError << "TPZGeoElT2d:SideToSideTransform sidefrom = " << sidefrom <<
            " sideto = " << sideto <<std::endl;
            return trans;
         }
         break;
      case 1:
      	if(sideto == 3) trans.Sum()(0,0) = 1.;
         else if(sideto == 4) trans.Sum()(0,0) = -1.;
         else if(sideto == 6) {
         	trans.Sum()(0,0) = 1.;
            trans.Sum()(1,0) = 0.;
         } else {
            PZError << "TPZGeoElT2d:SideToSideTransform sidefrom = " << sidefrom <<
            " sideto = " << sideto <<std::endl;
            return trans;
         }
         break;
      case 2:
      	if(sideto == 4) trans.Sum()(0,0) = 1.;
         else if(sideto == 5) trans.Sum()(0,0) = -1.;
         else if(sideto == 6) {
         	trans.Sum()(0,0) = 0.;
            trans.Sum()(1,0) = 1.;
         } else {
            PZError << "TPZGeoElT2d:SideToSideTransform sidefrom = " << sidefrom <<
            " sideto = " << sideto <<std::endl;
            return trans;
         }
         break;
   }
   return trans;
}

TPZCompEl *TPZGeoElT2d::CreateBCCompEl(int side,int bc,TPZCompMesh &cmesh) {
  if(side==6) {
    TPZManVector<int64_t> nodes(3);
	for (int ii = 0; ii < 3; ii++) nodes[ii] = fNodeIndexes[ii];
    TPZGeoElT2d *gel = CreateGeoEl(nodes,bc,*Mesh());
    TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,0));
    TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(this,1));
    TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(this,2));
    TPZGeoElSide(gel,3).SetConnectivity(TPZGeoElSide(this,3));
    TPZGeoElSide(gel,4).SetConnectivity(TPZGeoElSide(this,4));
    TPZGeoElSide(gel,5).SetConnectivity(TPZGeoElSide(this,5));
    TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(this,6));
	int64_t index;
    return gel->CreateCompEl(cmesh,index);
  }
  else if(side>-1 && side<3) {
    TPZBndCond *bcptr = (TPZBndCond *) cmesh.FindMaterial(bc);
    if(!bcptr) {
       PZError << "TPZGeoElT2d::CreateBCCompEl has no bc.\n";
       return 0;
    }
    TPZCompEl *cel = Reference();
    if(!cel) {
       PZError << "TPZGeoElT2d::CreateBCCompEl has no computational element\n";
       return 0;
    }
    // Philippe 11/8/2000 update the way nodal boundary condition are created
      TPZGeoElPoint *gel;
      TPZManVector<int64_t> node(1);
	  node[0] = fNodeIndexes[side];
	  gel = new TPZGeoElPoint(node,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,side));
      int64_t index;
      return gel->CreateCompEl(cmesh,index);
  }
  else if(side > 2 && side < 6) {
    TPZManVector<int64_t> nodes(2);
    nodes[0] = NodeIndex(side-3);
    nodes[1] = NodeIndex((side-2)%3);
    TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*Mesh());
    TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,side-3));
    TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(this,(side-2)%3));
    TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(this,side));
	int64_t index;
    return gel->CreateCompEl(cmesh,index);
  }
  else PZError << "TPZGeoElT2d::CreateBCCompEl. Side = " << side <<std::endl;
  return 0;
}

TPZGeoElSide TPZGeoElT2d::Father2(int side){//Augusto:09/01/01
	TPZGeoEl *fFather = TPZGeoEl::Father();

	if(side<0 || side>6 || fFather==0){
		PZError << "TPZGeoET2d::Father2 called error" <<std::endl;
	}
	if(fFather->SubElement(0)==this){
		if(side==0) return TPZGeoElSide(fFather,0);
		if(side==1 || side==3)	return TPZGeoElSide(fFather,3);
		if(side==4 || side==6)	return TPZGeoElSide(fFather,6);
		/*if(side==2 || side==5)*/	return TPZGeoElSide(fFather,5);
	}
	if(fFather->SubElement(1)==this){
		if(side==1) return TPZGeoElSide(fFather,1);
		if(side==0 || side==3)	return TPZGeoElSide(fFather,3);
		if(side==5 || side==6)	return TPZGeoElSide(fFather,6);
		/*if(side==2 || side==4)*/	return TPZGeoElSide(fFather,4);
	}
	if(fFather->SubElement(2)==this){
		if(side==2) return TPZGeoElSide(fFather,2);
		if(side==0 || side==5)	return TPZGeoElSide(fFather,5);
		if(side==3 || side==6)	return TPZGeoElSide(fFather,6);
		/*if(side==1 || side==4)*/	return TPZGeoElSide(fFather,4);
	}
	if(fFather->SubElement(3)==this){
		if(side==0) return TPZGeoElSide(fFather,4);
		if(side==1) return TPZGeoElSide(fFather,5);
		if(side==2) return TPZGeoElSide(fFather,3);
		return TPZGeoElSide(fFather,6);
	}
    return TPZGeoElSide();
}

static int subeldata[7][7][2] =
{
	{{0,0}},/*side=0 {isub0{0,1},isub1{0,1},isub2{0,1},...}*/
	{{1,1}},/*side=1*/
	{{2,2}},/*side=2*/
	{{0,3},{0,1},{1,3}},/*side=3*/
	{{1,4},{1,2},{2,4}},/*side=4*/
	{{2,5},{2,0},{0,5}},/*side=5*/
	{{0,6},{1,6},{2,6},{3,6},{2,3},{0,4},{1,5}}/*side=6*/
};

static int nsubeldata[7] = {1,1,1,3,3,3,7};


void TPZGeoElT2d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){//Augusto:09/01/01

  subel.Resize(0);
	if(side<0 || side>6 || !HasSubElement(side)){
		PZError << "TPZGeoElT2d::GetSubElements2 called error" <<std::endl;
		return;
	}
	int nsub = nsubeldata[side];
	int isub;
	for(isub=0; isub<nsub; isub++) subel.Push(TPZGeoElSide(fSubEl[subeldata[side][isub][0]],
		                                                          subeldata[side][isub][1]));
}

static REAL buildt[4][7][3][2] = {//por colunas
/*S0*/
     {/*0*/{{0.,0.},{0.,0.},{0.,0.}},
      /*1*/{{0.,0.},{0.,0.},{0.,0.}},
      /*2*/{{0.,0.},{0.,0.},{0.,0.}},
      /*3*/{{.5,0.},{0.,0.},{-.5,0.}},
      /*4*/{{-.25,.25},{0.,0.},{.25,.25}},
      /*5*/{{.5,0.},{0.,0.},{.5,0.}},
      /*6*/{{.5,0.},{0.,.5},{0.,0.}}},
/*S1*/
     {/*0*/{{0.,0.},{0.,0.},{0.,0.}},
      /*1*/{{0.,0.},{0.,0.},{0.,1.}},
      /*2*/{{0.,0.},{0.,0.},{0.,0.}},
      /*3*/{{.5,0.},{0.,0.},{.5,0.}},
      /*4*/{{.5,0.},{0.,0.},{-.5,0.}},
      /*5*/{{0.,-.25},{0.,0.},{.5,.25}},
      /*6*/{{.5,0.},{0.,.5},{.5,0.}}},
/*S2*/
     {/*0*/{{0.,0.},{0.,0.},{0.,0.}},
      /*1*/{{0.,0.},{0.,0.},{0.,0.}},
      /*2*/{{0.,0.},{0.,0.},{0.,1.}},
      /*3*/{{.25,0.},{0.,0.},{.25,.5}},
      /*4*/{{.5,0.},{0.,0.},{.5,0.}},
      /*5*/{{.5,0.},{0.,0.},{-.5,0.}},
      /*6*/{{.5,0.},{0.,.5},{0.,.5}}},
/*S3*/
     {/*0*/{{0.,0.},{0.,0.},{0.,0.}},
      /*1*/{{0.,0.},{0.,0.},{0.,0.}},
      /*2*/{{0.,0.},{0.,0.},{0.,0.}},
      /*3*/{{-.25,0.},{0.,0.},{.25,.5}},
      /*4*/{{.25,-.25},{0.,0.},{.25,.25}},
      /*5*/{{0.,.25},{0.,0.},{.5,.25}},
      /*6*/{{-.5,0.},{0.,-.5},{.5,.5}}}
};

TPZTransform<STATE> TPZGeoElT2d::BuildTransform2(int side, TPZGeoEl * /*father*/){//Augusto:09/01/01
	TPZGeoEl *fFather = TPZGeoEl::Father();

	if(side<0 || side>7 || !fFather){
  	PZError << "TPZGeoElT2d::BuildTransform2 side out of range or father null\n";
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

static REAL MidSideNode[7][3] = {
/*00*/{.0,0.},/*01*/{1.0,.0},/*02*/{0.,1.0},
/*03*/{.5,0.},/*04*/{0.5,.5},/*05*/{0.,0.5},
/*06*/{ 1./3.,1./3.} };

int TPZGeoElT2d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3), pss(3), pf(3), pfs(3);
  //point son, point side son, point father, point side father : elemento mestre
  pss[2] = 0.;//2d
  pfs[2] = 0.;
  pf[2] = 0.;
  for(sn=0;sn<4;sn++){
    TPZGeoEl *son = subs[sn];
    for(sd=0;sd<7;sd++){
      ps[0] = MidSideNode[sd][0];//element
      ps[1] = MidSideNode[sd][1];//master point
      ps[2] = MidSideNode[sd][2];// = 0
      if(son->WhichSide(ps) != sd) std::cout << "Lado nao bate\n";
      TPZTransform<STATE> telsd = pzshape::TPZShapeTriang::TransformElementToSide(sd);//2x2
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform<STATE> t;
	  son->BuildTransform2(sd, gel,t);
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = pzshape::TPZShapeTriang::TransformSideToElement(sdfat);//2x2
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(6).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao furada\n";
        PZError << "son    = " << (son->Id()) <<std::endl;
        PZError << "father = " << ((son->Father2(6).Element())->Id()) <<std::endl;
        PZError << "side   = " << sd <<std::endl <<std::endl;
        int ok;
		std::cin >> ok;
      } else {
		  std::cout << "Transformacao OK!\n";
		  std::cout << "Filho/lado : " << son->Id() << "/" << sd <<std::endl;
		  std::cout << "Pai : " << son->Father2(6).Element()->Id() <<std::endl <<std::endl;
      }
    }
  }
  return 1;
}

/** Compute the measure of the geometrical element - Jorge 17/7/99*/
REAL TPZGeoElT2d::Mesure(int dim) {
  if(dim!=2) return 0.;
  REAL fMesure = Volume();
  if(IsZero(fMesure)) {
    TPZGeoNode &nod1 = Mesh()->NodeVec()[fNodeIndexes[0]];
    REAL xx, yy, x0 = nod1.Coord(0), y0 = nod1.Coord(1);
    TPZGeoNode &nod2 = Mesh()->NodeVec()[fNodeIndexes[1]];
    xx = x0 - nod2.Coord(0);
    yy = nod2.Coord(1) - y0;
    TPZGeoNode &nod3 = Mesh()->NodeVec()[fNodeIndexes[2]];
    xx *= (nod3.Coord(1) - y0);
    yy *= (nod3.Coord(0) - x0);
    fMesure = 0.5 * fabs(xx + yy);
  }
  return fMesure;
}

void TPZGeoElT2d::Center(TPZVec<REAL> &center) {
  center[0] = center[1] = 1./3.;
}

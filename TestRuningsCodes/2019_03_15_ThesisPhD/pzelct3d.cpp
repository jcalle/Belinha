#include "pzelct3d.h"
#include "pzelgt3d.h"
#include "pzmatrix.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzerror.h"
#include "pzconnect.h"
#include "pzgraphel.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include <math.h>
#include <stdio.h>
#include "TPZMaterial.h"


//TShortMatrix TPZCompElT3d::gEqNumbers;
                                       //linha i � filho i, {a,b,c} = {lado do filho atual,irm�o vizinho,lado do vizinho}
int TPZCompElT3d::InNeigh[6][16][3] = { {{1,1,0},{2,4,3},{3,4,4},{5,5,6},{8,4,9},{9,4,12},{12,4,17},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
												     {{0,5,1},{2,2,1},{3,3,1},{6,5,10},{7,5,5},{9,5,9},{13,5,14},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                         {{0,0,2},{1,5,4},{3,3,2},{4,5,11},{7,4,7},{8,5,12},{11,5,16},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
                                         {{0,0,3},{1,4,1},{2,4,2},{4,4,10},{5,4,6},{6,4,11},{10,4,15},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
												     {{0,0,1},{1,5,0},{2,5,3},{3,5,2},{4,3,0},{5,1,7},{6,5,8},{7,5,7},{8,0,5},{9,0,8},{10,3,4},{11,3,6},{12,0,9},{13,5,13},{15,3,10},{17,0,12}},
                                         {{0,1,3},{1,4,0},{2,2,0},{3,2,3},{4,1,2},{5,4,5},{6,4,8},{7,2,7},{8,3,5},{9,1,9},{10,1,6},{11,2,4},{12,2,8},{13,4,13},{14,1,13},{16,2,11}} };

int TPZCompElT3d::CornerSons[6][5] = { {0,4,6,7,-1},{4,1,5,8,-1},{6,5,2,9,-1},{7,8,9,3,-1},{4,8,9,6,7},{8,4,6,9,5} };

int TPZCompElT3d::FaceSons[4][4] = { {0,1,2,5},{0,1,3,4},{1,2,3,5},{0,2,3,4} };

int TPZCompElT3d::MiddleFace[4] = {15,14,17,16 };//face do elemento associada a face do subelemento tetraedro contida nela

int TPZCompElT3d::FaceNodes[4][3]  = { {0,1,2},{0,1,3},{1,2,3},{0,2,3} };

int TPZCompElT3d::SideNodes[6][2]  = { {0,1},{1,2},{2,0},{0,3},{1,3},{2,3} };
//n�s medios dos lados como pertencentes a algum filho : (filho,n�)
int TPZCompElT3d::MidSideNodes[6][2]  = { {0,1},{1,2},{2,0},{3,0},{3,1},{3,2} };

int TPZCompElT3d::FaceSides[4][3] = { {4,5,6},{4,8,7},{5,9,8},{6,9,7} };

REAL TPZCompElT3d::MidCoord[6][3] = { {.5,0.,0.},{.5,.5,0.},{0.,.5,0.},{0.,0.,.5},{.5,0.,.5},{0.,.5,.5} };

REAL TPZCompElT3d::MasterCoord[4][3] = { {0.,0.,0.,},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.} };

int TPZCompElT3d::ShapeFaceId[4][3] = { {0,1,2},{0,1,3},{1,2,3},{0,2,3} };

int TPZCompElT3d::FaceConnectLocId[4][7] = { {0,1,2,4,5,6,10},{0,1,3,4,8,7,11},
									                  {1,2,3,5,9,8,12},{0,2,3,6,9,7,13} };

TPZCompElT3d::TPZCompElT3d(TPZCompMesh &mesh,TPZGeoElT3d *ref, int64_t &index) :
	TPZIntelGen<pzshape::TPZShapeTetra>(mesh,ref,index), fIntRule(2) {
  int i;                                     //dois pontos por eixo
  for(i=0;i<11;i++) {
    fSideOrder[i] = mesh.GetDefaultOrder();
    fPreferredSideOrder[i] = mesh.GetDefaultOrder();
  }
  ref->SetReference(this);
  for(i=0;i<4;i++) {//4 cantos
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
  }
  for(;i<15;i++) {//6 lados e 4 faces 1 centro
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
	 IdentifySideOrder(i);
  }
  TPZVec<int> order(3,2*fSideOrder[10]);//integra variavel ao quadrado em cada dire��o
  fIntRule.SetOrder(order);
}

TPZCompElT3d::TPZCompElT3d(TPZCompMesh &mesh,TPZGeoElT3d *ref, int64_t &index,int /*noconnects*/) :
	TPZIntelGen<pzshape::TPZShapeTetra>(mesh,ref,index), fIntRule(2) {
  int i;
  for(i=0;i<11;i++) {
    fSideOrder[i] = mesh.GetDefaultOrder();
    fPreferredSideOrder[i] = mesh.GetDefaultOrder();
  }
}

void TPZCompElT3d::Shape(TPZVec<REAL> &x, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &dphi) {
  TPZVec<int64_t> id(4);
  for(int i=0;i<4;i++) id[i] = Reference()->NodeIndex(i);
  TPZManVector<int> ord(11);
  for (int j = 0; j < 11; j++) ord[j] = fSideOrder[j];
  pzshape::TPZShapeTetra::Shape(x,id,ord,phi,dphi);
}

int TPZCompElT3d::NConnectShapeF(int side) {
   if(side<4) return 1;//0 a 3
   int s = side-4;
   if(side<10) return fSideOrder[s]-1;//4 a 9
   if(side<14) {//10 a 13
   	int sum = 0;
      for(int i=0;i<fSideOrder[s]-1;i++) sum += i;
   	return sum;
   }
   if(side==14) {
   	int totsum = 0,sum;
      for(int i=1;i<fSideOrder[s]-2;i++) {
         sum = i*(i+1) / 2;
         totsum += sum;
      }
      return totsum;
   }
   PZError << "TPZCompElT3d::NConnectShapeF, bad parameter side " << side <<std::endl;
   return 0;
}

int TPZCompElT3d::NSideConnects(int side) {
	if(side<0)   return -1;
	if(side<4)   return 1;//cantos : 0 a 3
   if(side<10)  return 3;//lados : 4 a 9
   if(side<14)	 return 7;//faces : 10 a 13
   if(side==14) return 15;//centro : 15
   return -1;
}

int TPZCompElT3d::SideConnectLocId(int node, int side) {
	if(side<0 || side>15) return -1;
   if(side<4) {
   	if(node==0) return side;
   } else
   if(side<7) {//4,5,6
	   int s = side-4;//0,1,2
   	if(!node) return s;//0,1,2
      if(node==1) return (s+1)%3;//1,2,0
      if(node==2) return side;//4,5,6
   } else
   if(side<10) {//7,8,9
   	int s = side-7;//0,1,2
   	if(!node) return s;//0,1,2
      if(node==1) return 3;//4,4,4
      if(node==2) return side;//7,8,9
   } else
   if(side<14) {//10 a 13
   	int s = side-10;
   	if(node<7) return FaceConnectLocId[s][node];
   } else
   if(side==14 && node<15){
   	return node;
   }
	PZError << "TPZCompElT3d::SideConnectLocId called for node = "
   	      << node << " and side = " << side << "\n";
   return -1;
}

int TPZCompElT3d::SideOrder(int side) {
	//0 <= side <= 15
  if(side<4 || side>15) return 0;//cantos ou side ruim
  return fSideOrder[side-4];//ordens dos lados e faces
                            //side = 4,5,6,7,8,9,10,11,12,13
}

void TPZCompElT3d::FaceOrder(int face,int &ord1,int &ord2) {

 if (face<0 || face>3) std::cout << "\nTPZCompElT3d::FaceOrder called whit face = " << face <<std::endl;
 //ordem minima para cada dire��o do plano da face
 ord2 = 0;
 if (face == 0 || face == 10) {//10
	 ord1 = (fSideOrder[0] <= fSideOrder[1]) ? fSideOrder[0] : fSideOrder[1];
    ord1 = (   ord1       <= fSideOrder[2]) ?     ord1      : fSideOrder[2];
 } else
 if (face == 1 || face == 11) {//11
	 ord1 = (fSideOrder[0] <= fSideOrder[3]) ? fSideOrder[0] : fSideOrder[3];
    ord1 = (   ord1       <= fSideOrder[4]) ?     ord1      : fSideOrder[4];
 } else
 if (face == 2 || face == 12) {//12
	 ord1 = (fSideOrder[1] <= fSideOrder[4]) ? fSideOrder[1] : fSideOrder[4];
   ord1 = (    ord1       <= fSideOrder[5]) ?     ord1      : fSideOrder[5];
 } else
 if (face == 3 || face == 13) {//13
	 ord1 = (fSideOrder[2] <= fSideOrder[3]) ? fSideOrder[2] : fSideOrder[3];
    ord1 = (   ord1       <= fSideOrder[5]) ?     ord1      : fSideOrder[5];
 }
}

/*void TPZCompElT3d::SetFaceOrder(int face, int order) {
	fFaceOrder[face] = order;//bota ordem de interpola��o de face
}  */

void TPZCompElT3d::Print(std::ostream &out) {
	TPZInterpolatedElement::Print(out);
/*   return;
	out << "TPZCompElT3d element:\n";
   if(fMaterial) {
		out << "material number = " << fMaterial->Id() << "\n";
   } else {
   	out << "no material defined\n";
   }
	out << "number of dof nodes = " << 15 << "\n"
		<< "nodes : ";
	for (int nod=0; nod<15; nod++) {
		out << "\nNode number = " << nod << " number of shape functions = "
				<< NConnectShapeF(nod) << "\n";
		Mesh()->ConnectVec()[fConnectIndexes[nod]].Print(*Mesh(),out);
	}
	out << "\nlocal equation matrix = \n";  */

	//EqNumber(gEqNumbers);
	//gEqNumbers.Print("",out);
}

// compare the level of two elements
short TPZCompElT3d::CompareLevel(TPZCompElT3d &element,TPZCompElT3d &neighbour)	{
		// if the level of elements is equal return 1, else return 0
    if (   ( (element.Reference()) -> Level() )	==
   		  	  ( (neighbour.Reference()) -> Level() )   )	return 1;
    else return 0;
}

void TPZCompElT3d::NormalVector(int side,TPZVec<REAL> &int_point,
										TPZVec<REAL> &normal, TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &norm_r3){
	((TPZGeoElT3d *)Reference())->NormalVector(side,int_point,normal,axes,norm_r3);
}

void TPZCompElT3d::EqNumber(TPZShortMatrix &/*mat*/) {

/*	mat.Redim(SideOrder(0),SideOrder(1));//ordem x e ordem y da face 0
	mat *= 0;
	mat -= 1;
	mat(0,0) = 0;
	mat(1,0) = 1;
	mat(1,1) = 2;
	mat(0,1) = 3;
	int ieq = 4;
	for(int is =0; is < 4; is++) {
		int OrSide = SideOrder(is);
		for(int i=2 ; i < OrSide; i++) {
			switch (is) {
				case 0 : mat(i,0) = ieq;
				break;
				case 1 : mat(1,i) = ieq;
				break;
				case 2 : mat(i,1) = ieq;
				break;
				case 3 : mat(0,i) = ieq;
				break;
				default :
				break;
			}
			ieq++;
		}
	}
	for(int i=2;i<fSideOrder[0];i++) {
		for(int j=2;j<fSideOrder[1];j++) {
			mat(i,j) = ieq++;
		}
	} */
}

void TPZCompElT3d::VarRange(int var,REAL &min,REAL &max) {
	PZError << "TPZCompElT3d::VarRange is not defined.\n";
   if(var>-1) max = min = 0.;
}

void TPZCompElT3d::Load() {
	PZError << "TPZCompElT3d::Load is called.\n";
}

void TPZCompElT3d::SetConnectIndex(int i,int connectindex) {
  if(i>-1 && i<15) fConnectIndexes[i] = connectindex;
  else {
    PZError << "TPZCompElT3d::SetConnectIndex. Bad parameter i.\n";
    PZError.flush();
  }
}

int TPZCompElT3d::ConnectIndex(int i) {
	if(i<0 || i>15) {
   	PZError << "TCompElT2d::ConnectIndex. Bad parameter i.\n";
      return -1;
   }
   return fConnectIndexes[i];
}

void TPZCompElT3d::SetInterpolationOrder(TPZVec<int> &ord) {
  if(ord.NElements()!=11)
  	PZError << "TPZCompElT3d::SetInterpolationOrder. ord has bad number of elements.\n";
  for(int i=0;i<11;i++) fPreferredSideOrder[i] = ord[i];
}

void TPZCompElT3d::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(11);
  for(int i=0;i<11;i++) ord[i] = fSideOrder[i];
}

TPZIntPoints *TPZCompElT3d::CreateSideIntegrationRule(int side) {
   if(side<0 || side>15) {
		PZError << "TPZCompElT3d::CreateSideIntegrationRule. bad side number.\n";
   	return new TPZInt1d(0);
   }
   //SideOrder corrige sides de 4 a 14 para 0 a 10
   if(side<4)   return new TPZInt1d(0);//cantos 0 a 3 : cria regra com um ponto
   if(side<10)  return new TPZInt1d(2*SideOrder(side));//lados 4 a 9
   if(side<14)  {//faces : 10 a 13
   	return new TPZIntTriang(2*SideOrder(side));
   }
   if(side==14) {//integra��o do elemento
   	return new TPZIntTetra3D(3*SideOrder(14));
   }
   return new TPZInt1d(0);//TPZInt1Point(1)
}

int TPZCompElT3d::PreferredSideOrder(int side) {
   if(side>-1 && side<4) return 0;//cantos
   if(side<15) return fPreferredSideOrder[side-4];//lados,faces e centro (ou interior)
   PZError << "TPZCompElT3d::PreferredSideOrder called for side = " << side << "\n";
   return 0;
}

void TPZCompElT3d::SetPreferredSideOrder(int side, int order) {
  if(side>-1 && side<4) return;//cantos
  if(side<15) fPreferredSideOrder[side-4] = order;//sides 4 a 14
  PZError << "TPZCompElT3d::SetPreferredSideOrder called for side = " << side << "\n";
  return;
}

void TPZCompElT3d::SetSideOrder(int side, int order) {
	if(side<0 || side>14 || order<1) PZError << "TPZCompElT3d::SetSideOrder. Bad paramenter side.\n";
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

void TPZCompElT3d::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {

   if(side<0 || side>15) PZError << "TPZCompElT3d::SideShapeFunction. Bad paramenter side.\n";
   else if(side==14) Shape(point,phi,dphi);
   else if(side<4) phi(0,0)=1.;
   else if(side<10) {//4 a 9
       TPZVec<int64_t> id(2);
	   TPZManVector<int> faceorder(1);
	   faceorder[0] = SideOrder(side);
	   int s = side-4;                 //n� local
       id[0] = Reference()->NodeIndex(SideNodes[s][0]);
       id[1] = Reference()->NodeIndex(SideNodes[s][1]);//n� global
       pzshape::TPZShapeLinear::Shape(point,id,faceorder,phi,dphi);
   }
   else if(side<14) {//faces 10,11,12,13
       TPZVec<int64_t> id(3);
   	 int face = side-10;
       id[0] = Reference()->NodeIndex(FaceNodes[face][0]);
       id[1] = Reference()->NodeIndex(FaceNodes[face][1]);//n� global
       id[2] = Reference()->NodeIndex(FaceNodes[face][2]);
	   TPZManVector<int> ord(4);
	   ord[0] = fSideOrder[FaceSides[face][0]];//lado
       ord[1] = fSideOrder[FaceSides[face][1]];
       ord[2] = fSideOrder[FaceSides[face][2]];
       ord[3] = fSideOrder[side];//face
       pzshape::TPZShapeTriang::Shape(point,id,ord,phi,dphi);
   }
}

void TPZCompElT3d::SetIntegrationRule(int order) {
	TPZIntTetra3D inttetra(order);
   SetIntegrationRule(inttetra);
}

void TPZCompElT3d::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
  //if(dimension == 3) new TPZGraphEl(this,&grmesh);
}


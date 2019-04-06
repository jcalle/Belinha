/**METHODS DEFINITION FOR CLASS ELEM1D*/

#include "pzelgpoint.h"
#include "pzgnode.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzelcpoint.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include <math.h>


static TPZCompEl *CreateEl(TPZGeoElPoint *gel, TPZCompMesh &mesh, int64_t &index) {
  return new TPZCompElPoint(mesh,gel,index);
}

TPZCompEl *(*TPZGeoElPoint::fp)(TPZGeoElPoint *,TPZCompMesh &, int64_t &) = CreateEl;

TPZGeoElPoint::TPZGeoElPoint(int id,TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh)
  : TPZGeoElRefLess<pzgeom::TPZGeoPoint>(id,nodeindices,matind,mesh) {

  int nnod = nodeindices.NElements();
  if(nnod != 1) {
    PZError << "TPZGeoElPoint Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }
  fNodeIndexes = nodeindices[0];
  fSubEl = 0;

}

TPZGeoElPoint::TPZGeoElPoint( TPZVec<int64_t> &nodeindices, int matind, TPZGeoMesh &mesh)
  : TPZGeoElRefLess<pzgeom::TPZGeoPoint>(nodeindices,matind,mesh) {

  int nnod = nodeindices.NElements();
  if(nnod != 1) {
    PZError << "TPZGeoElPoint Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }
  fNodeIndexes = nodeindices[0];
  fSubEl = 0;
}

TPZGeoElPoint *TPZGeoElPoint::CreateGeoEl(TPZVec<int64_t> &np, int matind, TPZGeoMesh &mesh) {
  return new TPZGeoElPoint(np,matind,mesh);
}

TPZGeoElPoint::~TPZGeoElPoint() {
}

int TPZGeoElPoint::NSubElements() {
	return 0;
}

void TPZGeoElPoint::Jacobian(TPZVec<REAL> &/*fl*/,TPZFMatrix<STATE> &result,TPZFMatrix<STATE> &axes,REAL &detjac,TPZFMatrix<STATE> &jacinv){

   result.Redim(1,1);
   jacinv.Redim(1,1);
   result(0,0) = 1.;
   detjac = 1.;
   jacinv(0,0) = 1.;
   axes.Zero();
   axes.Redim(3,3);
   axes.Zero();
   axes(0,0) = 1.;
   axes(1,1) = 1.;
   axes(2,2) = 1.;
}


void TPZGeoElPoint::X(TPZVec<REAL> & /*par*/, TPZVec<REAL> &result){

  result[0] = NodePtr(0)->Coord(0);
  result[1] = NodePtr(0)->Coord(1);
  result[2] = NodePtr(0)->Coord(2);
}

int TPZGeoElPoint::NSideNodes(int side) {
  if(side == 0) return 1;
  return 0;
}

int64_t TPZGeoElPoint::SideNodeIndex(int side, int /*node*/){
  switch(side) {
     case 0:
       return fNodeIndexes;
     default:
       PZError << "TPZGeoElPoint::SideNodeIndex. Bad parameter side.\n";
     }
  return -1;
}

void TPZGeoElPoint::MidSideNodeIndex(int side, int64_t &index) {
  switch(side) {
  case 0:
    index = fNodeIndexes;
    return;
  default:
    PZError << "TPZGeoElPoint::MidSideNode called for side " << side <<std::endl;
    PZError.flush();
    index = -1;
  }
}

int TPZGeoElPoint::SideDimension(int side) {
  switch(side) {
  case 0:
    return 0;
  default:
    PZError << "TPZGeoElPoint::SideDimension. Bad parameter side.\n";
    return -1;
  }
}

int64_t TPZGeoElPoint::NodeIndex(int node) {
  if(node == 0) return fNodeIndexes;
  return -1;
}

void TPZGeoElPoint::Divide(TPZVec<TPZGeoEl *> &pv) {

  fSubEl = 0;
  pv.Resize(0);

}

void TPZGeoElPoint::CreateNewNodes(int64_t *gnodindex) {

  //int numnod = 1;
	int64_t midsidenodeindex;
  MidSideNodeIndex(2,midsidenodeindex);
  TPZGeoNode *midsidenode;
  if(midsidenodeindex == -1) {
    //first look over all neighbours the middle node        // falta pesquisar nos vizinhos (1-Jorge)
    TPZGeoElSide neigh = Neighbour(2);
    while(neigh.Exists() && this!=neigh.Element()) {
      neigh.Element()->MidSideNodeIndex(neigh.Side(),midsidenodeindex);
      if(midsidenodeindex!=-1) break;
      neigh = neigh.Neighbour();
    }
    if(midsidenodeindex==-1) {
      TPZVec<REAL> gco(3);
      TPZVec<REAL> par(1);
      par[0] = 0.;
      X(par,gco);
      midsidenodeindex = Mesh()->NodeVec().AllocateNewElement();
      midsidenode = &Mesh()->NodeVec()[midsidenodeindex];
      midsidenode->Initialize(gco,*Mesh());
    }
  }
  //gnodindex[0] = NodeIndex(0);
  //gnodindex[2] = midsidenodeindex;
  gnodindex[0] = midsidenodeindex;
  //gnodindex[4] = NodeIndex(numnod-1);

}

void TPZGeoElPoint::NormalVector(int /*side*/, TPZVec<REAL> &/*loc*/, TPZVec<REAL> &normal,
			      TPZFMatrix<STATE> &axes, TPZFMatrix<STATE> &/*jac*/){
    axes(0,0) = 1.;
    axes(1,1) = 1.;
    axes(2,2) = 1.;
    normal[0] = 0.;
    normal[1] = 0.;
    normal[2] = 0.;
}

/**Accumulates the transformation of the jacobian which maps the current
master element space into the space of the master element of the father*/
void TPZGeoElPoint::BuildTransform(int /*side*/, TPZGeoEl * /*father*/, TPZTransform<STATE> &/*t*/) {
	std::cout <<  "\nTPZGeoElPoint::BuildTransform called\n";
}

TPZTransform<STATE> TPZGeoElPoint::SideToSideTransform(int sidefrom, int sideto) {
  if(sideto != 0 && sidefrom !=0) {
    PZError << "TPZGeoElPoint:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
  }
  TPZTransform<STATE> t(1,0);                                                          //   !!!   ???
  t.Sum()(0,0) = .0;
  return t;
}

int TPZGeoElPoint::NSideSubElements(int side) {
  if(side==0) return 0;
  return -1;
}

TPZGeoElSide TPZGeoElPoint::SideSubElement(int side,int /*subel*/) {
  if(side != 0) {
    PZError << "TGeoElQ2d::SideSubElement called for side " << side <<std::endl;
    return TPZGeoElSide(0,0);
  }
  if(!fSubEl) return TPZGeoElSide(0,0);
  return TPZGeoElSide(fSubEl,side);
}

/**Neigh é o lado de um elemento vizinho ao elemento atual El que esta sendo dividido rfnd guarda
   os dois nós globais deste lado no sentido antihorario do elemento El, posições 0 e 1 de rfnd*/
void TPZGeoElPoint::GetSubElement(int side,TPZVec<int> &/*rfndindex*/,TPZVec<TPZGeoElSide> &sub) {
  if(!HasSubElement(side)) {
    sub.Resize(0);
    return;
  }
  sub[0] = TPZGeoElSide();
}



TPZGeoElSide TPZGeoElPoint::Father(int /*side*/) {

  return TPZGeoElSide();
}

void TPZGeoElPoint::LowerDimensionSides(int side,TPZStack<TPZGeoElSide> &smallsides) {

    if(side!=0) return;
    smallsides.Push(TPZGeoElSide(this,0));

}

TPZGeoElSide TPZGeoElPoint::HigherDimensionSides(int /*side*/,int /*targetdimension*/) {
   return TPZGeoElSide();
}

/**Inicializa os coeficientes do par de nós do lado I do elemento de referencia*/
void TPZGeoElPoint::SideMasterCo(int /*side*/,TPZVec<REAL> &IVec,TPZVec<REAL> &JVec) {

  IVec.Fill(0.,0);
  JVec.Fill(0.,0);
}

void TPZGeoElPoint::SideMasterCo(int side,TPZFMatrix<STATE> &coord) {
  if(side == 0) coord(0,0) = 1.;

}

TPZCompEl *TPZGeoElPoint::CreateBCCompEl(int side, int bc, TPZCompMesh &cmesh) {

  if(side==0) {
	  TPZManVector<int64_t> nodeindexes(1);
    nodeindexes[0] = fNodeIndexes;
    TPZGeoElPoint *gel = CreateGeoEl(nodeindexes,bc,*Mesh());
    TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,0));
	int64_t index;
    return gel->CreateCompEl(cmesh,index);
  }
  else PZError << "TPZGeoElPoint::CreateBCCompEl. Side = " << side <<std::endl;
  return 0;
}

TPZGeoElSide TPZGeoElPoint::Father2(int /*side*/){//Augusto:09/01/01
	std::cout << "TPZGeoElPoint::Father2 to be implemented\n";
  return TPZGeoElSide();
}

TPZTransform<STATE> TPZGeoElPoint::BuildTransform2(int /*side*/, TPZGeoEl * /*father*/){//Augusto:09/01/01
	std::cout << "TPZGeoElPoint::BuildTransform2 to be implemented\n";
  return TPZTransform<STATE>(0,0);
}


void TPZGeoElPoint::GetSubElements2(int /*side*/, TPZStack<TPZGeoElSide> &/*subel*/){//Augusto:09/01/01
	std::cout << "TPZGeoElPoint::GetSubElements2 to be implemented\n";
}

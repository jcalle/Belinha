//METHODS DEFINITION FOR CLASS ELEM1D

#include "pzelg1d.h"
#include "pzelgpoint.h"
#include "pzgnode.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzelc1d.h"
#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzelctemp.h"
#include "pzvec.h"
#include <math.h>


static TPZCompEl *CreateEl(TPZGeoEl1d *gel, TPZCompMesh &mesh, int64_t &index) {
  return new TPZIntelGen<pzshape::TPZShapeLinear>(mesh,gel,index);
}

TPZCompEl *(*TPZGeoEl1d::fp)(TPZGeoEl1d *,TPZCompMesh &, int64_t &) = CreateEl;

TPZGeoEl1d::TPZGeoEl1d(int id,TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh,int refind)
  : TPZGeoElRefLess<pzgeom::TPZGeoLinear>(id,nodeindices,matind,mesh) {

  fYAxisIndex = refind;
  int nnod = nodeindices.NElements();
  if(nnod <2 || nnod > 3) {
    PZError << "TPZGeoEl1d Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }
  else fNodeIndexes[2] = -1;
  fSubEl[0] = 0;
  fSubEl[1] = 0;
}

TPZGeoEl1d::TPZGeoEl1d( TPZVec<int64_t> &nodeindices, int matind, TPZGeoMesh &mesh, int refind)
  : TPZGeoElRefLess<pzgeom::TPZGeoLinear>(nodeindices,matind,mesh) {

  fYAxisIndex = refind;
  int nnod = nodeindices.NElements();
  if(nnod <2 || nnod > 3) {
    PZError << "TPZGeoEl1d Constuctor, number of nodes : " << nnod <<std::endl;
    return;
  }
  fSubEl[0] = 0;
  fSubEl[1] = 0;
}

TPZGeoEl1d *TPZGeoEl1d::CreateGeoEl(TPZVec<int64_t> &np, int matind, TPZGeoMesh &mesh, int refind) {
  return new TPZGeoEl1d(np,matind,mesh,refind);
}

TPZGeoEl1d::~TPZGeoEl1d() {
}

int TPZGeoEl1d::NSubElements() {
	return 2;
}

void TPZGeoEl1d::Jacobian(TPZVec<REAL> &fl,TPZFMatrix<STATE> &result,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<STATE> &jacinv){

  REAL spacephi[9], spacedphi[18];
  int numnodes = 3;
  if(fNodeIndexes[2] == -1) numnodes = 2;
  TPZFMatrix<STATE> phi(numnodes,1,spacephi,9),dphi(1,numnodes,spacedphi,18);
  Shape1d(fl[0],numnodes,phi,dphi);

  TPZGeoNode *np;
  int ic;
  TPZVec<REAL> v1(3,0.), vi(3,0.), v2(3,0.), v3(3,0.);
  REAL mod1 = 0;
  REAL modi = 0;
  REAL mod2 = 0;

  for(int i=0; i < numnodes; i++) {
    np = NodePtr(i);
    for(ic = 0; ic < 3; ic++) {
      v1[ic] += np->Coord(ic)*dphi(0,i);
    }
  }

  for(ic=0; ic<3; ic++) {
    mod1 += v1[ic]*v1[ic];
  }
  mod1 = sqrt(mod1);
  result(0,0) = mod1;
  detjac = mod1;
  jacinv(0,0) = 1./detjac;

  if(fYAxisIndex == -1) {
    axes.Redim(3,3);
	axes.Identity();
    return;
  }
  TPZGeoNode &yax = Mesh()->NodeVec()[fYAxisIndex];
  for(ic=0; ic<3; ic++) {
    vi[ic] = yax.Coord(ic);
    modi += vi[ic]*vi[ic];
  }
  modi = sqrt(modi);

  for(ic=0; ic < 3; ic++) {
    v1[ic] = v1[ic]/mod1;
    vi[ic] = vi[ic]/modi;
  }

  REAL kk = 0;
  for(ic=0; ic<3; ic++) kk += vi[ic]*v1[ic];
  for(ic=0; ic<3; ic++) {
    v2[ic] = vi[ic]-kk*v1[ic];
  }
  for(ic=0; ic<3; ic++) {
    mod2 += v2[ic]*v2[ic];
  }
  mod2 = sqrt(mod2);
  for(ic=0; ic<3; ic++) {
    v2[ic] = v2[ic]/mod2;
  }
  v3[0] = v1[1]*v2[2]-v2[1]*v1[2];
  v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
  v3[2] = v1[0]*v2[1]-v1[1]*v2[0];

  for(ic=0; ic<3; ic++) {
    axes(0, ic) = v1[ic];
    axes(1, ic) = v2[ic];
    axes(2, ic) = v3[ic];
  }
}


void TPZGeoEl1d::X(TPZVec<REAL> & par, TPZVec<REAL> &result){

  int numnodes = NNodes();
  REAL spacephi[9], spacedphi[18];
  TPZFMatrix<STATE> phi(numnodes,1,spacephi,9),dphi(1,numnodes,spacedphi,18);
  Shape1d(par[0],numnodes,phi,dphi);
  int in;
  for(in=0; in<3; in++) result[in] = 0.;
  for(in = 0; in < numnodes; in++) {
    int ic;
    for(ic=0; ic<3 ; ic++) {
      result[ic] += NodePtr(in)->Coord(ic)*phi(in,0);
    }
  }
}

int TPZGeoEl1d::NSideNodes(int side) {
  if(side == 0 || side == 1) return 1;
  if(side == 2) return 2;
  return 0;
}

int TPZGeoEl1d::SideNodeIndex(int side, int node){
  switch(side) {
  case 0:
  case 1:
    return fNodeIndexes[side];
  case 2:
    if(node>2) return -1;
    return fNodeIndexes[node];
  default:
    PZError << "TPZGeoEl1d::SideNodeIndex. Bad parameter side.\n";
  }
  return -1;
}

void TPZGeoEl1d::MidSideNodeIndex(int side,int64_t &index) {
  switch(side) {
  case 0:
    index = fNodeIndexes[0];
    return;
  case 1:
    index = fNodeIndexes[1];
    return;
  case 2:
    index = fNodeIndexes[2];
    if(index != -1) return;
    if(HasSubElement(2))
      index = SubElement(1)->NodeIndex(0);
    return;
  default:
    PZError << "TPZGeoEl1d::MidSideNode called for side " << side <<std::endl;
    PZError.flush();
    index = -1;
  }
}

int TPZGeoEl1d::SideDimension(int side) {
  switch(side) {
  case 0:
  case 1:
    return 0;
  case 2:
    return 1;
  default:
    PZError << "TPZGeoEl1d::SideDimension. Bad parameter side.\n";
    return -1;
  }
}

int64_t TPZGeoEl1d::NodeIndex(int node) {
  if(0 <= node && node< 3) return fNodeIndexes[node];
  return -1;
}

void TPZGeoEl1d::Divide(TPZVec<TPZGeoEl *> &pv) {
	TPZGeoElRefLess<pzgeom::TPZGeoLinear>::Divide(pv);
	return;
/*
  int locnodid[2][3];   // = {{0,2,1},{2,4,3}};             // Apenas modificado aqui (2-Jorge)
	locnodid[0][0] = 0;
	locnodid[0][1] = 2;
	locnodid[0][2] = 1;
	locnodid[1][0] = 2;
	locnodid[1][1] = 4;
	locnodid[1][2] = 3;
  int gnod[5] = {-1,-1,-1,-1,-1};
  int numnod = NNodes();
  pv.Resize(2);
  
  // in the case the element is already divided, return the sons
  
  if(fSubEl[0]) {
    pv[0] = fSubEl[0];
    pv[1] = fSubEl[1];
    return;
  }
  // identify the new nodes of the small elements
  CreateNewNodes(gnod);
  TPZVec<int64_t> np(numnod);
  int in,isub;
  TPZGeoEl1d *subel[2];
  for(isub=0; isub<2; isub++) {
    for(in=0; in<numnod; in++) np[in] = gnod[locnodid[isub][in]];
    subel[isub] = CreateGeoEl(np,MaterialId(),*Mesh(),fYAxisIndex);
    subel[isub]->SetFather(this);
  }

  // update the connectivity of the new elements
  TPZGeoElSide neighbour;

  // set the external connectivity of the elements
  TPZVec<int64_t> refnodes(2);
  TPZManVector<TPZGeoElSide> subelem(2);  // Is necessary extel and subelem  ???
  int side;
  for(side=0; side<3; side++) {
    subelem.Resize(0);
    TPZGeoElSide thisside(this,side);
    neighbour = Neighbour(side);
    if(!neighbour.Exists()) {
		subelem[0] = TPZGeoElSide(0,0);
		subelem[1] = TPZGeoElSide(0,0);
    } else {
      refnodes[0] = SideNodeIndex(side,0);
      refnodes[1] = SideNodeIndex(side,1);
      while(this != neighbour.Element()) {
	if(neighbour.HasSubElement()) break;
	neighbour = neighbour.Neighbour();
      }
      neighbour.GetSubElement(refnodes,subelem);   // if neighbour.Element() == this then the
    }                                              // subelements are neighbours of them self ???
    int cap = subelem.NElements();
    TPZGeoElSide gels0(subel[0],0), gels1(subel[1],1);
    if(side == 2) {
      if(cap) {
		  subelem[0] = subelem[0];
		  subelem[1] = subelem[1];
	TPZVec<int64_t> cornernode(1);
	cornernode[0] = gnod[2];
	int cornerside = subelem[0].Element()->WhichSide(cornernode);   //Para que ???
	if(cornerside < 0) {                                          // ^
	  PZError << "TPZGeoEl1d::Divide I dont understand\n";        // ^
	}
	TPZGeoElSide gels(subelem[0].Element(),cornerside);             // ^
	gels0.SetSide(1);                                             // para realizar a co-
	gels0.SetConnectivity(gels);          //nectividade em um ponto? Deveria ser pelo
	gels1.SetSide(0);                     //lado 2 do subelemento tambem ???
	gels1.SetConnectivity(gels0);
      } else  {
	gels0.SetSide(1);
	gels1.SetSide(0);
	gels0.SetConnectivity(gels1);
      }
    }
    switch (side) {
    case 0:
      if(cap) {
	gels0.SetSide(0);
	gels0.SetConnectivity(subelem[0]);
      } else {
	subel[0]->SetSideDefined(side);
      }
      break;
    case 1:
      if(cap) {
	//TPZGeoElSide(subel[1],1).SetConnectivity(subelem[1]);//tinha esta linha ativa
	gels1.SetSide(1);
	gels1.SetConnectivity(subelem[0]);//Cedric
	//Cedric:falta ainda atualizar a conectividade 1 do irmao para ele (=subel[1])
				//	TPZGeoElSide(subel[0],1).SetConnectivity(TPZGeoElSide(subel[1],0));//Cedric
      } else {
	subel[1]->SetSideDefined(side);
      }
      break;
    case 2:
      if(cap) {
	gels0.SetSide(2);
	gels1.SetSide(2);
	
	gels0.SetConnectivity(subelem[0]);
	
	gels1.SetConnectivity(subelem[1]);
      } else {
	subel[0]->SetSideDefined(side);
	subel[1]->SetSideDefined(side);
      }
      break;
    }
  }
  // fill in the pointers of the new elements
  fSubEl[0] = subel[0];
  fSubEl[1] = subel[1];
  
  pv[0] = subel[0];
  pv[1] = subel[1];
  */
}

void TPZGeoEl1d::CreateNewNodes(int64_t *gnodindex) {
  
  int numnod = NNodes();
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
  gnodindex[0] = NodeIndex(0);
  gnodindex[2] = midsidenodeindex;
  gnodindex[4] = NodeIndex(numnod-1);
  if(numnod == 3) {
    // look for the midsidenodes of the sons of the neighbouring elements
    TPZVec<int64_t> nodid(2),subelside(2);
    int side = 2;
    TPZStack<TPZGeoElSide> subelements;
    nodid[0] = NodeIndex(0);
    nodid[1] = NodeIndex(1);
    TPZGeoElSide neighbour;
    neighbour = Neighbour(side);//this e o elemento atual que esta sendo dividido
    if(neighbour.Exists()) {
      // try to find the midsidenodes of the neighbouring elements
      while(this != neighbour.Element()) {
	if(neighbour.HasSubElement()) {
	  neighbour.GetSubElements2(subelements);
	  TPZGeoEl *gel = subelements[0].Element();
	  if(gnodindex[1]==-1 && gel) gel->MidSideNodeIndex(subelements[0].Side(),gnodindex[1]);
	  gel = subelements[1].Element();
	  if(gnodindex[3]==-1 && gel) gel->MidSideNodeIndex(subelements[1].Side(),gnodindex[3]);
	}
	neighbour = neighbour.Neighbour();
      }
    }
    if(gnodindex[1] != -1 && gnodindex[3] != -1) return;
    TPZVec<REAL> gco(3);
    TPZVec<REAL> par(1);
    if(gnodindex[1] == -1) {
      par[0] = -0.5;
      X(par,gco);
      gnodindex[1] = Mesh()->NodeVec().AllocateNewElement();
      TPZGeoNode *gnodptr = &Mesh()->NodeVec()[gnodindex[1]];
      gnodptr->Initialize(gco,*Mesh());
    }
    if(gnodindex[3] == -1) {
      par[0] = 0.5;
      X(par,gco);
      gnodindex[3] = Mesh()->NodeVec().AllocateNewElement();
      TPZGeoNode *gnodptr = &Mesh()->NodeVec()[gnodindex[3]];
      gnodptr->Initialize(gco,*Mesh());
    }
  }
}

void TPZGeoEl1d::NormalVector(int side, TPZVec<REAL> &loc, TPZVec<REAL> &normal,
			      TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &/*jac*/){
	int ic;
	for(ic=1;ic<normal.NElements();ic++)
		normal[ic] = 0.;

	if(!side) {
		normal[0] = -1;
		axes.Identity();
		return;
	}
	else if(side==1) {
		normal[0] = 1.;
		axes.Identity();
		return;
	}

  REAL spacephi[9], spacedphi[18];
  int numnodes = NNodes();
  TPZFMatrix<STATE> phi(numnodes,1,spacephi,9),dphi(1,numnodes,spacedphi,18);
  Shape1d(loc[0],numnodes,phi,dphi);                                     //Cuidado com Shape1d (orientacao) (3-Jorge)
  TPZGeoNode *np;

  TPZVec<REAL> yaxis(3,0.),v1(3,0.), vi(3,0.), v2(3,0.), v3(3,0.);
  REAL mod1 = 0;
  REAL modi = 0;
  REAL mod2 = 0;
  if(fYAxisIndex != -1) {
    TPZGeoNode &node = Mesh()->NodeVec()[fYAxisIndex];
    yaxis[0] = node.Coord(0);
    yaxis[1] = node.Coord(0);
    yaxis[2] = node.Coord(0);
  }
  
  for(int i=0; i < numnodes; i++) {
    np = NodePtr(i);
    for(ic = 0; ic < 3; ic++) {
      v1[ic] += np->Coord(ic)*dphi(0,i);
      vi[ic] = yaxis[i];
    }
  }

  for(ic=0; ic<3; ic++) {
    mod1 += v1[ic]*v1[ic];
    modi += vi[ic]*vi[ic];
  }
  
  mod1 = sqrt(mod1);
  modi = sqrt(modi);
  
  for(ic=0; ic < 3; ic++) {
    v1[ic] = v1[ic]/mod1;
    if(modi) vi[ic] = vi[ic]/modi;
  }

  REAL kk = 0;
  for(ic=0; ic<3; ic++) kk += vi[ic]*v1[ic];

  for(ic=0; ic<3; ic++) {
    v2[ic] = vi[ic]-kk*v1[ic];
  }
  
  for(ic=0; ic<3; ic++) {
    mod2 += v2[ic]*v2[ic];
  }
  
  mod2 = sqrt(mod2);

  if(mod2) {
    for(ic=0; ic<3; ic++) {
      v2[ic] = v2[ic]/mod2;
    }
  }
  
  v3[0] = v1[1]*v2[2]-v2[1]*v1[2];
  v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
  v3[2] = v1[0]*v2[1]-v1[1]*v2[0];
  
  for(ic=0; ic<3; ic++) {
    axes(0, ic) = v1[ic];
    axes(1, ic) = v2[ic];
    axes(2, ic) = v3[ic];
    normal[ic] = v1[ic];
  }
  switch(side) {
  case 1:
    normal[0] *= -1.;
    normal[1] *= -1.;
    normal[2] *= -1.;
    break;
  case 0:
  case 2:
  default:
    break;
  }
}

/*
void TPZGeoEl1d::SideSubElements(int side, VoidPtrVec &sub) {
        if(!fSubEl) {
        cout << "Error : TPZGeoEl1d::SideSubElements called for an element without subelements\n";
      Print();
        return;
   }
   if(side < 0 || side > 2) {
        cout << "Error : TPZGeoEl1d::SideSubElements called for side " << side <<std::endl;
      Print();
        return;
   }
   if(side < 2) {
   	sub.resize(1);
	   sub[0] = (*fSubEl)[side];
   } else {
   	sub.resize(2);
	   sub[0] = (*fSubEl)[0];
	   sub[1] = (*fSubEl)[1];
	}
}
*/

/**Accumulates the transformation of the jacobian which maps the current
master element space into the space of the master element of the father*/
void TPZGeoEl1d::BuildTransform(int side, TPZGeoEl *father, TPZTransform<STATE> &t) {
  if(side != 2) return;
  int isub = 0;
  TPZGeoElSide locfather = Father(side);
  if(!locfather) {
	  std::cout << "TPZGeoEl1d::BuildTransform could not identify the father element\n";
    return;
  }
  for(; isub<2; isub++) {
    if(locfather.HasSubElement()) break;
  }
  TPZTransform<STATE> tloc(1);
  REAL store[2];
  TPZFMatrix<STATE> mult(1,1,store,1);
  TPZFMatrix<STATE> sum(1,1,store+1,1);
  switch(isub) {
  case 0:
    mult(0,0) = 0.5;
    sum(0,0) = -0.5;
    break;
  case 1:
    mult(0,0) = 0.5;
    sum(0,0) = 0.5;
    break;
  default:
    mult(0,0) = 1.;
    sum(0,0) = 0.;
	std::cout << "TPZGeoEl1d::BuildTransform subelement not detected within the"
      " father element\n";
  }
  tloc.SetMatrix(mult,sum);
  t = tloc.Multiply(t);
  if(locfather.Element() != father) locfather.Element()->BuildTransform2(side,father,t);
}

TPZTransform<STATE> TPZGeoEl1d::SideToSideTransform(int sidefrom, int sideto) {
  if(sideto != 2 && (sidefrom >=1 || sidefrom < 0)) {
    PZError << "TPZGeoEl1d:SideToSideTransform sidefrom = " << sidefrom <<
      " sideto = " << sideto <<std::endl;
  }
  TPZTransform<STATE> t(1,0);                                                          //   !!!   ???
  switch(sidefrom) {
  case 0:
    t.Sum()(0,0) = -1.;
    break;
  case 1:
    t.Sum()(0,0) = 1.;
    break;
  }
  return t;
}

int TPZGeoEl1d::NSideSubElements(int side) {
  if(side==2) return 2;
  return 1;
}

TPZGeoElSide TPZGeoEl1d::SideSubElement(int side,int subel) {
  if(side < 0 || side > 2) {
    PZError << "TGeoElQ2d::SideSubElement called for side " << side <<std::endl;
    return TPZGeoElSide(0,0);
  }
  if(!fSubEl[0]) return TPZGeoElSide(0,0);

  if(side==2)
    return TPZGeoElSide(fSubEl[subel],2);
  return TPZGeoElSide(fSubEl[side],side);
}

/**Neigh é o lado de um elemento vizinho ao elemento atual El que esta sendo dividido rfnd guarda
   os dois nós globais deste lado no sentido antihorario do elemento El, posições 0 e 1 de rfnd*/
void TPZGeoEl1d::GetSubElement(int side,TPZVec<int> &rfndindex,TPZVec<TPZGeoElSide> &sub) {
  if(!HasSubElement(side)) {
    sub.Resize(0);
    return;
  }
  switch(side) {
  case 0:
  case 1:
    sub.Resize(1);
    sub[0] = TPZGeoElSide(SubElement(side),side);
    break;
  case 2:
    sub.Resize(2);
    if(SubElement(0)->NodeIndex(0) == rfndindex[0]) {
      sub[0] = TPZGeoElSide(SubElement(0),2);
      sub[1] = TPZGeoElSide(SubElement(1),2);
    } else {
      sub[0] = TPZGeoElSide(SubElement(1),2);
      sub[1] = TPZGeoElSide(SubElement(0),2);
    }
  }
}

/**Inicializa os coeficientes do par de nós do lado I do elemento de referencia*/
void TPZGeoEl1d::SideMasterCo(int side,TPZVec<REAL> &IVec,TPZVec<REAL> &JVec) {
  // modified Philippe 28/7/97
  // if this method is called from a 3d element, IVec and JVec must be zeroed

  IVec.Fill(0.,0);
  JVec.Fill(0.,0);
  if (side==2) {
    IVec[0]=-1.;
    JVec[0]=1.;
  } else if(side==0) {
    PZError << "TPZGeoEl1d::SideMasterCo called for side 0\n";
    IVec[0]=-1.;
  } else if(side==1) {
    PZError << "TPZGeoEl1d::SideMasterCo called for side 1\n";
    IVec[0]=1.;
  }
}

void TPZGeoEl1d::SideMasterCo(int side,TPZFMatrix<REAL> &coord) {
  int row = coord.Rows();
  if(side == 0 || side == 1) {
    coord.Redim(row,1);
    if(side == 0) coord(0,0) = -1.;
    else coord(0,0) = 1.;
  }
  else {
    coord.Redim(row,2);
    coord(0,0) = -1.;
    coord(0,1) = 1.;
  }
}

TPZGeoElSide TPZGeoEl1d::Father(int side) {
	TPZGeoEl *fFather = TPZGeoEl::Father();
	if (!fFather) return TPZGeoElSide();
	int is;
//  for(is=0; is<2; is++) if(Father(side).SubElement(is) == this) break;
  if(is> 1) {
	  std::cout << "TPZGeoEl1d::Father is fishy\n";
    return TPZGeoElSide();
  }
  if(is == side || side == 2) {
    return Father(side);
  }
  return TPZGeoElSide();
}

TPZCompEl *TPZGeoEl1d::CreateBCCompEl(int side, int bc, TPZCompMesh &cmesh) {
  if(side==2) {
    TPZManVector<int64_t> nodes(3);
    TPZGeoEl1d *gel = CreateGeoEl(nodes,bc,*Mesh(),fYAxisIndex);
    TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,0));
    TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(this,1));
    TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(this,2));
	int64_t index;
    return gel->CreateCompEl(cmesh,index);
  }
  else if(side==0 || side==1) {
      TPZMaterial *bcptr = cmesh.FindMaterial(bc);
      if(!bcptr) {
      PZError << "TPZGeoEl1d::CreateBCCompEl has no bc.\n";
      return 0;
      }
      // Philippe 1/12/99
      //      TPZCompEl *cel = Reference();
      //      if(!cel) {
      //      PZError << "TPZGeoEl1d::CreateBCCompEl has no computational element\n";
      //      return 0;
      //      }
      //Cedric 15/04/99
      TPZGeoElPoint *gel;
      //Jorge 8/5/99

      TPZManVector<int64_t> node(1);
	  node[0] = fNodeIndexes[side];
	  gel = new TPZGeoElPoint(node,bc,*Mesh());
      TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(this,side));
	  int64_t index;
      return gel->CreateCompEl(cmesh,index);

       /*int bcindex = cmesh.BCConnectVec().AllocateNewElement();
       TPZConnectBC &bcconnect = cmesh.BCConnectVec()[bcindex];

       TPZConnect *connect = &cel->Connect(side);
       bcconnect.fConnect = connect;
       bcconnect.fBC = (TPZBndCond *) bcptr;
       return 0;    */ //Cedric
  }
  else PZError << "TPZGeoEl1d::CreateBCCompEl. Side = " << side <<std::endl;
  return 0;
}

void TPZGeoEl1d::LowerDimensionSides(int side,TPZStack<TPZGeoElSide> &smallsides) {

  switch(side) {
  case 0:
  case 1:
    if(side==0 || side==1) return;
    break;
  case 2:
    smallsides.Push(TPZGeoElSide(this,0));
    smallsides.Push(TPZGeoElSide(this,1));
    break;
  }

}

TPZGeoElSide TPZGeoEl1d::HigherDimensionSides(int side,int targetdimension) {
  if(side > 1 && targetdimension < 1) {
    PZError << "TPZGeoEl1d::HigherDimensionSides called with side = " << side
	    << " targetdimension = " << targetdimension <<std::endl;
    return TPZGeoElSide();
  }
  switch(targetdimension) {
	  case 0 :
	    return TPZGeoElSide(this,side);
	  case 1 :
	    return TPZGeoElSide(this,2);
	  default :
	    return TPZGeoElSide();
  }
}

TPZGeoElSide TPZGeoEl1d::Father2(int side) {
	return TPZGeoElRefLess<pzgeom::TPZGeoLinear>::Father2(side);
/*	if(side>2 || side < 0 || !Father()){
		PZError << "TPZGeoEl1d::Father2 called error" <<std::endl;
	}
	if(Father()->SubElement(0)==this){
		if(side==0) return TPZGeoElSide(Father(),0);
		else return TPZGeoElSide(Father(),2);
	}
	if(Father()->SubElement(1)==this){
		if(side==1) return TPZGeoElSide(Father(),1);
		else return TPZGeoElSide(Father(),2);
	}
    return TPZGeoElSide();*/
}

static int subeldata[3][3][2] =
{
	{{0,0}},/*side=0 { isub0{0,1},isub1{0,1} }*/
	{{1,1}},/*side=1*/
	{{0,2},{0,1},{1,2}}/*side=2*/
};

static int nsubeldata[3] = {1,1,3};


void TPZGeoEl1d::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel){//Augusto:09/01/01

  subel.Resize(0);
	if(side<0 || side>2 || !HasSubElement(side)){
		PZError << "TPZGeoEl1d::GetSubElements2 called error" <<std::endl;
		return;
	}
	int nsub = nsubeldata[side];
	int isub;
	for(isub=0; isub<nsub; isub++) subel.Push(TPZGeoElSide(fSubEl[subeldata[side][isub][0]],
		                                                            subeldata[side][isub][1]));
}

static REAL buildt[2][3][2] = {//por colunas
/*S0*/
     {/*0*/{0.,-1.},
      /*1*/{0.,0.},
      /*2*/{.5,-.5}},
/*S1*/
     {/*0*/{0.,0.},
      /*1*/{0.,1.},
      /*2*/{.5,.5}}
};

TPZTransform<STATE> TPZGeoEl1d::BuildTransform2(int side, TPZGeoEl * /*father*/){//Augusto:09/01/01
	if(side<0 || side>2 || !Father(side).Element()){
  	PZError << "TPZGeoElT2d::BuildTransform2 side out of range or father null\n";
    return TPZTransform<STATE>(0,0);
  }
  TPZTransform<STATE> trans(1,1);
  int son = WhichSubel();
  trans.Mult()(0,0) = buildt[son][side][0];
  trans.Sum() (0,0) = buildt[son][side][1];

  return trans;
}

static REAL MidSideNode[3][3] = {
/*00*/{-1.},/*01*/{1.},/*02*/{0.} };

int TPZGeoEl1d::main(TPZGeoEl *gel){

  TPZVec<TPZGeoEl *> subs;
  gel->Divide(subs);
  int sn,sd;
  TPZVec<REAL> x1(3),x2(3);//x1 no filho deformado, x2 no pai deformado
  TPZManVector<REAL> ps(3),pss(3),pf(3),pfs(3);
                    //point son, point side son, point father, point side father : elemento mestre
  pss[1] = 0.; pss[2] = 0.;//1d
  pfs[1] = 0.; pfs[2] = 0.;
  pf[1] = 0.; pf[2] = 0.;
  for(sn=0;sn<2;sn++){
    TPZGeoEl *son = subs[sn];
    for(sd=0;sd<3;sd++){
      ps[0] = MidSideNode[sd][0];//element master point
      ps[1] = MidSideNode[sd][1];// = 0
      ps[2] = MidSideNode[sd][2];// = 0
      if(son->WhichSide(ps) != sd) std::cout << "Lado nao bate\n";
      TPZTransform<STATE> telsd = pzshape::TPZShapeLinear::TransformElementToSide(sd);//2x2
      telsd.Apply(ps,pss);//son element -> side
      son->X(ps,x1);//ponto deformado filho
	  TPZTransform<STATE> t;
		  son->BuildTransform2(sd,gel,t);
      t.Apply(pss,pfs);//son side -> fat side
      int sdfat = son->Father2(sd).Side();
      telsd = pzshape::TPZShapeLinear::TransformSideToElement(sdfat);//2x2
      telsd.Apply(pfs,pf);//lado do pai -> pai
      son->Father2(2).Element()->X(pf,x2);
      if( sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) ) > 1.e-10 ){
      	PZError << "\nTransformacao furada\n";
        PZError << "son    = " << (son->Id()) <<std::endl;
        PZError << "father = " << ((son->Father2(2).Element())->Id()) <<std::endl;
        PZError << "side   = " << sd <<std::endl <<std::endl;
        int ok;
		std::cin >> ok;
      } else {
		  std::cout << "Transformacao OK!\n";
		  std::cout << "Filho/lado : " << son->Id() << "/" << sd <<std::endl;
		  std::cout << "Pai : " << son->Father2(2).Element()->Id() <<std::endl <<std::endl;
      }
    }
  }
  return 1;
}

/** Jorge 17/7/99 */
/** Compute the measure of the geometrical element */
REAL TPZGeoEl1d::Mesure(int dim) {
  if(dim!=1) return 0.;
  REAL measure = 0;
  if(IsZero(Volume())) {
    TPZGeoNode &nod1 = Mesh()->NodeVec()[fNodeIndexes[0]];
    REAL x0 = nod1.Coord(0), y0 = nod1.Coord(1), z0 = nod1.Coord(2);
    TPZGeoNode &nod2 = Mesh()->NodeVec()[fNodeIndexes[1]];
    x0 -= nod2.Coord(0);
    y0 -= nod2.Coord(1);
    z0 -= nod2.Coord(2);
    for(int i=0;i<3;i++) measure = sqrt(x0 * x0 + y0 * y0 + z0 * z0);
  }
  if(IsZero(measure-Volume()))
	  return measure;
  return 0.;
}

void TPZGeoEl1d::Center(TPZVec<REAL> &center) {
  center[0] = 0.;
}



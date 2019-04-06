
#include "myheader.h"
#include "elcpointgd.h"
#include "interface.h"
#include "pzblock.h"
#include "pzcmesh.h"
#include "pzconnect.h"
#include "pzelgpoint.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZMaterial.h"

TPZCompEl *TCompElPointGD::CreateElDisc(TPZGeoElPoint *gel,TPZCompMesh &mesh, int64_t &index) {
  return new TCompElPointGD(mesh,gel,index);
}

TCompElPointGD::TCompElPointGD(TPZCompMesh &mesh, TPZGeoElPoint *ref, int64_t &index)
  : TPZCompElPoint(mesh,ref,index,0) {

  /**Determining the discontinuity upon the neighbours through side connect*/
    fConnectDisc = NeighbourDiscontinuous(0);
  /**Creating the new internal connect, it belongs only to the current discontinuous element*/
  int nvar = 1;
  if(Material()) nvar = Material()->NStateVariables();
  fConnect = mesh.AllocateNewConnect(1,nvar,1);
  TPZConnect &connect = mesh.ConnectVec()[fConnect];
  mesh.Block().Set(connect.SequenceNumber(),nvar*NConnectShapeF(1));
  connect.IncrementElConnected();

  /**Creating the connect over the sides of the current element*/
  SetConnectIndex(0,CreateMidSideConnect(0));
  mesh.ConnectVec()[ConnectIndex(0)].IncrementElConnected();

  ref->SetReference(this);
  /**Making discontinuous side connect of the element except when exist a continuous neighbour
  if(fConnectDisc)
    	//Can to be continuous if exist neighbour which is continuous element
    MakeConnectDiscontinuous(0);
  else
    MakeConnectContinuous(0); */
}

TCompElPointGD::~TCompElPointGD() {
  TPZConnect *con = &Mesh()->ConnectVec()[fConnect];
  Mesh()->Block().Set(con->SequenceNumber(),0);
  con->ResetElConnected();
  Mesh()->ExpandSolution();
}

int TCompElPointGD::NeighbourDiscontinuous(int /*icon*/) {
  TPZStack<TPZCompElSide> elvec;
  TPZCompElSide thisside(this,0);
  thisside.EqualLevelElementList(elvec,0,0);
  /**Asking if exist neighbour that it can not to be discontinuous*/
  if(elvec.NElements() && ((TPZInterpolatedElement *)(elvec[0].Element()))->ConnectIndex(elvec[0].Side()) < 0)
    return 0;
  return 1;
}

/**Make discontinuous the side connect of the discontinuous element.
   It is tranfered the values associated from side connect to fConnect*/
void TCompElPointGD::MakeConnectDiscontinuous(int /*icon*/) {
  if(fConnectDisc) return;

  TPZCompElSide thisside(this,0);

  /** Check whether there is a neighbouring element which cannot be discontinuous*/
  TPZStack<TPZCompElSide> elvec;
  thisside.EqualLevelElementList(elvec,0,0);
  int64_t j,n=elvec.NElements();
  /**Asking if exist neighbour that it can not to be discontinuous*/
//  for(j=0;j<n;j++) {
//    if(!((TPZInterpolatedElement *)(elvec[j].Element()))->CanBeDiscontinuous()) {
//      return;
//    }
//  }

  /**Removing restrains respect to large element and in small element connected*/
  TPZStack<TPZCompElSide> smallelvec;
  TPZGeoElSide gelside = thisside.Reference();
  gelside = gelside.Neighbour();
  if(!gelside.Exists()) PZError << "TCompElPointGD::MakeConnectDiscontinuous.. Not exist neighboard.";
  else gelside.HigherLevelCompElementList2(smallelvec,1,1);
//  thisside.HigherLevelElementList(smallelvec,1,1);
  thisside.ExpandConnected(smallelvec,1);
  n = smallelvec.NElements();
  /**Asking if exist neighbour that it can not to be discontinuous
  for(j=0;j<n;j++) {
    if(!((TPZInterpolatedElement *)(smallelvec[j].Element()))->CanBeDiscontinuous()) {
      return;
    }
  }*/
  TPZCompElSide large = thisside.LowerLevelElementList(1);
  /**Asking if neighbour (large) can not to be discontinuous
  if(!((TPZInterpolatedElement *)(large.Element()))->CanBeDiscontinuous()) {
    return;
  }*/
  /** Remove restricoes nos elementos vizinhos de maior nivel */
  for(j=0;j<n;j++) {
    TPZInterpolatedElement *el = (TPZInterpolatedElement *)smallelvec[j].Element();
    el->RemoveSideRestraintWithRespectTo(smallelvec[j].Side(),thisside);
//	elvec.Push(smallelvec[j]);
  }
  if(large.Exists()) RemoveSideRestraintWithRespectTo(0,large);

  /**Making discontinuous all the neighbour elements*/
  elvec.Push(thisside);
  n = elvec.NElements();
  for(j=0;j<n;j++) {
//    TPZInterpolatedElement *cel = (TPZInterpolatedElement *)(elvec[j].Element());
 //   cel->SetConnectDiscontinuous(elvec[j].Side());
  }

  /**Redimensioning the block of the icon connect*/
  TPZConnect &c = Connect(0);
  int seqnum = c.SequenceNumber();
  Mesh()->Block().Set(seqnum,0);
  Mesh()->ExpandSolution();

  /**Restablecing restraint into the large*/
  n = smallelvec.NElements();
  for(j=0;j<n;j++) {
    TPZInterpolatedElement *el = (TPZInterpolatedElement *)smallelvec[j].Element();
    el->RestrainSide(smallelvec[j].Side(),this,0);
  }
  if(large.Exists()) RestrainSide(0,(TPZInterpolatedElement *)large.Element(),large.Side());
}

/**Make continuous the icon connect of the discontinuous element
   It is tranfered the values associated from fConnect to icon connect*/
void TCompElPointGD::MakeConnectContinuous(int /*icon*/) {
  if(!fConnectDisc) return;
  TPZCompElSide thisside(this,0);

  /** Check whether there is a neighbouring element which cannot be discontinuous*/
  TPZStack<TPZCompElSide> elvec;
  thisside.EqualLevelElementList(elvec,0,0);

  TPZCompElSide large = thisside.LowerLevelElementList(1);
  if(large.Exists()) {
    RemoveSideRestraintWithRespectTo(0,large);
  }

  /**Redimensioning the block of the connect icon and copying the adequated
     values from discontinuous connect*/
  elvec.Push(thisside);
  int condition = StablizingConnectContinuous(0,elvec);
  if(condition) PZError << "TCompElPointGD::MakeConnectContinuous is not done.\n";

  /** Restablecing restraint into the large */
  if(large.Exists()) RestrainSide(0,(TPZInterpolatedElement *)large.Element(),large.Side());
}

/**To redimensioning the block of the icon connect that is being actived
   in continuous and copying the adequated values from discontinuous connect*/
int TCompElPointGD::StablizingConnectContinuous(int icon,TPZStack<TPZCompElSide> &elvec) {
#ifndef NOTDEBUG
  if(icon) PZError << "TCompElPointGD::StablizingConnectContinuous. Bad parameter icon (!=0)." << std::endl;
#endif
  int i, j;
  int n = elvec.NElements();
  /**Verificando se todos os elementos em elvec tem o connect continuo ou descontinuo*/
  for(i=0;i<n;i++)
	if(fConnectDisc == ((TPZInterpolatedElement *)elvec[i].Element())->ConnectIndex(elvec[i].Side()))
	  return 1;   //existe um vizinho com continuidade diferente
  if(!fConnectDisc) return 0;   //Ja eh continuo.

  TPZCompMesh *mesh = Mesh();
  TPZBlock<STATE> &block = mesh->Block();
  int nvar = Material()->NStateVariables();
  TPZConnect &c = Connect(0);
  int seqnum = c.SequenceNumber();
  int rows = TPZCompElPoint::NConnectShapeF(0)*nvar;
  int cols = mesh->Solution().Cols();

  TPZFMatrix<STATE> values(rows,cols,0.);

  /**Redimensioning the block of the first connect of the element*/
  block.Set(seqnum,rows);
  Mesh()->ExpandSolution();

  /**Making connect continuous over all the neighbour elements*/
  for(j=0;j<n;j++) {
    TPZInterpolatedElement *cel = (TPZInterpolatedElement *)(elvec[j].Element());
 //   cel->SetConnectContinuous(elvec[j].Side(),values);
  }

  /**Computing the mean values*/
  for(i=0;i<rows;i++)
	for(j=0;j<cols;j++)
	  values(i,j) /= ((double)n);

  /**Copying values into the block of the actual continuous connect*/
  for(i=0;i<rows;i++)
    for(j=0;j<cols;j++)
      block(seqnum,0,i,j) = values(i,j);

  return 0;
}
void TCompElPointGD::SetConnectContinuous(int icon,TPZFMatrix<STATE> &values) {
  #ifndef NOTDEBUG
  if(icon) {
    PZError << "TCompElPointGD::SetConnectContinuous. Bad parameter icon.\n";
    return;
  }
  #endif

  if(!fConnectDisc) return;

  /**Modify the data structure of the connect to be discontinuous.
     Interchange the values in icon connect block to fConnect block*/
  TPZBlock<STATE> &block = Mesh()->Block();
  int nvar = Material()->NStateVariables();
  int cols = Mesh()->Solution().Cols();
  int i, j, conseq = Connect(1).SequenceNumber();
  for(i=0;i<nvar;i++)
	for(j=0;j<cols;j++)
	  values(i,j) += block(conseq,0,i,j);

  /** Adapting the datastructure of the connect and redimensioning the block of fConnect*/
  fConnectDisc = 0;
  block.Set(conseq,0);
  /**Adequating solution vector*/
  Mesh()->ExpandSolution();
}

void TCompElPointGD::SetConnectIndex(int i, int64_t connectindex) {
 #ifndef NOTDEBUG
  if(i!=0 && i!=1) {
    PZError << "TCompEl1d::SetConnectIndex. Bad parameter i = " << i << " .\n";
    PZError.flush();
  }
 #endif
  if(!i) TPZCompElPoint::SetConnectIndex(0,connectindex);
  else fConnect = connectindex;
}

int64_t TCompElPointGD::ConnectIndex(int i) {
 #ifndef NOTDEBUG
  if(i!=0 && i!=1) {
    PZError << "TCompElPointGD::ConnectIndex. Bad parameter i = " << i << " .\n";
    PZError.flush();
    return -1;
  }
 #endif
  if(!i) return TPZCompElPoint::ConnectIndex(0);
  return fConnect;
}

int TCompElPointGD::NConnectShapeF(int connect) {
 #ifndef NOTDEBUG
  if(connect!=0 && connect!=1) {
    PZError << "TCompElPointGD::NConnectShapeF. Bad parameter connect.\n";
    return 0;
  }
 #endif
  if(connect==1) {
    if(fConnectDisc) return 1;
    return 0;
  }
  if(fConnectDisc) return 0;
  return 1;
}

int TCompElPointGD::NSideConnects(int side) {
  if(!side) return 2;
  PZError << "TCompElPointGD::NSideConnects. Bad parameter side = " << side << " .\n";
  return 0;
}

int TCompElPointGD::SideConnectLocId(int c,int side) {
 #ifndef NOTDEBUG
  if(c!=0 && c!=1) {
    PZError << "TPZCompEl1d::SideConnectLocId called with connect = " << c << std::endl;
    return 0;
  }
 #endif
  if(!side)
    return c;
  PZError << "TPZCompEl1d::SideConnectLocId called with side = " << side << std::endl;
  return 0;
}

void TCompElPointGD::SetConnectDiscontinuous(int icon) {
  #ifndef NOTDEBUG
  if(icon) {
    PZError << "TCompElPointGD::SetConnectDiscontinuous. Bad parameter icon.\n";
    return;
  }
  #endif

  if(fConnectDisc) return;

  /**Modify the data structure of the connect to be discontinuous.
     Interchange the values in icon connect block to fConnect block*/
  TPZBlock<STATE> &block = Mesh()->Block();
  /**Searching the position and size of the values to connect(side) that must to
     be into the block corresponding to internal connect*/
  int nvar = 1;
  if(Material()) nvar = Material()->NStateVariables();
  int rows = nvar;
  int cols = Mesh()->Solution().Cols();

  /** Adapt the datastructure of the connect */
  fConnectDisc = 1;
  /**Redimensioning the block of the fConnect connect*/
  int conseq = Connect(0).SequenceNumber();
  int conseqdis = Connect(1).SequenceNumber();
  block.Set(conseqdis,rows);
  /** Adequating solution vector*/
  Mesh()->ExpandSolution();

  /** Put the solution back into the fConnect block from the block of the first connect*/
  for(int i=0;i<rows;i++) {
	 for(int j=0;j<cols;j++) {
		block(conseqdis,0,i,j) = block(conseq,0,i,j);
	 }
  }
}

int TCompElPointGD::IsConnectContinuous(int icon) {
  #ifndef NOTDEBUG
  if(icon) {
    PZError << "TCompElPointGD::IsSideContinuous. Bad parameter side.\n";
    return 1;
  }
  #endif
  return !fConnectDisc;
}

void TCompElPointGD::Print(std::ostream &out) {
  out << std::endl << "Discontinuous element\n";
  TPZInterpolatedElement::Print(out);
}

/*******                                                             *******/
/***    Methods to discontinuous computational one-dimensional element   ***/

TPZCompEl *TCompElPointWI::CreateElDiscWI(TPZGeoElPoint *gel,TPZCompMesh &mesh, int64_t &index) {
  return new TCompElPointWI(mesh,gel,index);
}

TCompElPointWI::TCompElPointWI(TPZCompMesh &mesh, TPZGeoElPoint *ref, int64_t &index)
  : TCompElPointGD(mesh,ref,index) {

  fInterface = CreateInterface(0,mesh);
}

TCompElPointWI::~TCompElPointWI() {
  DeleteInterfaces();
}

void TCompElPointWI::DeleteInterfaces() {
  if(fInterface == -1) return;
  if(fInterface==-2) {
    TPZStack<TPZCompElSide> elvec;
    TPZCompElSide thisside(this,0);
    TPZGeoElSide gelside = thisside.Reference();
    gelside = gelside.Neighbour();
    if(!gelside.Exists()) PZError << "TCompElPointGD::DeleteInterfaces. Not exist neighboard.";
    gelside.HigherLevelCompElementList2(elvec,1,0);
//    thisside.HigherLevelElementList(elvec,1,0);
    TPZInterpolatedElement *intel = (TPZInterpolatedElement *)elvec[0].Element();
    int index = intel->Interface(elvec[0].Side());
    TInterfaceElement *inter_face = (TInterfaceElement *)fMesh->ElementVec()[index];
    if(!inter_face) return;
    delete inter_face;
    fMesh->ElementVec()[index] = 0;
    intel->SetInterface(elvec[0].Side(),-1);
    fInterface = -1;       // Jorge 17/5/2000
    return;
  }	
  /** When one interface element is deleted it delete respective fInterface and
      change into Right computational element by -1 the respective fInterface*/
  TInterfaceElement *inter_face = (TInterfaceElement *)Mesh()->ElementVec()[fInterface];
  TPZInterpolatedElement *intel = inter_face->LeftEl();
  int side;
  if(intel==this) {
    intel = inter_face->RightEl();
    side = inter_face->RightSide();
  }
  else 
    side = inter_face->LeftSide();
  delete inter_face;
  fMesh->ElementVec()[fInterface] = 0;
  fInterface = -1;       // Jorge 17/5/2000
  if(intel->Interface(side)!=-2)
    intel->SetInterface(side,-1);
}

/** If fInterface[side] == -2 means exist several neighboards of the lower level*/
int TCompElPointWI::CreateInterface(int side,TPZCompMesh &mesh) {
  int i, index = -1, indexneigh = -1, nneighs;

  TPZStack<TPZCompElSide> elvec;
/**Creating the interfaces over side with one neighboards discontinuous element*/
  /** Searching over the side neighboards with same level. Apenas eh permitida
	    a construcao de interface entre elementos do mesmo nivel */
  TPZCompElSide thisside(this,0);
  thisside.EqualLevelElementList(elvec,1,1);
  nneighs = elvec.NElements();
	if(!nneighs) {
	  PZError << "TCompElPoint::CreateInterface, have not neighboards\n";
		return index;
	}
  if(nneighs!=1)
	  PZError << "TCompElPoint::CreateInterface, has many neighboards\n";
  /**Verificando se todos os elementos computacionais podem ter interface*/
  for(i=0;i<nneighs;i++)
    if(!((TPZInterpolatedElement *)elvec[i].Element())->CanHaveInterface())
      return -1;

  /**The dimension of the thisside must to be equal to the dimension of the common
     side on the neighboard element in elvec*/
  int dim = thisside.Reference().Dimension();
  int dimmat = Material()->Dimension();
  if(dim != elvec[0].Reference().Dimension()) return index;
  int fluxtype = Material()->FluxType();
  /** Creating interface over this side boundary. The neighboard element cann't
      to has interface element over this common side boundary*/
  for(i=0;i<nneighs;i++) {
//    new TInterfaceElement(mesh,thisside,elvec[i],fluxtype,index);
  //  ((TPZInterpolatedElement *)elvec[i].Element())->SetInterface(elvec[i].Side(),index);
  }
  return index;
}

void TCompElPointWI::Print(std::ostream &out) {
  TCompElPointGD::Print(out);
  /** Information over interfaces*/
  out << "Interfaces information :\t";
  if(fInterface==-1) out << "no interface\n";
  else out << fInterface << std::endl;
}

int TCompElPointWI::Interface(int side) {
#ifndef NOTDEBUG
  if(side) {
    PZError << "TCompEl1dWI::Interface. Bad side.\n";
    return -1;
  }
#endif NOTDEBUG
  return fInterface;
}

void TCompElPointWI::SetInterface(int side, int64_t indexinterface) {
#ifndef NOTDEBUG
  if(side) {
    PZError << "TCompEl1dWI::SetInterface. Bad side.\n";
    return;
  }
#endif NOTDEBUG
  fInterface = indexinterface;
}



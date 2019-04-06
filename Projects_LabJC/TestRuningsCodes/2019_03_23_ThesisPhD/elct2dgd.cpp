#include "elct2dgd.h"
#include "elc1dgd.h"
#include "interface.h"
#include "myheader.h"

#include "pzblock.h"
#include "pzcmesh.h"
#include "pzconnect.h"
#include "pzquad.h"
#include "pzelgt2d.h"
#include "pztrigraphd.h"
#include "pzelmat.h"
#include "TPZMaterial.h"
#include "TPZInterfaceEl.h"

#include "pzshapelinear.h"

TPZCompEl *TCompElT2dGD::CreateElDisc(TPZGeoElT2d *gel,TPZCompMesh &mesh,int64_t &index) {
  return new TCompElT2dGD(mesh,gel,index);
}

TCompElT2dGD::TCompElT2dGD(TPZCompMesh &mesh, TPZGeoElT2d *ref, int64_t &index)
  : TPZCompElT2d(mesh,ref,index,0) {

  int icon;

  /**Determining the discontinuity upon the neighbours through icon connect*/
  for(icon=0;icon<7;icon++)
    fConnectDisc[icon] = NeighbourDiscontinuous(icon);
  /**Creating the new internal connect, it belongs only to the current discontinuous element*/
  int nvar = 1;
  if(Material())
    nvar = Material()->NStateVariables();
  fConnect = mesh.AllocateNewConnect(NShapeF(), nvar, GetPreferredOrder());
//  fConnect = mesh.AllocateNewConnect();
  TPZConnect &connect = mesh.ConnectVec()[fConnect];
  mesh.Block().Set(connect.SequenceNumber(),nvar*NConnectShapeF(7));
  connect.IncrementElConnected();
  SetPositions();

  /**Creating the connect over the sides of the current element*/
  for(icon=0;icon<7;icon++) {
    SetConnectIndex(icon,CreateMidSideConnect(icon));
    mesh.ConnectVec()[ConnectIndex(icon)].IncrementElConnected();
  }

  TPZVec<int> order(2,2*SideOrder(6));
  GetIntegrationRule().SetOrder(order);

  ref->SetReference(this);

  SetPositions();
}

TCompElT2dGD::TCompElT2dGD(TPZCompMesh &mesh, TPZGeoElT2d *ref, int64_t &index,TPZVec<int> &maskdisc)
  : TPZCompElT2d(mesh,ref,index,0) {
  int icon;

  /**Determining the discontinuity upon the neighbours through icon connect*/
  for(icon=0;icon<7;icon++)
    fConnectDisc[icon] = NeighbourDiscontinuous(icon);
  /**Creating the new internal connect, it belongs only to the current discontinuous element*/
  int nvar = 1;
  if(Material())
    nvar = Material()->NStateVariables();
  fConnect = mesh.AllocateNewConnect(NShapeF(), nvar, GetPreferredOrder());
  TPZConnect &connect = mesh.ConnectVec()[fConnect];
  mesh.Block().Set(connect.SequenceNumber(),nvar*NConnectShapeF(7));
  connect.IncrementElConnected();
  SetPositions();

  /**Creating the connect over the sides of the current element*/
  for(icon=0;icon<7;icon++) {
    SetConnectIndex(icon,CreateMidSideConnect(icon));
    mesh.ConnectVec()[ConnectIndex(icon)].IncrementElConnected();
  }

  TPZVec<int> order(2,2*SideOrder(6));
  GetIntegrationRule().SetOrder(order);

  ref->SetReference(this);
  /**Connects mask with maskdisc are doing discontinuous.
  Can to be continuous if exist a neighbour which is continuous element*/
  for(icon=0;icon<7;icon++) {
    if(maskdisc[icon]) {
      MakeConnectDiscontinuous(icon);
    }
    else {
      MakeConnectContinuous(icon);
    }
  }

  SetPositions();
}

TCompElT2dGD::~TCompElT2dGD() {
  TPZConnect *con = &Mesh()->ConnectVec()[fConnect];
  Mesh()->Block().Set(con->SequenceNumber(),0);
  con->ResetElConnected();
  Mesh()->ExpandSolution();
}

/**Return 1 if not exist neighbour or if neighbour is discontinuous*/
int TCompElT2dGD::NeighbourDiscontinuous(int icon) {
  TPZStack<TPZCompElSide> elvec;
  TPZCompElSide thisside(this,icon);
  thisside.EqualLevelElementList(elvec,0,0);
//  if(!SideOrder(icon)) return 1;
  /**Asking if exist neighbour that it can not to be discontinuous*/
//  if(elvec.NElements() && ((TPZInterpolatedElement *)(elvec[0].Element()))->IsConnectContinuous(elvec[0].Side()))
 //   return 0;
  return 1;
}

/**Make discontinuous the icon connect of the discontinuous element.
   It is tranfered the values associated from icon connect to fConnect*/
void TCompElT2dGD::MakeConnectDiscontinuous(int icon) {
  #ifndef NOTDEBUG
  if(icon<0 || icon>6) {
    PZError << "TCompElT2dGD::MakeConnectDiscontinuous. Bad parameter side.\n";
    return;
  }
  #endif
  if(fConnectDisc[icon]) return;

  TPZCompElSide thisside(this,icon);

  /** Check whether there is a neighbouring element which cannot be discontinuous*/
  TPZStack<TPZCompElSide> elvec;       // Porque mexe na fSolution??? !!!
  thisside.EqualLevelElementList(elvec,0,0);
  int j, n=elvec.NElements();
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
  if(!gelside.Exists()) PZError << "TCompElT2dGD::MakeConnectDiscontinuous.. Not exist neighboard.";
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
  if(large.Exists()) RemoveSideRestraintWithRespectTo(icon,large);

  /**Making discontinuous all the neighbour elements*/
  elvec.Push(thisside);
  n = elvec.NElements();
  for(j=0;j<n;j++) {
//    TPZInterpolatedElement *cel = (TPZInterpolatedElement *)(elvec[j].Element());
//    cel->SetConnectDiscontinuous(elvec[j].Side());
  }

  /**Redimensioning the block of the icon connect*/
  TPZConnect &c = Connect(icon);
  int seqnum = c.SequenceNumber();
  Mesh()->Block().Set(seqnum,0);
  Mesh()->ExpandSolution();

  /**Restablecing restraint into the large and short neighbour elements */
  n = smallelvec.NElements();
  for(j=0;j<n;j++) {
    TPZInterpolatedElement *el = (TPZInterpolatedElement *)smallelvec[j].Element();
    el->RestrainSide(smallelvec[j].Side(),this,icon);
  }
  if(large.Exists()) RestrainSide(icon,(TPZInterpolatedElement *)large.Element(),large.Side());
}

/**Make continuous the icon connect of the discontinuous element
   It is tranfered the values associated from fConnect to icon connect*/
void TCompElT2dGD::MakeConnectContinuous(int icon) {
  #ifndef NOTDEBUG
  if(icon<0 || icon>6) {
    PZError << "TCompElT2dGD::MakeConnectContinuous. Bad parameter side.\n";
    return;
  }
  #endif

  if(!fConnectDisc[icon]) return;

  TPZCompElSide thisside(this,icon);

  /** Check whether there is a neighbouring element which cannot be discontinuous*/
  TPZStack<TPZCompElSide> elvec;
  thisside.EqualLevelElementList(elvec,0,0);
  int i,n;

  /**Removing restrains respect to large element and in small element connected*/
  TPZStack<TPZCompElSide> smallelvec;
  TPZGeoElSide gelside = thisside.Reference();
  gelside = gelside.Neighbour();
  TPZStack<TPZGeoElSide> gsmallvec;
  if (!gelside.Exists()) PZError << "TCompEl1dGD::MakeConnectContinuous.. Not exist neighboard.";
  else gelside.GetSubElements2(gsmallvec);
  for (i = 0; i < gsmallvec.NElements(); i++)
	  smallelvec.Push(gsmallvec.Pop().Reference());

//  thisside.HigherLevelElementList(smallelvec,1,1);
  thisside.ExpandConnected(smallelvec,1);
  n = smallelvec.NElements();
  for(i=0;i<n;i++) {
    TPZInterpolatedElement *el = (TPZInterpolatedElement *)smallelvec[i].Element();
    el->RemoveSideRestraintWithRespectTo(smallelvec[i].Side(),thisside);
  }
//  for(i=0;i<n;i++) elvec.Push(smallelvec[i]);
  TPZCompElSide large = thisside.LowerLevelElementList(1);
  if(large.Exists()) {
    RemoveSideRestraintWithRespectTo(icon,large);
//	elvec.Push(large);
  }

  /**Redimensioning the block of the connect icon and copying the adequated
     values from discontinuous connect*/
  elvec.Push(thisside);
  int condition = StablizingConnectContinuous(icon,elvec);
  if(condition) PZError << "TCompElT2dGD::MakeConnectContinuous is not done.\n";

  /**Restablecing restraint into the large*/
  for(i=0;i<n;i++) {
    TPZInterpolatedElement *el = (TPZInterpolatedElement *)smallelvec[i].Element();
    el->RestrainSide(smallelvec[i].Side(),this,icon);
  }
  if(large.Exists()) RestrainSide(icon,(TPZInterpolatedElement *)large.Element(),large.Side());
}

/**To redimensioning the block of the icon connect that is being actived
   in continuous and copying the adequated values from discontinuous connect*/
int TCompElT2dGD::StablizingConnectContinuous(int icon,TPZStack<TPZCompElSide> &elvec) {
  int i, j;
  int n = elvec.NElements();
  /**Verificando se todos os elementos em elvec tem o connect continuo ou descontinuo*/
  for(i=0;i<n;i++)
//	if(fConnectDisc[icon] == ((TPZInterpolatedElement *)elvec[i].Element())->IsConnectContinuous(elvec[i].Side()))
	//  return 1;   //existe um vizinho com continuidade diferente
  if(!fConnectDisc[icon]) return 0;   //Ja eh continuo.

  TPZCompMesh *mesh = Mesh();
  TPZBlock<STATE> &block = mesh->Block();
  int nvar = Material()->NStateVariables();
  TPZConnect &c = Connect(icon);
  int seqnum = c.SequenceNumber();
  int rows = TPZCompElT2d::NConnectShapeF(icon)*nvar;
  int cols = mesh->Solution().Cols();

  TPZFMatrix<STATE> values(rows,cols,0.);

  /**Redimensioning the block of the actual continuous connect*/
  block.Set(seqnum,rows);
  mesh->ExpandSolution();

  /**Making continuous all the neighbour elements*/
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
void TCompElT2dGD::SetConnectContinuous(int icon,TPZFMatrix<STATE> &values) {
  #ifndef NOTDEBUG
  if(icon<0 || icon>6) {
    PZError << "TCompElT2dGD::SetConnectContinuous. Bad parameter side.\n";
    return;
  }
  #endif

  if(!fConnectDisc[icon]) return;

  /**Modify the data structure of the connect to be discontinuous.
     Interchange the values in icon connect block to fConnect block*/
  TPZBlock<STATE> &block = Mesh()->Block();
  /**Searching the position and size of the values to connect(side) that must to
     be into the block corresponding to internal connect*/
  int nvar = 1;
  if(Material()) nvar = Material()->NStateVariables();
  int nshape = NShapeF();
  int cols = Mesh()->Solution().Cols();
  TPZFMatrix<STATE> elsol(nshape*nvar,cols,0.);
  int nconnect = Reference()->NSides();
  int rows, con, conseq;
  int i,j,pos=0;
  for(con = 0; con < nconnect; con++) {
    if(!fConnectDisc[con]) {
      /**sequence numbers of the intervening connects and auxiliar stores of data to transfer*/
      conseq = Connect(con).SequenceNumber();
      rows = block.Size(conseq);
      for(i=0;i<rows;i++) {
        for(j=0;j<cols;j++) {
          elsol(pos+i,j) = block(conseq,0,i,j);
		}
      }
    }
    pos += TPZCompElT2d::NConnectShapeF(con)*nvar;
  }
  int discpos=0;
  pos = 0;
  conseq = Connect(7).SequenceNumber();
  for(con = 0; con < nconnect; con++) {
	/** Recuperando os valores em values do connect descontinuo que passara a ser continuo */
	if(con == icon) {
      rows = TPZCompElT2d::NConnectShapeF(icon)*nvar;
      for(i=0;i<rows;i++) {
        for(j=0;j<cols;j++) {
          values(i,j) += block(conseq,0,discpos+i,j);
		}
      }
      discpos += rows;
	  continue;
	}
    if(fConnectDisc[con]) {
      /**sequence numbers of the intervening connects and auxiliar stores of data to transfer*/
      rows = TPZCompElT2d::NConnectShapeF(con)*nvar;
      for(i=0;i<rows;i++) {
        for(j=0;j<cols;j++) {
          elsol(pos+i,j) = block(conseq,0,discpos+i,j);
		}
      }
      discpos += rows;
    }
    pos += TPZCompElT2d::NConnectShapeF(con)*nvar;
  }

  /** Adapting the datastructure of the connect */
  fConnectDisc[icon] = 0;

  /**Redimensioning the block of the fConnect connect*/
  int newblocksize = NConnectShapeF(7)*nvar;
//  conseq = Connect(7).SequenceNumber();
  block.Set(conseq,newblocksize);

  /**Adequating solution vector*/
  Mesh()->ExpandSolution();

  /** Put the solution back for the discontinuous connect*/
  discpos=0;
  pos = 0;
  for(con=0;con<nconnect;con++) {
    if(fConnectDisc[con] && con!=icon) {
      /**sequence numbers of the intervening connects and auxiliar stores of data to transfer*/
      rows = TPZCompElT2d::NConnectShapeF(con)*nvar;
      for(i=0;i<rows;i++) {
        for(j=0;j<cols;j++) {
          block(conseq,0,discpos+i,j) = elsol(pos+i,j);
		}
      }
      discpos += rows;
    }
    pos += TPZCompElT2d::NConnectShapeF(con)*nvar;
  }
  /**Ordening the relations of the connects with the shape functions*/
  SetPositions();
}

/**sets the interpolation order of side to order
This method only updates the datastructure of the element and
updates the blocksize of the associated connect object */
void TCompElT2dGD::SetSideOrder(int side, int order) {
  if(IsConnectContinuous(side))
	  TPZCompElT2d::SetSideOrder(side,order);
	else {
	  fSideOrder[side-3] = order;
		int seqnum = Connect(7).SequenceNumber();
		int nvar = 1;
		if(Material()) nvar = Material()->NStateVariables();
	  Mesh()->Block().Set(seqnum,NConnectShapeF(7)*nvar);
	}
	SetPositions();
}
		
void TCompElT2dGD::SetConnectIndex(int i,int connectindex) {
  if(i>-1 && i<7) TPZCompElT2d::SetConnectIndex(i,connectindex);
  else if(i==7) fConnect = connectindex;
  else {
    PZError << "TCompElT2d::SetConnectIndex. Bad parameter i = " << i << " .\n";
    PZError.flush();
  }
}

int TCompElT2dGD::ConnectIndex(int i) {
  if(i>-1 && i<7) return TPZCompElT2d::ConnectIndex(i);
  if(i==7) return fConnect;
  PZError << "TPZCompElT2d::ConnectIndex. Bad parameter i = " << i << " .\n";
  PZError.flush();
  return -1;
}

int TCompElT2dGD::NConnectShapeF(int connect) {
  if(connect<0 || connect>7) {
    PZError << "TCompElT2dGD::NConnectShapeF. Bad parameter connect.\n";
    return 0;
  }
  if(connect==7) {
    if(!SideOrder(6)) return 1;
    int i,num=0;
    for(i=0;i<7;i++)
      if(fConnectDisc[i]) num += TPZCompElT2d::NConnectShapeF(i);
    return num;
  }
  if(fConnectDisc[connect]) return 0;
  return TPZCompElT2d::NConnectShapeF(connect);
}

int TCompElT2dGD::NShapeF() {
  if(!SideOrder(6)) return 1;
  return TPZInterpolatedElement::NShapeF();
}

int TCompElT2dGD::NSideConnects(int side) {
#ifndef NOTDEBUG
  if(side<0 || side>6) {
    PZError << "TCompElT2dGD::NSideConnects. Bad parameter side.\n";
    return 0;
  }
#endif
  if(side<3) return 2;
  else if(side<6) return 4;
  else if(side==6) return 8;
  return 0;
}

int TCompElT2dGD::SideConnectLocId(int c,int side) {
  switch(side) {
  case 0:
  case 1:
  case 2:
    if(!c) return side;
    return 7;
  case 3:
  case 4:
  case 5:
    if(!c) return side-3;
    if(c==1) return (side-2)%3;
    if(c==2) return side;
    return 7;
  case 6:
    return c;
  default:
    PZError << "TPZCompElT2d::SideConnectLocId, connect = " << c << std::endl;
    return -1;
  }
}

/**Ordening the positions into matrix phi of the shape functions to connects*/
void TCompElT2dGD::SetPositions() {
  int i,n=0;
  /**To continuous connects is established the positions*/
  for(i=0;i<7;i++) {
    if(!fConnectDisc[i]) {
      fPositions[i] = n;
      n += TPZCompElT2d::NConnectShapeF(i);
    }
  }
  /**The first position to discontinuous connects is stored into the position
     of the internal connect*/
  fPositions[7] = n;
  /**To discontinuous connects is established the positions*/
  for(i=0;i<6;i++) {
    if(fConnectDisc[i]) {
      fPositions[i] = n;
      n += TPZCompElT2d::NConnectShapeF(i);
    }
  }
  if(fConnectDisc[6]) fPositions[6] = n;
}

void TCompElT2dGD::MakeDiscontinuous(TPZVec<int> &disc) {
  if(disc.NElements()!=7) {
    PZError << "TCompElT2dGD::SetDiscontinuities. Vector of discontinuities is incompatible.\n";
    return;
  }
  for(int icon=0;icon<7;icon++) {
    if(disc[icon]) MakeConnectDiscontinuous(icon);
    else MakeConnectContinuous(icon);
  }
}

void TCompElT2dGD::SetConnectDiscontinuous(int icon) {
#ifndef NOTDEBUG
  if(icon<0 || icon>6) {
    PZError << "TCompElT2dGD::SetConnectDiscontinuous. Bad parameter side.\n";
    return;
  }
#endif

  if(fConnectDisc[icon]) return;

  /**Modify the data structure of the connect to be discontinuous.
     Interchange the values in icon connect block to fConnect block*/
  TPZBlock<STATE> &block = Mesh()->Block();
  /**Searching the position and size of the values to connect(side) that must to
     be into the block corresponding to internal connect*/
  int nvar = 1;
  if(Material()) nvar = Material()->NStateVariables();
  int nshape = NShapeF();
  int cols = Mesh()->Solution().Cols();
  TPZFMatrix<STATE> elsol(nshape*nvar,cols);
  int nconnect = Reference()->NSides();
  int con;
  int i,j,pos=0;

  int conseq;

  for(con = 0; con < nconnect; con++) {
    if(!fConnectDisc[con]) {
      /**sequence numbers of the intervening connects and auxiliar stores of data to transfer*/
      conseq = Connect(con).SequenceNumber();
      int rows = block.Size(conseq);
      for(i=0;i<rows;i++) {
        for(j=0; j< cols; j++) {
          elsol(pos+i,j) = block(conseq,0,i,j);
		}
      }
    }
    pos += TPZCompElT2d::NConnectShapeF(con)*nvar;
  }
  int discpos=0;
  pos = 0;

  conseq = Connect(7).SequenceNumber();
  for(con = 0; con < nconnect; con++) {
    if(fConnectDisc[con]) {
      /**sequence numbers of the intervening connects and auxiliar stores of data to transfer*/
      int rows = TPZCompElT2d::NConnectShapeF(con)*nvar;
      for(i=0;i<rows;i++) {
        for(j=0; j< cols; j++) {
          elsol(pos+i,j) = block(conseq,0,discpos+i,j);
		}
      }
      discpos += rows;
    }
    pos += TPZCompElT2d::NConnectShapeF(con)*nvar;
  }

  /** Adapting the datastructure of the icon connect */
  fConnectDisc[icon] = 1;
  /**Redimensioning the block of the fConnect connect*/
  int newblocksize = NConnectShapeF(7)*nvar;
  conseq = Connect(7).SequenceNumber();
  block.Set(conseq,newblocksize);

  /**Adequating solution vector*/
  Mesh()->ExpandSolution();

  /** Put the solution back for the discontinuous connect*/
  discpos=0;
  pos = 0;
  for(con = 0; con < nconnect; con++) {
    if(fConnectDisc[con]) {
      /**sequence numbers of the intervening connects and auxiliar stores of data to transfer*/
      int rows = TPZCompElT2d::NConnectShapeF(con)*nvar;
      for(i=0;i<rows;i++) {
        for(j=0; j< cols; j++) {
          block(conseq,0,discpos+i,j) = elsol(pos+i,j);
		}
      }
      discpos += rows;
    }
    pos += TPZCompElT2d::NConnectShapeF(con)*nvar;
  }

  /**Ordening the relations of the connects with the shape functions*/
  SetPositions();
}

int TCompElT2dGD::IsConnectContinuous(int icon) {
#ifndef NOTDEBUG
  if(icon<0 || icon>6) {
    PZError << "TCompElT2dGD::IsSideContinuous. Bad parameter side.\n";
    return 1;
  }
#endif
  return !(fConnectDisc[icon]);
}

void TCompElT2dGD::Shape(TPZVec<REAL> &x,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {
  if(!SideOrder(6)) {
    phi(0,0) = 1.;
    dphi(0,0) = 0.;
    dphi(1,0) = 0.;
    return;
  }
  TPZCompElT2d::Shape(x,phi,dphi);
  OrdeningPhi(phi,dphi);
}

void TCompElT2dGD::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {
  #ifndef NOTDEBUG
  if(side<0 || side>6) PZError << "TCompElT2dGD::SideShapeFunction. Bad paramenter side.\n";
  #endif
  if(!SideOrder(6)) {
    phi(0,0) = 1.;
    if(side<3) dphi.Redim(1,0);
    else if(side<6) {
      dphi.Redim(1,1);
      dphi(0,0) = 0.;
    }
    else {
      dphi.Redim(2,1);
      dphi(0,0) = dphi(1,0) = 0.;
    }
    return;
  }
  phi.Zero();
  dphi.Zero();
  switch(side) {
    case 0:
      phi(0,0) = 1.;
      dphi.Redim(1,0);
      return;
    case 1:
    case 2:
    {
      int pos = fPositions[side]-fPositions[7];
      pos = (pos<0) ? 0 : pos;
      phi(pos,0) = 1.;
      dphi.Redim(1,0);
      return;
    }
    case 6:
      Shape(point,phi,dphi);
      return;
    case 3:
    case 4:
    case 5:
    {
      int posdisc = fPositions[7];
      int dif = fPositions[side] - posdisc;
      int pos = 0;
      TPZVec<int64_t> id(2);
	  TPZVec<int> ord(1);
	  ord[0] = SideOrder(side) + 1;
      id[0] = Reference()->NodeIndex(side-3);
      id[1] = Reference()->NodeIndex((side-2)%3);
      TPZFMatrix<STATE> auxphi(phi), auxdphi(dphi);
      pzshape::TPZShapeLinear::Shape(point,id,ord,auxphi,auxdphi);
      int nodloc = SideConnectLocId(0,side);
      int dif0 = fPositions[nodloc] - posdisc;
      int dif1 = fPositions[(nodloc+1)%3] - posdisc;
      if(dif0<0) {
        phi(pos,0) = auxphi(0,0);
        dphi(0,pos++) = auxdphi(0,0);
      }
      if(dif1<0) {
        phi(pos,0) = auxphi(1,0);
        dphi(0,pos++) = auxdphi(0,1);
      }
      int i,nshape2 = TPZCompElT2d::NConnectShapeF(side);
      if(dif<0) {
        for(i=0;i<nshape2;i++) {
          dphi(0,pos) = auxdphi(0,2+i);
          phi(pos++,0) = auxphi(2+i,0);
        }
        if(pos==nshape2+2) return;
      }
      if(!(dif0<0)) {
        dif0 += pos;
        phi(dif0,0) = auxphi(0,0);
        dphi(0,dif0) = auxdphi(0,0);
      }
      if(!(dif1<0)) {
        dif1 += pos;
        phi(dif1,0) = auxphi(1,0);
        dphi(0,dif1) = auxdphi(0,1);
      }
      if(!(dif<0)) {
        dif += pos;
        for(i=0;i<nshape2;i++) {
          dphi(0,dif) = auxdphi(0,2+i);
          phi(dif++,0) = auxphi(2+i,0);
        }
      }
      return;
    }
    default :
      PZError << "TCompElT2dGD::SideShapeFunction bad parameter side.\n";
  }
}

void TCompElT2dGD::OrdeningPhi(TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {

  TPZFMatrix<STATE> auxphi = phi;
  TPZFMatrix<STATE> auxdphi = dphi;
  int i,j,n,m=0,nsideshapes;

  /**fSideDisc[i] store the first position of the shape functions associated with connect i */
  for(i=0;i<7;i++) {
    n = fPositions[i];
    nsideshapes = TPZCompElT2d::NConnectShapeF(i);
    for(j=0;j<nsideshapes;j++) {
      phi(n+j,0) = auxphi(m,0);
      dphi(0,n+j) = auxdphi(0,m);
      dphi(1,n+j) = auxdphi(1,m++);
    }
  }
}

void TCompElT2dGD::Print(std::ostream &out) {
  out << std::endl << "Discontinuous element\n";
  TPZInterpolatedElement::Print(out);
}

void TCompElT2dGD::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
	  if(dimension == 2) new TPZGraphElTd(this,&grmesh);
}

void TCompElT2dGD::EvaluateError(void (*fp)(TPZVec<REAL>&loc,TPZVec<REAL> &val,TPZFMatrix<STATE> &deriv),
        REAL &true_error,REAL &L2_error,TPZBlock<STATE>* /*flux*/,REAL &estimate) {

  true_error=0.; L2_error=0.; estimate=0.;
  if(Material() == NULL){
     PZError << "TCompElT2dGD::EvaluateError : no material for this element\n";
     Print(PZError);
     return;
  }
// Adjust the order of the integration rule
  TPZIntPoints &intrule = GetIntegrationRule();
  int dimension = Dimension();

  TPZManVector<int> prevorder(dimension), order(dimension);
  intrule.GetOrder(prevorder);

  TPZManVector<int> interpolation(8);
  GetInterpolationOrder(interpolation);

// compute the interpolation order of the shapefunctions squared
// This is wrong!
  int maxorder = interpolation[0];
  int dim;
  for(dim=0; dim<interpolation.NElements(); dim++) {
    maxorder = maxorder < interpolation[dim] ? interpolation[dim] : maxorder;
  }
  for(dim=0; dim<dimension; dim++) {
	  order[dim] = 2*maxorder+2;
  }
  intrule.SetOrder(order);


   int ndof = Material()->NStateVariables();
   int nflux = Material()->NFluxes();
   TPZFMatrix<STATE> phi(NShapeF(),1);
   TPZFMatrix<STATE> dphi(dimension,NShapeF());
   REAL jacobianstore[9],axesstore[9],xstore[3];
   TPZFMatrix<STATE> jacobian(dimension,dimension,jacobianstore,9);
   TPZFMatrix<REAL> axes(3,3,axesstore,9);
   TPZManVector<REAL> x(3);
   TPZVec<REAL> u_exact(ndof,0.);
   TPZFMatrix<STATE> du_exact(dimension,ndof,0.);
   TPZVec<REAL> intpoint(dimension,0.);
   TPZVec<REAL> values(3,0.);
   REAL detjac,weight;
   TPZManVector<REAL> u(ndof,0.);
   TPZFMatrix<STATE> dudx(dimension,ndof,0.);
   TPZVec<REAL> flux_el(nflux,0.);

   int ncon = NConnects();
   int nshape = NShapeF();
   TPZBlock<STATE> &block = Mesh()->Block();
   TPZFMatrix<STATE> dphix(dimension,nshape);
   TPZFMatrix<STATE> jacinv(dimension,dimension);
   int ieq;

   for(int nint=0; nint<GetIntegrationRule().NPoints(); nint++) {

      GetIntegrationRule().Point(nint,intpoint,weight);
	  Reference()->Jacobian( intpoint , jacobian, axes, detjac , jacinv);
      Shape(intpoint,phi,dphi);
	  Reference()->X( intpoint , x);
      weight *= fabs(detjac);
      switch(dimension) {
         case 0:
            break;
         case 1:
          dphix = dphi*(1./detjac);
          break;
         case 2:
            for(ieq = 0; ieq < nshape; ieq++) {
               dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
               dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
            }
            break;
         case 3:
         default:
          PZError << "pzintel.c please implement the " << dimension << "d Jacobian and inverse\n";
          PZError.flush();
      }
      int iv=0,in,jn,d;
      TPZConnect *df;
      u.Fill(0.);
      dudx.Zero();
      for(in=0; in<ncon; in++) {
          df = &Connect(in);
          int dfseq = df->SequenceNumber();
          int dfvar = block.Size(dfseq);
          for(jn=0; jn<dfvar; jn++) {
            u[iv%ndof] += phi(iv/ndof,0)*block(dfseq,0,jn,0);
            for(d=0; d<dimension; d++){
               dudx(d,iv%ndof) += dphix(d,iv/ndof)*block(dfseq,0,jn,0);
            }
          iv++;
          }
      }
      TPZVec<REAL> sol(ndof,0.);
	  TPZMaterialData data;
	  data.sol[0] = u;
	  data.dsol[0] = dudx;
	  data.axes = axes;
	  Material()->Solution(data,0,sol);
	  Material()->Flux(x,u,dudx,axes,flux_el);
      //contribuções dos erros
      if(fp) {
          fp(x,u_exact,du_exact);
		  Material()->Errors(x,u,dudx,axes,flux_el,u_exact,du_exact,values);
          true_error += values[0]*weight;
          L2_error += values[1]*weight;
          estimate = (estimate<values[2])?values[2]:estimate;
      }
   }//fim for : integration rule
   //Norma sobre o elemento
   L2_error = sqrt(L2_error);
   intrule.SetOrder(prevorder);
}

/*******                                                             *******/
/***    Methods to discontinuous computational one-dimensional element   ***/

TPZCompEl *TCompElT2dWI::CreateElDiscWI(TPZGeoElT2d *gel,TPZCompMesh &mesh,int64_t &index) {
  return new TCompElT2dWI(mesh,gel,index);
}

TCompElT2dWI::TCompElT2dWI(TPZCompMesh &mesh, TPZGeoElT2d *ref, int64_t &index)
  : TCompElT2dGD(mesh,ref,index) {

  int side;
  /**We do not let to create interface over pontoal sides */
  for(side=0;side<3;side++) fInterface[side] = -1;
  for(side=3;side<7;side++)
    fInterface[side] = CreateInterface(side,mesh);
}

TCompElT2dWI::TCompElT2dWI(TPZCompMesh &mesh, TPZGeoElT2d *ref, int64_t &index,TPZVec<int> &continuous)
  : TCompElT2dGD(mesh,ref,index,continuous) {

  int side;
  for(side=0;side<3;side++) fInterface[side] = -1;
  for(side=3;side<7;side++)
    fInterface[side] = CreateInterface(side,mesh);
}

TCompElT2dWI::~TCompElT2dWI() {
  DeleteInterfaces();
}
void TCompElT2dWI::Divide(int64_t index,TPZVec<int64_t> &pv,int interpolatesolution) {
	DeleteInterfaces();
	TPZInterpolatedElement::Divide(index,pv,interpolatesolution);
}
void TCompElT2dWI::DeleteInterfaces() {
  TInterfaceElement *inter_face;
//	cout << "Delete interface to element Id " << Index() << endl;
  for(int side=0;side<7;side++) {
    int index = fInterface[side];
    fInterface[side] = -1;         // Jorge 17/5/2000
    if(index==-1) continue;
    else if(index==-2) {
      /**Entao este elemento tem elementos de nivel maior conectados ao longo do lado side*/
      TPZStack<TPZCompElSide> elvec;
      TPZCompElSide thisside(this,side);
      thisside.HigherLevelElementList(elvec,1,0);
/*      TPZGeoElSide gelside = thisside.Reference();
      /**Os elementos comp com interface apenas podem ter um vizinho do mesmo
         nivel em cada lado
      gelside = gelside.Neighbour();
      if(!gelside.Exists()) PZError << "TCompElT2dWI::DeleteInterfaces. Not exist neighboard.\n";
      gelside.GetCompSubElements(elvec,1,0);*/
      int neighs = elvec.NElements();
      for(int i=0;i<neighs;i++) {
        TPZInterpolatedElement *neighcel = dynamic_cast<TPZInterpolatedElement *> (elvec[i].Element());
        index = neighcel->Interface(elvec[i].Side());
		inter_face = dynamic_cast<TInterfaceElement *>(fMesh->ElementVec()[index]);
        if(!inter_face) continue;
//				cout << "Deleting interface with index " << index << endl;
        delete inter_face;
        fMesh->ElementVec()[index] = 0;
        neighcel->SetInterface(elvec[i].Side(),-1);
      }
    }
    else {
    /** When one interface element is deleted it delete respective fInterface and
        change into Right computational element by -1 the respective fInterface*/
		inter_face = dynamic_cast<TInterfaceElement *>(fMesh->ElementVec()[index]);
      if(!inter_face) continue;
      if(this==inter_face->LeftEl()) {
        if(inter_face->RightEl()->Interface(inter_face->RightSide())!=-2)
			inter_face->RightEl()->SetInterface(inter_face->RightSide(),-1);
      }
      else {
        if(inter_face->LeftEl()->Interface(inter_face->LeftSide())!= -2)
			inter_face->LeftEl()->SetInterface(inter_face->LeftSide(),-1);
      }
//			cout << "Deleting interface with index " << index << endl;
      delete inter_face;
      fMesh->ElementVec()[index] = 0;
    }
  }
}

/** If fInterface[side] == -2 means exist several neighboards of the lower level*/
int TCompElT2dWI::CreateInterface(int side,TPZCompMesh &mesh) {
  int i,  nneighs;
  int64_t index = -1, indexneigh = -1;

  TPZStack<TPZCompElSide> elvec;
/**Creating the interfaces over side with one neighboards discontinuous element*/
  /** Searching over the side neighboards with same level. Apenas eh permitida
	    a construcao de interface entre elementos do mesmo nivel */
  TPZCompElSide thisside(this,side);
  thisside.EqualLevelElementList(elvec,1,1);
  nneighs = elvec.NElements();
  switch(nneighs) {
  case 0: {
    TPZCompElSide large = thisside.LowerLevelElementList(1);
    if(large.Exists()) {
      if(!((TPZInterpolatedElement *)large.Element())->CanHaveInterface()) return index;
      elvec.Push(large);
      indexneigh = -2;
      nneighs = 1;
    }
    else {
      TPZGeoElSide gelside = thisside.Reference();
//      thisside.HigherLevelElementList(elvec,1,0);
      /**Os elementos comp com interface apenas podem ter um vizinho do mesmo
         nivel em cada lado*/
      gelside = gelside.Neighbour();
      if(!gelside.Exists()) return -1;
      gelside.HigherLevelCompElementList2(elvec,1,0);
      nneighs = elvec.NElements();
      if(!nneighs) return -1;
      for(i=0;i<nneighs;i++)
        if(!((TPZInterpolatedElement *)elvec[i].Element())->CanHaveInterface())
	        return index;
      index = -2;
    }
    break;
  }
  case 1:
    if(!((TPZInterpolatedElement *)elvec[0].Element())->CanHaveInterface()) return -1;
    break;
  default:
  /**Asking if exist only one neighbour and whether it can to has interface*/
    return -1;
  }

  /**The dimension of the thisside must to be equal to the dimension of the common
     side on the neighboard element in elvec*/
  int dim = thisside.Reference().Dimension();
  int dimmat = Material()->Dimension();
  if(dim != elvec[0].Reference().Dimension() || dim < dimmat-1) return -1;
	TInterfaceElement *interf;
  int fluxtype = Material()->FluxType();
  /** Creating interface over this side boundary. The neighboard element cann't
      to has interface element over this common side boundary*/
//  TPZInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *geo, int64_t &index, TPZCompElSide & left, TPZCompElSide &right);

  for(i=0;i<nneighs;i++) {
    if(index==-2) interf = new TInterfaceElement(mesh,elvec[i],thisside,fluxtype,indexneigh);
    else {
      interf = new TInterfaceElement(mesh,thisside,elvec[i],fluxtype,index);
      if(indexneigh!=-2) indexneigh = index;
    }
//		cout << "Creating interface with index " << interf->Index() << endl;
    ((TPZInterpolatedElement *)elvec[i].Element())->SetInterface(elvec[i].Side(),indexneigh);
  }
  return index;
}

void TCompElT2dWI::Print(std::ostream &out) {
  TCompElT2dGD::Print(out);
  /** Information over interfaces*/
  int i;
  out << "Interfaces information :\n";
  for(i=0;i<7;i++) {
    out << "   Side " << i << " : ";
    if(fInterface[i]==-1) out << "no interface\n";
    else out << "interface " << fInterface[i] << std::endl;
  }
}

int TCompElT2dWI::Interface(int side) {
#ifndef NOTDEBUG
  if(side<0 || side>6) {
    PZError << "TCompElT2dWI::Interface. Bad side.\n";
    return -1;
  }
#endif NOTDEBUG
  return fInterface[side];
}

void TCompElT2dWI::SetInterface(int side,int64_t indexinterface) {
#ifndef NOTDEBUG
  if(side<0 || side>6) {
    PZError << "TCompElT2dWI::SetInterface. Bad side.\n";
    return;
  }
#endif NOTDEBUG
  fInterface[side] = indexinterface;
}

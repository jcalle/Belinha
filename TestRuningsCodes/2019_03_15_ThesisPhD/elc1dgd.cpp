#include "elc1dgd.h"
#include "interface.h"
#include "pzelmat.h"
#include "pzblock.h"
#include "pzcmesh.h"
#include "pzconnect.h"
#include "pzquad.h"
#include "pzelg1d.h"
#include "pzgraphel1dd.h"
#include "myheader.h"
#include "TPZMaterial.h"


TPZCompEl *TCompEl1dGD::CreateElDisc(TPZGeoEl1d *gel,TPZCompMesh &mesh,int64_t &index) {
  return new TCompEl1dGD(mesh,gel,index);
}

TCompEl1dGD::TCompEl1dGD(TPZCompMesh &mesh, TPZGeoEl1d *ref, int64_t &index)
  : TPZCompEl1d(mesh,ref,index,0) {

  int icon;

  /**Determining the discontinuity upon the neighbours through icon connect*/
  for(icon=0;icon<3;icon++)
    fConnectDisc[icon] = NeighbourDiscontinuous(icon);
  /**Creating the new internal connect, it belongs only to the current discontinuous element*/
  int nvar = 1;
  if(Material())
    nvar = Material()->NStateVariables();
  fConnect = mesh.AllocateNewConnect(NShapeF(),nvar,GetPreferredOrder());
  TPZConnect &connect = mesh.ConnectVec()[fConnect];
  mesh.Block().Set(connect.SequenceNumber(),nvar*NConnectShapeF(3));
  connect.IncrementElConnected();
  SetPositions();

  /**Creating the connect over the sides of the current element*/
  for(icon=0;icon<3;icon++) {
    SetConnectIndex(icon,CreateMidSideConnect(icon));
    mesh.ConnectVec()[ConnectIndex(icon)].IncrementElConnected();
  }

  TPZVec<int> order(1,2*SideOrder(2));
  GetIntegrationRule().SetOrder(order);

  ref->SetReference(this);

  SetPositions();
}

TCompEl1dGD::TCompEl1dGD(TPZCompMesh &mesh, TPZGeoEl1d *ref, int64_t &index,TPZVec<int> &maskdisc)
  : TPZCompEl1d(mesh,ref,index,0) {

  int icon;

  /**Determining the discontinuity upon the neighbours through icon connect*/
  for(icon=0;icon<3;icon++)
    fConnectDisc[icon] = NeighbourDiscontinuous(icon);
  /**Creating the new internal connect, it belongs only to the current discontinuous element*/
  int nvar = 1;
  if(Material())
    nvar = Material()->NStateVariables();
  fConnect = mesh.AllocateNewConnect(NShapeF(), nvar, GetPreferredOrder());
  TPZConnect &connect = mesh.ConnectVec()[fConnect];
  mesh.Block().Set(connect.SequenceNumber(),nvar*NConnectShapeF(3));
  connect.IncrementElConnected();
  SetPositions();

  /**Creating the connect over the sides of the current element*/
  for(icon=0;icon<3;icon++) {
    SetConnectIndex(icon,CreateMidSideConnect(icon));
    mesh.ConnectVec()[ConnectIndex(icon)].IncrementElConnected();
  }

  TPZVec<int> order(1,2*SideOrder(2));
  GetIntegrationRule().SetOrder(order);

  ref->SetReference(this);
  /**Connects mask with maskdisc are doing discontinuous.
  Can to be continuous if exist a neighbour which is continuous element*/
  for(icon=0;icon<3;icon++) {
    if(maskdisc[icon]) {
      MakeConnectDiscontinuous(icon);
    }
    else {
      MakeConnectContinuous(icon);
    }
  }

  SetPositions();
}

TCompEl1dGD::~TCompEl1dGD() {
  TPZConnect *con = &Mesh()->ConnectVec()[fConnect];
  Mesh()->Block().Set(con->SequenceNumber(),0);
  con->ResetElConnected();
  Mesh()->ExpandSolution();
// if(!con->HasDependency())  Mesh()->ConnectVec().SetFree(fConnect);
}

int TCompEl1dGD::NeighbourDiscontinuous(int icon) {
  TPZStack<TPZCompElSide> elvec;
  TPZCompElSide thisside(this,icon);
  thisside.EqualLevelElementList(elvec,0,0);
  /**Asking if exist neighbour that it can not to be discontinuous*/
//  if(!SideOrder(icon)) return 1;
  if(elvec.NElements() && ((TPZInterpolatedElement *)(elvec[0].Element()))->ConnectIndex(elvec[0].Side()) < 0)
    return 0;
  return 1;
}

/**Make discontinuous the icon connect of the discontinuous element.
   It is tranfered the values associated from icon connect to fConnect*/
void TCompEl1dGD::MakeConnectDiscontinuous(int icon) {
  #ifndef NOTDEBUG
  if(icon<0 || icon>2) {
    PZError << "TCompEl1dGD::MakeConnectDiscontinuous. Bad parameter side.\n";
    return;
  }
  #endif
  if(fConnectDisc[icon]) return;

  TPZCompElSide thisside(this,icon);
  /** Check whether there is a neighbouring element which cannot be discontinuous*/
  TPZStack<TPZCompElSide> elvec;       // Porque mexe na fSolution??? !!!
  thisside.EqualLevelElementList(elvec,0,0);
  int j, n=elvec.NElements();
  /**Asking if exist neighbour that it can not to be discontinuous*
  for(j=0;j<n;j++) {
    if(!((TPZCompElDisc *)(elvec[j].Element()))->CanBeDiscontinuous()) {
      return;
    }
  }
  */
  /**Removing restrains respect to large element and in small element connected*/
  TPZStack<TPZCompElSide> smallelvec;
  TPZGeoElSide gelside = thisside.Reference();
  gelside = gelside.Neighbour();
  if(!gelside.Exists()) PZError << "TCompEl1dGD::MakeConnectDiscontinuous.. Not exist neighboard.";
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
//    elvec.Push(smallelvec[j]);
  }
  if(large.Exists()) RemoveSideRestraintWithRespectTo(icon,large);

  /**Making discontinuous all the neighbour elements*/
  elvec.Push(thisside);
  n = elvec.NElements();
  for(j=0;j<n;j++) {
    TPZInterpolatedElement *cel = (TPZInterpolatedElement *)(elvec[j].Element());
//    cel->SetConnectDiscontinuous(elvec[j].Side());
  }

  /**Redimensioning the block of the icon connect*/
  TPZConnect &c = Connect(icon);
  int seqnum = c.SequenceNumber();
  Mesh()->Block().Set(seqnum,0);
  Mesh()->ExpandSolution();

  /**Restablecing restraint into the large*/
  n = smallelvec.NElements();
  for(j=0;j<n;j++) {
    TPZInterpolatedElement *el = (TPZInterpolatedElement *)smallelvec[j].Element();
    el->RestrainSide(smallelvec[j].Side(),this,icon);
  }
  if(large.Exists()) RestrainSide(icon,(TPZInterpolatedElement *)large.Element(),large.Side());
}

/**Make continuous the icon connect of the discontinuous element
   It is tranfered the values associated from fConnect to icon connect*/
void TCompEl1dGD::MakeConnectContinuous(int icon) {
  #ifndef NOTDEBUG
  if(icon<0 || icon>2) {
    PZError << "TCompEl1dGD::MakeConnectContinuous. Bad parameter side.\n";
    return;
  }
  #endif
  if(!fConnectDisc[icon]) return;

  TPZCompElSide thisside(this,icon);

  /** Find neighboard elements in the same level of the current element */
  TPZStack<TPZCompElSide> elvec;
  thisside.EqualLevelElementList(elvec,0,0);
  int i,n;

  /**Removing restrains respect to large element and in small element connected*/
  TPZStack<TPZCompElSide> smallelvec;
  TPZGeoElSide gelside = thisside.Reference();
  gelside = gelside.Neighbour();
  TPZStack<TPZGeoElSide> gsmallvec;
  if(!gelside.Exists()) PZError << "TCompEl1dGD::MakeConnectContinuous.. Not exist neighboard.";
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
  if(condition) PZError << "TCompEl1dGD::MakeConnectContinuous is not done.\n";

  /**Restablecing restraint into the large*/
  for(i=0;i<n;i++) {
    TPZInterpolatedElement *el = (TPZInterpolatedElement *)smallelvec[i].Element();
    el->RestrainSide(smallelvec[i].Side(),this,icon);
  }
  if(large.Exists()) RestrainSide(icon,(TPZInterpolatedElement *)large.Element(),large.Side());
}

/**To redimensioning the block of the icon connect that is being actived
   in continuous and copying the adequated values from discontinuous connect*/
int TCompEl1dGD::StablizingConnectContinuous(int icon,TPZStack<TPZCompElSide> &elvec) {
  int i, j;
  int n = elvec.NElements();
  /**Verificando se todos os elementos em elvec tem o connect continuo ou descontinuo*/
  for(i=0;i<n;i++)
//	if(fConnectDisc[icon] == ((TPZInterpolatedElement *)elvec[i].Element())->IsConnectContinuous(elvec[i].Side()))
	  return 1;   //existe um vizinho com continuidade diferente
  if(!fConnectDisc[icon]) return 0;   //Ja eh continuo.

  TPZCompMesh *mesh = Mesh();
  TPZBlock<STATE> &block = mesh->Block();
  int nvar = 1;
//	if(fMaterial) fMaterial->NStateVariables();
  TPZConnect &c = Connect(icon);
  int seqnum = c.SequenceNumber();
  int rows = TPZCompEl1d::NConnectShapeF(icon)*nvar;
  int cols = mesh->Solution().Cols();

  TPZFMatrix<STATE> values(rows,cols,0.);

  /**Redimensioning the block of the actual continuous connect*/
  block.Set(seqnum,rows);
  mesh->ExpandSolution();

  /**Making continuous all the neighbour elements*/
  for(j=0;j<n;j++) {
    TPZInterpolatedElement *cel = (TPZInterpolatedElement *)(elvec[j].Element());
//    cel->SetConnectContinuous(elvec[j].Side(),values);
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
void TCompEl1dGD::SetConnectContinuous(int icon,TPZFMatrix<STATE> &values) {
  #ifndef NOTDEBUG
  if(icon<0 || icon>2) {
    PZError << "TCompEl1dGD::SetConnectContinuous. Bad parameter icon.\n";
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
  int cols=Mesh()->Solution().Cols();
  TPZFMatrix<STATE> elsol(nshape*nvar,cols);
  int nconnect = Reference()->NSides();
  int rows, con, conseq;
  int i,j,pos=0;
  for(con = 0; con < nconnect; con++) {
    if(!fConnectDisc[con]) {
      /**sequence numbers of the intervening connects and auxiliar stores of data to transfer*/
      conseq = Connect(con).SequenceNumber();
      rows = block.Size(conseq);
      for(i=0;i<rows;i++) {
        for(j=0; j< cols; j++) {
          elsol(pos+i,j) = block(conseq,0,i,j);
		}
      }
    }
    pos += TPZCompEl1d::NConnectShapeF(con)*nvar;
  }
  int discpos=0;
  pos = 0;
  conseq = Connect(3).SequenceNumber();
  for(con = 0; con < nconnect; con++) {
	/** Recuperando os valores em values do connect descontinuo que passara a ser continuo */
	if(con == icon) {
      rows = TPZCompEl1d::NConnectShapeF(icon)*nvar;
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
      rows = TPZCompEl1d::NConnectShapeF(con)*nvar;
      for(i=0;i<rows;i++) {
        for(j=0; j< cols; j++) {
          elsol(pos+i,j) = block(conseq,0,discpos+i,j);
		}
      }
      discpos += rows;
    }
    pos += TPZCompEl1d::NConnectShapeF(con)*nvar;
  }

  /** Adapting the datastructure of the connect */
  fConnectDisc[icon] = 0;

  /**Redimensioning the block of the fConnect connect*/
  int newblocksize = NConnectShapeF(3)*nvar;
//  int conseq = Connect(3).SequenceNumber();
  block.Set(conseq,newblocksize);
  
  /**Adequating solution vector*/
  Mesh()->ExpandSolution();
  
  /** Put the solution back for the continuous connect*/
  discpos=0;
  pos = 0;
  for(con = 0; con < nconnect; con++) {
    if(fConnectDisc[con] && con!=icon) {
      /**sequence numbers of the intervening connects and auxiliar stores of data to transfer*/
      rows = TPZCompEl1d::NConnectShapeF(con)*nvar;
      for(i=0;i<rows;i++) {
        for(j=0; j< cols; j++) {
          block(conseq,0,discpos+i,j) = elsol(pos+i,j);
		}
      }
      discpos += rows;
    }
    pos += TPZCompEl1d::NConnectShapeF(con)*nvar;
  }
  /**Ordening the relations of the connects with the shape functions*/
  SetPositions();
}

/**sets the interpolation order of side to order
This method only updates the datastructure of the element and
updates the blocksize of the associated connect object */
void TCompEl1dGD::SetSideOrder(int side, int order) {
  if(IsConnectContinuous(side))
	  TPZCompEl1d::SetSideOrder(side,order);
	else {
	  fSideOrder = order;
		int seqnum = Connect(3).SequenceNumber();
		int nvar = 1;
		if(Material()) nvar = Material()->NStateVariables();
	  Mesh()->Block().Set(seqnum,NConnectShapeF(3)*nvar);
	}
	SetPositions();
}
		
void TCompEl1dGD::SetConnectIndex(int i, int64_t connectindex) {
  if(i>-1 && i<3) TPZCompEl1d::SetConnectIndex(i,connectindex);
  else if(i==3) fConnect = connectindex;
  else {
    PZError << "TCompEl1d::SetConnectIndex. Bad parameter i = " << i << " .\n";
    PZError.flush();
  }
}

int64_t TCompEl1dGD::ConnectIndex(int i) {
  if(i>-1 && i<3) return TPZCompEl1d::ConnectIndex(i);
  if(i==3) return fConnect;
  PZError << "TCompEl1dGD::ConnectIndex. Bad parameter i = " << i << " .\n";
  PZError.flush();
  return -1;
}

int TCompEl1dGD::NConnectShapeF(int connect) {
  if(connect<0 || connect>3) {
    PZError << "TCompEl1dGD::NConnectShapeF. Bad parameter connect.\n";
    return 0;
  }
  if(connect==3) {
    if(!SideOrder(2)) return 1;
    int i,num=0;
    for(i=0;i<3;i++)
      if(fConnectDisc[i]) num += TPZCompEl1d::NConnectShapeF(i);
    return num;
  }
  if(fConnectDisc[connect]) return 0;
  return TPZCompEl1d::NConnectShapeF(connect);
}

int TCompEl1dGD::NShapeF() {
  if(!SideOrder(2)) return 1;
  return TPZInterpolatedElement::NShapeF();
}

int TCompEl1dGD::NSideConnects(int side) {
  if(side==0 || side==1) return 2;
  else if(side==2) return 4;
  PZError << "TCompEl1dGD::NSideConnects. Bad parameter side = " << side << " .\n";
  return 0;
}

int TCompEl1dGD::SideConnectLocId(int c,int side) {
  switch(side) {
  case 0:
  case 1:
    if(c)
      return 3;
    return side;
  case 2:
    return c;
  default:
    PZError << "TPZCompEl1d::SideConnectLocId called with side = " << side << std::endl;
    return -1;
  }
}

/**Ordening the positions into matrix phi of the shape functions to connects*/
void TCompEl1dGD::SetPositions() {
  int i,n=0;
  /**To continuous connects is established the positions*/
  for(i=0;i<3;i++) {
    if(!fConnectDisc[i]) {
      fPositions[i] = n;
      n += TPZCompEl1d::NConnectShapeF(i);
    }
  }
  /**Store the first position to discontinuous connects as the position of the internal connect*/
  fPositions[3] = n;
  /**To discontinuous connects is established the positions*/
  for(i=0;i<2;i++) {
    if(fConnectDisc[i]) {
      fPositions[i] = n;
      n += TPZCompEl1d::NConnectShapeF(i);
    }
  }
  if(fConnectDisc[2]) fPositions[2] = n;
}

void TCompEl1dGD::MakeDiscontinuous(TPZVec<int> &disc) {
  if(disc.NElements()!=3) {
    PZError << "TCompEl1dGD::SetDiscontinuities. Vector of discontinuities is incompatible.\n";
    return;
  }
  for(int icon=0;icon<3;icon++) {
    if(disc[icon]) MakeConnectDiscontinuous(icon);
    else MakeConnectContinuous(icon);
  }
}

void TCompEl1dGD::SetConnectDiscontinuous(int icon) {
#ifndef NOTDEBUG
  if(icon<0 || icon>2) {
    PZError << "TCompEl1dGD::SetConnectDiscontinuous. Bad parameter icon.\n";
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
    pos += TPZCompEl1d::NConnectShapeF(con)*nvar;
  }
  int discpos=0;
  pos = 0;

  conseq = Connect(3).SequenceNumber();
  for(con = 0; con < nconnect; con++) {
    if(fConnectDisc[con]) {
      /**sequence numbers of the intervening connects and auxiliar stores of data to transfer*/
      int rows = TPZCompEl1d::NConnectShapeF(con)*nvar;
      for(i=0;i<rows;i++) {
        for(j=0; j< cols; j++) {
          elsol(pos+i,j) = block(conseq,0,discpos+i,j);
		}
      }
      discpos += rows;
    }
    pos += TPZCompEl1d::NConnectShapeF(con)*nvar;
  }

  /** Adapting the datastructure of the icon connect */
  fConnectDisc[icon] = 1;
  /**Redimensioning the block of the fConnect connect*/
  int newblocksize = NConnectShapeF(3)*nvar;
  conseq = Connect(3).SequenceNumber();
  block.Set(conseq,newblocksize);

  /**Adequating solution vector*/
  Mesh()->ExpandSolution();

  /** Put the solution back for the discontinuous connect*/
  discpos=0;
  pos = 0;
  for(con = 0; con < nconnect; con++) {
    if(fConnectDisc[con]) {
      /**sequence numbers of the intervening connects and auxiliar stores of data to transfer*/
      int rows = TPZCompEl1d::NConnectShapeF(con)*nvar;
      for(i=0;i<rows;i++) {
        for(j=0; j< cols; j++) {
          block(conseq,0,discpos+i,j) = elsol(pos+i,j);
		}
      }
      discpos += rows;
    }
    pos += TPZCompEl1d::NConnectShapeF(con)*nvar;
  }

  /**Ordening the relations of the connects with the shape functions*/
  SetPositions();
}

int TCompEl1dGD::IsConnectContinuous(int icon) {
#ifndef NOTDEBUG
  if(icon<0 || icon>2) {
    PZError << "TCompEl1dGD::IsSideContinuous. Bad parameter side.\n";
    return 1;
  }
#endif
  return !(fConnectDisc[icon]);
}

void TCompEl1dGD::Shape(TPZVec<REAL> &x,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {
  if(!SideOrder(2)) {
    phi(0,0) = 1.;
    dphi(0,0) = 0.;
    return;
  }
  TPZCompEl1d::Shape(x,phi,dphi);
  OrdeningPhi(phi,dphi);
}

/** !!! Revistar TPZCompEl1d::SideShapeFunction, o valor de dphi */
void TCompEl1dGD::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi) {
  phi.Zero();
  dphi.Zero();
  switch(side) {
  case 0:
    phi(0,0) = 1.;
    return;
  case 1:
    if(fConnectDisc[0] && fConnectDisc[1] && SideOrder(2)) phi(1,0) = 1.;
    else phi(0,0) = 1.;
    return;
  case 2:
    Shape(point,phi,dphi);
    return;
  default:
    PZError << "TCompEl1dGD::SideShapeFunction. Bad parameter side.\n";
  }
}

void TCompEl1dGD::OrdeningPhi(TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphix) {

  if(phi.Rows()==1) return;
  TPZFMatrix<STATE> auxphi(phi);
  TPZFMatrix<STATE> auxdphix(dphix);
  int i,j,n,m=0,nsideshapes;

  /**fSideDisc[i] store the first position of the shape functions associated with connect i */
  for(i=0;i<3;i++) {
    n = fPositions[i];
    nsideshapes = TPZCompEl1d::NConnectShapeF(i);
    for(j=0;j<nsideshapes;j++) {
      phi(n+j,0) = auxphi(m,0);
      dphix(0,n+j) = auxdphix(0,m++);
    }
  }
}

void TCompEl1dGD::Print(std::ostream &out) {
  out << std::endl << "Discontinuous element\n";
  TPZInterpolatedElement::Print(out);
}

void TCompEl1dGD::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) {
  if(dimension == 1) new TPZGraphEl1dd(this,&grmesh);
}

void TCompEl1dGD::EvaluateError(void (*fp)(TPZVec<REAL>&loc,TPZVec<REAL> &val,TPZFMatrix<STATE> &deriv),
        REAL &true_error,REAL &L2_error,TPZBlock<STATE>* /*flux*/,REAL &estimate) {

  true_error=0.; L2_error=0.; estimate=0.;
  if(Material() == NULL){
     PZError << "TPZInterpolatedElement::EvaluateError : no material for this element\n";
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
	  data.axes = axes;
	  data.sol[0] = u;
	  data.dsol[0] = dudx;
	  Material()->Solution(data,0,sol);
	  Material()->Flux(x,u,dudx,axes,flux_el);
      //contribuções dos erros
      if(fp) {
          fp(x,u_exact,du_exact);
 //         TPZMaterial *matp = (TPZMaterial *) fMaterial;
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

/** Return total mass contained into the element
REAL TCompEl1dGD::MeanSolution(int var) {
#ifndef NOTDEBUG
  if(var<0 || var>Material()->NStateVariables()) {
    PZError << "TCompEl1dGD::MeanSolution. var is out range.\n";
    return 0.;
  }
#endif
  if(fMeanSolution[var]==-1.) return TPZInterpolatedElement::MeanSolution(var);
  return fMeanSolution[var];
}
*/
/*******                                                             *******/
/***    Methods to discontinuous computational one-dimensional element   ***/

TPZCompEl *TCompEl1dWI::CreateElDiscWI(TPZGeoEl1d *gel,TPZCompMesh &mesh, int64_t &index) {
  return new TCompEl1dWI(mesh,gel,index);
}

TCompEl1dWI::TCompEl1dWI(TPZCompMesh &mesh, TPZGeoEl1d *ref, int64_t &index)
  : TCompEl1dGD(mesh,ref,index) {

  int side;
  for(side=0;side<3;side++)
    fInterface[side] = CreateInterface(side,mesh);
}

TCompEl1dWI::TCompEl1dWI(TPZCompMesh &mesh, TPZGeoEl1d *ref, int64_t &index,TPZVec<int> &continuous)
  : TCompEl1dGD(mesh,ref,index,continuous) {

  int side;
  for(side=0;side<3;side++)
    fInterface[side] = CreateInterface(side,mesh);
}

TCompEl1dWI::~TCompEl1dWI() {
  DeleteInterfaces();
}

void TCompEl1dWI::Divide(int64_t index,TPZVec<int64_t> &pv,int interpolatesolution) {
	DeleteInterfaces();
	TPZInterpolatedElement::Divide(index,pv,interpolatesolution);
}
void TCompEl1dWI::DeleteInterfaces() {
	TInterfaceElement *inter_face;
//	cout << "Delete interface to element Id " << Index() << endl;
  for(int side=0;side<3;side++) {
	  int64_t index = fInterface[side];
    fInterface[side] = -1;          //Jorge 17/5/2000
    if(index==-1) continue;
    else if(index==-2) {
      /**Entao este elemento tem elementos de nivel maior conectados ao longo do lado side*/
      TPZStack<TPZCompElSide> elvec;
      TPZCompElSide thisside(this,side);
//      thisside.HigherLevelElementList(elvec,1,0);
      TPZGeoElSide gelside = thisside.Reference();
      /**Os elementos comp com interface apenas podem ter um vizinho do mesmo
         nivel em cada lado*/
	  int i;
      gelside = gelside.Neighbour();
      if(!gelside.Exists()) PZError << "TCompEl1dWI::DeleteInterfaces. Not exist neighboard.\n";
	  TPZStack<TPZGeoElSide> gsmallvec;
	  gelside.GetSubElements2(gsmallvec);
	  for (i = 0; i < gsmallvec.NElements(); i++)
		  elvec.Push(gsmallvec.Pop().Reference());

      int neighs = elvec.NElements();
      for(i=0;i<neighs;i++) {
        TPZInterpolatedElement *neighcel = (TPZInterpolatedElement *)elvec[i].Element();
        index = neighcel->Interface(elvec[i].Side());
        inter_face = (TInterfaceElement *)fMesh->ElementVec()[index];
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
      inter_face = (TInterfaceElement *)fMesh->ElementVec()[index];
      if(!inter_face) continue;
//      if(this==inter_face->LeftElement()) {
  //      if(inter_face->RightEl()->Interface(inter_face->RightSide())!=-2)
    //      inter_face->RightEl()->SetInterface(inter_face->RightSide(),-1);
      //}
//      else {
  //      if(inter_face->LeftElement()->Interface(inter_face->LeftElementSide())!= -2)
    //      inter_face->LeftElement()->SetInterface(inter_face->LeftElementSide(),-1);
      //}
//			cout << "Deleting interface with index " << index << endl;
      delete inter_face;
      fMesh->ElementVec()[index] = 0;
    }
  }
}

/** If fInterface[side] == -2 means exist several neighboards of the lower level*/
int64_t TCompEl1dWI::CreateInterface(int side,TPZCompMesh &mesh) {
  int i,  nneighs;
  int64_t index = -1;
  int64_t  indexneigh = -1;
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
      /**Os elementos comp com interface apenas podem ter um vizinho do mesmo
         nivel em cada lado*/
      gelside = gelside.Neighbour();
      if(!gelside.Exists()) return -1;
      gelside.HigherLevelCompElementList2(elvec,1,0);
//      thisside.HigherLevelElementList(elvec,1,0);
      nneighs = elvec.NElements();
      if(!nneighs) return -1;
      /**Verificando se todos os elementos computacionais podem ter interface*/
      for(i=0;i<nneighs;i++)
        if(!((TPZInterpolatedElement *)elvec[i].Element())->CanHaveInterface())
          return -1;
      index = -2;
    }
    break;
  }
  case 1:
    if(!((TPZInterpolatedElement *)elvec[0].Element())->CanHaveInterface()) return index;
    break;
  default:
  /**Asking if exist only one neighbour and whether it can to has interface*/
    return -1;
  }

  /**The dimension of the thisside must to be equal to the dimension of the common
     side on the neighboard element in elvec*/
  int dim = thisside.Reference().Dimension();
  int dimmat = Material()->Dimension();
  TInterfaceElement *interf;
  if(dim != elvec[0].Reference().Dimension() || dim < dimmat-1) return -1;
  int fluxtype = Material()->FluxType();
  /** Creating interface over this side boundary. The neighboard element cann't
      to has interface element over this common side boundary*/
  for(i=0;i<nneighs;i++) {
    if(index==-2) 
		interf = new TInterfaceElement(mesh,elvec[i],thisside,fluxtype,indexneigh);
    else {
      interf = new TInterfaceElement(mesh,thisside,elvec[i],fluxtype,index);
      if(indexneigh!=-2) indexneigh = index;
    }
//		cout << "Creating interface with index " << interf->Index() << endl;
    ((TPZInterpolatedElement *)elvec[i].Element())->SetInterface(elvec[i].Side(),indexneigh);
  }
  return index;
}

void TCompEl1dWI::Print(std::ostream &out) {
  TCompEl1dGD::Print(out);
  /** Information over interfaces*/
  int i;
  out << "Interfaces information :\n";
  for(i=0;i<3;i++) {
    out << "   Side " << i << " : ";
    if(fInterface[i]==-1) out << "no interface\n";
    else out << fInterface[i] << std::endl;
  }
}

int TCompEl1dWI::Interface(int side) {
#ifndef NOTDEBUG
  if(side<0 || side>2) {
    PZError << "TCompEl1dWI::Interface. Bad side.\n";
    return -1;
  }
#endif NOTDEBUG
  return fInterface[side];
}

void TCompEl1dWI::SetInterface(int side, int64_t indexinterface) {
#ifndef NOTDEBUG
  if(side<0 || side>2) {
    PZError << "TCompEl1dWI::SetInterface. Bad side.\n";
    return;
  }
#endif NOTDEBUG
  fInterface[side] = indexinterface;
}

/**Compute the contribution to stiffness matrix and load vector on the element*/
/*void TCompEl1dGD::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef) {

  if(fMaterial == NULL){
    PZError << "TPZInterpolatedElement::CalcStiff : no material for this element\n";
    Print(PZError);
    return;
  }
  if(fMaterial->Id()<0) {
    delete ef.fMat;
    ef.fMat = 0;
    return;
  }

  int i, numdof = fMaterial->NStateVariables();
  int ncon = NConnects(), dim = fMaterial->Dimension();
  int nshape = NShapeF();
  int numeq = nshape*numdof;
  
  TPZBlock &block = Mesh()->Block();
  // clean ek and ef
  if(!ek.fMat) ek.fMat = new TPZFMatrix();
  if(!ef.fMat) ef.fMat = new TPZFMatrix();
  if(!ek.fBlock) ek.fBlock = new TPZBlock(ek.fMat);
  if(!ef.fBlock) ef.fBlock = new TPZBlock(ef.fMat);
#ifndef NOTDEBUG
  if( !ek.fMat || !ef.fMat || !ek.fBlock || !ef.fBlock){
    cout << "TPZInterpolatedElement.calc_stiff : not enough storage for local stifness"
      " matrix \n";
    Print(cout);
    return;
  }
#endif

  ek.fMat->Redim(numeq,numeq);
  ef.fMat->Redim(numeq,1);
  ek.fBlock->SetNBlocks(ncon);
  ef.fBlock->SetNBlocks(ncon);
  TPZVec<REAL> sol(numdof);
  TPZFMatrix dsol(dim,numdof);
  for (i = 0; i < ncon ; i++)	{
    ek.fBlock->Set(i,NConnectShapeF(i)*numdof);
    ef.fBlock->Set(i,NConnectShapeF(i)*numdof);
  }

  int in,jn;
  TPZConnect *df;

  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);

  for(i=0; i<ncon; ++i){
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }

  TPZFMatrix phi(nshape,1);
  TPZFMatrix dphi(dim,nshape);
  TPZFMatrix dphix(dim,nshape);
  TPZFMatrix axes(3,3,0.);

  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL weight = 0.;
*/
  /**Trabalhando com o centro do elemento para formar Beta para difussividade */
/*  int dfseq, dfvar, iv;
  REAL coef;

  for(int int_ind = 0; int_ind < GetIntegrationRule().NPoints(); ++int_ind){

    GetIntegrationRule().Point(int_ind,intpoint,weight);
    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    fReference->X(intpoint, x);
    weight *= fabs(detjac);
    Shape(intpoint,phi,dphi);

    iv=0;
    int ieq, l;
    switch(dim) {
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
      PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
      PZError.flush();
    }
    for(in=0;in<numdof;in++) sol[in] = 0.;
    dsol.Zero();
    for(in=0; in<ncon; in++) {
      df = &Connect(in);
      dfseq = df->SequenceNumber();
      dfvar = block.Size(dfseq);
      for(jn=0; jn<dfvar; jn++) {
        coef = block(dfseq,0,jn,0);
        sol[iv%numdof] += phi(iv/numdof,0)*coef;
        for(l=0;l<dim;l++)                                 // Jorge 7/9/99
          dsol(l,iv%numdof) += dphix(l,iv/numdof)*coef;
        iv++;
      }
    }
    fMaterial->Contribute(x,jacinv,sol,dsol,weight,faxes,phi,dphix,*ek.fMat,*ef.fMat);
    if(SideOrder(2)) {
      if(IsZero(fCFLDiffussion)) continue; 
      fMaterial->IncrementDiffusion(x,jacinv,sol,dsol,weight,axes,phi,dphix,*ek.fMat,*ef.fMat,fCFLDiffussion);
    }
  }
}
*/
/**Compute the contribution to right hand vector on the element
   This method just used when the scheme is explicit*/
/*void TCompEl1dGD::CalcRhs(TPZElementMatrix &ef) {

  if(fMaterial == NULL){
    PZError << "TPZCompEl1d::CalcRhs : no material for this element\n";
    return;
  }

  if(fMaterial->Id()<0) {
    delete ef.fMat;
    ef.fMat = 0;
    return;
  }
  int i, numdof = fMaterial->NStateVariables();
  int ncon = NConnects();
  int dim = Dimension();
  int nshape = NShapeF();
  TPZBlock &block = Mesh()->Block();
  // clean ek and ef
  if(!ef.fMat) ef.fMat = new TPZFMatrix();
  if(!ef.fBlock) ef.fBlock = new TPZBlock(ef.fMat);

  int numeq = nshape*numdof;
  ef.fMat->Redim(numeq,1);
  ef.fBlock->SetNBlocks(ncon);
  TPZVec<REAL> sol(numdof);
  TPZFMatrix dsol(dim,numdof);

  for(i=0;i<ncon;i++)	{
    ef.fBlock->Set(i,NConnectShapeF(i)*numdof);
  }

  int in,jn;
  TPZConnect *df;

  ef.fConnect.Resize(ncon);

  for(i=0; i<ncon; ++i){
    (ef.fConnect)[i] = ConnectIndex(i);
  }

  TPZFMatrix phi(nshape,1);
  TPZFMatrix dphi(dim,nshape);
  TPZFMatrix dphix(dim,nshape);
  TPZFMatrix axes(3,3,0.);

  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  REAL detjac;
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL weight = 0.;
*/
  /** Encontrando a solucao no centro do elemento para difussao */
/*  int dfseq, dfvar, iv = 0;
  REAL coef;

  for(int int_ind = 0; int_ind < GetIntegrationRule().NPoints(); ++int_ind){
    GetIntegrationRule().Point(int_ind,intpoint,weight);
    fReference->Jacobian(intpoint,jacobian,axes,detjac,jacinv);
    fReference->X(intpoint,x);

    weight *= fabs(detjac);
    Shape(intpoint,phi,dphi);

    iv=0;
    int ieq, l;
    switch(dim) {
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
      PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
      PZError.flush();
    }
    for(in=0;in<numdof;in++) sol[in] = 0.;
    dsol.Zero();
    for(in=0; in<ncon; in++) {
      df = &Connect(in);
      dfseq = df->SequenceNumber();
      dfvar = block.Size(dfseq);
      for(jn=0; jn<dfvar; jn++) {
        coef = block(dfseq,0,jn,0);
        sol[iv%numdof] += phi(iv/numdof,0)*coef;
        for(l=0;l<dim;l++)                                 // Jorge 7/9/99
          dsol(l,iv%numdof) += dphix(l,iv/numdof)*coef;
        iv++;
      }
    }
    fMaterial->ContributeRhs(x,sol,weight,axes,phi,dphix,*ef.fMat);
    if(SideOrder(2)) {
      if(IsZero(fCFLDiffussion)) continue;
      TPZFMatrix temp(numeq,numeq);
      fMaterial->IncrementDiffusion(x,jacinv,sol,dsol,weight,axes,phi,dphix,temp,*ef.fMat,fCFLDiffussion);
    }
  }
}

*/

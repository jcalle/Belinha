/*******       File : gmesh.c

This file contains the function definitions for the class
TGeoMesh.
This file is based on pzgrid.c

Programmer : Jorge Lizardo Diaz Calle

*******/

#include "gmesh.h"
#include "pzelg1d.h"
#include "pzelgt2d.h"
#include "pzelgq2d.h"
//#include "elgp2d.h"
#include "pzerror.h"
#include "pzmatrix.h"
#include "pzcompel.h"
#include "myheader.h"


// CONSTRUTORS
TGeoMesh::TGeoMesh() : TPZGeoMesh() {
  if(Name()[0]=='\0') SetName("Geometrical mesh to work with no connect dependents");
}

/*To create new nodes, when we have a vector of your coordinates and bc type
  //Return bc nodes in bcnodes, and nodal bc id in idbc, to computational mesh
  void TGeoMesh::AddNodes(TPZVec<REAL> &Coord,TPZVec<int> &refbc,TPZVec<TPZGeoNode> &bcnodes,TPZVec<int> &idbc) {
  int i;
  int nnodes=refbc.NElements();
  if(!Coord.NElements()%nnodes) {
  PZError << "\nBad dimension of nodes.\n";
  return;
  }
  i=NoZeroNumbers(refbc);
  bcnodes.Resize(i);
  idbc.Resize(i);
  int nbc=0,id=CreateUniqueNodeId();

  for(i=0;i<nnodes;i++) {
  TPZGeoNode *newnode = new TPZGeoNode(id,Coord,*this);
  NodeVec()[id++] = *newnode;
  if(refbc[i]>0) { //Fill fBCNodeMap to the nodes with boundary condition associated
  bcnodes[nbc] = *newnode;
  idbc[nbc++]=refbc[i];
  }
  }
  } */

void TGeoMesh::FindTwoBoundaryNodes(TPZGeoEl *gel,int &side,int &node,int &lastnode) {
  int i,nsides,nelems=NElements();
  for(i=0;i<nelems;i++) {
    gel = ElementVec()[i];
    nsides = gel->NSides()-1;   //Last side is the self element
    //Identify any boundary edge.
    for(side=0;side<nsides;side++) {
      if(!gel->Neighbour(side).Element()) {
	node = gel->SideNodeIndex(side,0);
	lastnode = gel->SideNodeIndex(side,1);
	return;
      }
    }
  }
  node = -1;
}

int TGeoMesh::NSides1d() {
  int i,j,nsides=0,nelems=NElements();
  TPZGeoEl *gel;
  for(i=0;i<nelems;i++) {
    gel = ElementVec()[i];
    for(j=0;j<gel->NSides();j++)
      if(gel->SideDimension(j)==1) nsides++;
  }
  return nsides;
}
int TGeoMesh::MaxNSidesByEl() {
  int i,max=0,nsides,nelems=NElements();
  TPZGeoEl *gel;
  for(i=0;i<nelems;i++) {
    gel = ElementVec()[i];
    nsides = gel->NSides();
    max = (nsides<max) ? max : nsides;
  }
  return max;
}

void TGeoMesh::MiddleSideNode2d(TPZGeoEl *gel,int side,int &index) {
  //Node=Transformation1d(zero). We consider edge as one-dimensional cell.
  index = gel->SideNodeIndex(side,2);
  if(index < 0) {
    TPZVec<REAL> coord(3,0.);
    for(int i=0;i<3;i++)
      coord[i]=.5*(gel->SideNodePtr(side,0)->Coord(i)+gel->SideNodePtr(side,1)->Coord(i));
    index = NodeVec().AllocateNewElement();
    NodeVec()[index].Initialize(coord,*this);
  }
}

void TGeoMesh::IndexElement(TPZGeoEl *gel,int nelems,int &index) {
  int i;
  for(i=0;i<nelems;i++)
    if(gel==ElementVec()[i]) {
      index = i;
      return;
    }
  PZError << "TGeoMesh::IndexElement. This geometrical element is not belongs to mesh.\n";
  index = -1;
}

/**Create a new geometrical element with node ids, and material id. Return index gel*/
void TGeoMesh::CreateGeoEl2d(TPZVec<int64_t> &elnodes,int mat,int &indexgel) {
  //Construct the poligonal geometrical elements
  if(elnodes.NElements()==4) new TPZGeoElQ2d(elnodes,mat,*this);
  else if(elnodes.NElements()==3) new TPZGeoElT2d(elnodes,mat,*this);
//  else new TGeoElP2d(elnodes,mat,*this);
  indexgel = this->NElements()-1;
}

//Create new central nodes for each cell. Return centernodeid in idcenters
//Don't create bc nodes.
void TGeoMesh::AddCenterNodes(TPZVec<int> &idcenters,int eldim,int /*CenterType*/) {
  int ncells=NElements();
  if(idcenters.NElements()!=ncells) idcenters.Resize(ncells);

//  TPZGeoNode *center;
  TPZVec<REAL> coord(3);
//  int n,k=0;
  if(eldim==1) {
//    TPZGeoEl1d *cell;
    /*      Pix i=fCellMap.first();
	    while(i) {
	    cell=(TGeoCell1d *) fCellMap.contents(i);
	    fCellMap.next(i);
	    cell->Center(coord);
	    center=new TNode(id,coord);
	    fNodeMap[id]=center;
	    idcenters[k++]=id++;
	    }
	    }
	    else if(celldim==2) {
	    TGeoCell *cell;
	    Pix i=fCellMap.first();
	    while(i) {
	    cell=(TGeoCell *)fCellMap.contents(i);
	    fCellMap.next(i);
	    n=cell->NEdges();
	    if(n==3) (TGeoCellT *) cell;
	    else if(n==4) (TGeoCellQ *) cell;
	    cell->Center(coord,CenterType);
	    center=new TNode(id,coord);
	    fNodeMap[id]=center;
	    idcenters[k++]=id++;
	    }*/
  }
  //   else PZError << "\n TGeoMesh::CenterNodes not implemented for dimension >= 3\n";
}
//Create new nodes. It does the middle points of all triangle edges a new geometrical node
//We suppose the triangular cells are initials into the fCellMap. Only to two-dimensional cells
//Return bc nodes in bcnodes, and nodal bc id in idbc, to computational mesh
//and resize this vectors at end. The initial size of these must to be sufficient
void TGeoMesh::AddMiddle1dSideNodes(TPZVec<int> &indexmidnodes) {

  int indexgel,indexneigh,newnode,nsides,i,j,nelems=NElements();
  TPZGeoEl *gel;
  TPZGeoElSide neigh;
  int maxnsides = MaxNSidesByEl();
  indexmidnodes.Resize(maxnsides*nelems,-1);

  for(indexgel=0;indexgel<nelems;indexgel++) {
    gel = ElementVec()[indexgel];
    nsides = gel->NSides();
    j = indexgel*maxnsides;
    for(i=0;i<nsides;i++) {
      if(gel->SideDimension(i)!=1 && indexmidnodes[j+i]!=-1) continue;
      neigh=gel->Neighbour(i);
      MiddleSideNode2d(gel,i,newnode);
      indexmidnodes[j+i] = newnode;
      if(neigh.Element()) {
	IndexElement(neigh.Element(),nelems,indexneigh);
	indexmidnodes[indexneigh*maxnsides+neigh.Side()] = newnode;
      }
    }
  }
}

/*Return bc nodes in bcnodes, and nodal bc id in idbc, to computational mesh
void TGeoMesh::AddMiddleEdgeNodes(TPZVec<TPZGeoNode> &midnodes,TPZVec<int> &refbc,TPZVec<TPZGeoNode> &bcnodes,
   		TPZVec<int> &idbc,int n) {
	if(fDim==1) {       //n indicate edges number of each cell, must to be uniform
   	cout << "\nThis function was implemented only to two-dimensional cells.\n";
   }
   //Constructing auxiliar vector to know whether the middle edge node exist already
   Long id=NCells();
   TIntVec existnode(n*id,0);

   id=CreateNodeId();
   TNode *newnode;

   if(bcnodes.capacity()!=midnodes.capacity())
   	cout << "ERROR TGeoMesh::AddMiddleEdgeNode. Bad length in vectors.\n";
   int num=0; // number of bc nodes finding

   TGeoCell *cell;
   int idcell,edgecell;
   TGeoCell *neigh;
   int idneigh,edgeneigh;

   Pix i=fCellMap.first();
   int j=0;
	while(i) {
   	cell=(TGeoCell *)fCellMap.contents(i);
      fCellMap.next(i);

      idcell=cell->Id();
      for(edgecell=0;edgecell<n;edgecell++) {
      	if(existnode[n*idcell+edgecell]!=0) continue;
         existnode[n*idcell+edgecell]=1;
         neigh=(TGeoCell *)cell->Neighbour(edgecell);
      	newnode=cell->NewMiddleEdgeNode(id,edgecell);
         fNodeMap[id++]=newnode;
         if(!neigh) {
         	bcnodes[num]=newnode;
            idbc[num++]=RefNodes[(cell->EdgeNodeId(edgecell,0))-1];
         }
         else {
         	idneigh=neigh->Id();
         	edgeneigh=cell->NeighbourEdge(edgecell);
            existnode[n*idneigh+edgeneigh]=1;
         }
         midnodes[j]=newnode;
			j++;
      }
   }
   bcnodes.resize(num);
   idbc.Resize(num);
	if(cell->Dim()!=2)
   	cout << "TGeoMesh::AddMiddleEdgeNode only defined to two-dimensional cell.\n";
   if(j!=midnodes.capacity())
   	cout<<"Warning : MiddleEdgesNodes. Number of middles different of the fNEdges.\n";
}

//Create new elements. nnodes is array of the number of nodes of each new cell
//and nedges is array of the number of edges of each new cell
void TGeoMesh::AddGCells(TPZVec<int> &nnodes,TPZVec<int> &nedges,TPZVec<int> &nodes,int celldim) {
	if(celldim<1 || celldim>2) {
      cout << "\nTGeoMesh::AddGCells not implemented for this cell dimension.\n";
      return;
   }

	int i,j,m=0;
   int ncells=nnodes.Capacity();
   int n=nnodes.Sum();
	if(n!=nodes.capacity()) {
   	cerr << "\nTGeoMesh::AddGCells bad number of nodes.\n";
      return;
   }
   if(nedges.Capacity()!=ncells) {
		cerr << "\nTGeoMesh::AddGCells bad number of cells.\n";
      return;
   }

   Long id=CreateCellId();
   for(i=0;i<ncells;i++) {
   	n=nnodes[i];
     	LongVec Node(n);
		for(j=0;j<n;j++) Node[j] = nodes[m+j];
      m+=n;
		TGeoCell *newcell;
      if(celldim==1) (TGeoCell1d *)newcell=new TGeoCell1d(id,Node);
      else {
         if(nedges[i]<3) cout<<"ERROR : This geometrical element isn't two-dimensional.\n";
         else if(nedges[i]==3) (TGeoCellT *)newcell = new TGeoCellT(id,Node);
         else if(nedges[i]==4) (TGeoCellQ *)newcell = new TGeoCellQ(id,Node);
         else newcell = new TGeoCell(id,Node,nedges[i]);
      }
		fCellMap[id++] = newcell;
	}
}
void TGeoMesh::AddGCells(int nnodes,int nedges,TPZVec<int> &nodes,int celldim) {
	if(celldim<1 || celldim>2) {
      cout << "\nTGeoMesh::CenterNodes not implemented for this cell dimension.\n";
      return;
   }

	int i,j;
   int ncells=nodes.NElements()/nnodes;

	if(nedges*ncells!= nodes.NElements())
   	cerr << "\nTGeoMesh::AddGCells Bad number of nodes.";

   Long id=CreateCellId();
	LongVec Node(nnodes);
	for(i=0;i<ncells;i++) {
		for(j=0;j<nnodes;j++) Node[j] = nodes[i*nnodes+j];
		TGeoCell *newcell;
      if(celldim==1) (TGeoCell1d *)newcell=new TGeoCell1d(id,Node);
      else {
         if(nedges<3) cout<<"ERROR : This geometrical element isn't two-dimensional.\n";
         else if(nedges==3) (TGeoCellT *)newcell = new TGeoCellT(id,Node);
         else if(nedges==4) (TGeoCellQ *)newcell = new TGeoCellQ(id,Node);
         else newcell = new TGeoCell(id,Node,nedges);
      }
		fCellMap[id++] = newcell;
	}
} */


/***
    void TGeoGrid::BuildConnectivity() {

    TLongVec NeighNode(NumNodes(),-1),SideNum(NumNodes(),-1);
    TIntAVLMap Ids(-1);
    Pix inod = fNodeMap.first();
    int iseq = 0;
    while(inod) {
    TGeoNod *gn = (TGeoNod *) fNodeMap.contents(inod);
    Ids[gn->Id()] = iseq++;
    fNodeMap.next(inod);
    }
    Pix iel = fElementMap.first();

    long numsearch =1;
    // if there are no elements, do nothing
    iel = fElementMap.first();
    while(iel) {
    TGeoEl *el = (TGeoEl *) fElementMap.contents(iel);

    int numsides = el->NumSides();
    int side;
    for(side = 0;side<numsides;side++) {

    // check whether all entries in NeighNode are equal

    int equalnode = 1;
    int numsidenodes = el->NumSideNodes(side);
    long sidenode = el->SideNodeId(side,0);
    int locnod = Ids[sidenode];
    long neigh = NeighNode[locnod];
    int sidenumber = SideNum[locnod];
    for(int sn = 0;sn < numsidenodes; sn++) {
    sidenode = el->SideNodeId(side,sn);
    locnod = Ids[sidenode];
    if (neigh != NeighNode[locnod]){
    equalnode=0;
    break;
    }
    }

    if(equalnode && neigh == -1) {
    if(!el->Neighbour(side) && !el->Bc(side) && el->NeighbourSide(side)) {
    int elloaded = 0;
    for(int in=0; in<el->NumberOfNodes(); in++) {
    if(NeighNode[Ids[el->NodePtr(in)->Id()]] == el->Id()) elloaded = 1;
    }
    // this element is not loaded and its side is undefined

    // load the element side in the NeighNode vector
    for(int sn=0;!elloaded && sn < numsidenodes; sn++) {
    sidenode = el->SideNodeId(side,sn);
    int locnod = Ids[sidenode];
    NeighNode[locnod] = el->Id();
    SideNum[locnod] = side;
    }
    numsearch++;
    }
    } else if(equalnode && side == sidenumber && neigh == el->Id()) {
    // unload the element side
    for(int sn=0;sn < numsidenodes; sn++) {
    sidenode = el->SideNodeId(side,sn);
    int locnod = Ids[sidenode];
    NeighNode[locnod] = -1;
    SideNum[locnod] = -1;
    }
    // if no neighbouring element was detected during the loop
    //    define the element side as undefined
    if(!el->Neighbour(side)) el->SetSide(side,0);
    numsearch++;
    } else if(equalnode && neigh != el->Id()) {
    // we found a neigbour
    TGeoEl *neighbour = FindElem(neigh);
    TLongVec SideNodes(numsidenodes);
    // detect which side of the neigbour is loaded witin NeighNode
    if(!neighbour) {
    cout << "TGeoGrid::BuildConnectivity element neigh = " << neigh <<
    " not found\n";
    exit(-1);
    }
    for(int sn=0;sn < numsidenodes; sn++) {
    sidenode = el->SideNodeId(side,sn);
    SideNodes[sn] = sidenode;
    }
    // WhichSide will tell the side number which contains the vector
    //    of node ids SideNodes
    int neighside = neighbour->WhichSide(SideNodes);
    if(neighside != -1 && !el->NeighbourExists(neighbour,side)){
    el->SetConnectivity(side,neighbour,neighside);
    }
    }
} // loop over the sides
fElementMap.next(iel);
if(!iel && numsearch) {
numsearch = 0;
iel = fElementMap.first();
}
}
}
TGeoMesh::TGeoMesh(TGeoMesh &gr) : fCellMap(gr.fCellMap),
  fNodeMap(gr.fNodeMap), fCosysMap(gr.fCosysMap) {

  fChecked = gr.fChecked;
  SetName(gr.fName);
  fDim  = gr.fDim;
  gCurrent = this;
  fReference = 0;
  }

  // DELETE ITEMS (CoSys, cells and nodes) in lists
  void TGeoMesh::CleanUp() {
  Pix i = fCellMap.first();
  while(i) {
  delete ((TGeoCell *) fCellMap.contents(i));
  fCellMap.next(i);
  }
  fCellMap.clear();
  i = fNodeMap.first();
  while(i) {
  delete ((TNode *) fNodeMap.contents(i));
  fNodeMap.next(i);
  }
  fNodeMap.clear();
  i = fCosysMap.first();
  while(i) {
  delete ((TCosys *) fCosysMap.contents(i));
  fCosysMap.next(i);
  }
  fCosysMap.clear();
  }

  // Create an id for a new geometrical element of a list of pointers (map).
  Long TGeoMesh::CreateNodeId() {
  Pix i = fNodeMap.last();
  if(!i) return 0;
  TNode *node = (TNode *) fNodeMap.contents(i);
  return (node->Id())+1;
  }
  Long TGeoMesh::CreateCellId() {
  Pix i = fCellMap.last();
  if(!i) return 0;
  TGeoCell *gcell = (TGeoCell *) fCellMap.contents(i);
  return (gcell->Id())+1;
  }
  // Acoording to order in the CellMap return the center node pointers in centers
  // Can't exist bc nodes
  void TGeoMesh::AddCenterNodes(TPZVec<TPZGeoNode> &centers,int celldim,int CenterType) {
  int i,index;

  TPZGeoNode center;
  TPZVec<REAL> coord(3,0.);

  int nelem = NElements();
  if(nelem != centers.NElements()) {
  PZError << "TGeoMesh::AddCenterNodes. Bad dimension of the centers vector.\n";
  centers.Resize(nelem);
  }
  TPZGeoEl *gel;

  for(i=0;i<nelem;i++) {




  gel = ElementVec()
  if(celldim==1) {
  TGeoCell1d *cell;
  Pix i=fCellMap.first();
  while(i) {
  cell=(TGeoCell1d *) fCellMap.contents(i);
  fCellMap.next(i);
  cell->Center(coord);
  center=new TNode(id,coord);
  fNodeMap[id++]=center;
  centers[k++]=center;
  }
  }
  else if(celldim==2) {  //here must to be celldim = 2
  TGeoCell *cell;
  Pix i=fCellMap.first();
  while(i) {
  cell=(TGeoCell *)fCellMap.contents(i);
  fCellMap.next(i);
  n=cell->NEdges();
  if(n==3) (TGeoCellT *) cell;
  else if(n==4) (TGeoCellQ *) cell;
  cell->Center(coord,CenterType);
  center=new TNode(id,coord);
  fNodeMap[id++]=center;
  centers[k++]=center;
  }
  }
  else PZError << "\n TGeoMesh::CenterNodes not implemented for dimension >= 3\n";
  }
*/

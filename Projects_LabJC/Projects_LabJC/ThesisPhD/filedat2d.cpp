/*******       File : modwrm.c

This file contains the method definitions for class TDatafile2d.
One object of this class construct a geometrical and a computational
mesh from data file generated from Modulef, Aranha and/or Triangle
packets.

*******       Programmer : Jorge Lizardo Diaz Calle       *******/
#include "conslaw.h"

#include "filedat2d.h"
#include "gmesh.h"
#include "pzcmesh.h"
#include "pzerror.h"
#include "elct2dgd.h"
#include "elcpointgd.h"
#include "elc1dgd.h"
#include "pzelgpoint.h"
#include "pzadmchunk.h"
#include "pzelgt2d.h"
#include "pzelgq2d.h"

#include "pzgeoelbc.h"

#include "pzelg1d.h"

#include "pzelct2d.h"
#include "pzelcq2d.h"

#include "myheader.h"
#include "hadaptive.h"

#include <stdlib.h>
#include <string.h>

TDatafile2d::TDatafile2d(std::istream &input,char *name) : RefNodes(0), fDataMesh() {
  fFileError = 1;
	fData = &(input);
  fCellType=1;
  fCenterType=1;
  fFilename[0]='\0';
  GetDataCommented(input,fFilename,256);
  if(!strncmp(fFilename,"MESH DATA",9)) {
    fFilename[0]='\0';
    strncpy(fFilename,name,64);
	  fFileType = 1;
    fFileError = 0;
    return;
  }
  fDataMesh.open(fFilename, std::ios::in);
  while(!fDataMesh.rdbuf()->is_open()) {
    PZError << "\nERROR(TPZDatafile::const)-> Cannot open " << fFilename << "\n";
    PZError << "\nBad file name. Data file name (\"new\" create new file; \"quit\" exit) : ";
    fFilename[0]='\0';
	std::cin >> fFilename;
    if(!strncmp(fFilename,"quit",4)) exit(1);
    fDataMesh.close();
    fDataMesh.open(fFilename, std::ios::in);
  }
	fFileType = 0;
  fFileError=0;
	fDeleted = 0;
}
TDatafile2d::TDatafile2d(char *file) : RefNodes(0), fDataMesh() {
  fFileError = 1;
  strcpy(fFilename,file);
  if(!strncmp(fFilename,"quit",4)) exit(1);
  fDataMesh.close();
  fDataMesh.open(fFilename, std::ios::in);
  while(!fDataMesh.rdbuf()->is_open()) {
	  std::cerr << "\nERROR(TPZDatafile::const)-> Cannot open " << fFilename << "\n";
	  std::cerr << "\nBad file name. Data file name (\"new\" create new file; \"quit\" exit) : ";
    fFilename[0]='\0';
	std::cin >> fFilename;
    if(!strncmp(fFilename,"quit",4)) exit(1);
    fDataMesh.close();
    fDataMesh.open(fFilename, std::ios::in);
  }
  fFileType = 0;
  fFileError=0;
}

void TDatafile2d::RequireCellType() {
  int type;
  std::cout << "\nFour types of mesh can to be constructed. With\n\t1. Triangular Cell";
  std::cout << "\n\t2. Cell with boundary between two initial triangle centers";
  std::cout << "\n\t3. Cell with boundary join middle point side and center of initial triangle";
  std::cout << "\n\t4. Cell edge type. In each triangle is joined the center with the three vertexs";
  std::cout << "Choice Cell Type : ";
  std::cin >> type;
  if(type<1 || type>4) {
	  std::cout << "\nBad choice. Default Type = 1";
    type=1;
  }
  fCellType=type;
}
void TDatafile2d::RequireCenterType() {
  int type;
  std::cout << "\nYou can to choice between four center points of the triangle.\n";
  std::cout << "\t1. Centroid \t2. Circuncenter";
  std::cout << "\n\t3. Incenter \t4. Orthocenter";
  std::cout << "Choice center point : ";
  std::cin >> type;
  if(type<1 || type>4) {
    PZError << "TDatafile2d::RequireCenterType. Bad choice. Default Center = 1";
    type=1;
  }
  fCenterType=type;
}

void TDatafile2d::SetCellType(int type) {
  if(type<1 || type>4) {
	  std::cout << "TDatafile2d::SetCellType. Bad choice. Default Type = 1";
    fCellType=1;
  }
  else fCellType=type;
}

/**To initialize triangular geometrical mesh from modulef file data*/
void TDatafile2d::InitializeMesh(TPZCompMesh *c, TConservationLaw *mat) {
  TPZStack<int> type;
  TPZGeoMesh *g = c->Reference();

  if(!fFileType) Read(g,type);     //Read data of the file and put into vectors.
  else ComputeData(g,type);       //Compute data mesh from short information.

  /**Determining the computational cell type and center type*/
  GetDataCommented(*fData,2,fCenterType);      //center types
  if(fCenterType>4 || fCenterType<1) fCenterType=1;
  GetDataCommented(*fData,2,fCellType);     //Dual cell type
  if(fCellType>4 || fCellType<1) fCellType=1;

  /**Determining the geometrical connectivity between the geometrical elements*/
  g->BuildConnectivity();
  /**Initialize material and boundary conditions into the computational mesh*/
  ReadMaterialBC(c,mat,type);   //Use type

  /**Constructs a dual geometrical grid depending of the fCellType whether
     fCellType==1 the geometrical mesh is the same original mesh g=c->Reference*/
  if(fCellType==1) c->AutoBuild();     // ((TCompMesh *)c)->AutoBuildDiscontinuous();
  else {
    /**To initialize a geo. grid with type center-center cell*/
    if(fCellType==2) DualMeshCC(c);
    //To initialize a comp. grid with type center-middle_point cell
    else if(fCellType==3) DualMeshCM(c);
    //To initialize a comp. grid with type edge cell
    else if(fCellType==4) DualMeshEdge(c);
    else if(fCellType!=1)
      PZError << "TDatafile2d::InitializeCMesh. Bad value of the fCellType.\n";
  }
}

void TDatafile2d::ReadMaterialBC(TPZCompMesh *c, TConservationLaw *mate,TPZStack<int> &type) {
  /**We must to have applied BuildConnectivity of the geometrical mesh*/
  TPZGeoMesh *g = c->Reference();
  int nnodes=g->NNodes(),nelems=g->NElements();
  int i,j,k,n,nn,id,nvar,refnode=0,nbcs;
  TPZVec<int> MaskNodes(nnodes,-1);
  std::ofstream output("filebc.dat");
  output << "Boundary element - boundary side - start node " << std::endl;
  for(i=0;i<type.NElements();i++) {
	if(i>0) mate = (TConservationLaw *)mate->NewMaterial();  //Is necessary to enable us to create same type material
    c->InsertMaterialObject(mate);
    mate->SetId(type[i]);
	mate->SetData(*fData);
  }
  for(i=0;i<nelems;i++) {
    TPZGeoEl *gel = g->ElementVec()[i];
		if(!gel) continue;
    n = gel->NSides()-1;
    nn = gel->NNodes();
    for(j=nn;j<n;j++) {
      if(!gel->Neighbour(j).Element()) {
        double x[4];
        double xmedio[2];
        TPZGeoNode *node1 = gel->NodePtr(j-3);
        TPZGeoNode *node2 = gel->NodePtr((j-2)%3);
        x[0] = node1->Coord(0);
        x[1] = node1->Coord(1);
        x[2] = node2->Coord(0);
        x[3] = node2->Coord(1);
        output << gel->Id() << "\t" << j << "\t" << x[0] << " " << x[1];
        output << "\t" << x[2] << " " << x[3] << std::endl;
        xmedio[0] = 0.5*(x[0] + x[2]);
        xmedio[1] = 0.5*(x[1] + x[3]);
        k = mate->IdBC(xmedio);
        if(k<0 && k!= -4) TPZGeoElBC(gel,j,k);
      }
    }
  }
  GetDataCommented(*fData,nbcs);
  nvar = mate->NStateVariables();
  GetCommentary(*fData,2);
  TPZFMatrix<STATE> value1(nvar,nvar),value2(nvar,1);
  int mask = 0;
  for(i=0;i<nbcs;i++) {
    int typebc;
    (*fData) >> id;       // id of boundary condition
    (*fData) >> typebc;   // type of boundary condition
    if(typebc==5 && !mask) {
      mask = 1;
	  std::cout << "Create the periodic data file : bcperiod.dat.\n";
	  std::cout << "Type any key whether the bcperiod.dat is created.\n";
      Pause(typebc);
	  std::ifstream input("bcperiod.dat");
      BCPeriodic(g,input);
      input.close();
    }
    for(j=0;j<nvar;j++) for(k=0;k<nvar;k++) (*fData) >> value1(j,k);
    for(j=0;j<nvar;j++) (*fData) >> value2(j,0);
    c->InsertMaterialObject((TPZMaterial *)(mate->CreateBC(id,typebc,value1,value2)));
  }
  output.close();
}

void TDatafile2d::BCPeriodic(TPZGeoMesh *gmesh, std::istream &input) {
  int i,numberpar,indexgel,side;
  TPZGeoEl *gel, *neigh;
  TPZGeoElSide gelsidenull(0,-1);
  /** Apagando as connectividades nos geoelside fronteira */
  input >> numberpar;
  for(i=0;i<numberpar;i++) {
    input >>  indexgel >> side;
    gel = gmesh->ElementVec()[indexgel];
		if(!gel) continue;
    TPZGeoElSide gelside(gel,side);
    gelside.SetNeighbour(gelsidenull);
  }
  /** Atribuindo a nova connectividade nos elementos fronteira com condicoes periodicas */
  input >> numberpar;
  for(i=0;i<numberpar;i++) {
    input >>  indexgel >> side;
    gel = gmesh->ElementVec()[indexgel];
		if(!gel) continue;
    TPZGeoElSide gelside(gel,side);
    input >> indexgel >> side;
    neigh = gmesh->ElementVec()[indexgel];
    TPZGeoElSide neighside(neigh,side);
    gelside.SetConnectivity(neighside);
  }
}

void TDatafile2d::ComputeData(TPZGeoMesh *g,TPZStack<int> &mattypes) {
	/** Short mesh information 1, large mesh information 2 */
	GetDataCommented(*fData,fFileType);
	if(!fFileType) {
	  PZError << "TDatafile2d::ComputeData.ComputeData is called with type 0.\n";
		exit(1);
  }
	if(fFileType!=1) {
    PZError << "TDatafile2d::ComputeData. You need to implement filetype no 1.\n";
		return;
	}

  TPZVec<REAL> Coord1(3,0.);
  TPZVec<REAL> Coord2(3,0.);
  TPZVec<REAL> dx(3,0.);
  int i,dim=3,eltype,nelemx,nelemy;

  /**The coordinates of the first and last nodes*/
  GetDataCommented(*fData,Coord1);   // Coordinates of the first node
  GetDataCommented(*fData,Coord2);   // Coordinates of the last node

  /**Number of subdivisions-X and subdivisions-Y*/
  GetDataCommented(*fData,nelemx);   // Subdivision-x, number of intervales
  GetDataCommented(*fData,nelemy);   // Subdivision-y, number of intervales

	/** Read parameter to deleted. If fDeleted then destroy all geometrical elements with xcenter > fxZero and ycenter < fyZero */
	GetDataCommented(*fData,fDeleted);
	GetDataCommented(*fData,fxZero);
	GetDataCommented(*fData,fyZero);

  /**Computing dx, dy and dz*/
  dx[0] = (Coord2[0]-Coord1[0])/nelemx;
  dx[1] = (Coord2[1]-Coord1[1])/nelemy;
  dx[2] = (Coord2[2]-Coord1[2])/nelemy;

  /**To construct the geometrical element */
  int nmat,idmat;
  GetDataCommented(*fData,eltype);   // type of element : triangular(3) or quadrilateral(4)
  if(eltype < 3 || eltype > 4) {
    PZError << "TDatafile2d::ComputeData. Bad parameter eltype.\n"; return;
  }

  GetDataCommented(*fData,nmat);   //Number of materials
  GetDataCommented(*fData,idmat);   //Ids of the materials
  mattypes.Push(idmat);
  for(i=1;i<nmat;i++) {
    (*fData) >> idmat;
    mattypes.Push(idmat);
  }

  WriteCoord(g,dim,eltype,nelemx,nelemy,Coord1,dx);
  WriteElements(g,eltype,nelemx,nelemy,mattypes);
	
	/** Deleting geometrical elements if fDeleted */
	nelemx = g->NElements();
	TPZVec<REAL> center(3,0.);
	TPZVec<REAL> centerx(3,0.);
	TPZGeoEl *gel;
	for(i=0;i<nelemx;i++) {
		gel = g->ElementVec()[i];
		if(!gel) continue;
		gel->CenterPoint(gel->NSides()-1,center);
		gel->X(center,centerx);
		if(centerx[0] > fxZero && centerx[1] < fyZero) {
			delete gel;
			g->ElementVec()[i] = 0;
		}
	}
}

/**Compute coordinates of the geometrical nodes*/
void TDatafile2d::WriteCoord(TPZGeoMesh *g,int dim,int /*eltype*/,int nelemx,int nelemy,
			     TPZVec<REAL> &FirstCoord,TPZVec<REAL> &dx) {

  TPZVec<REAL> Coord(FirstCoord);
  nelemx++;
  nelemy++;
  int i,j,index;

  if(dim==2) {
    for(i=0;i<nelemy;i++) {
      Coord[0] = FirstCoord[0];
      for(j=0;j<nelemx;j++) {
        index = g->NodeVec().AllocateNewElement();
        g->NodeVec()[index].Initialize(Coord,*g);
				if(i<nelemy-1 && j<nelemx-1) {
	        Coord[0] += (.5*dx[0]);
					Coord[1] += (.5*dx[1]);
					index = g->NodeVec().AllocateNewElement();
					g->NodeVec()[index].Initialize(Coord,*g);
					Coord[0] -= (.5*dx[0]);
					Coord[1] -= (.5*dx[1]);
				}
				Coord[0] += dx[0];
      }
      Coord[1] += dx[1];
    }
  }
  else if(dim==3) {                                 //   !!!   ???  Revisar !
    for(i=0;i<nelemy;i++) {
      Coord[0] = FirstCoord[0];
      Coord[2] = FirstCoord[2];
      for(j=0;j<nelemx;j++) {
        index = g->NodeVec().AllocateNewElement();
        g->NodeVec()[index].Initialize(Coord,*g);
				if(i<nelemy-1 && j<nelemx-1) {
	        Coord[0] += (.5*dx[0]);
					Coord[1] += (.5*dx[1]);
					Coord[2] += (.5*dx[2]);
					index = g->NodeVec().AllocateNewElement();
					g->NodeVec()[index].Initialize(Coord,*g);
					Coord[0] -= (.5*dx[0]);
					Coord[1] -= (.5*dx[1]);
					Coord[2] -= (.5*dx[2]);
				}
        Coord[0] += dx[0];
        Coord[2] += dx[2];
      }
      Coord[1] += dx[1];
    }
  }
  else PZError << "TDatafile2d::WriteCoord. Bad parameter dimension.\n";

  RefNodes.Resize(g->NNodes()+1,0);
}
/**Compute node indexes for geometrical elements*/
void TDatafile2d::WriteElements(TPZGeoMesh *g,int eltype,int nelemx,int nelemy,TPZStack<int> &mattypes) {
  int i,j,k,index;
  TPZVec<int64_t> indexnodes(eltype);

  if(eltype==4) {
    for(i=0;i<nelemy;i++) {
      for(k=0;k<2;k++) indexnodes[k] = i*(nelemx+1)+k;
      indexnodes[2] = (i+1)*(nelemx+1)+1;
      indexnodes[3] = (i+1)*(nelemx+1);
      for(j=0;j<nelemx;j++) {
				new TPZGeoElQ2d(indexnodes,mattypes[0],*g);
				for(k=0;k<4;k++) indexnodes[k]++;
      }
    }
    return;
  }
  int lastindex, zeroindex;
  for(i=0;i<nelemy-1;i++) {
		zeroindex = i*(2*nelemx + 1);
		lastindex = (i+1)*(2*nelemx + 1);
		for(j=0;j<nelemx;j++) {
			indexnodes[0] = zeroindex + 2*j;
			indexnodes[1] = zeroindex + 2*(j+1);
			indexnodes[2] = zeroindex + 2*j+1;
      new TPZGeoElT2d(indexnodes,mattypes[0],*g);
      indexnodes[0] = indexnodes[2];
      indexnodes[2] = lastindex + 2*(j+1);
      new TPZGeoElT2d(indexnodes,mattypes[0],*g);
			indexnodes[1] = indexnodes[2];
      indexnodes[2] = lastindex + 2*j;
      new TPZGeoElT2d(indexnodes,mattypes[0],*g);
			indexnodes[1] = indexnodes[2];
      indexnodes[2] = zeroindex + 2*j;
      new TPZGeoElT2d(indexnodes,mattypes[0],*g);
    }
	}
	zeroindex = (nelemy-1)*(2*nelemx + 1);
	lastindex = nelemy*(2*nelemx + 1);
	for(j=0;j<nelemx;j++) {
		indexnodes[0] = zeroindex + 2*j;
		indexnodes[1] = zeroindex + 2*(j+1);
		indexnodes[2] = zeroindex + 2*j+1;
    new TPZGeoElT2d(indexnodes,mattypes[0],*g);
    indexnodes[0] = indexnodes[2];
    indexnodes[2] = lastindex + j+1;
    new TPZGeoElT2d(indexnodes,mattypes[0],*g);
		indexnodes[1] = indexnodes[2];
    indexnodes[2] = lastindex + j;
    new TPZGeoElT2d(indexnodes,mattypes[0],*g);
		indexnodes[1] = indexnodes[2];
    indexnodes[2] = zeroindex + 2*j;
    new TPZGeoElT2d(indexnodes,mattypes[0],*g);
  }
	/*
    for(k=0;k<2;k++) indexnodes[k] = i*(nelemx+1)+k;
    indexnodes[2] = (i+1)*(nelemx+1);

    for(j=0;j<nelemx-1;j++) {
      new TPZGeoElT2d(indexnodes,mattypes[0],*g);
      indexnodes[0] = indexnodes[2];
      indexnodes[2]++;
      new TPZGeoElT2d(indexnodes,mattypes[0],*g);
			indexnodes[0] = indexnodes[1];
      indexnodes[1]++;
    }
    if(i<nelemy/2) {
      new TPZGeoElT2d(indexnodes,mattypes[0],*g);
      indexnodes[0] = indexnodes[2];
      indexnodes[2]++;
      new TPZGeoElT2d(indexnodes,mattypes[0],*g);
    }
    else {
      indexnodes[2]++;
      new TPZGeoElT2d(indexnodes,mattypes[0],*g);
			indexnodes[1] = indexnodes[0];
      indexnodes[0] = indexnodes[2]-1;
      new TPZGeoElT2d(indexnodes,mattypes[0],*g);
    }
  }*/
  int nmat = mattypes.NElements();
  k = g->NElements();
  if(nmat!=1) {
	  std::cout << "\nGeometrical mesh : Number of elements = " << g->NElements();
	  std::cout << "Material Id = " << mattypes[0] << ". Last element = ";
	std::cin >> index;
    if(index > k) index = k;
    for(i=1;i<nmat;i++) {
		std::cout << "Material Id = " << mattypes[i] << ". Last element = ";
		std::cin >> lastindex;
      if(lastindex > k) lastindex = k;
      for(j=index;j<lastindex;j++)
	g->ElementVec()[j]->SetMaterialId(mattypes[i]);
    }
  }
}

//Initialize center-center type geometrical dual mesh from triangular mesh
//Remenber: The node ids are numered from one. The element ids are numered from zero
//because of the experiment with BuildConnectivity, the node id iniciate from 1
//but the element id can to iniciate from zero. Is necessary to analise this problem!!!
//Construct one geometrical poligonal element for each vertex of the triangularization
void TDatafile2d::DualMeshCC(TPZCompMesh *c) {
  TGeoMesh *g = (TGeoMesh *)(c->Reference());
  int i,j,mat,matneigh=0,firstnewel=-1;

  //We store a)the node pointers for each new poligonal elements into nodes, b)the number of
  //nodes of each new elements, c)the associated node for each new computational element
  //Warning : for simplicity we assume the material is unique, for now.
  //The sum of sizes of all nodes is 3*(fNElements+fNBCNodes);
  int side,nnodes=g->NNodes(),nelems=g->NElements();  //4 is the minimum number of edges to poligonal boundary elements
  TPZVec<int> mask(nnodes,0);  //mask[id] identify whether (id+1) node has poligonal element
  TPZVec<int> centers(nelems);   //To store the indexes of the center nodes
  ((TGeoMesh *)g)->AddCenterNodes(centers,2,fCenterType); //Into fNodeMap have fNNodes vertexs plus fNTriangle baricenters

  TPZAdmChunkVector<TPZGeoEl *> *elmap=&(g->ElementVec());
  elmap->CompactDataStructure();

  if(!nnodes) PZError << "\nERROR : TDatafile2d::DualMeshBB. Geometric mesh is empty";

  TPZGeoEl *gel = 0;
  TPZGeoElSide neigh;

  /**Identifying two consecutive boundary nodes to start the construction of the boundary elements*/
  int node,lastnode;
  ((TGeoMesh *)g)->FindTwoBoundaryNodes(gel,side,node,lastnode);
  if(node < 0) {
    PZError << "TDatafile2d::DualMeshCC. Don't exist any boundary side?\n";
    return;
  }

  //We are working with the others bc nodes (clock).
  int newnode,newnode1,firstnewnode,indexgel,k=1;                //indexes
  //To new bc nodes in boundary middle edge
  TPZStack<int> bcnodes;   //To store index boundary node and bc id
  TPZStack<int64_t> elnodes;

  mat = gel->MaterialId();
  g->MiddleSideNode2d(gel,side,newnode);
  if(newnode<0) {
    if(RefNodes.NElements()<newnode+1) RefNodes.Resize(newnode+1,0);
    RefNodes[newnode]=RefNodes[node];
  }
  g->IndexElement(gel,nelems,indexgel);
  firstnewnode = newnode;

  while(k) {
    elnodes[0]=node;
    elnodes[1]=newnode;
    g->IndexElement(gel,nelems,indexgel);
    elnodes[2]=centers[indexgel];

    bcnodes.Push(newnode);
    bcnodes.Push(RefNodes[node]);

    side=(side+2)%3;
    neigh = gel->Neighbour(side);
    //Filling the nodes
    while(neigh.Element()) {  //loop over elements connected with nodefirst
      matneigh = neigh.Element()->MaterialId();
      if(mat!=matneigh) {
        matneigh = -1;
        break;
      }
      side = neigh.Side();
      gel = neigh.Element();
      g->IndexElement(gel,nelems,indexgel);
      elnodes.Push(centers[indexgel]);
      side = (side+2)%3;
      neigh = gel->Neighbour(side);
    }
    if(lastnode==node && matneigh!=-1) {
      elnodes.Push(firstnewnode);
      k=0; //Don't exist more boundary cells. Out to loop
    }
    else {
      if(matneigh != -1) node = gel->SideNodeIndex(side,0);
      g->MiddleSideNode2d(gel,side,newnode);
      if(newnode<0) {
        if(RefNodes.NElements()<newnode+1) RefNodes.Resize(newnode+1,0);
        if(matneigh != -1) RefNodes[newnode]=RefNodes[node];
      }
      elnodes.Push(newnode);
    }

    g->CreateGeoEl2d(elnodes,mat,indexgel);
    if(firstnewel<0) firstnewel=indexgel;

    /**Fill the boundary condition number over the boundary sides
       newgel->SetEdge(0,-fRefNodes[(int)((TNode *)nodes[0])->Id()-1]);
       if(k!=0) newgel->SetEdge(n0-1,-fRefNodes[(int)node->Id()-1]);
       else newgel->SetEdge(n0-1,-fRefNodes[(int)gel->EdgeNodeId(edge,0)-1]);
    */
    if(matneigh != -1) mask[node]=1;         //mask the node with associated element
    else {
      gel = neigh.Element();
      side = neigh.Side();
    }
    mat = gel->MaterialId();
    elnodes.Resize(3);
    matneigh = 0;
  }

  //Construct the internal poligonal elements. We use the information in fNodesT
  for(i=0;i<nnodes;i++) {
    elnodes.Resize(0);
    if(mask[i]!=0) continue;
    node = -1;
    //else we create a new internal poligonal element associated with j+1 node
    //find element id with j+1 node as first vertex of a edge
    while(1) {
      j=0;
      gel = g->ElementVec()[j++];
      for(k=0;k<gel->NNodes();k++)
	if(i==gel->NodeIndex(k)) { node=i; break; }
      if(k<0) break;
    }
    int gelside;
    for(k=0;k<gel->NSides();k++) {
      if(i==gel->SideNodeIndex(k,0)) {
        gelside=k;
        break;
      }
    }

    g->IndexElement(gel,nelems,indexgel);
    elnodes.Push(centers[indexgel]);
    side=(gelside+2)%3;
    mat = gel->MaterialId();
    neigh=gel->Neighbour(side);
    side=neigh.Side();  //neigh exist
    TPZGeoEl *neighel = neigh.Element();
    newnode1=0;
    while(gel!=neighel) {
      /**If gel and neigh have different materials*/
      matneigh = neigh.Element()->MaterialId();
      if(mat!=matneigh) {
	if(!newnode1) {
	  g->MiddleSideNode2d(neighel,side,newnode);
	  if(newnode<0 && RefNodes.NElements()<newnode+1) RefNodes.Resize(newnode+1,0);
	  elnodes.Push(newnode);
	  elnodes.Push(node);
	  g->MiddleSideNode2d(gel,gelside,newnode1);
	  if(newnode1<0 && RefNodes.NElements()<newnode1+1) RefNodes.Resize(newnode1+1,0);
	  elnodes.Push(newnode1);
	}
	else {
	  int newnode2;
	  g->MiddleSideNode2d(neighel,side,newnode2);
	  if(newnode2<0 && RefNodes.NElements()<newnode2+1) RefNodes.Resize(newnode2+1,0);
	  elnodes.Push(newnode2);
	  elnodes.Push(node);
	  elnodes.Push(newnode);
	  newnode = newnode2;
	}

	g->CreateGeoEl2d(elnodes,mat,indexgel);
	elnodes.Resize(0);
	mat = matneigh;
      }

      g->IndexElement(neighel,nelems,indexgel);
      elnodes.Push(centers[indexgel]);
      side = (side+2)%3;
      neigh = neighel->Neighbour(side);
      side = neigh.Side();
      neighel = neigh.Element();
    }

    if(newnode1) {
      elnodes.Push(newnode1);
      elnodes.Push(node);
      elnodes.Push(newnode);
    }
    g->CreateGeoEl2d(elnodes,mat,indexgel);
  }

  /**Constructing a computational mesh*/
  nelems = g->NElements();
  TPZAdmChunkVector<TPZGeoEl *> newelmap(nelems);
  for(i=0;i<nelems;i++) newelmap[i] = (*elmap)[i];
  for(i=0;i<nelems-firstnewel;i++) (*elmap)[i] = (*elmap)[firstnewel++];
  for(;i<nelems;i++) (*elmap)[i] = 0;
  elmap->Resize(nelems-firstnewel);
  g->BuildConnectivity();
  elmap->Resize(nelems);
  for(i=0;i<firstnewel;i++) (*elmap)[nelems-firstnewel+i] = newelmap[i];
  c->AutoBuild();
}
void TDatafile2d::DualMeshCM(TPZCompMesh *c) {
	TGeoMesh *g = (TGeoMesh *)(c->Reference());
  int i,j,mat,matneigh=0,firstnewel=-1;

  //We store a)the node pointers for each new poligonal elements into nodes, b)the number of
  //nodes of each new elements, c)the associated node for each new computational element
  //Warning : for simplicity we assume the material is unique, for now.
  //The sum of sizes of all nodes is 3*(fNElements+fNBCNodes);
  int side,nnodes=g->NNodes(),nelems=g->NElements();
  if(!nnodes) { PZError << "\nERROR : TDatafile2d::DualMeshBB. Geometric mesh is empty"; return; }
  TPZVec<int> mask(nnodes,0);  //mask[id] identify whether (id+1) node has poligonal element

  TPZAdmChunkVector<TPZGeoEl *> *elmap=&(g->ElementVec());
  elmap->CompactDataStructure();

  TPZVec<int> centers(nelems);   //To store the indexes of the center nodes
  ((TGeoMesh *)g)->AddCenterNodes(centers,2,fCenterType); //Into fNodeMap have fNNodes vertexs plus fNTriangle baricenters
  TPZVec<int> midsides(0);
  g->AddMiddle1dSideNodes(midsides);

  TPZGeoEl *gel = 0;
  TPZGeoElSide neigh;

  /**Identifying two consecutive boundary nodes to start the construction of the boundary elements*/
  int node,lastnode;
  g->FindTwoBoundaryNodes(gel,side,node,lastnode);
  if(node < 0) {
    PZError << "TDatafile2d::DualMeshCC. Don't exist any boundary side?\n";
    return;
  }






  //We are working with the others bc nodes (clock).
  int newnode,newnode1,firstnewnode,indexgel,k=1;                //indexes
  //To new bc nodes in boundary middle edge
  TPZStack<int> bcnodes;   //To store index boundary node and bc id
  TPZStack<int64_t> elnodes;

  mat = gel->MaterialId();
  g->MiddleSideNode2d(gel,side,newnode);
  if(newnode<0) {
    if(RefNodes.NElements()<newnode+1) RefNodes.Resize(newnode+1,0);
    RefNodes[newnode]=RefNodes[node];
  }
  g->IndexElement(gel,nelems,indexgel);
  firstnewnode = newnode;

  while(k) {
    elnodes[0]=node;
    elnodes[1]=newnode;
    g->IndexElement(gel,nelems,indexgel);
    elnodes[2]=centers[indexgel];

    bcnodes.Push(newnode);
    bcnodes.Push(RefNodes[node]);

    side=(side+2)%3;
    neigh = gel->Neighbour(side);
    //Filling the nodes
    while(neigh.Element()) {  //loop over elements connected with nodefirst
      matneigh = neigh.Element()->MaterialId();
      if(mat!=matneigh) {
	matneigh = -1;
	break;
      }
      side = neigh.Side();
      gel = neigh.Element();
      g->IndexElement(gel,nelems,indexgel);
      elnodes.Push(centers[indexgel]);
      side = (side+2)%3;
      neigh = gel->Neighbour(side);
    }
    if(lastnode==node && matneigh!=-1) {
      elnodes.Push(firstnewnode);
      k=0; //Don't exist more boundary cells. Out to loop
    }
    else {
      if(matneigh != -1) node = gel->SideNodeIndex(side,0);
      g->MiddleSideNode2d(gel,side,newnode);
      if(newnode<0) {
	if(RefNodes.NElements()<newnode+1) RefNodes.Resize(newnode+1,0);
	if(matneigh != -1) RefNodes[newnode]=RefNodes[node];
      }
      elnodes.Push(newnode);
    }

    g->CreateGeoEl2d(elnodes,mat,indexgel);
    if(firstnewel<0) firstnewel=indexgel;

    /**Fill the boundary condition number over the boundary sides
       newgel->SetEdge(0,-fRefNodes[(int)((TNode *)nodes[0])->Id()-1]);
       if(k!=0) newgel->SetEdge(n0-1,-fRefNodes[(int)node->Id()-1]);
       else newgel->SetEdge(n0-1,-fRefNodes[(int)gel->EdgeNodeId(edge,0)-1]);
    */
    if(matneigh != -1) mask[node]=1;         //mask the node with associated element
    else {
      gel = neigh.Element();
      side = neigh.Side();
    }
    mat = gel->MaterialId();
    elnodes.Resize(3);
    matneigh = 0;
  }

  //Construct the internal poligonal elements. We use the information in fNodesT
  for(i=0;i<nnodes;i++) {
    elnodes.Resize(0);
    if(mask[i]!=0) continue;
    node = -1;
    //else we create a new internal poligonal element associated with j+1 node
    //find element id with j+1 node as first vertex of a edge
    while(1) {
      j=0;
      gel = g->ElementVec()[j++];
      for(k=0;k<gel->NNodes();k++)
	if(i==gel->NodeIndex(k)) { node=i; break; }
      if(k<0) break;
    }
    int gelside;
    for(k=0;k<gel->NSides();k++) {
      if(i==gel->SideNodeIndex(k,0)) {
	gelside=k;
//	k=-1;
	break;
      }
    }

    g->IndexElement(gel,nelems,indexgel);
    elnodes.Push(centers[indexgel]);
    side=(gelside+2)%3;
    mat = gel->MaterialId();
    neigh=gel->Neighbour(side);
    side=neigh.Side();  //neigh exist
    TPZGeoEl *neighel = neigh.Element();
    newnode1=0;
    while(gel!=neighel) {
      /**If gel and neigh have different materials*/
      matneigh = neigh.Element()->MaterialId();
      if(mat!=matneigh) {
	if(!newnode1) {
	  g->MiddleSideNode2d(neighel,side,newnode);
	  if(newnode<0 && RefNodes.NElements()<newnode+1) RefNodes.Resize(newnode+1,0);
	  elnodes.Push(newnode);
	  elnodes.Push(node);
	  g->MiddleSideNode2d(gel,gelside,newnode1);
	  if(newnode1<0 && RefNodes.NElements()<newnode1+1) RefNodes.Resize(newnode1+1,0);
	  elnodes.Push(newnode1);
	}
	else {
	  int newnode2;
	  g->MiddleSideNode2d(neighel,side,newnode2);
	  if(newnode2<0 && RefNodes.NElements()<newnode2+1) RefNodes.Resize(newnode2+1,0);
	  elnodes.Push(newnode2);
	  elnodes.Push(node);
	  elnodes.Push(newnode);
	  newnode = newnode2;
	}

	g->CreateGeoEl2d(elnodes,mat,indexgel);
	elnodes.Resize(0);
	mat = matneigh;
      }

      g->IndexElement(neighel,nelems,indexgel);
      elnodes.Push(centers[indexgel]);
      side = (side+2)%3;
      neigh = neighel->Neighbour(side);
      side = neigh.Side();
      neighel = neigh.Element();
    }

    if(newnode1) {
      elnodes.Push(newnode1);
      elnodes.Push(node);
      elnodes.Push(newnode);
    }
    g->CreateGeoEl2d(elnodes,mat,indexgel);
  }

  /**Constructing a computational mesh*/
  nelems = g->NElements();
  TPZAdmChunkVector<TPZGeoEl *> newelmap(nelems);
  for(i=0;i<nelems;i++) newelmap[i] = (*elmap)[i];
  for(i=0;i<nelems-firstnewel;i++) (*elmap)[i] = (*elmap)[firstnewel++];
  for(;i<nelems;i++) (*elmap)[i] = 0;
  elmap->Resize(nelems-firstnewel);
  g->BuildConnectivity();
  elmap->Resize(nelems);
  for(i=0;i<firstnewel;i++) (*elmap)[nelems-firstnewel+i] = newelmap[i];
  TPZGeoEl1d::SetCreateFunction(&(TCompEl1dWI::CreateElDiscWI));
  TPZGeoElPoint::SetCreateFunction(&(TCompElPointWI::CreateElDiscWI));
  TPZGeoElT2d::SetCreateFunction(&(TCompElT2dWI::CreateElDiscWI));
  c->AutoBuild();






  /*	TGeoMesh *g = c->Reference();
	if(fNNodes!=g->NNodes()) { //because fNodeMap must to contain all vertexs
   	cout << "\nERROR : TDatafile2d::DualMeshBM. Geometric triangular mesh is not initial.";
	return;
	}
	int j,edge,n0=4;  //4 is the minimum number of edges to poligonal elements
	int *mask=new int[fNNodes];  //mask[id] identify whether (id+1) node has poligonal element
	for(j=0;j<fNNodes;j++) mask[j]=0;

	VoidPtrVec centers(fNElements);
	g->AddCenterNodes(centers,2,fCenterType); //Create new nodes, with all baricenters.
	//To new bc nodes
	VoidPtrVec bcnodes(fNEdges);
	TIntVec idbc(fNEdges,0); // we can improve the length of the vector
	VoidPtrVec midedges(fNEdges);
	g->AddMiddleEdgeNodes(midedges,fRefNodes,bcnodes,idbc,3);//Create new node for each middle edge point

	int idnewel=(int)g->CreateCellId();
	if(idnewel!=fNElements) cout<<"Warning DualMeshBM. idnewel != fNElements"; //Is idnewel=g->NCells() ?

	TNode *node;
	TNode *lastnode;
	VoidPtrVec nodes(n0);
	int iANodes=0;
	TVoidPtrMap gmap(0);
	VoidPtrVec AssocNodes(fNNodes); //associated node for each gel

	if(!g->NNodes()) cout << "\nERROR : TDatafile2d::DualMeshBM. Geometric mesh is empty";
	if(g->NCells()!=fNElements)
   	cout << "\nERROR : TDatafile2d::DualMeshBM. Geometrical mesh isn't triangular.";

	TGeoCellT *gel;  //gel=g->ElBCMap()[0]->fCell;
	TGeoCellT *neigh;
	TGeoCell *newgel;

	//Construct boundary cells
	Pix i = g->CellMap().first();
	while(i) {
	gel =(TGeoCellT *)g->CellMap().contents(i);
	g->CellMap().next(i);
	//Identify any boundary edge.
	for(edge=0;edge<3;edge++) {
	if(!gel->Neighbour(edge)) {
	node=gel->EdgeNode(edge,0);
	lastnode=gel->EdgeNode(edge,1);
	edge=5;
	break;
	}
	}
	if(edge>4) break;
	}
	int idnode,idgel=(int)gel->Id();
	//We are working with the others bc nodes (clock).
	int k=1;
	while(k) {
   	n0=4; nodes.resize(n0);
	idnode=(int)node->Id();
	nodes[0]=node;
	nodes[1]=midedges[3*idgel+edge];
	nodes[2]=centers[idgel];
	edge=(edge+2)%3;
	neigh=(TGeoCellT *)gel->Neighbour(edge);
	//Filling the nodes
	while(neigh) {  //loop over elements connected with nodefirst
      	n0+=2;
	nodes.resize(n0);
      	edge=gel->NeighbourEdge(edge);
      	gel=neigh;
	idgel=(int)gel->Id();
	nodes[n0-3]=midedges[3*idgel+edge];
      	nodes[n0-2]=centers[idgel];
      	edge=(edge+2)%3;
	neigh=(TGeoCellT *)gel->Neighbour(edge);
	}
	nodes[n0-1]=midedges[3*idgel+edge];
	AssocNodes[iANodes++]=node;
	if(lastnode==node) k=0;
	node=gel->EdgeNode(edge,0);
	//Construct the poligonal geometrical elements
	if(n0==3) newgel=new TGeoCellT(idnewel,nodes);  //falta material
	else if(n0==4) newgel=new TGeoCellQ(idnewel,nodes);  //falta material
	else newgel=new TGeoCell(idnewel,nodes,n0);  //falta material
	gmap[idnewel]=newgel;
   	g->CellMap()[idnewel++]=newgel;

	//Fill the boundary condition number over the boundary sides
	newgel->SetEdge(0,-fRefNodes[(int)((TNode *)nodes[0])->Id()-1]);
	newgel->SetEdge(n0-1,-fRefNodes[(int)node->Id()-1]);

	mask[idnode-1]=1;         //mask the node with associated element
  	}

	//Construct the internal poligonal elements. We use the information in fNodesT
	for(j=0;j<fNNodes;j++) {
	if(mask[j]!=0) continue;
	n0=2; nodes.resize(n0);
	//else we create a new internal poligonal element associated with j+1 node
	idnode=j+1;  //not necessary
	idgel=-1;
	k=0;
	//find element id with j+1 node as first vertex of a edge
	while(idgel<0) {
      	if(idnode==fNodesT[k]) {
	edge=k%3;
	idgel=k/3;
	break;
	}
	k++;
	}
	gel=(TGeoCellT *)g->FindCell(idgel);
	node=gel->EdgeNode(edge,0);
	if(idnode!=node->Id()) cout << "ERROR : TDatafile2d::DualMeshBM. Incompatible idnode";
	edge=(edge+2)%3;
	neigh=(TGeoCellT *)gel->Neighbour(edge);
	edge=gel->NeighbourEdge(edge);  //neigh exist
	while(gel!=neigh) {
      	nodes[n0-2]=centers[idgel];
	idgel=(int)neigh->Id();
	nodes[n0-1]=midedges[3*idgel+edge];
      	n0+=2;
	nodes.resize(n0);
	edge=(edge+2)%3;
	k=neigh->NeighbourEdge(edge);
	neigh=(TGeoCellT *)neigh->Neighbour(edge);
	edge=k;
	}
	nodes[n0-2]=centers[idgel];
	nodes[n0-1]=midedges[3*((int)gel->Id())+edge];
	AssocNodes[iANodes++]=node;
	if(n0==3) newgel=new TGeoCellT(idnewel,nodes);  //falta material
	else if(n0==4) newgel=new TGeoCellQ(idnewel,nodes);  //falta material
	else newgel=new TGeoCell(idnewel,nodes,n0);  //falta material
	gmap[idnewel]=newgel;
   	g->CellMap()[idnewel++]=newgel;
  	}*/
}
void TDatafile2d::DualMeshEdge(TPZCompMesh *c) {
  TPZGeoMesh *g = c->Reference();
  TPZGeoMesh gdual(*g);

//  int edge,n0=4,k;  //3 is the minimum number of edges to poligonal elements
  /*   VoidPtrVec centers(fNElements);
       g->AddCenterNodes(centers,2,fCenterType); //Create new nodes, with all baricenters.

       k=0;

       int idnewel; //idfirstnewel is useful to the connectivity of the dual mesh
       idnewel=(int)g->CreateCellId();
       if(idnewel!=fNElements) cout<<"Warning DualMeshBV. idnewel != fNElements"; //Is idnewel=g->NCells()
       VoidPtrVec nodes(n0);

       if(!g->NNodes()) cout << "\nERROR : TDatafile2d::DualMeshBV. Geometric mesh is empty";
       if(g->NCells()!=fNElements)
       cout << "\nERROR : TDatafile2d::DualMeshBV. Geometrical mesh isn't triangular.";

       int idgel;
       TGeoCellT *gel;
       TGeoCellT *neigh;
       TGeoCellT *gelt;
       TGeoCellQ *gelq;

       TIntVec mask(3*fNElements,0);
       //Constructing the geometrical poligonal elements associated with each triangle edges
       //If edge is boundary then the new element is triangular otherwise is quadrilateral
       Pix i = g->CellMap().last();
       TGeoCellT *lastgel=(TGeoCellT *)g->CellMap().contents(i);
       i=g->CellMap().first();

       while(k<fNEdges) {
       gel = (TGeoCellT *) g->CellMap().contents(i);
       g->CellMap().next(i);
       idgel=(int)gel->Id();
       for(edge=0;edge<3;edge++) {
       if(mask[3*idgel+edge]) continue;
       neigh=(TGeoCellT *)gel->Neighbour(edge);
       if(!neigh) {
       n0=3;
       nodes.resize(n0);
       nodes[0]=gel->EdgeNode(edge,0);
       nodes[1]=gel->EdgeNode(edge,1);
       nodes[2]=centers[idgel];

       gelt=new TGeoCellT(idnewel,nodes);  //falta material
       gdual.CellMap()[idnewel]=gelt;
       g->CellMap()[idnewel++]=gelt;

       gelt->SetEdge(0,-fRefNodes[((int)((TNode *)nodes[0])->Id())-1]);
       k++;
       n0=4;
       nodes.resize(n0);
       }
       else {
       nodes[0]=gel->EdgeNode(edge,0);
       nodes[1]=centers[(int)neigh->Id()];
       nodes[2]=gel->EdgeNode(edge,1);
       nodes[3]=centers[idgel];

       gelq=new TGeoCellQ(idnewel,nodes);  //falta material
       gdual.CellMap()[idnewel]=gelq;
       g->CellMap()[idnewel++]=gelq;

       k++;
       mask[3*idgel+edge]=mask[3*((int)neigh->Id())+(gel->NeighbourEdge(edge))]=1;
       }
       }
       if(gel==lastgel) {
       if(k!=fNEdges)
       cout<<"Warning: DualMeshBV. Computational elements constructed less that expected.\n";
       break;
       }
       }*/
}

//Identify extension of a modulef file data
void TDatafile2d::FileExtension(char *ext) {
  int i=0,j;
  if(fFilename[0]=='\0') {
    PZError << "\n ERROR : Give filename with mesh data. Again \n\tFilename : ";
	std::cin >> fFilename;
    FileExtension(ext);
    i++;
    if(i==3) return;
  }
  i=1;
  while(fFilename[i]!='\0' && i<256) i++;
  if(i==256) {
    error("TDatafile2d::Extension","Filename is very long.");
    return;
  }
  j=i-1;
  while(fFilename[i]!='.' && i>0) i--;
  if(i==0) {
    error("TDatafile2d::Extension","Filename unknow, have not extension.");
    return;
  }
  strncpy(ext,fFilename+i+1,j-i);
  //No caso de utilizar strcpy ou strncpy observar que nao coloca '\0' no final
  //da copia, dai e conveniente acrescentar um caracter null no final da copia
}

//To read fFilename extension and identify filemodulef type: .am_fmt,.mesh,.msh,.amdba
void TDatafile2d::Read(TPZGeoMesh *g,TPZStack<int> &materialtypes) {
  if(!g) error("TDatafile2d::Read","Geometrical mesh is null");
  char ext[250];
  FileExtension(ext);
	
  //Verifying if file is Modulef file type.
  if(!strncmp(ext,"am_fmt",6)) Readam_fmt(g,materialtypes);
  else if(!strncmp(ext,"amdba",5)) Readamdba(g,materialtypes);
  else if(!strncmp(ext,"msh",3)) Readmsh(g,materialtypes);
  else if(!strncmp(ext,"mesh",4)) Readmesh(g,materialtypes);
  else if(!strncmp(ext,"fem",3)) Readara(g,materialtypes);
  else if(!strncmp(ext,"gid",3)) Readgid(g,materialtypes);
  else error("TDatafile2d::Read","Extension filename is bad");
}

/**The files from the modulef are the several types : .am_fmt, .amdba, .mesh, .msh
In file .am_fmt ->
		NNodes , NTriangles , caracters
		Array giving the three ids of the nodes to each triangles (counterclock).
		Array giving the two coordinates to each nodes.
		Array giving the reference number of the triangules(identify material number).
		Array giving the reference number of the nodes */
int TDatafile2d::Readam_fmt(TPZGeoMesh *g,TPZStack<int> &materialtypes) {
  if(fFileError) return -1;
  char linebuf[256];
  int index,nnodes,nelems;
  fDataMesh >> nnodes >> nelems;
  RefNodes.Resize(nnodes,0);


  TPZVec<REAL> coord(3,0.);
  int nnodesbyel = 3;
  TPZVec<int64_t> nodeindexes(nnodesbyel);
  fDataMesh.getline(linebuf,256);

  int i,j,firstindexel,firstindexnode;
  for(i=0;i<nelems;i++) {
    for(j=0;j<nnodesbyel;j++) {
		  fDataMesh >> index;
      nodeindexes[j] = index-1;
    }
    ((TGeoMesh *)g)->CreateGeoEl2d(nodeindexes,-1,index);
    if(i==0) firstindexel = index;
  }
  for(i=0;i<nnodes;i++) {
    fDataMesh >> coord[0] >> coord[1];
    index = g->NodeVec().AllocateNewElement();
    g->NodeVec()[index].Initialize(coord,*g);
    if(i==0)
      firstindexnode = index;
  }
	int indexmat;
  for(i=0;i<nelems;i++) {
    fDataMesh >> index;
		if(!i) indexmat = index+1;
		if(indexmat!=index) {
		  materialtypes.Push(index);
			indexmat = index;
		}
    g->ElementVec()[firstindexel++]->SetMaterialId(index);
  }

  TPZVec<int> RefNodes(2*nnodes);
  for(i=0;i<nnodes;i++) {
    RefNodes[2*i] = firstindexnode++;
    fDataMesh >> RefNodes[2*i+1];
  }
  return 1;
}
// In file .amdba ->
//		NNodes, NTriangles, caracters
//		One line for each node with : id of node, XfCoord, YfCoord, reference number of node
//		One line for each triangle with : id of triangle, ids of nodes triangle, material
int TDatafile2d::Readamdba(TPZGeoMesh *g,TPZStack<int> &materialtypes) {
  if(fFileError) return -1;
  int index,i,j,mat,n,nnodes,nelems;
  fDataMesh >> nnodes >> nelems;
	GetCommentary(fDataMesh);

  TPZVec<int> RefNodes(2*nnodes,-1);
  TPZVec<REAL> coord(3,0.);
  int nnodesbyel = 3;
  TPZVec<int64_t> nodeindexes(nnodesbyel);

  for(i=0;i<nnodes;i++) {
    fDataMesh >> n;
    fDataMesh >> coord[0] >> coord[1];
    index = g->NodeVec().AllocateNewElement();
    g->NodeVec()[index].Initialize(coord,*g);
    RefNodes[2*i] = index;
    fDataMesh >> RefNodes[2*i+1];
  }
  int indexmat;
  for(i=0;i<nelems;i++) {
    fDataMesh >> n;
    for(j=0;j<nnodesbyel;j++) {
      fDataMesh >> index;
      nodeindexes[j] = index-1;
    }
    fDataMesh >> mat;
		if(!i) indexmat = mat+1;
		if(indexmat!=mat) {
		  materialtypes.Push(mat);
			indexmat = mat;
		}
    if(nnodesbyel==3 || nnodesbyel==6)
      new TPZGeoElT2d(nodeindexes,mat,*g);
    else if(nnodesbyel==4 || nnodesbyel==8)
      new TPZGeoElQ2d(nodeindexes,mat,*g);
    else {
      PZError << "TDatafile2d::Readam_fmt. Bad parameter nnodesbyel.\n";
      return index;
    }
  }
  return 1;
}
// In file .msh ->
//		NNodes, NTriangles,
//		One line for each node with : XfCoord, YfCoord, reference number of node
//		One line for each triangle with : ids of nodes triangle, id of subdomain
int TDatafile2d::Readmsh(TPZGeoMesh *g,TPZStack<int> &materialtypes) {
  if(fFileError) return -1;
  int i, j, mat, nelembcs;
  int64_t index, nnodes, nelems;
  fDataMesh >> nnodes >> nelems >> nelembcs;

  TPZVec<int64_t> RefNodes(2*nnodes,-1);
  TPZVec<REAL> coord(3,0.);
  int nnodesbyel = 3;
  TPZVec<int64_t> nodeindexes(nnodesbyel);

  for(i=0;i<nnodes;i++) {
    fDataMesh >>  coord[0] >> coord[1];
    index = g->NodeVec().AllocateNewElement();
    g->NodeVec()[index].Initialize(coord,*g);
    RefNodes[2*i] = index;
    fDataMesh >> RefNodes[2*i+1];
  }
	int indexmat;
  for(i=0;i<nelems;i++) {
    for(j=0;j<nnodesbyel;j++) {
      fDataMesh >> index;
      nodeindexes[j] = index-1;
    }
    fDataMesh >> mat;
		if(!i) indexmat = mat+1;
		if(indexmat!=mat) {
		  materialtypes.Push(mat);
			indexmat = mat;
		}
    if(nnodesbyel==3 || nnodesbyel==6)
      new TPZGeoElT2d(nodeindexes,mat,*g);
    else if(nnodesbyel==4 || nnodesbyel==8)
      new TPZGeoElQ2d(nodeindexes,mat,*g);
    else {
      PZError << "TDatafile2d::Readmsh. Bad parameter nnodesbyel.\n";
      return 1;
    }
  }
  return 1;
}
// In file .mesh ->
//		NTriangles, NNodes, NEdgesBC, NSubDomains, caracters
//		Data nodes : XfCoord, YfCoord, ZfCoord, IdNode
//		Data edges : For each boundary edge, Id initial node and Id end node
//		Data triangle : three ids nodes of triangle, id of material
//		Data triangle : reference number of subdomain of the triangle
//		Data subdomain : reference of subdomains
//		Data nodes : reference number of the nodes over the boundary component type
//		Data edges : reference number of the edges over the boundary component type
//		Others especifications.
int TDatafile2d::Readmesh(TPZGeoMesh *g,TPZStack<int> &materialtypes) {
  if(fFileError) return -1;
  int index,nnodes,nelems,i,j,NEdgesBC,NSubDomains,mat;
  char linebuf[256];
  char ch;
  fDataMesh >> nelems >> nnodes >> NEdgesBC >> NSubDomains;
  fDataMesh.getline(linebuf,256);
	fDataMesh.getline(linebuf,256);

  TPZVec<REAL> coord(3,0.);
  int nnodesbyel = 3;
  TPZVec<int64_t> nodeindexes(nnodesbyel);
  TPZVec<int> RefNodes(2*nnodes);

  int DataEdge[3]; //For each edge : initial node, last node and
  int firstindexel,firstindexnode;
  for(i=0;i<nnodes;i++) {
    fDataMesh >>  coord[0] >> coord[1] >> coord[2] >> i;
    index = g->NodeVec().AllocateNewElement();
    g->NodeVec()[index].Initialize(coord,*g);
    if(i==0) firstindexnode = index;
  }

  fDataMesh.get(ch);
  fDataMesh.getline(linebuf,256);
  for(i=0;i<NEdgesBC;i++) fDataMesh >> DataEdge[0] >> DataEdge[1];

  fDataMesh.get(ch);
  fDataMesh.getline(linebuf,256);
  for(i=0;i<nelems;i++) {
    for(j=0;j<nnodesbyel;j++) {
      fDataMesh >> index;
      nodeindexes[j] = index-1;
    }
    fDataMesh >> mat;
    if(nnodesbyel==3 || nnodesbyel==6)
      new TPZGeoElT2d(nodeindexes,mat,*g);
    else if(nnodesbyel==4 || nnodesbyel==8)
      new TPZGeoElQ2d(nodeindexes,mat,*g);
    else {
      PZError << "TDatafile2d::Readamesh. Bad parameter nnodesbyel.\n";
      return 1;
    }
    if(!i) firstindexel = g->NElements()-1;
  }

  fDataMesh.get(ch);
  fDataMesh.getline(linebuf,256);
  int indexmat;
  for(i=0;i<nelems;i++) {
    fDataMesh >> mat;
		if(!i) indexmat = mat+1;
		if(indexmat!=mat) {
		  materialtypes.Push(mat);
			indexmat = mat;
	  }
    g->ElementVec()[firstindexel++]->SetMaterialId(mat);
  }
  fDataMesh.get(ch);
  fDataMesh.getline(linebuf,256);
  for(i=0;i<NSubDomains;i++) for(j=0;j<4;j++) fDataMesh.getline(linebuf,256);
  fDataMesh.getline(linebuf,256);
  for(i=0;i<nnodes;i++) {
    RefNodes[2*i] = firstindexnode++;
    fDataMesh >> RefNodes[2*i+1]; //reference number of nodes over the boundary component
  }
  fDataMesh.get(ch);
  fDataMesh.getline(linebuf,256);
  //Pointer of the edges over the boundary
  for(i=0;i<NEdgesBC;i++) fDataMesh >> DataEdge[2];
  //other references of boundary component type
  return 1;
}

// In file .ara -> This file is from Safe (aranha).
int TDatafile2d::Readara(TPZGeoMesh *g,TPZStack<int> &materialtypes) {
  if(fFileError) return -1;
  char lbuf[256];
  int i,j,k,n,m,nnodT=0;

  int nnodes,nelems=0,index;
  fDataMesh.getline(lbuf,256);    // Read *TITLE
  do {
    fDataMesh.getline(lbuf,256);
  }while(strncmp(lbuf,"*COORDINATES",12));

  fDataMesh >> nnodes;

  TPZVec<REAL> coord(3,0.);
  int nnodesbyel = 3,firstindexnode;
  TPZVec<int64_t> nodeindexes(nnodesbyel);

  for(i=0;i<nnodes;i++) {
    fDataMesh >> n >> coord[0] >> coord[1];
    index = g->NodeVec().AllocateNewElement();
    g->NodeVec()[index].Initialize(coord,*g);
    if(i==0) firstindexnode = index;
  }

  TPZVec<int> RefNodes(2*nnodes,0);

  fDataMesh.get(lbuf[0]);
  fDataMesh.getline(lbuf,256);    // Read *ELEMENT_GROUPS
  fDataMesh >> m;                 // m = number of areas with different material
  TPZVec<int> NElem(m);
  TPZVec<int> NodesByEl(m);
  for(j=0;j<m;j++) {
    fDataMesh >> i >> NElem[j];
//    nelems += NElem[j];
    for(i=0;i<4;i++) fDataMesh.get(lbuf[0]);   //Read TRI_
    fDataMesh >> NodesByEl[j];
    nnodT += (NodesByEl[j]*NElem[j]);
  }
  if(nnodT!=nnodes) PZError << "TDatafile2d::Readara. Bad number of nodes.\n";

  fDataMesh.get(lbuf[0]);
  fDataMesh.getline(lbuf,256);   // Read *INCIDENCES
  for(j=0;j<m;j++) {
    nnodesbyel = NodesByEl[j];
    for(i=0;i<NElem[j];i++) {
      fDataMesh >> n;
      nodeindexes.Resize(nnodesbyel);
      for(k=0;k<NodesByEl[j];k++)
	fDataMesh >> nodeindexes[k];
//      index = g->ElementVec().AllocateNewElement();
      if(nnodesbyel==3 || nnodesbyel==6)
        new TPZGeoElT2d(nodeindexes,-1,*g);
      else if(nnodesbyel==4 || nnodesbyel==8)
        new TPZGeoElQ2d(nodeindexes,-1,*g);
      else {
        PZError << "TDatafile2d::Readara. Bad parameter nnodesbyel.\n";
        return 1;
      }
    }
  }

  fDataMesh.get(lbuf[0]);
  fDataMesh.getline(lbuf,256);   // Read *LINES -> Boundary lines
  fDataMesh >> m;                // number of different boundary curves
  for(j=0;j<m;j++) {
    fDataMesh >> n >> nnodT;
    for(i=0;i<nnodT;i++) {
      fDataMesh >> k;
      RefNodes[2*(k-1)] = firstindexnode+k-1;
      RefNodes[2*k-1]=n;
    }
  }
  fDataMesh.get(lbuf[0]);
  fDataMesh.getline(lbuf,256);   // Read *KEYPOINTS

  return 1;
}
// In file .gid -> This file is from Gid.
int TDatafile2d::Readgid(TPZGeoMesh* /*g*/,TPZStack<int> &/*materialtypes*/) {
  if(fFileError) return -1;
  int nnodesbyel = 3;
  TPZVec<int> nodeindexes(nnodesbyel);
  return 1;
}

void TDatafile2d::error(char *function,char *message) {
  PZError << "\nERROR : " << function << ". " << message << std::endl;
}

/*
  switch(er) {
  case 1:
  cout << "Nodes number different of the fRefNodes capacity.";
  return;
  case 2:
  cout << "fCoord capacity incompatible.";
  return;
  case 3:
  cout << "Nodes number different of the fMaterial capacity.";
  return;
  case 4:
  cout << "fNodesT capacity incompatible.";
  return;
  case 10:
  cout << "Vectors capacity in GenerateGCells incompatible.";
  return;
  case 11:
  cout << "?.";
  return;
  default:
  cout << "Error non identified";
  }
  void TDatafile2d::ConstructTriangularGCell(TGeoMesh *g,VoidPtrVec &bcnodes,TIntVec &idbc) {
  //Verifying the measure of the vectors
  if(fNNodes!=fRefNodes.Capacity()) error(1,"ConstructTriangularCell");
  if((2*fNNodes)!=fCoord.capacity()) error(2,"ConstructTriangularCell");
  if(fNElements!=fMaterial.Capacity()) error(3,"ConstructTriangularCell");
  if((3*fNElements)!=fNodesT.capacity()) error(4,"ConstructTriangularCell");

  // create the geometric nodes and fill fNodeMap, and return in bcnodes the bc nodes
  g->AddNodes(fCoord,fRefNodes,bcnodes,idbc);
  // create the geometric elements (triangular) and fill fCellMap
  g->AddGCells(3,3,fNodesT,2);

  g->BuildConnectivity();
  }

TDatafile2d::TDatafile2d() : TPZDatafile("new") {
  while(fFileError) {
    cout << "\n TDatafile2d :: Data file name (\"quit\") = ";
    fFilename[0]='\0';
    cin >> fFilename;
    if(!strncmp(fFilename,"quit",4)) exit(1);
    if(fBuffer.open(fFilename,ios::in | ios::nocreate)) fFileError=0;
    else cout << "\nBad file name. Again.\n";
  }
  fNNodes=0;
  fNElements=0;

  fNNodesByElement = 3;
}
//To construct triangular geometrical grid.
//Return bc nodes in bcnodes and nodal bc id in idbc
//To initialize boundary conditions and material - In computational grid
void TDatafile2d::InitializeBCMat(TPZCompMesh *c,istream &input,TPZMaterial *material) {
  //n = number of different bc nodes, m = number of different materials
  int n = fRefNodes.DifferentNumbers() - 1; //Quit ref zero
  int m = fMaterial.DifferentNumbers();

  //Read bc and material data, with checking the number of bcids and materials
  TBCMatFile BCMat(input,material);
  BCMat.ReadFile(c,n,m);
}

*/

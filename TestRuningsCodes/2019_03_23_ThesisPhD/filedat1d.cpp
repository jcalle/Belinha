/*******       File : data1d.c

Gera um arquivo de dados com as coordenadas, os materiais e os dados proprios
para criar uma malha geometrica e uma computacional no caso uni-dimensional.

*******/

#include <iostream>
#include <fstream>

#include "filedat1d.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "pzcmesh.h"
#include "pzerror.h"
#include "myheader.h"
#include "hadaptive.h"
#include "TPZMaterial.h"
#include "conslaw.h"
#include "pzfmatrix.h"

TDatafile1d::TDatafile1d(std::istream &input,char *name) : fData() {
  fFileError = 1;
  fFilename[0]='\0';
  GetDataCommented(input,fFilename,256);
  while(fFileError) {
    if(!strncmp(fFilename,"new",3)) break;
    else if(!strncmp(fFilename,"MESH DATA",9)) {
      fFilename[0]='\0';
//      fData.close();
      fData.open(name, std::ios::in);
//	  GetDataCommented(((istream)fData),fFilename,256);
      while(strncmp(fFilename,"MESH DATA",9)) fData.getline(fFilename,256);
      fFileType = 1;
      fFileError=0;
      return;
    }
    else if(!strncmp(fFilename,"quit",4)) exit(1);
    fData.close();
    fData.open(fFilename, std::ios::in);
    if(!fData.rdbuf()->is_open()) {
      PZError << "\nERROR(TPZDatafile::const)-> Cannot open " << fFilename << "\n";
      PZError << "\nBad file name. Data file name (\"new\" create new file; \"quit\" exit) : ";
      fFilename[0]='\0';
	  std::cin >> fFilename;
    }
    else fFileError=0;
		fFileType = 0;
  }

  if(!fFileError) {
	  std::cout << "\n New TDatafile1d :: Data file name (\"quit\" exit)= ";
    fFilename[0]='\0';
	std::cin.getline(fFilename,256);
    if(!strncmp(fFilename,"quit",4)) exit(1);
	std::ofstream Data(fFilename);
    WriteDataFile(Data);
    Data.close();
    fFileError=0;
    fData.open(fFilename);
  }
}
TDatafile1d::TDatafile1d(char *file) : fData() {
  fFileError = 1;
	fFileType = 1;
  strcpy(fFilename,file);
  while(fFileError) {
    if(!strncmp(fFilename,"new",3)) break;
    else if(!strncmp(fFilename,"quit",4)) exit(1);
    fData.close();
    fData.open(fFilename, std::ios::in);
		fFileType = 0;
    if (!fData.rdbuf()->is_open()) {
		std::cerr << "\nERROR(TPZDatafile::const)-> Cannot open " << fFilename << "\n";
		std::cerr << "\nBad file name. Data file name (\"new\" create new file; \"quit\" exit) : ";
      fFilename[0]='\0';
	  std::cin >> fFilename;
    }
    else fFileError=0;
  }

  if(fFileError) {
	  std::cout << "\n New TDatafile1d :: Data file name (\"quit\" exit)= ";
    fFilename[0]='\0';
	std::cin >> fFilename;
    if(!strncmp(fFilename,"quit",4)) exit(1);
	std::ofstream Data(fFilename);
    WriteDataFile(Data);
    Data.close();
    fFileError=0;
    fData.open(fFilename, std::ios::in);
  }
}
TDatafile1d::~TDatafile1d() {
  fData.close();
}

/** To generate one-dimensional geometrical and computational grid */
void TDatafile1d::InitializeMesh(TPZCompMesh *c,TConservationLaw *mat) {

  /**First verify if file have data */
  if(fData.eof()) { //Return 1, whether is eof of the datafile
	  std::ofstream Data(fFilename);
    WriteDataFile(Data);
    Data.close();
    fData.open(fFilename, std::ios::in);
  }

  /**Reading geometrical information */
  TPZStack<int> type;
  TPZGeoMesh *g = c->Reference();
  /**Initialize geometrical mesh*/
  if(!fFileType) Read(g,type);
  else ComputeData(g,type);

  /**Determining the geometrical connectivity between the geometrical elements*/
  g->BuildConnectivity();
  /**Initialize material and boundary conditions*/
  ReadMaterialBC(c,mat,type);
  c->AutoBuild();
}

void TDatafile1d::ComputeData(TPZGeoMesh *g,TPZStack<int> &type) {
	/** Short mesh information 1, large mesh information 2 */
	GetDataCommented(fData,fFileType);
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
  int i,dim=3,nelem,regular,nnodes;

  /**The coordinates of the first and last nodes*/
  GetDataCommented(fData,Coord1);   // Coordinates of the first node
  GetDataCommented(fData,Coord2);   // Coordinates of the last node

  /**Number of subdivisions-X and subdivisions-Y*/
  GetDataCommented(fData,nelem);   // Subdivision-x, number of intervales
  GetDataCommented(fData,regular); // Is regular(1) or irregular(0) mesh

  /**Computing dx, dy and dz*/
  dx[0] = (Coord2[0]-Coord1[0])/nelem;
  dx[1] = (Coord2[1]-Coord1[1])/nelem;
  dx[2] = (Coord2[2]-Coord1[2])/nelem;

  /**To construct the geometrical element */
  int nodesbyel;
  GetDataCommented(fData,nodesbyel);
  if(nodesbyel!=2 && nodesbyel!=3) {
    PZError << "TDatafile1d::CalculeData. Bad number of nodes by element.\n";
    return;
  }
  nnodes = (nodesbyel-1)*nelem+1;

  if(regular) WriteCoord(g,nnodes,Coord1,0,dx,dx);
  else {
    TPZVec<REAL> alfa(3,0.);
    int NIrregCells;
    GetCommentary(fData);   // Factors to modify deltas
    for(i=0;i<dim;i++) fData >> alfa[i]; //const maior que 1, modifica o dx
    GetDataCommented(fData,NIrregCells); //3*N numero de celulas irregulares

    WriteCoord(g,nnodes,Coord1,NIrregCells,dx,alfa);
  }

  /**Over the element material*/
  GetDataCommented(fData,dim);   // Number of materials
  TPZVec<int> lastel(dim);
  //Material type and number of last element with this material type
  GetCommentary(fData);
  for(i=0;i<dim;i++) {
    fData >> regular;
    type.Push(regular);
    if(dim==1) lastel[0] = nelem;
    else fData >> lastel[i];
    if(i>0 && lastel[i]<lastel[i-1]) {
      int temp;
      PZError << "TDatafile1d::CalculeData. Bad number of last element with material information.\n";
      temp = lastel[i];
      lastel[i] = lastel[i-1];
      lastel[i-1] = temp;
      temp = type[i];
      type[i] = type[i-1];
      type[i-1] = temp;
    }
    if(!(lastel[i]<nelem)) {
      lastel[i] = nelem;
      if(i+1<dim) {
        lastel.Resize(i+1);
        break;
      }
    }
  }
  if(i==dim && lastel[i-1]<nelem)
    lastel[i-1]=nelem;
  WriteElements(g,type,lastel,nodesbyel);
}

int TDatafile1d::Read(TPZGeoMesh *g,TPZStack<int> &type) {
  int nelem,dim;
  GetDataCommented(fData,nelem);   //Number of elements
  GetDataCommented(fData,dim);     //Dimension of the spatial domain

  TPZVec<REAL> Coord(3,0.);
  int mat, i, j;
  int64_t index;

  GetCommentary(fData);   //Coordinates of the geometrical nodes
  for(i=0;i<nelem+1;i++) {
    for(j=0;j<dim;j++) fData >> Coord[j];
    index = g->NodeVec().AllocateNewElement();
    g->NodeVec()[index].Initialize(Coord,*g);
  }

  int indexnode = index-nelem;

  GetDataCommented(fData,dim);   //Number of nodes by element
  type.Resize(10);

  TPZVec<int64_t> indexnodes(dim);
  GetCommentary(fData);   //node ids for elements
  for(i=0;i<nelem;i++) {
    /**Indentifying differents material types*/
    int k=0,j;
    fData >> mat;
    for(j=0;j<k;j++) {
      if(type[j]==mat) break;
      if(j==k-1) j=2*k;
    }
    if(j==2*k) type[k++] = mat;

    for(j=0;j<dim;j++) indexnodes[j] = indexnode++;
    indexnode--;
	g->CreateGeoElement(EOned,indexnodes, 1, index);
	g->Element(index)->SetMaterialId(mat);
  }

  return 0;
}

void TDatafile1d::ReadMaterialBC(TPZCompMesh *c,TConservationLaw *mat,TPZStack<int> &type) {
  /**We must to have applied BuildConnectivity of the geometrical mesh*/
  TPZGeoMesh *g = c->Reference();
  int nnodes=g->NNodes(),nelems=g->NElements();
  int i,j,k,n,id,order,refnode=0,nbcs;
  TPZVec<int> MaskNodes(nnodes,-1);
  int PointPeriodic[4] = { -1, -1, -1, -1 };
  /**Materials and boundary conditions to computational mesh*/
  for(i=0;i<type.NElements();i++) {
    if(i>0) mat = (TConservationLaw *)mat->NewMaterial();
	mat->SetId(type[i]);
	mat->SetData(fData);
	c->InsertMaterialObject(mat);
  }
  for(i=0;i<nelems;i++) {
    TPZGeoEl *gel = g->ElementVec()[i];
    n = gel->NSides()-1;
    for(j=0;j<n;j++) {
      if(gel->Neighbour(j).Element()->Index() == gel->Index() && gel->Neighbour(j).Side() == j) {
        double x[2];
        double xmedio;
        TPZGeoNode *node1 = gel->NodePtr(j);
        PointPeriodic[refnode++] = gel->Id();
        PointPeriodic[refnode++] = j;
        x[0] = node1->Coord(0);
        x[1] = node1->Coord(1);
        xmedio = x[0];
        /**Material devolve o id da condicao fronteira para o elemento. Retorna -4 para ciclico?? */
        k = mat->IdBC(&xmedio);
        if(k<0 && k!=-4) TPZGeoElBC(gel,j,k);
      }
    }
  }
  GetDataCommented(fData,nbcs);
  order = mat->NStateVariables();
  GetCommentary(fData,2);
  TPZFMatrix<STATE> value1(order,order),value2(order,1);
  for(i=0;i<nbcs;i++) {
    int typebc;
    fData >> id;       // id of boundary condition
    fData >> typebc;   // type of boundary condition

	if(typebc==5) BCPeriodic(g,PointPeriodic);
    
	for(j=0;j<order;j++) 
		for (k = 0; k < order; k++) {
			fData >> value1(j, k);
		}
    for(j=0;j<order;j++) fData >> value2(j,0);
    c->InsertMaterialObject((TPZMaterial *)mat->CreateBC(id,typebc,value1,value2));
  }
}

void TDatafile1d::BCPeriodic(TPZGeoMesh *gmesh,int *periodicids) {
  TPZGeoEl *gel, *neigh;
  gel = gmesh->ElementVec()[periodicids[0]];
  TPZGeoElSide gelside(gel,periodicids[1]);
  neigh = gmesh->ElementVec()[periodicids[2]];
  TPZGeoElSide neighside(neigh,periodicids[3]);
  gelside.SetConnectivity(neighside);
//  neighside.SetConnectivity(gelside);
}

void TDatafile1d::Print(TPZCompMesh *c, std::ostream &out) {
  c->Reference()->Print(out);
  c->Print(out);
}

//To write datafile with coordinates and material of the nodes of the one-dimensional grid
void TDatafile1d::WriteCoord(TPZGeoMesh *g,int nnodes,TPZVec<REAL> &Coord,int NIrregCells,TPZVec<REAL> &dx,TPZVec<REAL> &alfa) {

  TPZVec<REAL> newdx(dx);
  int i,j,P,Mod;
  if(NIrregCells>nnodes) NIrregCells=nnodes; //The grid is all irregular.
  P=NIrregCells/3;
  Mod=NIrregCells%3;

  int index = g->NodeVec().AllocateNewElement();
  g->NodeVec()[index].Initialize(Coord,*g);

  for(i=0;i<P;i++) {
    if(newdx[0]>2*dx[0] || newdx[0]<.5*dx[0]) newdx=dx;
    else for(j=0;j<3;j++) newdx[j]=alfa[j]*dx[j];
    for(j=0;j<3;j++) Coord[j]+=dx[j];
    index = g->NodeVec().AllocateNewElement();
    g->NodeVec()[index].Initialize(Coord,*g);
    for(j=0;j<3;j++) Coord[j]+=newdx[j];
    index = g->NodeVec().AllocateNewElement();
    g->NodeVec()[index].Initialize(Coord,*g);
    for(j=0;j<3;j++) Coord[j]+=((1./alfa[j])*dx[j]);
    index = g->NodeVec().AllocateNewElement();
    g->NodeVec()[index].Initialize(Coord,*g);
  }
  for(i=0;i<Mod;i++) {
    for(j=0;j<3;j++) Coord[j]+=newdx[j];
    index = g->NodeVec().AllocateNewElement();
    g->NodeVec()[index].Initialize(Coord,*g);
  }

  //If NIrregCells=0 then the grid is regular.
  //Regular grid part
  for(i=NIrregCells;i<nnodes-1;i++) {
    for(j=0;j<3;j++) Coord[j]+=dx[j];
    index = g->NodeVec().AllocateNewElement();
    g->NodeVec()[index].Initialize(Coord,*g);
  }
}

void TDatafile1d::WriteElements(TPZGeoMesh *g,TPZStack<int> &type,TPZVec<int> &lastelement,	int nodesbyel) {
  int nmat = type.NElements();
  if(nmat!=lastelement.NElements())
    PZError << "TDatafile1d::WriteElement. Incompatible length of type and lastelement vectors.\n";
  int i,k,j=0,nn=2;
  int64_t index;
  TPZVec<int64_t> indexnode(3);
  indexnode[0] = 0;
  indexnode[1] = 2;
  indexnode[2] = 1;
  if(nodesbyel==2) { nn=1; indexnode[1]=1; }

  for(i=0;i<nmat;i++) {
    for(;j<lastelement[i];j++) {
		g->CreateGeoElement(EOned, indexnode, type[i],index);
      for(k=0;k<nodesbyel;k++) indexnode[k] += nn;
    }
  }
}

//To write data mesh into fData file
void TDatafile1d::WriteDataFile(std::ostream &Data) {
  int i,j,val,filetype=0,regular=0,num;
  REAL value;
  char c;
  std::cout << "\nLarge information? (y/n) ";
  std::cin >> c;
  if(c=='y') filetype=1;
  Data << "Type of information form contained in this file: short(0), and large(1) information.";
  Data << std::endl << filetype << std::endl;

  if(!filetype) {
	  std::cout << "\nRegular mesh? (y/n) ";
	  std::cin >> c;
    if(c=='y') regular=1;
    Data << "Regular or irregular mesh. 1 to regular, 0 otherwise.\n";
    Data << regular << std::endl;

    //Coordinates of nodes for short information to regular or irregular mesh
	std::cout << "\nFilling data into the one-dimensional mesh data file";
	std::cout << "\n\t filetype = 0. (short mesh information)";
	std::cout << "\nNumber of cells = ";
	std::cin >> val;
    Data << "GEOMETRICAL MESH INFORMATION : Number of elements" << std::endl;
    Data << val << std::endl;
	std::cout << "Dimension of the spatial domain = ";
	std::cin >> val;
    Data << "Dimension of the spatial domain" << std::endl << val << std::endl;
	std::cout << "Coordinate of the first node, x0 ";
    if(val==2) std::cout << ": y0";
    else if(val==3) std::cout << ": y0 : z0";
	std::cout << std::endl;
    Data << "Coordinates of the first node" << std::endl;
    for(i=0;i<val;i++) {
		std::cin >> value;
      Data << value;
      if(i<val-1) Data << "\t";
      else Data << std::endl;
    }
    int nodesbyel;
	std::cout << "Number of nodes by element ";
	std::cin >> nodesbyel;
    Data << "Number of nodes by element" << std::endl << nodesbyel << std::endl;
	std::cout << "Deltas : deltaX";
    if(val==2) std::cout << "  deltaY";
    else if(val==3) std::cout << "  deltaY  deltaZ";
	std::cout << std::endl;
    Data << "Deltas on the cartesian coordinates" << std::endl;
    for(i=0;i<val;i++) {
		std::cin >> value;
      Data << value;
      if(i<val-1) Data << "\t";
      else Data << std::endl;
    }

    if(!regular) {
		std::cout << "Modificator of the deltas : alfaX";
		std::cin >> value;
      if(val==2) std::cout << "  alfaY";
      else if(val==3) std::cout << "  alfaY  alfaZ";
	  std::cout << std::endl;
      Data << "Modificator of the deltas : alfas" << std::endl;
      for(i=0;i<val;i++) {
		  std::cin >> value;
	Data << value;
	if(i<val-1) Data << "\t";
	else Data << std::endl;
      }

	  std::cout << "Number of irregular elements = ";
	  std::cin >> val;
      Data << "Number of irregular elements" << std::endl << val << std::endl;
    }

	std::cout << "Number of material types ";
	std::cin >> val;
    Data << "Number of material types" << std::endl << val << std::endl;
    Data << "Material types ids and number of the last element with this material" << std::endl;
    for(i=0;i<val;i++) {
		std::cout << (i+1) << ". Material type Id = ";
		std::cin >> num;
      Data << num << "\t";
	  std::cout << "   Number of the last element with this material = ";
	  std::cin >> num;
      Data << num << std::endl;
    }
	std::cout << "Boundary condition to two boundary geometrical nodes" << std::endl;
	std::cin >> num;
    Data << "Boundary condition to two boundary geometrical nodes" << std::endl << num << "\t";
	std::cin >> num;
    Data << num << std::endl;
  }
  else {
    int dim;
    //Coordinates of nodes for short information to regular or irregular mesh
	std::cout << "\nFilling data into the one-dimensional mesh data file";
	std::cout << "\n\t filetype = 1. (large mesh information)";
	std::cout << "\nNumber of cells = ";
	std::cin >> val;
    Data << "GEOMETRICAL INFORMATION : Number of elements" << std::endl << val << std::endl;
	std::cout << "Dimension of the spatial domain = ";
	std::cin >> dim;
    Data << "Dimension of the spatial domain" << std::endl << dim << std::endl;
    Data << "Coordinates of the geometrical nodes" << std::endl;
	std::cout << "Coordinates of the nodes,\n";
    for(i=0;i<val+1;i++) {
      for(j=0;j<dim;j++) {
		  std::cin >> value;
	Data << "\t" << value;
      }
      Data << std::endl;
    }

	std::cout << "Number of nodes by element" << std::endl;
	std::cin >> num;
    Data << "Number of nodes by element" << std::endl << num << std::endl;
	std::cout << "To each element, identifier of material\n";
    Data << "Number of material for each element" << std::endl;
    for(i=0;i<val;i++) {
      if(!(i%10)) Data << std::endl;
	  std::cin >> num;
      Data << "\t" << num;
    }
    Data << std::endl;
	std::cout << "Boundary condition (Ids) for two boundary geometrical nodes\n";
	std::cin >> num >> val;
    Data << "Boundary condition (Ids) for two boundary geometrical nodes" << std::endl;
    Data << num << val << std::endl;
  }
  std::cout << "Number of boundary conditions = ";
  std::cin >> num;
  for(i=0;i<num;i++) {
	  std::cout << "Matrix order to boundary condition " << (i+1) << " = ";
	std::cin >> val;
    Data << "Matrix order to " << (i+1) << "-bc" << std::endl << val << std::endl << "Matrices data" << std::endl;
	std::cout << "Matrices data(by line), to val1("<<val<<"x"<<val << ")" << std::endl;
    for(j=0;j<val*val;j++) {
		std::cin >> value;
      Data << "\t" << value;
    }
    Data << std::endl;
	std::cout << "to val2("<<val<<"x1)"<< std::endl;
    for(j=0;j<val;j++) {
		std::cin >> value;
      Data << "\t" << value;
    }
    Data << std::endl;
  }
}

/*

** To fill a file with the material index

for(i=0;i<N1;i++) {
if(i%10) fprintf(out,"\t%d",j);
else fprintf(out,"\n%d",j);
}
if(0<N1) j++;
for(;i<N2;i++) {
  if(i%10) fprintf(out,"\t%d",j);
  else fprintf(out,"\n%d",j);
}
if(N1<N2) j++;
for(;i<N3;i++) {
  if(i%10) fprintf(out,"\t%d",j);
  else fprintf(out,"\n%d",j);
}
if(N3<N4) j++;
for(;i<N4;i++) {
  if(i%10) fprintf(out,"\t%d",j);
  else fprintf(out,"\n%d",j);
}
}

***/

/*
  //create (boundary condition) TBC and TMaterial for the computational elements
  input.getline(buf,256);
  TBCMatFile BCMat(input);
  input.get(buf[0]);
  BCMat.ReadFile(c);
*/

/*
  // create the geometric nodes
  int k,id,idbc;
  id = g->CreateNodeId();
  for(k=0;k<fNCells+1;k++) {
  TNode *newnode=new TNode(id,1,fCoord+k);
  g->NodeMap()[id++]=newnode;
  }

  //create (boundary condition) TBC and TMaterial for the computational elements
  input.getline(buf,256);
  TBCMatFile BCMat(input);
  input.get(buf[0]);
  BCMat.ReadFile(c);

	// create the geometric elements
	VoidPtrVec Node(2);
	TBC *bc;
	TBCNodal *nodal;
	Pix i = g->NodeMap().first();
	Node[0] = g->NodeMap().contents(i);
	g->NodeMap().next(i);
	TGeoCell1d *newel, *oldel=NULL;
	id = g->CreateCellId();
	for(k=0;k<fNCells;k++) {
	Node[1] = g->NodeMap().contents(i);
	g->NodeMap().next(i);

	newel = new TGeoCell1d(id,Node,fMat[k]);
	g->CellMap()[id++] = newel;
	if(k==0) {
      	bc = (TBC *)c->BCMap().contents(c->BCMap().first());   // To first cell and bc
	idbc = c->CreateBCNodalId();
      	nodal = new TBCNodal(idbc,(TNode *)Node[0],bc,c->Equation()->Order());
	c->BCNodalMap()[idbc++]=nodal;
      	newel->SetUnilateralConnectivity(0,oldel,-(bc->Id()));
	}
	else if(k==fNCells-1) {
	bc = (TBC *)c->BCMap().contents(c->BCMap().last());   // To last cell and bc
	nodal = new TBCNodal(idbc,(TNode *)Node[1],bc,c->Equation()->Order());
	c->BCNodalMap()[idbc++]=nodal;
	newel->SetUnilateralConnectivity(1,NULL,-(bc->Id()));
	newel->SetConnectivity(0,oldel,1);
	}
	else newel->SetConnectivity(0,oldel,1);
	oldel = newel;
	Node[0]=Node[1];
	}

	void TDatafile1d::GenerateCCells(TPZCompMesh *c) {
	TCompEl1d::gOrder=3;
	Pix i;
	TGeoMesh *g = (TGeoMesh *)c->Reference();
	i=g->CellMap().first();
	TGeoCell1d *gel;
	TCompEl1d *cel;
	long id;
	while(i) {
	gel = (TGeoCell1d *)(g->CellMap()).contents(i);
	(g->CellMap()).next(i);
	cel = (TCompEl1d *)gel->CreateCompEl();
	id = cel->Id();
	c->CellMap()[id] = cel;
	// the boundary nodal dofs are created into the construction of comp. cell
	}
	c->CreateTEdge();
	}
	void TDatafile1d::GenerateCCells(TPZCompMesh *c) {
	Pix i;
	char buf[256];
	int ord;
	TGeoMesh *g = (TGeoMesh *)c->Reference();
	i=g->CellMap().first();
	TGeoCell1d *gel;
	TCompCell1d *cel;
	fData.get(buf[0]);
	fData.getline(buf,256);   // integration order
	fData >> ord;
	fData.get(buf[0]);
	long id=c->CreateCellId();
	while(i) {
	gel = (TGeoCell1d *)(g->CellMap()).contents(i);
	(g->CellMap()).next(i);
	cel = (TCompCell1d *)gel->CreateCompCell(id);
	c->CellMap()[id++] = cel;
	cel->SetIntegrationOrder(ord);
	}
	c->CreateTEdge();
	}

*/

/*******   File : modwrm.h   -   Jorge Lizardo Diaz Calle       ***/

//A classe TModulefWRM,
//permite inicializar uma malha geometrica e uma computacional a partir de um
//arquivo tipo Modulef.
//Passo 1 - Inicializa a malha geom. com os triangulos ou quadril. segundo arquivo
//Passo 2 - Inicializa a malha comp. com o tipo de celula pedida a partir dos
//   elem. triangulares.
//Passo 3 - Cada elemento comp. se diferente do triangular geo. cria o elemento geo.
//   correspondiente e adiciona na lista de elementos geom.
//O tipo de material no elemento geometrico determina a que superelemento pertence


#ifndef TINPUTFILE
#define TINPUTFILE

#include <fstream>

#include "pzreal.h"
#include "pzvec.h"

class TPZCompMesh;
class TPZGeoEl;
class TPZGeoNode;


class TDatafile2d {
	std::ifstream  fDataMesh;
	std::istream  *fData;
  char  fFilename[256];
  int	fFileError;
	int fFileType;

  /** Cell type to computational grid */
  int fCellType; //1- triang. cell, 2- CellCC, 3- CellCM, 4- CellEdge ->Default 1
  /** Center type to construct dual mesh from triangulization */
  int fCenterType; //1- centroid, 2-circuncenter, 3- incenter, 4-orthocenter -> Default 1
  /** To store a reference number of the geometrical nodes*/
  TPZVec<int64_t> RefNodes;

	/** To deleted geometrical elements to deleted sub-domains */
	int fDeleted;
	REAL fxZero;
	REAL fyZero;

 public:
  TDatafile2d(std::istream &input,char *name=0);
  TDatafile2d(char *file);
  ~TDatafile2d() { fDataMesh.close(); }

  void RequireCellType();
  void RequireCenterType();
  void SetCellType(int type);

  /** To inicialize TGeoGrid and TCompGrid from several datafiles */
  void InitializeMesh(TPZCompMesh *c, TConservationLaw *material);

  /** To initialize a computational grid from triangular geometric grid */
  void DualMeshCC(TPZCompMesh *c); //Dual mesh : center-center type -> joined the triangle centers
  void DualMeshCM(TPZCompMesh *c); //Dual mesh : center-middle type -> joined the center to edge middle point
  void DualMeshEdge(TPZCompMesh *c); //Dual mesh : edge type -> joined the center with the vertexs

  /** Select extension of the fFilename */
  void FileExtension(char *ext);
  /** Verify if the data file is modulef type and call the right function to read mesh data */
  void Read(TPZGeoMesh *g,TPZStack<int> &materialtypes);
  /** Depending of type of the modulef file put mesh data in vectors */
  int Readam_fmt(TPZGeoMesh *g,TPZStack<int> &materialtypes);
  int Readamdba(TPZGeoMesh *g,TPZStack<int> &materialtypes);
  int Readmsh(TPZGeoMesh *g,TPZStack<int> &materialtypes);
  int Readmesh(TPZGeoMesh *g,TPZStack<int> &materialtypes);
  int Readara(TPZGeoMesh *g,TPZStack<int> &materialtypes);
  int Readgid(TPZGeoMesh *g,TPZStack<int> &materialtypes);

  /**To case regular. It is needed short information*/
  void ComputeData(TPZGeoMesh *g,TPZStack<int> &materialtypes);
  /**Compute coordinates of the nodes for geometrical triangular or quadrilateral mesh*/
  void WriteCoord(TPZGeoMesh *g,int dim,int eltype,int nelemx,int nelemy,TPZVec<REAL> &Coord,
		  TPZVec<REAL> &dx);
  /**Construct geometrical elements. Triangles or quadrilaterals when domain is a square*/
  void WriteElements(TPZGeoMesh *g,int eltype,int nelemx,int nelemy,TPZStack<int> &mattypes);

  /**To print error messages*/
  void error(char *function,char *message);

  void ReadMaterialBC(TPZCompMesh *c, TConservationLaw *mat,TPZStack<int> &type);
  void BCPeriodic(TPZGeoMesh *c, std::istream &inputbc= std::cin);
  
};

#endif


/*******   File : data1d.h

TDatafile1d class facilitate the opening of the input file
with a one-dimensional mesh data.

*******/

#ifndef TDATAFILE1D
#define TDATAFILE1D

#include <fstream>

#include "pzreal.h"
#include "pzstack.h"

class TPZCompMesh;
class TPZGeoMesh;
class TConservationLaw;


//***       Jorge Lizardo Diaz Calle       ***/

class TDatafile1d {

	std::ifstream  fData;
  char  fFilename[256];
  int	fFileError;
	int fFileType;

 public:
  TDatafile1d(std::istream &input,char *name=0);
  TDatafile1d(char *file);
  ~TDatafile1d();

  /**Read the data file and construct the geometrical and computational mesh with data*/
  int Read(TPZGeoMesh *mesh,TPZStack<int> &materialtypes);

  /**To inicialize one-dimensional TPZGeoMesh and TPZCompMesh*/
  void InitializeMesh(TPZCompMesh *c,TConservationLaw *material);
  void Print(TPZCompMesh *c, std::ostream &out);

  /**To write datafile to read the coordinates and material of the nodes*/
  void ComputeData(TPZGeoMesh *g,TPZStack<int> &materialtypes);
  void WriteCoord(TPZGeoMesh *g,int nelem,TPZVec<REAL> &coord,int NIrregCells,TPZVec<REAL> &dx,TPZVec<REAL> &alfa);
  void WriteElements(TPZGeoMesh *g,TPZStack<int> &type,TPZVec<int> &lastelement,int nodesbyel);

  void WriteDataFile(std::ostream &data);
  void BCPeriodic(TPZGeoMesh *gmesh, int *periodicids);

  void ReadMaterialBC(TPZCompMesh *c,TConservationLaw *mat,TPZStack<int> &materialtypes);

};


#endif


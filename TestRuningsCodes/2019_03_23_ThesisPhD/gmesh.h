/*******       File : gmesh.h

Header file for class TGeoMesh. TGeoMesh defines a geometrical mesh and
contains several lists of the nodes, cells, coordinate systems of the
geometrical mesh.

The nodes will have the same dimension as the geometrical mesh dimension.
The cells can to be one- two- or three-dimensional still the mesh dimension
is three. This is, the cell dimension is independent to mesh dimension.

Programmer : Jorge Lizardo Diaz Calle

*******/

#ifndef GEOMESHTODISCONTINUOUS
#define GEOMESHTODISCONTINUOUS

#include "pzgmesh.h"

class TGeoMesh : public TPZGeoMesh {

 public:
  /** constructor*/
  TGeoMesh();

  int NSides1d();
  int MaxNSidesByEl();
  // To initialize or increment cells and nodes into the list
  //   void AddNodes(TPZVec<REAL> &Coord,TPZVec<int> &refbc,TPZVec<TPZGeoNode> &bcnodes,TPZVec<int> &idbc);
  void AddCenterNodes(TPZVec<int> &indexcenters,int eldim,int CenterType=1);

  void AddMiddle1dSideNodes(TPZVec<int> &indexmidnodes);
  //	void AddGCells(TPZVec<int> &nnodes,TPZVec<int> &nedges,TPZVec<int> &nodes,int celldim);
  //	void AddGCells(int nnodes,int nedges,TPZVec<int> &nodes,int celldim);


  void FindTwoBoundaryNodes(TPZGeoEl *gel,int &side,int &node,int &lastnode);
  void MiddleSideNode2d(TPZGeoEl *gel,int side,int &index);
  void IndexElement(TPZGeoEl *gel,int nelems,int &index);

  void CreateGeoEl2d(TPZVec<int64_t> &elnodes,int mat,int &indexgel);
};


#endif



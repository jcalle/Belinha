//HEADER FILE FOR CLASS ELEMC3D

#ifndef ELGC3DHPP
#define ELGC3DHPP

#include <iostream>
#include "pzgeoel.h"
#include "TPZGeoCube.h"
#include "pzgeoelrefless.h"

/** Header file for the three dimensional hexaedral element class*/

class TPZGeoElC3d : public TPZGeoElRefLess<pzgeom::TPZGeoCube> {

  /** pointer to function which creates a computational element*/
  static TPZCompEl *(*fp)(TPZGeoElC3d *geoel,TPZCompMesh &mesh, int64_t &index);

 protected:
  /**Indices of the nodes which compose the element*/
  int64_t fNodeIndexes[8];
  /**Pointers to the neighbouring elements with corresponding sides*/
  TPZGeoElSide fNeighbours[27];
  /**Pointers to subelements. If has less than 4 sub-elements the pointer is NULL*/
  TPZGeoEl *fSubEl[8];

  /** implementation of two-dimensional bilinear interpolation*/
  void Shape(TPZVec<REAL> &x,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi);

 public:
  /**Constructors. Parameters: id - element id, nodeindexes - vector containing node indexes,
     matind - material index, refind - index of the node which indicates the reference direction*/
  TPZGeoElC3d(int id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh);
  TPZGeoElC3d(TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh);

  TPZGeoElC3d *CreateGeoEl(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh);
  /**Destructor*/
  ~TPZGeoElC3d();

  TPZCompEl *CreateCompEl(TPZCompMesh &mesh, int64_t &index);

  static void SetCreateFunction(TPZCompEl *(*f)(TPZGeoElC3d *el,TPZCompMesh &mesh, int64_t &index));

  /**Return the number of nodes of the element*/
  int NNodes();

  /**Return the index within the geometric mesh of node*/
  int64_t NodeIndex(int node);
  /**Initializes the node i of the element*/
  void SetNodeIndex(int i, int64_t nodeindex) { fNodeIndexes[i]=nodeindex; }

  /**Return the number of connectivities of the element*/
  int NSides() { return 27;}//uma para cada conectividade

  /**Return the number of nodes for a particular side*/
  int NSideNodes(int side);

  int NCornerNodes() { return 8;}

  int64_t SideNodeIndex(int side,int node);

  /**Returns the midsidenode along a side of the element*/
  void MidSideNodeIndex(int side, int64_t &index);

  /**Whether the middle side node not exist it will be created*/
  void NewMidSideNode(int side, int64_t &index);

  /**Determine the coordinate of the center of the element*/
  //int CenterIndex();//é determinado por NewMidSideNode(side=)

  /**return the dimension of transformation*/
  int SideDimension(int side);

  int Dimension() { return 3; }

  /**returns the face node ids */
  void NodeFaceIds(TPZVec<int> &ids,int face);

  /** */
  TPZGeoElSide HigherDimensionSides(int side,int targetdimension);

  void LowerDimensionSides(int side,TPZStack<TPZGeoElSide> &smallsides);

  int BestDimensionSideOfTwoFaces(int face1,int face2);

  /**returns 1 if the side has not been defined by buildconnectivity
     After construction the side is undefined. The buildconnectivity method
     loops over all elements and tries to identify neighbours along their
     uninitialized sides*/
  int SideIsUndefined(int side) { return fNeighbours[side].Side() == -1;}

  /**flags the side as defined, this means no neighbouring element was found*/
  void SetSideDefined(int side) { fNeighbours[side] = TPZGeoElSide(0,0); }

  /**returns a pointer to the neighbour and the neighbourside along side of the current element*/
  TPZGeoElSide Neighbour(int side) { return fNeighbours[side]; }

  /**fill in the data structure for the neighbouring information*/
  void SetNeighbour(int side,const TPZGeoElSide &neighbour) { fNeighbours[side]=neighbour; }

  /**Returns the coordinates in master element coordinate space of the coordinates of side*/
  void SideMasterCo(int side,TPZVec<REAL> &IVec,TPZVec<REAL> &JVec);

  /**Computes the jacobian*/
  void Jacobian(TPZVec<REAL>& par, TPZFMatrix<STATE> &jacobian, TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<STATE> &jacinv);

  /**Computes the geometric location*/
  void X(TPZVec<REAL>& par, TPZVec<REAL> &result);

  void NormalVector(int side, TPZVec<REAL> &loc, TPZVec<REAL> &normal,
		    TPZFMatrix<REAL> &axes, TPZFMatrix<STATE> &jac);

  /**Divide the elements and puts the result in the point vector*/
  void Divide(TPZVec<TPZGeoEl *> &SubElVec);
  
  /**return 1 if the element has subelements along side*/
  int HasSubElement(int side) {return fSubEl[0]!=0;}

  /**return the number of subelements of the element*/
  int NSubElements();

  /**returns a pointer to the subelement is*/
  TPZGeoEl *SubElement(int is) { return fSubEl[is];}

  /** returns the number of the subelements along side*/
  int NSideSubElements(int side);
  /** return the subcell at the side and the side of the subelement over side*/
  TPZGeoElSide SideSubElement(int side,int position);

  /** returns the subelements along side in sub*/
  void SideSubElements(int side,TPZVec<TPZGeoEl *> &sub);

  /**Returns the sub elements and their sides along a side*/
  void GetSubElement(int side,TPZVec<int> &refnodes,TPZVec<TPZGeoElSide> &sub);

  /**Returns the father and corresponding side for this element,
     if the element is not connected to the father along this side return null*/
  TPZGeoElSide Father(int side);

  /**accumulates the transformation of the jacobian which maps the current
     master element space into the space of the master element of the father*/
  void BuildTransform(int side,TPZGeoEl *father,TPZTransform<STATE> &t);

  TPZTransform<STATE> SideToSideTransform(int sidefrom, int sideto);

  /**method which creates a computational boundary condition element
     based on the current geometric element, a side and a boundary condition number*/
  TPZCompEl *CreateBCCompEl(int side, int bc, TPZCompMesh &cmesh);

  /** Methods which will implement the declaration of a refinemnt topology
  **/

  /**
  * returns the father/side of the father which contains the side of the
  * sub element
  **/
  virtual TPZGeoElSide Father2(int side);

  /**
  * returns the transformation which maps the parameter side of the element/side
  * into the parameter space of the father element/side
  **/
  virtual TPZTransform<STATE> BuildTransform2(int side, TPZGeoEl *father);

  /**
   * This method will return a partition of the side of the current element
   * as the union of sub elements/side which are put in the stack
  **/
  virtual void GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel);

  static int main(TPZGeoEl *gel);
  //static REAL MidSideNode[27][3];
  
};

inline TPZCompEl *TPZGeoElC3d::CreateCompEl(TPZCompMesh &mesh,int64_t &index){
  return fp(this,mesh,index);
}

inline void TPZGeoElC3d::SetCreateFunction(TPZCompEl *(*f)(TPZGeoElC3d *el,TPZCompMesh &mesh,int64_t &index)){
  fp = f;
}

#endif




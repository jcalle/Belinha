/**File : pzelg1d.h
   Header file for class TPZGeoEl, the one dimensional element class.
*/

#ifndef ELG1DHPP
#define ELG1DHPP

#include <iostream>
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeoelrefless.h"




/**TGeoEl1d defines a one dimensional geometric element. The geometric elements
   define the shape of the finite elements, not its interpolation function*/

class TPZGeoEl1d : public TPZGeoElRefLess<pzgeom::TPZGeoLinear> {

  /** pointer to function which creates a computational element*/
  static TPZCompEl *(*fp)(TPZGeoEl1d *geoel,TPZCompMesh &mesh, int64_t &index);

 protected:
  /**Index of a geometric node which indicates the direction of the reference axis*/
  int fYAxisIndex;

  /**Indices of the nodes which compose the element*/
  int64_t fNodeIndexes[3];
  /**Pointers to the neighbouring elements with corresponding sides*/
  TPZGeoElSide fNeighbours[3];

  /**Pointers to subelements*/
  TPZGeoEl *fSubEl[2];

  /**this method creates the new nodes for the subelements*/
  virtual void CreateNewNodes(int64_t *gnod);

 public:
  /**Constructors. Parameters: id - element id, nodeindexes - vector containing node indexes,
     matind - material index, refind - index of the node which indicates the reference direction*/
  TPZGeoEl1d(int id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh,int refind = -1);
  TPZGeoEl1d(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh,int refind = -1);

  virtual TPZGeoEl1d *CreateGeoEl(TPZVec<int64_t> &nodeindices,int matind,TPZGeoMesh &mesh,int yaxisindex);
  /**Destructor*/
  ~TPZGeoEl1d();

  virtual TPZCompEl *CreateCompEl(TPZCompMesh &mesh,int64_t &index);

  static void SetCreateFunction(TPZCompEl *(*f)(TPZGeoEl1d *el,TPZCompMesh &mesh,int64_t &index));

  /**Return the number of nodes of the element*/
  virtual int NNodes() { return (fNodeIndexes[2]==-1)?2:3; }

  /**Return the index within the geometric mesh of node*/
  virtual int64_t NodeIndex(int node);

  /**Return the number of connectivities of the element*/
  virtual int NSides(){ return 3;}

  /**Return the number of nodes for a particular side*/
  virtual int NSideNodes(int side);

  virtual int NCornerNodes() { return 2;}

  virtual int SideNodeIndex(int side, int node);

  /**return the dimension of transformation*/
  virtual int SideDimension(int side);
  int Dimension() { return 1; }

  /** */
  virtual TPZGeoElSide HigherDimensionSides(int side,int targetdimension);

  virtual void LowerDimensionSides(int side,TPZStack<TPZGeoElSide> &smallsides);

  /**Returns the midsidenode along a side of the element*/
  virtual void MidSideNodeIndex(int side,int64_t &index);

  /**Initializes the node i of the element*/
  virtual void SetNodeIndex(int i,int64_t nodeindex) { fNodeIndexes[i]=nodeindex; }

  /**returns 1 if the side has not been defined by buildconnectivity
     After construction the side is undefined. The buildconnectivity method
     loops over all elements and tries to identify neighbours along their
     uninitialized sides*/
  virtual int SideIsUndefined(int side){ return fNeighbours[side].Side() == -1;}

  /**flags the side as defined, this means no neighbouring element was found*/
  virtual void SetSideDefined(int side){ fNeighbours[side] = TPZGeoElSide(0,0);}

  /**returns a pointer to the neighbour and the neighbourside along side of the current element*/
  virtual TPZGeoElSide Neighbour(int side) {
    return fNeighbours[side];
  }

  /**fill in the data structure for the neighbouring information*/
  virtual void SetNeighbour(int side, const TPZGeoElSide &neighbour){
    fNeighbours[side] = neighbour;
  }

  /**return the coordinates in the master element space of the corner nodes of side
     the coordinates of each corner node is stored in a column of the matrix*/
  virtual void SideMasterCo(int side,TPZFMatrix<REAL> &coordinates);

  /**Returns the coordinates in master element coordinate space of the coordinates of side*/
  virtual void SideMasterCo(int side,TPZVec<REAL> &IVec,TPZVec<REAL> &JVec);

  /**Computes the jacobian*/
  virtual void Jacobian(TPZVec<REAL>& par, TPZFMatrix<STATE> &jacobian, TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<STATE> &jacinv);

  /**Computes the geometric location*/
  virtual void X(TPZVec<REAL>& par, TPZVec<REAL> &result);

  virtual void NormalVector(int side, TPZVec<REAL> &loc, TPZVec<REAL> &normal,
			    TPZFMatrix<REAL> &axes, TPZFMatrix<STATE> &jac);

  /**Divide the elements and puts the result in the point vector*/
  virtual void Divide(TPZVec<TPZGeoEl *> &pv);

  /**return 1 if the element has subelements along side*/
  int HasSubElement(int side) {return fSubEl[0]!=0;}

  /**return the number of subelements of the element*/
  int NSubElements();

  /**return the number of subelements through the side*/
  int NSideSubElements(int side);

  /**returns a pointer to the subelement is*/
  TPZGeoEl *SubElement(int is) { return fSubEl[is];}
  /**return a pointer and a side of the subelement of the element at the side
     and the indicated position. position = 0 indicate first subelement, ...*/
  TPZGeoElSide SideSubElement(int side,int position);

  /**Returns the sub elements and their sides along a side*/
  virtual void GetSubElement(int side, TPZVec<int> &refnodes, TPZVec<TPZGeoElSide> &sub);

  /**Returns the father and corresponding side for this element,
     if the element is not connected to the father along this side return null*/
  virtual TPZGeoElSide Father(int side);

  /**accumulates the transformation of the jacobian which maps the current
     master element space into the space of the master element of the father*/
  virtual void BuildTransform(int side, TPZGeoEl *father, TPZTransform<STATE> &t);

  TPZTransform<STATE> SideToSideTransform(int sidefrom, int sideto);

  /**method which creates a computational boundary condition element
     based on the current geometric element, a side and a boundary condition number*/
  virtual TPZCompEl *CreateBCCompEl(int side, int bc, TPZCompMesh &cmesh);

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
  //static REAL MidSideNode[3][3];

  /**Jorge 17/7/99*/
  /** Return the measure of the geometrical element - Area */
  REAL Mesure(int dim);
	/** 
	 * Return into the center a especial point of the geometrical element
	 * If 1-d => middle point,
	 * If 2-d => barycenter, orthocenter, etc
	 */
  void Center(TPZVec<REAL> &center);
};

inline TPZCompEl *TPZGeoEl1d::CreateCompEl(TPZCompMesh &mesh,int64_t &index){
  return fp(this,mesh,index);
}

inline void TPZGeoEl1d::SetCreateFunction(TPZCompEl *(*f)(TPZGeoEl1d *el,TPZCompMesh &mesh,int64_t &index)){
  fp = f;
}

#endif


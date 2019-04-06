#ifndef INTERFACECOMPELEMENTH
#define INTERFACECOMPELEMENTH

#include "pzcompel.h"
#include "pzerror.h"
#include "pzstack.h"
#include "TPZInterfaceEl.h"

#define MAXCONNECTS 5

struct TPZElementMatrix;
class TPZCompMesh;
class TPZInterpolatedElement;
class TPZConnect;
class TPZGeoEl;
class TPZGeoNode;


class TInterfaceElement : public TPZInterfaceElement {

  /** To store the computational elements and respective side of them related
      with this interface element.*/
  TPZCompElSide fElVec[2];
  /** Material of the interface. It can to be boundary condition*/
  TPZMaterial *fMaterial;

 protected:
  /** Whether different of zero, it compute the contribution in CalcStiff.
      The solution has discontinuity over the interface*/
  int fContribution;
  /** Dimension of the interface element. It is the dimension of the side of the
      connected elements from elvec. If dim = 0 (point_interface),
      dim =1 (linear_interface), dim = 2 (surface interface)*/
  int fDimension;
  /** Gets the connect related with this interface, respective to elside*/
  //  TPZConnect &GetConnect(int icon,TPZCompElSide &elside);
	/** Flux type number, to use at numerical flux */
	int fFluxType;

 public:
  /** Constructor and destructor*/
  TInterfaceElement(TPZCompMesh &mesh,TPZCompElSide &elsideleft,TPZCompElSide &elsideright,
                    int fluxtype, int64_t &index);
  ~TInterfaceElement();

  //ACCESS TO PRIVATE DATA
  /**return the type of the element
     The types are : ENoType, EOned, ETriangle, EQuadrilateral, EInterface, ESubstructure*/
  virtual MElementType Type();
	int IsInterface() { return 1; }
	
  virtual int Dimension() { return fDimension; }
  TPZMaterial *Material() { return fMaterial; }
  void SetMaterial(TPZMaterial* mat) { fMaterial = mat; }
	void SetFluxType(int fluxtype) { fFluxType = fluxtype; }

  TPZInterpolatedElement *RightEl() { return (TPZInterpolatedElement *)fElVec[1].Element(); }
  int RightSide() { return fElVec[1].Side(); }
  TPZInterpolatedElement *LeftEl() { return (TPZInterpolatedElement *)fElVec[0].Element(); }
  int LeftSide() { return fElVec[0].Side(); }

  /**The connects of the interface element are the side connects of the left element
     plus the side connects of the right element, independently whether some connects
     are duplicates*/
  int NConnects();
  /**Return number of side connects over one element of the fElVec*/
  int NConnects(int elnumber);
  /**return a reference to the ith connect*/
  TPZConnect &Connect(int i);
  /**Return a reference to the ith connect of the elnumber element of fElVec*/
  TPZConnect &Connect(int i,int elnumber);
	
	void SetConnectIndex(int i, int64_t connectindex);
	
  /**return the index of the ith connectivity of the element*/
	int64_t ConnectIndex(int i);

  /**active the contribution to the calcstiff*/
  void ActiveContribution() { fContribution = 1; }
  /**desactive the contribution to the calcstiff*/
  void NonContribution() { fContribution = 0; }

  void Print(std::ostream & out = std::cout);

  /**compute the integral over the interface if fContribution is not zero*/
  void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
  /**compute the integral over the interface if fContribution is not zero, only to Rhs*/
  void CalcRhs(TPZElementMatrix &ef);
  virtual void CalcResidual(TPZElementMatrix &ef);

};

inline void TInterfaceElement::SetConnectIndex(int i, int64_t connectindex) {
	std::cout << "TInterface::SetConnectIndex is called." << std::endl;
}

#endif



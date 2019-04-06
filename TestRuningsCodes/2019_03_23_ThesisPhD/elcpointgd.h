//HEADER FILE FOR CLASS TCompElPointGD

#ifndef COMPELPOINTGDH
#define COMPELPOINTGDH

#include <iostream>
#include <string.h>
#include "pzelcpoint.h"

struct 	TPZElementMatrix;
class  TPZConnect;
class  TPZBndCond;
class  TPZMat2dLin;
class  TPZGeoElPoint;
class  TPZCompMesh;
class  TPZGraphGrid;


class TCompElPointGD : public TPZCompElPoint {

  /** To identify whether the connect side is continuous(0) or discontinuous(1)*/
  int fConnectDisc;
  /** Index to internal connect, used if the element is discontinuous*/
  int fConnect;

 protected:
  /** Set fConnectDisc = 1 to indicate the side connect is discontinuous and
      set the correct positions to point to shape functions*/
  void SetConnectDiscontinuous(int icon);
  void SetConnectContinuous(int icon,TPZFMatrix<STATE> &values);
  /**To redimensioning the block of the icon connect that is being actived
     in continuous and copying the adequated values from discontinuous connect*/
  int StablizingConnectContinuous(int icon,TPZStack<TPZCompElSide> &elvec);

 public:

  static TPZCompEl *CreateElDisc(TPZGeoElPoint *gel,TPZCompMesh &mesh, int64_t &index);

  TCompElPointGD(TPZCompMesh &mesh,TPZGeoElPoint *ref, int64_t &index);
  ~TCompElPointGD();

  int CanBeDiscontinuous() { return 1; }

  /**To set and return the index of the connect i until to internal connect */
  void SetConnectIndex(int i, int64_t connectindex);
  int64_t ConnectIndex(int i);

  /**Returns the number of shapefunctions associated with a connect*/
  virtual int NConnectShapeF(int iconnect);

  /**return the number of connects along side iside*/
  virtual int NSideConnects(int iside);
  /**returns the local connect number of c along side*/
  virtual int SideConnectLocId(int c, int side);

  int NConnects() { return 2; }

  /** Set a connect icon continuous*/
  void MakeConnectContinuous(int icon=0);
  /** Set a connect icon discontinuous*/
  void MakeConnectDiscontinuous(int icon=0);
  /** Return 1 whether not exist neighbours or the connect in neighbours are discontinuous*/
  int NeighbourDiscontinuous(int disc=0);

  int IsConnectContinuous(const int i);

  void Print(std::ostream &out= std::cout);

};

/** The class TCompEl1dWI has the abstraction of the computational discontinuous
one-dimensional element enable to has interface over its common side boundary
with other (only one) discontinuous computational element with interface*/
class TCompElPointWI : public TCompElPointGD {

 public:
  static TPZCompEl *CreateElDiscWI(TPZGeoElPoint *gel,TPZCompMesh &mesh,int64_t &index);

  /** Constructors and destructor*/
  TCompElPointWI(TPZCompMesh &mesh,TPZGeoElPoint *ref, int64_t &index);
  ~TCompElPointWI();

  int CanToHaveInterface() { return 1; }

  /** Returns index of the interface over the side*/
  int Interface(int side=0);
  /** Set index of the interface over side */
  void SetInterface(int side, int64_t index);
  /** Creates one interface element over side of the element if it is possible*/
  int CreateInterface(int side,TPZCompMesh &mesh);

  void Print(std::ostream &out= std::cout);
  void DeleteInterfaces();
	
 protected:
  /** Vector of interface indexes */
	 int64_t fInterface;
};

#endif

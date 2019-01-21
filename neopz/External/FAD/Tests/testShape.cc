#include "pzreal.h"
#include "fadType.h"
#include "pzshapecube.h"
#include <iostream>
#include "error.h"

void error(char * string)
{
  if(string)cerr << endl << string << endl;

}

int main(){


  TPZVec<REAL> point(3,0.);
  TPZVec<int> id(8);
  int i;

  for(i = 0; i< 8; i ++)
  {
  	id[i] = i;
  }

  TPZVec<int> order(19);
  for(i = 0; i< 19; i ++)
  {
        order[i] = 3;
  }

  TPZVec<FADREAL> phi(64);
  TPZFMatrix OldPhi(64,1), OldDPhi(3,64);
  TPZFMatrix DiffPhi(64,1), DiffDPhi(3,64);

  TPZShapeCube::ShapeCube(point, id, order, phi);
  TPZShapeCube::ShapeCube(point, id, order, OldPhi, OldDPhi);

  cout << "Calculated by Fad" << phi;
  cout << "Old derivative method (phi)\n" << OldPhi;
  cout << "Old derivative method (dPhi)\n" << OldDPhi;

  shapeFAD::ExplodeDerivatives(phi, DiffPhi, DiffDPhi);
  DiffPhi-=OldPhi;
  DiffDPhi-=OldDPhi;
  cout << "FAD derivative method (phi)\n" << /*TPZFMatrix (OldPhi -*/ DiffPhi;
  cout << "FAD derivative method (dPhi)\n" <</* TPZFMatrix (OldDPhi -*/ DiffDPhi;
  return 0;
}

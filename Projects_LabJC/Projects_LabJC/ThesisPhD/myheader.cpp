#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "hadaptive.h"
#include "myheader.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzmanvector.h"
#include "pzelg1d.h"
#include "pzelgt2d.h"
#include "pzelgpoint.h"
#include "elcpointgd.h"
#include "interface.h"
#include "elc1dgd.h"
#include "elct2dgd.h"
#include "pzadmchunk.h"
#include "pzadmchunk.h"
#include "ttimeanalysis.h"
#include "pzgraphmesh.h"
#include "pzdxmesh.h"
#include "pzbndcond.h"

/** Global counter of the pause at the program */
int npausas = 0;

std::string plotfile("plot.vtk");

void GetCommentary(std::istream &input,int nlines) {
  char lixo[256];
  for(int i=0;i<nlines;i++) {
    input >> lixo[0];
    input.getline(lixo,256);
  }
}
void GetDataCommented(std::istream &input,int &value) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  input >> value;
}
void GetDataCommented(std::istream &input,int const nlines,int &value) {
  char lixo[256];
  input >> lixo[0];
  for(int j=0;j<nlines;j++) input.getline(lixo,256);
  input >> value;
}
void GetDataCommented(std::istream &input,REAL &value) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  input >> value;
}
void GetDataCommented(std::istream &input,TPZVec<REAL> &vector) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  int i, n = vector.NElements();
  for(i=0;i<n;i++) input >> vector[i];
}
void GetDataCommented(std::istream &input,char *string,int size) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  input >> string[0];
  string[1] = '\0';
  input.getline(lixo,256);
  lixo[size - 1] = '\0';
  strncat(string,lixo,size-1);
}
void GetDataCommented(std::istream &input, std::string &str, int size) {
	char lixo[256];
	input >> lixo[0];
	input.getline(lixo, 256);
	input >> lixo[0];
	lixo[1] = '\0';
	str = lixo;
	input.getline(lixo, 256);
	lixo[size-1] = '\0';
	str += lixo;
}


char Itoa(int num) {
  switch(num) {
  case 0: return '0';
  case 1: return '1';
  case 2: return '2';
  case 3: return '3';
  case 4: return '4';
  case 5: return '5';
  case 6: return '6';
  case 7: return '7';
  case 8: return '8';
  case 9: return '9';
  default: return 'u';
  }
}

/** Imprime o numero de elementos interpolados, interfaces e interpolados de 
    dimensao menor a dimensao do lei*/
void CountElements(TPZCompMesh &cmesh, std::ostream &out) {
  long nelem, i, dim, numel = 0, numeldimmenor = 0, numinterfaces = 0;
  /**Depending on the material dimension choose the maxime number of sub-elements*/
  TPZMaterial *mat = 0;
  std::map<int, TPZMaterial * >::const_iterator it;
  for (it = cmesh.MaterialVec().begin(); it != cmesh.MaterialVec().end(); it++) {
	  if (it->second->Id() > -1) break;
  }
  mat = it->second;
  if(!mat) {
    PZError << "CountElements. Not found materials.\n";
    exit(1);
  }
	dim = mat->Dimension();
  /**Inicia coarsing*/
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  TPZCompEl *cel;
  nelem = cmesh.NElements();
  for(i=0;i<nelem;i++) {
    cel = elvec[i];
		if(!cel) continue;
		if(!cel->IsInterface()) {
		  if(cel->Dimension()==dim) numel++;
			else numeldimmenor++;
		}
		else numinterfaces++;
	}
	out << numel << "   " << numeldimmenor << "   " << numinterfaces << std::endl;
}

/**To refine the complete computational mesh*/
void Refinements(TPZCompMesh &cmesh,int ncycles,int MaxLevelRefinements) {
  if(!ncycles) return;
  int j,cycle;
  int nelem,ninterp;
  TPZManVector<int64_t> subindex;
  TPZCompEl *el;
  TPZVec<int> vec;
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  for(cycle=0;cycle<ncycles;cycle++) {
    ninterp = 0;
    nelem = cmesh.NElements();
    vec.Resize(nelem);
    for(j=0;j<nelem;j++) {
      el = elvec[j];
      if(!el || el->IsInterface() || !el->Dimension()) continue;
      if(el->Reference()->Level() >= MaxLevelRefinements) continue;
      vec[ninterp++] = j;
    }
    for(j=0;j<ninterp;j++)
      cmesh.Divide(vec[j],subindex,1);
  }
//  cmesh.CleanInterfaces();
  cmesh.InitializeBlock();
}

void Refinements(TPZCompMesh &cmesh,int ncycles,TPZVec<int> &typecel,
        int interpolated,int MaxLevelRefinements) {
  int j,cycle;
	if(!ncycles) return;
  int nelem,ninterp;
  TPZManVector<int64_t> subindex;
  TPZCompEl *el;
  TPZVec<int> vec;
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  for(cycle=0;cycle<ncycles;cycle++) {
    if(!cycle) {
      switch(typecel.NElements()) {
      case 3:
        if(typecel[2]==1) TPZGeoElT2d::SetCreateFunction(&(TCompElT2dGD::CreateElDisc));
        else if(typecel[2]==2) TPZGeoElT2d::SetCreateFunction(&(TCompElT2dWI::CreateElDiscWI));
      case 2:
        if(typecel[0]==1) TPZGeoElPoint::SetCreateFunction(&(TCompElPointGD::CreateElDisc));
        else if(typecel[0]==2) TPZGeoElPoint::SetCreateFunction(&(TCompElPointWI::CreateElDiscWI));
        if(typecel[1]==1) TPZGeoEl1d::SetCreateFunction(&(TCompEl1dGD::CreateElDisc));
        else if(typecel[1]==2) TPZGeoEl1d::SetCreateFunction(&(TCompEl1dWI::CreateElDiscWI));
        break;
      default:
        PZError << "Refinements. Elements continuous.\n";
        break;
      }
    }
    ninterp = 0;
    nelem = cmesh.NElements();
    if(!nelem) continue;
    vec.Resize(nelem);
    for(j=0;j<nelem;j++) {
	  vec[j] = -1;
      el = elvec[j];
      if(!el || el->IsInterface() || !el->Dimension()) continue;
      if(el->Reference()->Level() >= MaxLevelRefinements) continue;
      vec[ninterp++] = j;
    }
    for(j=0;j<ninterp;j++) {
      cmesh.Divide(vec[j],subindex,interpolated);
	elvec[vec[j]] = 0;
	}
//    cmesh.CleanInterfaces();		
    cmesh.InitializeBlock();
  }
}

void InvJacob2d(TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &jacinv) {
	REAL det = 1./(jacinv(0,0)*jacinv(1,1));
	if(jacinv(0,0)<0.) det *= -1.;
	REAL coef1 = 1./(fabs(jacinv(1,1))*det), coef2 = -jacinv(0,1);

	jacinv(0,0) = coef1*axes(1,1)+coef2*axes(0,1);
	jacinv(0,1) = -(coef1*axes(1,0)+coef2*axes(0,0));
	jacinv(1,0) = -(jacinv(1,1)*axes(0,1));
	jacinv(1,1) = jacinv(1,1)*axes(0,0);
}

REAL MeanSolutionFromSubElements(TPZGeoEl *gel,int var) {
#ifndef NOTDEBUG
  if(!gel || !gel->HasSubElement())
    PZError << "MeanSolution(TPZGeoEl *).  Bad geometrical element.\n";
#endif
  TPZGeoEl *gelsub;
  int i, nsub = gel->NSubElements();
  REAL mean = 0.;
  for(i=0;i<nsub;i++) {
    gelsub = gel->SubElement(i);
    if(!gelsub->Reference()) mean += MeanSolutionFromSubElements(gelsub,var);
    else mean += ((TPZInterpolatedElement *)gelsub->Reference())->MeanSolution(var);
  }
  return mean/nsub;
}

REAL LinearApproximation(TPZVec<REAL> &qsi,TPZInterpolatedElement *el,int var) {
  TPZBlock<STATE> &block = el->Mesh()->Block();
  int nshape = el->NShapeF();
  int dim = el->Dimension();
  TPZFMatrix<STATE> phi(nshape,1), dphi(dim,nshape);
  TPZConnect *con;
  int ncorners = el->NCornerConnects();
  int ncon = el->NConnects();
  int numdof = el->Material()->NStateVariables();

  int i, iv = 0, mask = 0;
  int seqnum;
  REAL sol = 0.;
  el->Shape(qsi,phi,dphi);
  for(i=0;i<ncorners;i++) {
	  TPZCompElSide cels(el, i);
	  int64_t indexconn = cels.ConnectIndex();
    if(indexconn < 0) {
      mask++;
      continue;
    }
    con = &(el->Connect(i));
    seqnum = con->SequenceNumber();
    sol += (phi(iv,0)*block(seqnum,0,var,0));
    iv++;
  }
  if(!mask) return sol;
  con = &(el->Connect(ncon-1));
  seqnum = con->SequenceNumber();
  for(i=0;i<mask;i++) {
    sol += (phi(iv,0)*block(seqnum,0,i*numdof+var,0));
    iv++;
  }
  return sol;
}

void SetCreateFunctionToElement(int dim,int type) {
  switch(dim) {
  case 0: {
    if(type==1)
      TPZGeoElPoint::SetCreateFunction(&(TCompElPointGD::CreateElDisc));
    else if(type==2)
      TPZGeoElPoint::SetCreateFunction(&(TCompElPointWI::CreateElDiscWI));
    break;
  }
  case 1: {
    if(type==1)
      TPZGeoEl1d::SetCreateFunction(&(TCompEl1dGD::CreateElDisc));
    else if(type==2)
      TPZGeoEl1d::SetCreateFunction(&(TCompEl1dWI::CreateElDiscWI));
    break;
  }
  case 2: {
    if(type==1)
      TPZGeoElT2d::SetCreateFunction(&(TCompElT2dGD::CreateElDisc));
    else if(type==2)
      TPZGeoElT2d::SetCreateFunction(&(TCompElT2dWI::CreateElDiscWI));
    break;
  }
  default:
    PZError << "SetCreateFunctionToElement has bad dimension information.\n";
    break;
  }
}

/**Return band width to stiffness matrix to discontinuous elements*/
int BandWidth(TPZCompMesh *mesh) {
 int i,j, bw = 0;
 int nelem, nshape;
 TPZAdmChunkVector<TPZCompEl *> &elvec = mesh->ElementVec();
 nelem = elvec.NElements();
 for(i=0;i<nelem;i++) {
   TPZCompEl *el = elvec[i];
   if(!el || el->IsInterface()) continue;
   TPZInterpolatedElement *cel = (TPZInterpolatedElement *)el;
   int ncon = cel->Reference()->NSides();
   for (j = 0; j < ncon; j++) {
	   TPZCompElSide cels(cel, j);
	   int64_t connind = cels.ConnectIndex();
	   if (!(connind < 0))
		   return mesh->BandWidth();
   }
   nshape = cel->NShapeF()-1;
   if(bw<nshape) bw = nshape;
 }
 return bw;
}

   /**BCConnects - To apply bundary conditions*/
/*
   int index = cmesh1->BCConnectVec().AllocateNewElement();
   TPZConnectBC bccon(&(cmesh1->ElementVec()[0]->Connect(0)),bc);
   cmesh1->BCConnectVec()[index] = bccon;
   TPZConnectBC bccon1(&(cmesh1->ElementVec()[1]->Connect(1)),bc);
   index = cmesh1->BCConnectVec().AllocateNewElement();
   cmesh1->BCConnectVec()[index] = bccon1;
*/
   /**Turn the bcconnects as continuous connects*/
//  TCompEl1dWI *celgd = (TCompEl1dWI *) cmesh1->ElementVec()[0];
//  celgd->MakeConnectContinuous(0);
//  celgd = (TCompEl1dWI *) cmesh1->ElementVec()[1];
//  celgd->MakeConnectContinuous(0);

   /**Testing subdivision
  TPZManVector<int> subindex0(2);
  TPZManVector<int> subindex1(2);
  subindex0[0] = 0;
  subindex0[1] = 1;
  for(int k=0;k<30;k++) {
    cmesh1->Divide(subindex0[0],subindex);
    int sub2 = subindex[1];
    cmesh1->Divide(subindex[0],subindex);
  cmesh1->Divide(sub2,subindex);
  celgd->MakeConnectContinuous(1);
  celgd = (TCompEl1dGD *) cmesh1->ElementVec()[subindex[1]];
  celgd->MakeConnectContinuous(0);

  celgd->MakeConnectDiscontinuous(1);
  cmesh1->Divide(2,subindex);
  celgd = (TCompEl1dGD*) cmesh1->ElementVec()[subindex[1]];
  cmesh1->Divide(1,subindex);
  sub2 = subindex[1];
  cmesh1->Divide(subindex[0],subindex);
  cmesh1->Divide(sub2,subindex);
  cmesh1->Divide(subindex[0],subindex);
  cmesh1->Divide(subindex[0],subindex);
  celgd->MakeConnectContinuous(0);

  DivideAny(*cmesh1);

   /**Testing coarsen
   CoarsenAny(*cmesh1);
  /** Printing computational information
//   cmesh1->Print(outm1);
//   outm1.flush();
*/

void Pause(char *data,int nelem) {
  int n;
 std::cout << "Pause = " << npausas++;
 std::cout << "\tData ";
  if(nelem)
   std::cout << data << "\n";
 std::cout.flush();
 std::cin >> n;
  if(!n) exit(1);
}

int Pause(int data) {
  int n;
 std::cout << "Pause = " << npausas++;
 std::cout << "\tData ";
 std::cout << data << "\n";
 std::cout.flush();
 std::cin >> n;
  return n;
}

void Pause(float data) {
  int n;
 std::cout << "Pause = " << npausas++;
 std::cout << "\tData ";
 std::cout << data << "\n";
 std::cout.flush();
 std::cin >> n;
  if(!n) exit(1);
}

void Pause(double data) {
  int n;
 std::cout << "Pause = " << npausas++;
 std::cout << "\tData ";
 std::cout << data << "\n";
 std::cout.flush();
 std::cin >> n;
  if(!n) exit(1);
}

void Multiply(TPZMatrix<STATE> &A,TPZMatrix<STATE> &B,TPZMatrix<STATE> &product) {
  int rows = A.Rows(), cols = B.Cols();
  int n = A.Cols();
#ifndef NOTDEBUG
  if(n!=B.Rows()) {
    PZError << "myheader::Multiply. A.Cols and B.Rows are incompatibles.\n";
    return;
  }
#endif
  int i,j,k;
  for(i=0;i<rows;i++) {
    for(j=0;j<cols;j++) {
      product(i,j)=0.;
      for(k=0;k<n;k++)
        product(i,j) += A(i,k)*B(k,j);
    }
  }
}

/**Print information of the interfaces into comp. elements and the
   computational elements related with each interface*/
void PrintInterfaces(TPZCompMesh *cmesh, std::ostream &out) {
  int i, j, nelem = cmesh->NElements();
	int countcel = 0, countinterf = 0;
  TPZCompEl *cel;
  TInterfaceElement *inter_face;
  for(i=0;i<nelem;i++) {
    cel = cmesh->ElementVec()[i];
    if(!cel) continue;
    if(cel->IsInterface()) continue;
    out << "COMPEL " << i << "  Interf: ";
		countcel++;
    for(j=0;j<cel->Reference()->NSides();j++) out << ((TPZInterpolatedElement *)cel)->Interface(j) << "  ";
    out << "\n";
  }
	out << "Total COMPELs " << countcel << std::endl;
  for(i=0;i<nelem;i++) {
    cel = cmesh->ElementVec()[i];
    if(!cel || !cel->IsInterface()) continue;
    inter_face = (TInterfaceElement *)cel;
		countinterf++;
    out << "INTERF " << i << "  Left " << inter_face->LeftEl()->Index();
    out << "  Right " << inter_face->RightEl()->Index() << std::endl;
  }
	out << "Total interfaces " << countinterf << std::endl;
}

/**Print vector elements*/
void PrintVector(TPZVec<int> &vec, std::ostream &out,char *vecname) {
  int i, nelem = vec.NElements();
  out <<std::endl << vecname <<std::endl;
  for(i=0;i<nelem;i++) {
    out << vec[i];
    if((i+1)%5) out << "\t";
    else out <<std::endl;
  }
	out <<std::endl;
}
/** Print the evaluated erros for several computations with different runs */
void PrintErrors(TPZVec<REAL> &errvec,int ngrouped,char *name) {
  char filename[32];
  strncpy(filename,name,9);
  filename[9] = '\0';
  strncat(filename,".dat",5);
  std::ofstream outerr(filename);
  char *nameerr[] = {"L1Error","L2Error","NormInfinity"};
  int i, j, nerros = errvec.NElements()/3;
  double value;
  outerr << "Numero de Erros = " << nerros << std::endl;
  outerr << "Numero de Erros agrupados = " << ngrouped <<std::endl;
  for(i=0;i<3;i++) {
    outerr <<std::endl << nameerr[i] << " = {";
    for(j=0;j<nerros;j++) {
      value = errvec[3*j+i];
      if(!(j%ngrouped)) outerr << "{";
      if(value < 1.) outerr << value;
      else outerr << 1.;
      if((j+1)%ngrouped) outerr << ",";
      else {
        outerr << "}";
        if(j!=nerros-1) outerr << ",";
				outerr <<std::endl;
      }
    }
    for(;j%ngrouped;j++) {
      if((j+1)%ngrouped) outerr << 1. << ",";
      else outerr << 1. << "}";
    }
    outerr << "}" <<std::endl;
  }
	time_t t = time(0);
  outerr <<std::endl << "Current time " << ctime(&t) <<std::endl;
  outerr.close();
}
/** Print the evaluated erros for several computations with different runs */
void PrintErrorsDX(TPZVec<REAL> &errvec,int ngrouped,std::string &name) {
  char filename[32];
  strncpy(filename,name.c_str(),10);
  filename[10] = '\0';
  strncat(filename,".dax",5);
  std::ofstream outerr(filename);
  int i, j, nerros = errvec.NElements()/3;
  double value;
  outerr << nerros <<std::endl;
  outerr << ngrouped <<std::endl;
  for(i=0;i<3;i++) {
    for(j=0;j<nerros;j++) {
      value = errvec[3*j+i];
      if(value < 1.) outerr << value;
      else outerr << 1.;
      if((j+1)%ngrouped) outerr << "   ";
      else outerr <<std::endl;
    }
    for(;j%ngrouped;j++) {
      if((j+1)%ngrouped) outerr << 1. << "   ";
      else outerr << 1. <<std::endl;
    }
    outerr <<std::endl;
  }
	time_t t = time(0);
  outerr <<std::endl << "Current time " << ctime(&t) <<std::endl;
  outerr.close();
}
/** Recupera erros computados em rodadas anteriores*/
void GetOldErros(TPZVec<REAL> &errovec,int &erroini,char *ex_out) {
	std::ifstream input(ex_out);
	int i, n = -1, m = 0, nchars = 3;
	if(!input.rdbuf()->is_open()) return;
	char lixo[128];
	input.getline(lixo,128);
	while(!input.eof()) {
	  for(i=0;i<nchars;i++) input >> lixo[i];
  	while(strncmp(lixo,"Run",nchars)) {
		  input.getline(lixo,128);
			for(i=0;i<nchars;i++) input >> lixo[i];
	    if(input.eof()) {
			  input.close();
			  erroini = n+1;
			  return;
			}
		}
	  input >> n;
		input.getline(lixo,128);
		while(strncmp(lixo,"Estimated error",12))
		  input.getline(lixo,128);
		for(i=0;i<3;i++) {
      input.getline(lixo,128,'=');
			input >> errovec[m++];
		}
		input.getline(lixo,128);
	}
	erroini = n+1;
	input.close();
}

/**Imprime os valores da diagonal de matrix para verificacao desde
   matrix(start,start) ateh matrix(end-1,end-1) */
int PrintDiagonal(TPZMatrix<STATE> *matrix, int mask) {
  int valpause = Pause(mask);
  if(valpause) {
    int cols = matrix->Cols();
	int start = 10, end = 50;
//   std::cout << setprecision(3);
    for(int r=start;r<end;r++) {
	  if(cols>r)std::cout << (*matrix)(r,r);
	  else std::cout << (*matrix)(r,0) << "*";
	  if(r%10)std::cout << "\t";
	  else std::cout <<std::endl;
	}
   std::cout <<std::endl;
//   std::cout << setprecision(8);
  }
  return valpause;
}

/** To asign name to file for post-processing - name contains of the current law name */
void PlotFileName(char *name,int num,char *plotfile) {
  int index = 0, aux = num%100;
  char nameaux[256];
  nameaux[index++] = Itoa(num/100);
  nameaux[index++] = Itoa(aux/10);
  nameaux[index++] = Itoa(aux%10);
  nameaux[index++] = '.';
  nameaux[index++] = 'v';
  nameaux[index++] = 't';
  nameaux[index++] = 'k';
  nameaux[index++] = '\0';
  strncpy(plotfile,name,16);
  strncat(plotfile,nameaux,9);
}
void PlotFileName(std::string name, int num, std::string plotfile) {
	int index = 0, aux = num % 100;
	char nameaux[256];
	nameaux[index++] = Itoa(num / 100);
	nameaux[index++] = Itoa(aux / 10);
	nameaux[index++] = Itoa(aux % 10);
	nameaux[index++] = '.';
	nameaux[index++] = 'v';
	nameaux[index++] = 't';
	nameaux[index++] = 'k';
	nameaux[index++] = '\0';
	plotfile = name;
	plotfile += nameaux;
}

/** Determine whether the x point is into, at left or at right of the hyperplane
    a*x + b*y + c*z + d = 0   */
REAL Term_c = 0.;
REAL SecondTerm_c = -1.e12;
//REAL SecondTerm_c = .6;
REAL Coef[3] = {1., 0., 0.};
REAL SecondCoef[3] = {1., 0., 0.};
REAL u0[16] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0.,0.,0.,0.,0.,0.};
void PutDiscontinuousHyperplane(REAL term_c, REAL coef_x, REAL coef_y, REAL coef_z) {
  Coef[0] = coef_x;
  Coef[1] = coef_y;
  Coef[2] = coef_z;
  Term_c = term_c;
}
void PutSecondDiscontinuousHyperplane(REAL term_c, REAL coef_x, REAL coef_y, REAL coef_z) {
  SecondCoef[0] = coef_x;
  SecondCoef[1] = coef_y;
  SecondCoef[2] = coef_z;
  SecondTerm_c = term_c;
}
int DiscontinuousHyperplane(TPZVec<REAL> &x) {
  REAL r = Term_c;
  int dim = x.NElements();
#ifndef NOTDEBUG
  if(dim<1 || dim>3) {
    PZError << "myheader::DiscontinuousLine. Bad dimension of the point x.\n";
    return 0;
  }
#endif
  for(int i=0;i<dim;i++)
    r += Coef[i] * x[i];
  if(IsZero(r)) return 0;
  else if(r<0) return -1;
  return 1;
}
int SecondDiscontinuousHyperplane(TPZVec<REAL> &x) {
  REAL r = SecondTerm_c;
  int dim = x.NElements();
  for(int i=0;i<dim;i++)
    r += SecondCoef[i] * x[i];
  if(IsZero(r)) return 0;
  else if(r<0) return -1;
  return 1;
}

void UZero(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel) {
#ifndef NOTDEBUG
  if(x.NElements()!=3) PZError << "UZero. Bad dimension of the point x.\n";
#endif
  int ncon, disc = DiscontinuousHyperplane(x);
  /** If disc == 0 then x is into the discontinuous line */
  if(!disc) {
    ncon = cel->NCornerConnects();
    TPZVec<REAL> mastercenter(3,0.);
    TPZVec<REAL> center(3,0.);
    if(ncon == 3) mastercenter[0] = mastercenter[1] = 1./3.;
    cel->Reference()->X(mastercenter,center);
    disc = DiscontinuousHyperplane(center);
  }
  if(disc > 0) {
    disc = SecondDiscontinuousHyperplane(x);
    if(!disc) {
      ncon = cel->NCornerConnects();
      TPZVec<REAL> mastercenter(3,0.);
      TPZVec<REAL> center(3,0.);
      if(ncon==3) mastercenter[0] = mastercenter[1] = 1./3.;
      cel->Reference()->X(mastercenter,center);
      disc = SecondDiscontinuousHyperplane(center);
    }
    if(disc>0) UZeroRight_2(x,u);
    else UZeroRight(x,u);
    return;
  }
  else UZeroLeft(x,u);
}
    
void UZero(TPZVec<REAL> &x,TPZVec<REAL> &u) {
  int disc = DiscontinuousHyperplane(x);
  /** If disc == 0 then x is into the discontinuous line */
  if(!disc) PZError << "UZero. Point x is on the first discontinuous line.\n";
  if(disc > 0) {
    disc = SecondDiscontinuousHyperplane(x);
    if(!disc) PZError << "UZero. Point x is on the second discontinuous line.\n";
    if(disc>0) UZeroRight_2(x,u);
    else UZeroRight(x,u);
    return;
  }
  else UZeroLeft(x,u);
}

void UZeroLeft(TPZVec<REAL> &x,TPZVec<REAL> &u) {
  int i, nvar = u.NElements();
  for(i=0;i<nvar;i++)
    u[i] = u0[i];
}
void UZeroRight(TPZVec<REAL> &x,TPZVec<REAL> &u) {
  int i, nvar = u.NElements();
  for(i=0;i<nvar;i++)
    u[i] = u0[nvar+i];
}
void UZeroRight_2(TPZVec<REAL> &x,TPZVec<REAL> &u) {
  int i, nvar = u.NElements();
  for(i=0;i<nvar;i++)
    u[i] = u0[2*nvar+i];
}

void UZeroNull(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel) {
  u[0] = 0.;
}
void UZeroUnity(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel) {
  u[0] = 1.;
}

void UZeroSin(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel) {
  u[0] = sin(2*PI_Value*x[0]);
}
void UZeroTwo(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel) {
  if(x[0]<.25) u[0] = 1.;
	else if(x[0] > .5) u[0] = .5;
	else u[0] = -2.*x[0] + 1.5;
}

void UZeroSin2d(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel) {
  u[0] = sin(2*PI_Value*(x[0]+x[1]));
}

void UZeroLinearCircular(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel) {
  if(x[0]<0. && x[1]>0.)
    u[0] = sin(-2.*PI_Value*x[0])*sin(2.*PI_Value*x[1]);
  else
    u[0] = 0.;
}
void UZeroCone(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZInterpolatedElement *cel){
//  double xfocus = 2., yfocus = 1.5, radio = 1.5, h = 2.5;
  double xfocus = 0., yfocus = 2.5, radio = 1.5, h = 2.5;
//  double xfocus = -1.25*sqrt(2), yfocus = 1.25*sqrt(2), radio = 1.5, h = 2.5;
  double r1 = x[0]-xfocus, r2 = x[1]-yfocus;
  double distance = sqrt(r1*r1 + r2*r2);
  if(distance < radio)
    u[0] = h*(1 + cos(PI_Value*distance/radio));
  else u[0] = 0.;
}

int PutUZero(int nuzero,TTimeAnalysis *an, std::istream &input,int dim) {
  if(!nuzero) {
    PZError << "UZero(x) is undefined. The program is killed.\n";
    exit(0);
  }
  /**Applying initial function in T = 0 */
  int i, discontinuity;
  GetDataCommented(input,2,discontinuity);   // Is initial function discontinuous?
  if(discontinuity) {   // Fill the coefficients of the discontinuous hyperplanes.
    /**First discontinuous hyperplane*/
    input >> Term_c;
    for(i=0;i<dim;i++)
      input >> Coef[i];
    int order = an->OrderLaw();
    for(i=0;i<order;i++) input >> u0[order+i];
    if(discontinuity>1) {
      input >> SecondTerm_c;
      for(i=0;i<dim;i++)
        input >> SecondCoef[i];
      for(i=0;i<order;i++) input >> u0[2*order+i];
    }
    an->SetUZero(UZero);
    return 1;
  }
  if(nuzero==1) input >> nuzero;
  switch(nuzero) {
  case 0:
   std::cout << "UZero(x) := 0" <<std::endl;
    an->SetUZero(UZeroNull);
    break;
  case 1:
   std::cout << "UZero(x) := sin(2*pi*x).\n";
    an->SetUZero(UZeroSin);
    break;
  case 2:
   std::cout << "UZero(x) is continuous with two constant values.\n";
    an->SetUZero(UZeroTwo);
    break;
  case 3:
   std::cout << "UZero(x,y) := sin(2*pi*(x+y)).\n";
    an->SetUZero(UZeroSin2d);
    break;
  case 4:
   std::cout << "UZero(x,y) := sin(2*pi*x)*sin(2*pi*y). Displacement circular.\n";
    an->SetUZero(UZeroLinearCircular);
    break;
  case 5:
   std::cout << "UZero(x,y) := 2.5(1+cos(PI*r/1.5)) if r <= 1.5  caso contrario 0\n";
    an->SetUZero(UZeroCone);
    break;
  case 10:
   std::cout << "UZero(x) := 0" <<std::endl;
    an->SetUZero(UZeroNull);
    break;
  case 11:
   std::cout << "UZero(x,y) := 1" <<std::endl;
    an->SetUZero(UZeroUnity);
    break;
  default:
    PZError << "UZero(x) is undefined. Write nuzero (Quit -> 0) = ";
   std::cin >> nuzero;
    PutUZero(nuzero,an,input);
    return 1;
  }
  return 1;
}

/**Count which different Dirichlet conditions has the boundary conditions */
int NumberOfDirichletConditions(TPZCompMesh &cmesh){
	int n = 0;
	TPZBndCond *bc;
	std::map<int, TPZMaterial * >::const_iterator it;
	for (it = cmesh.MaterialVec().begin(); it != cmesh.MaterialVec().end(); it++) {
		if (it->second->Id() > -1) continue;
		bc = (TPZBndCond *)it->second;
		if (!bc->Type()) n++;
	}

  return n;
}

/**Detect the Dirichlet values from boundary conditions into the cmesh */
REAL DetectDirichletValue(TPZCompMesh *cmesh, int nvar) {
  TPZBndCond *bc;
  std::map<int, TPZMaterial * >::const_iterator it;
  for (it = cmesh->MaterialVec().begin(); it != cmesh->MaterialVec().end(); it++) {
	  if (it->second->Id() > -1) continue;
	  bc = (TPZBndCond *)it->second;
	  if (bc->Type()) continue;
	  return bc->Val2()(nvar, 0);
  }

  PZError << "DetectDirichletValue. Not found Dirichlet condition to var " << nvar << ".\n";
  return 1.;
}

/**Turn continuous all the connects into the cmesh, excepts the connects from
boundary elements */
int AllConnectContinuous(TPZCompMesh *cmesh) {
  int nelem = cmesh->NElements();
  int i, j, dim, icon, k=0;
  std::map<int, TPZMaterial * >::const_iterator it;
  for (it = cmesh->MaterialVec().begin(); it != cmesh->MaterialVec().end(); it++) {
	  if (it->second->Id() < 0) continue;
	  dim = it->second->Dimension();
	  break;
  }

  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh->ElementVec();
  TPZCompEl *el;
  TPZInterpolatedElement *cel;
  //Apagando os elementos pontoais se estamos em dimensao 2 ou superior
  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el) continue;
//	if(!el->IsInterpolated()) {
//		if(((TInterfaceElement *)el)->LeftEl()->Dimension() == ((TInterfaceElement *)el)->RightEl()->Dimension())
//			((TInterfaceElement *)el)->NonContribution();
//	}
	if(el->Dimension()==0 && dim == 2) {
		delete el;
		elvec[i] = 0;
	}
  }
  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el || el->IsInterface()|| el->Dimension()!=dim) continue;
    cel = (TPZInterpolatedElement *)el;
    icon = cel->NCornerConnects();
    for(j=0;j<icon;j++) {
//      cel->MakeConnectContinuous(j);
      k++;
    }
  }
  cmesh->InitializeBlock();

	for(i=0;i<nelem;i++) {
	  el = elvec[i];
		if(!el || el->IsInterface()) continue;
		if(el->Dimension() == dim - 1) {
			cel = (TPZInterpolatedElement *)el;
			icon = cel->NCornerConnects();
			for(j=0;j<icon;j++) {
//			  cel->MakeConnectDiscontinuous(j);
			  k--;
			}
		}
	}
  cmesh->InitializeBlock();
  return k;
}

/**Turn discontinuous all the connects into the cmesh, excepts the connects from
boundary elements */
int AllConnectDiscontinuous(TPZCompMesh *cmesh) {
  int nelem = cmesh->NElements();
  int i, j, dim, icon, k=0;
  std::map<int, TPZMaterial * >::const_iterator it;
  for (it = cmesh->MaterialVec().begin(); it != cmesh->MaterialVec().end(); it++) {
	  if (it->second->Id() < 0) continue;
	  dim = it->second->Dimension();
	  break;
  }

  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh->ElementVec();
  TPZCompEl *el;
  TPZInterpolatedElement *cel;

  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el || el->IsInterface()) continue;
    cel = (TPZInterpolatedElement *)el;
    icon = cel->Reference()->NSides();
    for(j=0;j<icon;j++) {
//      cel->MakeConnectDiscontinuous(j);
      k++;
    }
  }
  return k;
}

/**Create a *.xls file with excell format, to data visualization */
void CreateXLS(char *name,int nobjects,int ncoords) {
  if(ncoords<1 || ncoords>3) ncoords = 2;
  std::ifstream input(name);
  int npoints;
  REAL value;
  int i,j=0;
  for(i=0;i<16;i++)
    if(name[i]=='.') {
      name[i+1] = '\0';
      break;
    }
  strcat(name,"xls");
  std::ofstream output(name);
  char lixo[256];
  while(!j) {
    input >> lixo;
    if(!strncmp(lixo,"items",5)) j=1;
  }
  j=0;
  input >> npoints;
  input.getline(lixo,256);
  output << "Coordinates" <<std::endl;
  switch(ncoords) {
  case 1:
  {
    REAL y,z;
    for(i=0;i<npoints;i++) {
      input >> value >> y >> z;
      output << value << '\t';
    }
    break;
  }
  case 2:
  {
    REAL *y = new REAL[npoints];
    REAL z;
    for(i=0;i<npoints;i++) {
      input >> value >> y[i] >> z;
      output << value << '\t';
    }
    output <<std::endl <<std::endl;
    for(i=0;i<npoints;i++) output << y[i] << '\t';
    delete[] y;
  }
  case 3:
  {
    REAL *y = new REAL[npoints];
    REAL *z = new REAL[npoints];
    for(i=0;i<npoints;i++) {
      input >> value >> y[i] >> z[i];
      output << value << '\t';
    }
    output <<std::endl <<std::endl;
    for(i=0;i<npoints;i++) output << y[i] << '\t';
    output <<std::endl <<std::endl;
    for(i=0;i<npoints;i++) output << z[i] << '\t';
    delete[] y;
    delete[] z;
  }
  break;
  }
  output <<std::endl <<std::endl;
  GetCommentary(input,3);
  for(i=0;i<nobjects;i++) {
    output << "object " << i <<std::endl;
    while(!j) {
      input.getline(lixo,256);
      if(!strncmp(lixo,"#",1)) {
        GetCommentary(input);
        j = 1;
      }
    }
    for(j=0;j<npoints;j++) {
      input >> value;
      output << value << '\t';
    }
    output <<std::endl <<std::endl;
    j = 0;
    GetCommentary(input,2);
    while(!j) {
      input.getline(lixo,256);
      if(!strncmp(lixo,"end",3)) {
        input.close();
        output.close();
        return;
      }
      if(!strncmp(lixo,"#",1)) {
        GetCommentary(input);
        j = 1;
      }
    }
  }
  input.close();
  output.close();
}

int NoZeroNumbers(TPZVec<int> &vec) {
  int i,nnozeros=0,nelem = vec.NElements();
  for(i=0;i<nelem;i++) if(vec[i]!=0) nnozeros++;
  return nnozeros;
}

int DifferentNumbers(TPZVec<int> &vec) {
  int i,j,ndif,nelem = vec.NElements();
  ndif = nelem;
  for(i=0;i<nelem;i++) {
    for(j=i+1;j<nelem;j++) {
      if(vec[i]==vec[j]) {
	ndif--;
	break;
      }
    }
  }
  return ndif;
}


/*******       File : NLinLaw1.c       *******/

#include "buck1d.h"

/*******       TBuckleyLeverett       *******/

TBuckleyLeverett::TBuckleyLeverett(std::istream &input) : TConservationLaw(input) {
	//if(fDim!=1) std::cout << "ERROR : TBuckleyLeverett, bad spacial domain dimension." << std::endl;
   //if(fOrder!=1) std::cout << "ERROR : TBuckleyLeverett, bad system order." << std::endl;
	fSonicPoint =0.;
	std::cout << "\n Coefficient, fa = ";
	std::cin >> fa;
}
/*TBuckleyLeverett::TBuckleyLeverett(char *file,double Time):TConservationLaw(1,file,Time,1) {
	fSonicPoint =0.;
	std::cout << "\n Coefficient, fa = ";
	std::cin >> fa;
}
TBuckleyLeverett::TBuckleyLeverett(char *file) : TConservationLaw(file) {
	if(fDim!=1) std::cout << "ERROR : TBuckleyLeverett, bad spacial domain dimension." << std::endl;
   if(fOrder!=1) std::cout << "ERROR : TBuckleyLeverett, bad system order." << std::endl;
   std::cout << "\n Coefficient, fa = ";
   std::cin >> fa;
	fSonicPoint =0.;
}
TBuckleyLeverett::TBuckleyLeverett(TCompMesh *c,double Time):TConservationLaw(c,Time,1) {
	if(fDim!=1) cout << "ERROR : TBuckleyLeverett, bad spacial domain dimension." << endl;
	cout << "\n Coefficient, fa = ";
	cin >> fa;
	fSonicPoint =0.;
}*/
TBuckleyLeverett::TBuckleyLeverett() : TConservationLaw(1) {
	std::cout << "\n Coefficient, fa = ";
	std::cin >> fa;
//   if(fOrder!=1) std::cout << "ERROR : TBuckleyLeverett, bad system order." << std::endl;
	fSonicPoint =0.;
}
void TBuckleyLeverett::Flux(double *Ui,double *funcao) {
	double a=fa;
	double U=Ui[0];
	a*=((1.-U)*(1.-U));
	U*=U;
	a+=U;
	funcao[0]=U/a;
}
void TBuckleyLeverett::JacobFlux(double *Ui,double *jacob) {
	double a = fa;
	double U=Ui[0];
	a+=(U*U);
	U*=(1.-U);
	U*=(2*a);
	jacob[0]=U/(a*a);
}
void TBuckleyLeverett::ValJacobFlux(double *Ui,double *valjacob) {
	JacobFlux(Ui,valjacob);
   valjacob[0]=abs(valjacob[0]);
}



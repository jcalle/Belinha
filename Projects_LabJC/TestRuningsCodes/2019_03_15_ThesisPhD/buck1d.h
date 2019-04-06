
#ifndef BUCK1DH
#define BUCK1DH

#include "conslaw.h"

/*******       TBuckleyLeverett       *******/

class TBuckleyLeverett : public TConservationLaw {   //funcao nao convexa

	double fSonicPoint;
	double fa;

   public :
	TBuckleyLeverett(std::istream &input);
	//TBuckleyLeverett(char *file,double Time);
//	TBuckleyLeverett(TCompMesh *c,double Time);
//	TBuckleyLeverett(char *file);
	TBuckleyLeverett();
	~TBuckleyLeverett() {
	}

	/** @brief Returns the integrable dimension of the material */
	virtual int Dimension() const { return 1; }

	void Flux(double *Ui,double *funcao);
	void JacobFlux(double *Ui,double *jacob);
	void ValJacobFlux(double *Ui,double *valjacob);
   double MaxEigJacob() { }

};


#endif

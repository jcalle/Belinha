/**
 * @file
 * @brief Contains the implementation of the TPZDohrSubstruct methods. 
 */

#include "tpzdohrsubstruct.h"
#include <iostream>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("substruct.dohrsubstruct"));
#endif

using namespace std;

template<class TVar>
typename TPZDohrSubstruct<TVar>::EWeightType TPZDohrSubstruct<TVar>::fWeightType = TPZDohrSubstruct<TVar>::CorrectWeight;


template<class TVar>
TPZDohrSubstruct<TVar>::TPZDohrSubstruct()
{
	//Inicializacao
}

template<class TVar>
TPZDohrSubstruct<TVar>::~TPZDohrSubstruct()
{
	//Limpesa
}

template<class TVar>
void TPZDohrSubstruct<TVar>::ContributeTestV1(TPZFMatrix<TVar> &testV1, int NumCoarse) {
	/* temp1 will be equal to W(i)*Phi(i) */
	TPZFMatrix<TVar> temp1(fNEquations,fCoarseIndex.NElements());
	int i,j;
	for(i=0;i<fPhiC.Rows();i++) {
		for(j=0;j<fPhiC.Cols();j++) {
			temp1(i,j) = fPhiC(i,j)*fWeights[i];
		}
	}
	/* And now temp2 will be equal to temp1*R(ci) */
	TPZFMatrix<TVar> col(fNEquations,1);
	TPZFMatrix<TVar> temp2(fNEquations,NumCoarse,0.);
	for(i=0;i<fCoarseIndex.NElements();i++) {
		temp1.GetSub(0,i,temp1.Rows(),1,col);
		temp2.PutSub(0,fCoarseIndex[i],col);
	}
	/* And testV1 will be equal to R(i)_trans*temp2 */
	col.Resize(1,NumCoarse);
	for(i=0;i<fNEquations;i++) {
		temp2.GetSub(i,0,1,temp2.Cols(),col);
		//testV1.PutSub(fGlobalIndex[i],0,col);
		for(j=0;j<NumCoarse;j++) {
			testV1(fGlobalIndex[i],j) += col(0,j);
		}
	}
}

template<class TVar>
void TPZDohrSubstruct<TVar>::ContributeResidual(TPZFMatrix<TVar> &u, TPZFMatrix<TVar> &r){
	TPZFMatrix<TVar> ulocal(fNEquations,1,0.);
	TPZFMatrix<TVar> reslocal(fNEquations,1);
	int i;
	int neqs = fGlobalEqs.NElements();
	for (i=0;i<neqs;i++) 
	{
		std::pair<int,int> ind = fGlobalEqs[i];
		ulocal(ind.first,0) = u(ind.second,0);
	}
	fStiffness->Residual(ulocal, fLocalLoad, reslocal);
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		r(ind.second,0) += reslocal(ind.first,0);
	}
}

template<class TVar>
void TPZDohrSubstruct<TVar>::LoadWeightedResidual(const TPZFMatrix<TVar> &r_global){
	int i;
	fLocalWeightedResidual.Resize(fNEquations,1);
	int neqs = fGlobalEqs.NElements();
	for (i=0;i<neqs;i++) 
	{
		std::pair<int,int> ind = fGlobalEqs[i];
		fLocalWeightedResidual(ind.first,0) = fWeights[ind.first]*r_global.GetVal(ind.second,0);
	}
}

template<class TVar>
void TPZDohrSubstruct<TVar>::Contribute_rc(TPZFMatrix<TVar> &rc) {
	int i;
	TPZFMatrix<TVar> temp;
	fPhiC.Multiply(fLocalWeightedResidual, temp, 1);
	for (i=0;i<fCoarseIndex.NElements();i++) {
		rc(fCoarseIndex[i],0) += temp(i,0);
	}
}

/**
 * It computes the local contribution to r(c).
 * The method LoadWeightedResidual must be called before this one.
 */
template<class TVar>
void TPZDohrSubstruct<TVar>::Contribute_rc_local(TPZFMatrix<TVar> &residual_local, TPZFMatrix<TVar> &rc_local) const
{
	fPhiC_Weighted_Condensed.Multiply(residual_local, rc_local, 1);
}


template<class TVar>
void TPZDohrSubstruct<TVar>::Contribute_Kc(TPZMatrix<TVar> &Kc, TPZVec<int> &coaseindex) {
	int i;
	int j;
	for (i=0;i<fCoarseIndex.NElements();i++) {
		for (j=0;j<fCoarseIndex.NElements();j++) {
			Kc(fCoarseIndex[i],fCoarseIndex[j]) += fKCi(i,j);
		}
	}
}

template<class TVar>
void TPZDohrSubstruct<TVar>::Contribute_v1(TPZFMatrix<TVar> &v1, TPZFMatrix<TVar> &invKc_rc) {
	int i;
	//  int j;
	TPZFMatrix<TVar> temp(fCoarseIndex.NElements(), 1);
	TPZFMatrix<TVar> temp2;
	//temp = R(ci)*K(c)_inverted*r(c)
	for (i=0;i<fCoarseIndex.NElements();i++) {
		temp(i, 0) = invKc_rc(fCoarseIndex[i],0);
	}
	//temp2 = Phi*temp
	fPhiC.Multiply(temp, temp2, 0);
	//v1 += R(i)_transposta*W(i)*temp2
#ifdef ZERO_INTERNAL_RESIDU
	for(i=0; i<fInternalEqs.NElements(); i++)
	{
		temp2(fInternalEqs[i],0) = 0.;
	}
#endif
	int neqs = fGlobalEqs.NElements();
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		v1(ind.second,0) += fWeights[ind.first]*temp2(ind.first,0);
	}
}

template<class TVar>
void TPZDohrSubstruct<TVar>::Contribute_v1_local(TPZFMatrix<TVar> &v1_local, TPZFMatrix<TVar> &invKc_rc_local) const {
	int neqs = fGlobalEqs.NElements();
	v1_local.Resize(neqs, 1);
	fPhiC_Weighted_Condensed.Multiply(invKc_rc_local,v1_local);
}

template<class TVar>
void TPZDohrSubstruct<TVar>::Contribute_v2(TPZFMatrix<TVar> &v2) {
	int i;
	SolveSystemZi();
#ifdef ZERO_INTERNAL_RESIDU
	for(i=0; i<fInternalEqs.NElements(); i++)
	{
		fzi(fInternalEqs[i],0) = 0.;
	}
#endif
	int neqs = fGlobalEqs.NElements();
	for (i=0;i<neqs;i++) 
	{
		std::pair<int,int> ind = fGlobalEqs[i];
		v2(ind.second,0) += fWeights[ind.first] * fzi(ind.first,0);
	}
}

/*
 * It computes the local contribution to v2.
 */
template<class TVar>
void TPZDohrSubstruct<TVar>::Contribute_v2_local(TPZFMatrix<TVar> &residual_local, TPZFMatrix<TVar> &v2_local)
{
    int ncols = residual_local.Cols();
	TPZFMatrix<TVar> LocalWeightedResidual(fNEquations,ncols,0.);
	int neqs = fGlobalEqs.NElements();
	int i;
    for (int ic=0; ic<ncols; ic++) 
    {
        for (i=0;i<neqs;i++) 
        {
            std::pair<int,int> ind = fGlobalEqs[i];
            LocalWeightedResidual(ind.first,ic) += fWeights[ind.first] * residual_local(i,ic);
        }
    }
	int ncoarse = fCoarseIndex.NElements();
	// size of the kernel
	int nnull = fNullPivots.Rows();
	// number of global indices
	int nglob = fNEquations;
	/** Solving the system for zi */
	//Constructing I star is the same I star for Phi
	//C star is the same C star for Phi
	//Constructing I_lambda
	TPZFMatrix<TVar> I_lambda(ncoarse+nnull,ncols);
	I_lambda.Zero();
	//K_star_inv*C_star_trans is the same used for Phi
	/* Computing K_star_inv*W(i)*R(i)*r = K_star_inv*fLocalWeightedResidual */
	TPZFMatrix<TVar> KWeightedResidual(nglob,ncols);
	fInvertedStiffness.Solve(LocalWeightedResidual,KWeightedResidual);
	//Obtaining lambda_star
	TPZFMatrix<TVar> Lambda_star(ncoarse+nnull,1);
	TPZFMatrix<TVar> CstarKW(ncoarse+nnull,ncols);
	fC_star.MultAdd(KWeightedResidual,KWeightedResidual,CstarKW,-1,0,0);
	TPZFMatrix<TVar> temp2(ncoarse+nnull,ncols);
	I_lambda.Add(CstarKW,temp2);
	finv.Solve(temp2, Lambda_star);
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		Lambda_star.Print("Lambda_star ",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	//Obtaining z(i)
	TPZFMatrix<TVar> zi(nglob,ncols);
	temp2.Resize(nglob,ncols);
	fKeC_star.Multiply(Lambda_star,temp2);
	temp2 *= -1.;
	temp2.Add(KWeightedResidual,zi);
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		zi.Print("zi ",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
#ifdef ZERO_INTERNAL_RESIDU
	for(i=0; i<fInternalEqs.NElements(); i++)
	{
		zi(fInternalEqs[i],0) = 0.;
	}
#endif
	v2_local.Resize(neqs, ncols);
    for (int ic=0; ic<ncols; ic++) 
    {
        for (i=0;i<neqs;i++) 
        {
            std::pair<int,int> ind = fGlobalEqs[i];
            v2_local(i,ic) = fWeights[ind.first] * zi(ind.first,ic);
        }
    }	
	
	
}

template<class TVar>
void TPZDohrSubstruct<TVar>::Contribute_v3(TPZFMatrix<TVar> &v3, const TPZFMatrix<TVar> &r, TPZFMatrix<TVar> &v1Plusv2) const {
	TPZFNMatrix<100,TVar> vec_t(fNEquations,1,0.);
	TPZFNMatrix<100,TVar> vec_t2(fNEquations,1,0.);
	TPZFNMatrix<100,TVar> vec_t3(fInternalEqs.NElements(),1,0.);
	TPZFNMatrix<100,TVar> inv_sys(fInternalEqs.NElements(),1,0.);
	int i;
	int neqs = fGlobalEqs.NElements();
	//vec_t=R(i)(v1+v2)
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		vec_t(ind.first,0) = v1Plusv2(ind.second,0);
	}
	//vec_t2=K(i)*vec_t
	fStiffness->Multiply(vec_t, vec_t2, 0);
	//vec_t=R(i)r
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		vec_t(ind.first,0) = r.GetVal(ind.second,0);
	}
	//vec_t3=R(Ii)*(vec_t - vec_t2)
	for (i=0;i<fInternalEqs.NElements();i++) {
		vec_t3(i,0) = vec_t(fInternalEqs[i],0) - vec_t2(fInternalEqs[i],0);
	}
	//inv_sys = temp_inverted * vec_t3
	fInvertedInternalStiffness.Solve(vec_t3, inv_sys);
	//vec_t=R(Ii)_transposed * inv_sys
	vec_t.Zero();
	for (i=0;i<fInternalEqs.NElements();i++) {
		vec_t(fInternalEqs[i],0) = inv_sys(i,0);
	}
#ifndef MAKEINTERNAL
	//v3=v3+R(i)_transposed*vec_t
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		v3(ind.second,0) += vec_t(ind.first,0);
	}
#else
	v3 = vec_t;
#endif
	
}

template<class TVar>
void TPZDohrSubstruct<TVar>::Contribute_v3_local(TPZFMatrix<TVar> &v1Plusv2, TPZFMatrix<TVar> &v3) const {
	TPZFNMatrix<100,TVar> vec_t2(fNEquations,1,0.);
	TPZFNMatrix<100,TVar> vec_t3(fInternalEqs.NElements(),1,0.);
	TPZFNMatrix<100,TVar> inv_sys(fInternalEqs.NElements(),1,0.);
	int i;
	//vec_t2=K(i)*vec_t
	fStiffness->Multiply(v1Plusv2, vec_t2);
	for (i=0;i<fInternalEqs.NElements();i++) {
		vec_t3(i,0) = vec_t2(fInternalEqs[i],0);
	}
	//inv_sys = temp_inverted * vec_t3
	fInvertedInternalStiffness.Solve(vec_t3, inv_sys);
	//vec_t=R(Ii)_transposed * inv_sys
	v3.Redim(fNEquations,1);
	for (i=0;i<fInternalEqs.NElements();i++) {
		v3(fInternalEqs[i],0) = -inv_sys(i,0);
	}	
}

template<class TVar>
void TPZDohrSubstruct<TVar>::Print(std::ostream &out) const
{
	out << "++++++++++++++++++++++++++++++++++++" << std::endl;
	out << "Coarse Index " << fCoarseIndex << std::endl;
	out << "Global Index " << fGlobalIndex << std::endl;
	out << "Internal Nodes " << fInternalEqs << std::endl;
	out << "Coarse Nodes " << fCoarseNodes << std::endl;
	//fEigenVectors.Print("Eigen vectors",out);
	fC_star.Print("fC_star",out);
	fKeC_star.Print("fKeC_star", out);
	fNullPivots.Print("fNullPivots", out);
	out << "fNEquations = " << fNEquations << std::endl;
	out << "fCoarseNodes " << fCoarseNodes << std::endl;
	
	out << "fCoarseIndex " << fCoarseIndex << std::endl;
	
	out << "fGlobalEqs " << fGlobalEqs << std::endl;
	out << "fInternalEqs " << fInternalEqs << std::endl;
	out << "fBoundaryEqs " << fBoundaryEqs << std::endl;
	fLocalLoad.Print("fLocalLoad",out);
	fLocalWeightedResidual.Print("fLocalWeightedResidual",out);
	fzi.Print("fzi", out);
	fAdjustSolution.Print("fAdjustSolution", out);
	
	
	fKCi.Print("Coarse Matrix",out);
	fStiffness->Print("Stiffness Matrix",out,EMathematicaInput);
	fInvertedStiffness.Matrix()->Print("Matrix Ke",out,EMathematicaInput);
	fC.Print("Matrix Ci data structure fC",out,EMathematicaInput);
	//fEigenVectors.Print("Eigenvectors",out,EMathematicaInput);
	fPhiC.Print("fPhiC = ",out,EMathematicaInput);
	out << "fWeights = " << fWeights  << endl;
	fPhiC_Weighted_Condensed.Print("fPhiC_Weighted_Condensed = ", out);
}

template<class TVar>
void TPZDohrSubstruct<TVar>::SolveSystemPhi() {
	int ncoarse = fCoarseIndex.NElements();
	int i;
	//  I_star.Print("fIStar = ",out,EMathematicaInput);
	//  std::cout << "Nci: ";
	//  std::cout << ncoarse;// << endl;
	//  std::cout << endl;
	//Constructing I_lambda
	TPZFMatrix<TVar> I_lambda(ncoarse+fNullPivots.Rows(),ncoarse);
	I_lambda.Zero();
	for (i=0;i<ncoarse;i++) {
		I_lambda(i,i)=1;
	}
	//  I_lambda.Print("ILambda = ",out,EMathematicaInput);
	
	//Obtaining lambda_star
	TPZFMatrix<TVar> Lambda_star(ncoarse+fNullPivots.Rows(),ncoarse);
	//  temp1->Print("temp1 = ",out,EMathematicaInput);
	//  out.flush();
	finv.Solve(I_lambda, Lambda_star);
	
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		Lambda_star.Print("matrix lambda star",sout );
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	//Obtaining Phi
	/** Acho q posso comentar essas duas linhas abaixo, pq naum tou vendo aplicacao pra temp2 */
	/*
	 TPZFMatrix<TVar> temp2(fNEquations,ncoarse);
	 C_star.MultAdd(Lambda_star,Lambda_star,temp2,-1,0,1,1);
	 */
	fPhiC.Resize(fNEquations,ncoarse);
	fKeC_star.Multiply(Lambda_star,fPhiC);
	fPhiC *= -1.;
	ComputeCoarseStiffness();
	//  fPhiC.Print("PhiC = ",out,EMathematicaInput);
	//  fC.Print("Ci = ",std::cout,EMathematicaInput);
}

template<class TVar>
void TPZDohrSubstruct<TVar>::SolveSystemZi() {
	int ncoarse = fCoarseIndex.NElements();
	// size of the kernel
	int nnull = fNullPivots.Rows();
	// number of global indices
	int nglob = fNEquations;
	/* Solving the system for zi */
	//Constructing I star is the same I star for Phi
	//C star is the same C star for Phi
	//Constructing I_lambda
	TPZFMatrix<TVar> I_lambda(ncoarse+nnull,1);
	I_lambda.Zero();
	//K_star_inv*C_star_trans is the same used for Phi
	/* Computing K_star_inv*W(i)*R(i)*r = K_star_inv*fLocalWeightedResidual */
	TPZFMatrix<TVar> KWeightedResidual(nglob,1);
	fInvertedStiffness.Solve(fLocalWeightedResidual,KWeightedResidual);
	//Obtaining lambda_star
	TPZFMatrix<TVar> Lambda_star(ncoarse+nnull,1);
	TPZFMatrix<TVar> CstarKW(ncoarse+nnull,1);
	fC_star.MultAdd(KWeightedResidual,KWeightedResidual,CstarKW,-1,0,0);
	TPZFMatrix<TVar> temp2(ncoarse+nnull,1);
	I_lambda.Add(CstarKW,temp2);
	finv.Solve(temp2, Lambda_star);
	//Obtaining z(i)
	fzi.Resize(nglob,1);
	temp2.Resize(nglob,1);
	fKeC_star.Multiply(Lambda_star,temp2);
	temp2 *= -1.;
	temp2.Add(KWeightedResidual,fzi);
}

template<class TVar>
void TPZDohrSubstruct<TVar>::ComputeCoarseStiffness() {
	TPZFMatrix<TVar> temp1(fNEquations,fCoarseIndex.NElements());
	fKCi.Resize(fCoarseIndex.NElements(),fCoarseIndex.NElements());
	fStiffness->Multiply(fPhiC,temp1);
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		temp1.Print("stiffness matrix times phi",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	fPhiC.MultAdd(temp1,temp1,fKCi,1,0,1);
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		fKCi.Print("Coarse stiffness matrix",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

template<class TVar>
void TPZDohrSubstruct<TVar>::ContributeGlobalDiagonal(TPZFMatrix<TVar> &StiffnessDiag) {
	int i;
	fWeights.Resize(fNEquations);
	//fWeights = diag(kci) if u(i) is on coarse or diag(Ki)
	for (i=0;i<fNEquations;i++) {
		fWeights[i] = fStiffness->operator()(i,i);
	}
	if(fWeightType == CorrectWeight)
	{
		for (i=0;i<fNEquations;i++) {
			fWeights[i] = fStiffness->operator()(i,i);
		}
		for (i=0;i<fCoarseNodes.NElements();i++) {
			fWeights[fCoarseNodes[i]] = fKCi(i,i);
		}
	} else {
		for (i=0;i<fNEquations;i++) {
			fWeights[i] = 1.;
		}
		for (i=0;i<fCoarseNodes.NElements();i++) {
			fWeights[fCoarseNodes[i]] = 1.;
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Weight used for assembly" << fWeights;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	//
	int neqs = fGlobalEqs.NElements();
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		StiffnessDiag(ind.second,0) += fWeights[ind.first];
	}
}

template<class TVar>
void TPZDohrSubstruct<TVar>::ContributeDiagonalLocal(TPZFMatrix<TVar> &StiffnessDiagLocal) {
	int i;
	fWeights.Resize(fNEquations);
	//fWeights = diag(kci) if u(i) is on coarse or diag(Ki)
	for (i=0;i<fNEquations;i++) {
		fWeights[i] = fStiffness->operator()(i,i);
	}
	if(fWeightType == CorrectWeight)
	{
		for (i=0;i<fNEquations;i++) {
			fWeights[i] = fStiffness->operator()(i,i);
		}
		for (i=0;i<fCoarseNodes.NElements();i++) {
			fWeights[fCoarseNodes[i]] = fKCi(i,i);
		}
	} else {
		for (i=0;i<fNEquations;i++) {
			fWeights[i] = 1.;
		}
		for (i=0;i<fCoarseNodes.NElements();i++) {
			fWeights[fCoarseNodes[i]] = 1.;
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Weight used for assembly" << fWeights;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	int neqs = fGlobalEqs.NElements();
	StiffnessDiagLocal.Resize(neqs,1);
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		StiffnessDiagLocal(i,0) = fWeights[ind.first];
	}
}

template<class TVar>
void TPZDohrSubstruct<TVar>::ComputeWeights(TPZFMatrix<TVar> &StiffnessDiag) {
	int i;
	//fWeights.Fill(1.);
	int neqs = fGlobalEqs.NElements();
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		fWeights[ind.first] = fWeights[ind.first] / StiffnessDiag(ind.second,0);
	}
	for(i=0; i<fInternalEqs.NElements(); i++)
	{
		fWeights[fInternalEqs[i]] = 1.;
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Weights = " <<  fWeights;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

template<class TVar>
void TPZDohrSubstruct<TVar>::ComputeWeightsLocal(TPZFMatrix<TVar> &StiffnessDiagLocal) {
	int i;
	int neqs = fGlobalEqs.NElements();
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		fWeights[ind.first] = fWeights[ind.first] / StiffnessDiagLocal(i,0);
	}
	for(i=0; i<fInternalEqs.NElements(); i++)
	{
		fWeights[fInternalEqs[i]] = 1.;
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Weights = " <<  fWeights;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	int c,nc = fPhiC.Cols();
#ifdef ZERO_INTERNAL_RESIDU
	// eliminate the values of v2 for internal equations
	int nint = fInternalEqs.NElements();
	for(i=0; i<nint; i++)
	{
		for(c=0; c<nc; c++)
		{
			fPhiC(fInternalEqs[i],c) = 0.;
		}
	}
#endif
	fPhiC_Weighted_Condensed.Resize(neqs,fPhiC.Cols());
	for (i=0;i<neqs;i++) 
	{
		for(c=0; c<nc; c++)
		{
			std::pair<int,int> ind = fGlobalEqs[i];
			fPhiC_Weighted_Condensed(i,c) = fPhiC(ind.first,c)*fWeights[ind.first];
		}
	}
}

template<class TVar>
void TPZDohrSubstruct<TVar>::ContributeKU(const TVar alpha, const TPZFMatrix<TVar> &uglobal, TPZFMatrix<TVar> &z) const
{
	int i,j;
	int nglob = fNEquations;
	int nglobglob = uglobal.Rows();
	int ncols = uglobal.Cols();
	TPZFMatrix<TVar> uloc(nglob,ncols);
	TPZFNMatrix<100,TVar> temp1(nglob,ncols),v3(nglob,ncols);
	TPZFNMatrix<100,TVar> temp1glob(nglobglob,ncols,0.),zero(nglobglob,ncols,0.),v3glob(nglobglob,ncols,0.);
	TPZFMatrix<TVar> kuloc;
	uloc.Zero();
	int neqs = fGlobalEqs.NElements();
	/* Performing R(i)*u=uloc */
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		for(j=0; j<ncols; j++)
		{	
			uloc(ind.first,j) = uglobal.Get(ind.second,j);
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Value of the local solution = ";
		uloc.Print("uloc " ,sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
#ifdef ZERO_INTERNAL_RESIDU
	// convert the local temp1 matrix para global
	for (i=0;i<neqs;i++) {
		for(j=0; j<ncols; j++)
		{
			std::pair<int,int> ind = fGlobalEqs[i];
			temp1glob(ind.second,j) = uloc(ind.first,j);
		}
	}
	// zero out the internal residu
	this->Contribute_v3(v3glob,zero,temp1glob);
#ifndef MAKEINTERNAL
	// convert the local temp1 matrix para global
	for (i=0;i<neqs;i++) {
		for(j=0; j<ncols; j++)
		{
			std::pair<int,int> ind = fGlobalEqs[i];
			v3(ind.first,j) = v3glob(ind.second,j);
		}
	}
#else
	v3 = v3glob;
#endif
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Contribution of v3";
		v3.Print("v3",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	uloc += v3;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "New value of uloc";
		uloc.Print("uloc",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
#endif
	/* Computing K(i)*uloc and storing in temp1 */
	temp1.Resize(nglob,ncols);
	fStiffness->Multiply(uloc,temp1);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "After multiplying the solution temp1 = ";
		temp1.Print("temp1 " ,sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	int zcols = z.Cols();
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		/* Sum row "i" of temp1 with row "fGlobalNodes[i]" of z */
		for (j=0;j<zcols;j++) {
			z(ind.second,j) += alpha*temp1(ind.first,j);
		}
	}
}

template<class TVar>
void TPZDohrSubstruct<TVar>::ContributeKULocal(const TVar alpha, const TPZFMatrix<TVar> &u, TPZFMatrix<TVar> &z) const
{
	int i,j;
	int nglob = fNEquations;
	int ncols = u.Cols();
	TPZFMatrix<TVar> uloc(nglob,ncols,0.);
	TPZFNMatrix<100,TVar> temp1(nglob,ncols),v3(nglob,ncols);
	TPZFMatrix<TVar> kuloc;
	uloc.Zero();
	int neqs = u.Rows();
	/* Performing R(i)*u=uloc */
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		for(j=0; j<ncols; j++)
		{	
			uloc(ind.first,j) = u.Get(i,j);
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Value of the local solution = ";
		uloc.Print("uloc " ,sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	// *****************************************************************
#ifdef ZERO_INTERNAL_RESIDU
	// zero out the internal residu
	this->Contribute_v3_local(uloc,v3);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Contribution of v3";
		v3.Print("v3",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	uloc += v3;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "New value of uloc";
		uloc.Print("uloc",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
#endif
	/* Computing K(i)*uloc and storing in temp1 */
	temp1.Resize(nglob,ncols);
	fStiffness->Multiply(uloc,temp1);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "After multiplying the solution temp1 = ";
		temp1.Print("temp1 " ,sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	int zcols = z.Cols();
	for (i=0;i<neqs;i++) {
		std::pair<int,int> ind = fGlobalEqs[i];
		/* Sum row "i" of temp1 with row "fGlobalNodes[i]" of z */
		for (j=0;j<zcols;j++) {
			z(i,j) += alpha*temp1(ind.first,j);
		}
	}
}

template<class TVar>
void TPZDohrSubstruct<TVar>::Initialize() {
	PrepareSystems();
	SolveSystemPhi();
	ComputeCoarseStiffness();
	
}

template<class TVar>
void TPZDohrSubstruct<TVar>::PrepareSystems() {
	int ncoarse = fCoarseIndex.NElements();
	//TPZFMatrix<TVar> C_transposed;
	/****************************************
	 * TEMPORARIO!!!!!!!!!!!!!!
	 */
	//  int nc = fC.Cols();
	//  fC.Redim(1,nc);
	//  fC(0,nc-1) = 1.;
	
	/*
	 *******************************************/
	/** Tou mantendo isso abaixo só pra fatoracao de K ser TVarizada e se obter os pivos nulos (caso haja) */
	
	//  ofstream out("saida.nb");
	//  this->fStiffness->Print("fK = ",out,EMathematicaInput);
	//  fC.Print("fCi = ",out,EMathematicaInput);
	TPZFMatrix<TVar> C_transposed(fC.Cols(),fC.Rows());
	fC.Transpose(&C_transposed);
	//Ke_inv*C(i)_trans
	TPZFMatrix<TVar> KeC(fNEquations,ncoarse);
	fInvertedStiffness.Solve(C_transposed,KeC);
	/** */
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Number of zero pivots of the substructure : " << fInvertedStiffness.Singular().size();
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
	//  std::cout << KeC;
	//  int neq = fInvertedStiffness.Matrix()->Rows();
	//Constructing C barra = NullPivots
	fNullPivots.Resize(fInvertedStiffness.Singular().size(),fInvertedStiffness.Matrix()->Cols());
	fNullPivots.Zero();
	//  fEigenVectors.Redim(fInvertedStiffness.Matrix()->Rows(),fInvertedStiffness.Singular().size());
	if(fNullPivots.Rows())
	{
		//    out << "NullPivots = {";
		int count = 0;
		std::list<int64_t>::iterator it = fInvertedStiffness.Singular().begin();
		for(;it != fInvertedStiffness.Singular().end(); it++,count++)
		{
			//      out << *it;
			//      if(count < fInvertedStiffness.Singular().size()-1) out << ",";
			fNullPivots(count,*it) = 1.;
		}
		//    out << "};\n";
		//    fInvertedStiffness.Solve(fEigenVectors,fEigenVectors);
	}
	//  NullPivots.Print("fNullPivots = ",out,EMathematicaInput);
	
	//Constructing I star
	TPZFMatrix<TVar> I_star(ncoarse+fNullPivots.Rows(),ncoarse+fNullPivots.Rows());
	I_star.Zero();
	int i;
	for (i=ncoarse;i<ncoarse+fNullPivots.Rows();i++) {
		I_star(i,i)=1;
	}
	//Constructing C_star and storing in fC_star
	fC_star.Resize(ncoarse+fNullPivots.Rows(),fNEquations);
	fC_star.PutSub(0,0,fC);
	fC_star.PutSub(ncoarse,0,fNullPivots);
	//  C_star.Print("CStar = ",out,EMathematicaInput);
	//constructing K_star_inv*C_star_trans and storing it in fKeC_star
	fKeC_star.Resize(fNEquations,ncoarse+fNullPivots.Rows());
	TPZFMatrix<TVar> C_star_trans(fNEquations,ncoarse+fNullPivots.Rows());
	fC_star.Transpose(&C_star_trans);
	fInvertedStiffness.Solve(C_star_trans,fKeC_star);
	//Constructing StepSolver finv
	TPZFMatrix<TVar> *temp1 = new TPZFMatrix<TVar>(ncoarse+fNullPivots.Rows(),ncoarse+fNullPivots.Rows());
	fC_star.MultAdd(fKeC_star,I_star,*temp1,-1,1,0);
	finv.SetMatrix(temp1);
	finv.SetDirect(ELU);
}

/**
 * Adjust the residual to reflect a static condensation
 * The residual corresponding to the internal nodes will be zeroed
 */
template<class TVar>
void TPZDohrSubstruct<TVar>::AdjustResidual(TPZFMatrix<TVar> &r_global)
{
	int nglob = fNEquations;
	int nint = this->fInternalEqs.NElements();
	TPZFMatrix<TVar> rint(nint,1,0.),radapt(nglob,1,0.);
	int i;
	for(i=0; i<nint; i++)
	{
		if(fGlobalIndex[fInternalEqs[i]] >= 0)
		{
			rint(i,0) = r_global(fGlobalIndex[fInternalEqs[i]]);
		}
	}
	
	// initializar radapt com zero
	int neqs = fGlobalEqs.NElements();
	for(i=0; i<neqs; i++)
	{
		std::pair<int,int> ind = fGlobalEqs[i];
		radapt(ind.first,0) = r_global(ind.second,0);
	}
	fInvertedInternalStiffness.Solve(rint,fAdjustSolution);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Adjusted internal solution";
		fAdjustSolution.Print("fAdjustSolution",sout);
		rint.Print("Internal residual",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	TPZFMatrix<TVar> rloc(nglob,1,0.),radjust(nglob,1,0.);
	for(i=0; i<nint; i++)
	{
		rloc(fInternalEqs[i],0) = fAdjustSolution(i,0);
	}
	fStiffness->MultAdd(rloc,radapt,radjust,-1.,1.);
	for(i=0; i<nint; i++)
	{
		if(fabs(radjust(fInternalEqs[i],0)) >= 1.e-10)
		{
			std::cout << "Internal node " << i << " was not zeroed " << radjust(fInternalEqs[i],0) << std::endl;
		}
		radjust(fInternalEqs[i],0) = 0.;
	}
	
	// somar a contribuicao no r_global (tarefa)
	for(i=0; i<neqs; i++)
	{
		std::pair<int,int> ind = fGlobalEqs[i];
		r_global(ind.second,0) = radjust(ind.first,0);
	}
}

/**
 * Add the internal solution to the final result
 */
template<class TVar>
void TPZDohrSubstruct<TVar>::AddInternalSolution(TPZFMatrix<TVar> &sol)
{
	int nglob = fNEquations;
	int nint = this->fInternalEqs.NElements();
	TPZFMatrix<TVar> rint(nint,1),radapt(nglob,1,0.);
	int i,c,nc = sol.Cols();
	for(i=0; i<nint; i++) for(c=0; c<nc; c++)
	{
		if(fGlobalIndex[fInternalEqs[i]] >=0) sol(fGlobalIndex[fInternalEqs[i]],c) += this->fAdjustSolution(i,c);
	}
	
}

template class TPZDohrSubstruct<float>;
template class TPZDohrSubstruct<double>;
template class TPZDohrSubstruct<long double>;

//#ifdef STATE_COMPLEX
template class TPZDohrSubstruct<std::complex<float> >;
template class TPZDohrSubstruct<std::complex<double> >;
template class TPZDohrSubstruct<std::complex<long double> >;
//#endif
/**
 * @file
 * @brief Contains the implementation of the TPZSequenceSolver methods.
 */

#include "pzseqsolver.h"
#include "TPZPersistenceManager.h"

using namespace std;

template<class TVar>
TPZSequenceSolver<TVar>::TPZSequenceSolver(TPZMatrix<TVar> *refmat) : TPZRegisterClassId(&TPZSequenceSolver::ClassId),
TPZMatrixSolver<TVar>(refmat), fSolvers() {
}

template<class TVar>
TPZSequenceSolver<TVar>::TPZSequenceSolver(const TPZSequenceSolver<TVar> & copy): TPZRegisterClassId(&TPZSequenceSolver::ClassId),
TPZMatrixSolver<TVar>(copy) {
    int nums = copy.fSolvers.NElements();
    int s;
    for(s=0; s<nums; s++) AppendSolver(*copy.fSolvers[s]);
}

template <class TVar>
TPZSolver<TVar> * TPZSequenceSolver<TVar>::Clone() const {
    return new TPZSequenceSolver(*this);
}

template<class TVar>
void TPZSequenceSolver<TVar>::AppendSolver(TPZMatrixSolver<TVar> & solve){
    fSolvers.Push((TPZMatrixSolver<TVar> *) solve.Clone());
}

template <class TVar>
void TPZSequenceSolver<TVar>::ResetSolver() {
    int nums = fSolvers.NElements();
    int s;
    for(s=0; s<nums; s++) delete fSolvers.Pop();
}

template<class TVar>
void TPZSequenceSolver<TVar>::Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual){
	if(!this->Matrix()) {
		cout << "TPZSequenceSolver::Solve called without a matrix pointer\n";
	}
	TPZAutoPointer<TPZMatrix<TVar> > mat = this->Matrix();
	if(result.Rows() != mat->Rows() || result.Cols() != F.Cols()) {
		result.Redim(mat->Rows(),F.Cols());
	}
	
	this->fScratch = F;
	TPZFMatrix<TVar> delu(result.Rows(),result.Cols(),0.);
	TPZFMatrix<TVar> resloc(F.Rows(),F.Cols(),0.);
	result.Zero();
	if(this->fScratch.Rows() != result.Rows() || this->fScratch.Cols() != result.Cols()) {
		this->fScratch.Redim(result.Rows(),result.Cols());
	}
    int nums = fSolvers.NElements();
    int s;
    for(s=0; s<nums; s++) {
        fSolvers[s]->Solve(this->fScratch,delu,&resloc);
        result += delu;
        mat->Residual(result,F,this->fScratch);
    }
    if(residual) *residual = this->fScratch;
}


/**
 This method will reset the matrix associated with the solver
 This is useful when the matrix needs to be recomputed in a non linear problem
 */
template<class TVar>
void TPZSequenceSolver<TVar>::ResetMatrix()
{
    int nums = fSolvers.NElements();
    int s;
    for(s=0; s<nums; s++) {
        fSolvers[s]->ResetMatrix();
    }
    TPZMatrixSolver<TVar>::ResetMatrix();
}

/**
 Updates the values of the preconditioner based on the values of the matrix
 */
template <class TVar>
void TPZSequenceSolver<TVar>::UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > matrix)
{
    int nums = fSolvers.NElements();
    int s;
    for(s=0; s<nums; s++) {
        fSolvers[s]->UpdateFrom(matrix);
    }
    TPZMatrixSolver<TVar>::UpdateFrom(matrix);
}

template<class TVar>
void TPZSequenceSolver<TVar>::Write(TPZStream &buf, int withclassid) const
{
	TPZMatrixSolver<TVar>::Write(buf, withclassid);
	int StackSz = fSolvers.NElements();
	buf.Write(&StackSz, 1);
	int i = 0;
	for(i = 0; i < StackSz; i++)
	{
            TPZPersistenceManager::WritePointer(fSolvers[i], &buf);
	}
	
}
template <class TVar>
void TPZSequenceSolver<TVar>::Read(TPZStream &buf, void *context)
{
	TPZMatrixSolver<TVar>::Read(buf, context);
	int StackSz = 0;
	buf.Read(&StackSz, 1);
	fSolvers.Resize(StackSz);
	for(int i = 0; i< StackSz; i++)
	{
            fSolvers[i] = dynamic_cast<TPZMatrixSolver<TVar> *>(TPZPersistenceManager::GetInstance(&buf));
	}
}

template class TPZSequenceSolver<float>;
template class TPZSequenceSolver<std::complex<float> >;
template class TPZSequenceSolver<double>;
template class TPZSequenceSolver<std::complex<double> >;
template class TPZSequenceSolver<long double>;
template class TPZSequenceSolver<std::complex<long double> >;

#ifndef BORLAND
template class TPZRestoreClass< TPZSequenceSolver<float>>;
template class TPZRestoreClass< TPZSequenceSolver<std::complex<float>>>;
template class TPZRestoreClass< TPZSequenceSolver<double>>;
template class TPZRestoreClass< TPZSequenceSolver<std::complex<double>>>;
template class TPZRestoreClass< TPZSequenceSolver<long double>>;
template class TPZRestoreClass< TPZSequenceSolver<std::complex<long double>>>;
#endif

/**
 * @file
 * @brief Contains the implementation of the TPZParFrontMatrix methods.
 */

#include "TPZParFrontMatrix.h"
#include "TPZFrontMatrix.h"
#include "pzsfulmat.h"
#include "TPZFront.h"
#include "pzstack.h"
#include "pzreal.h"
#include <math.h>
#include "pz_pthread.h"

#include "tpzeqnarray.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

#include "pzlog.h"

#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("pz.strmatrix.frontstructmatrix"));
static LoggerPtr loggerfw(Logger::getLogger("pz.frontal.frontmatrix.fw"));

#endif

/** @brief Initializing semaphore */
pthread_mutex_t mutex_write = PTHREAD_MUTEX_INITIALIZER;
/** @brief Initializing condition */
pthread_cond_t conda_write = PTHREAD_COND_INITIALIZER;

// At the class constructor creates a thread
// this thread will be active while ParFrontMatrix is active
// It will check if a stack contains some equations to be writen to the disk

template<class TVar, class store, class front>
TPZParFrontMatrix<TVar, store, front>::TPZParFrontMatrix():
TPZRegisterClassId(&TPZParFrontMatrix::ClassId), fFinish(0)
{
	fEqnStack.Resize(0);
	pthread_mutex_t mlocal = PTHREAD_MUTEX_INITIALIZER;
	fwritelock = mlocal;
	pthread_cond_t clocal = PTHREAD_COND_INITIALIZER;
	fwritecond = clocal;
	/*	fFront.Reset();
	 fStorage.Reset();
	 fNumElConnected.Resize(0);
	 fLastDecomposed = -1;
	 fNumEq=0;
	 */
}

template<class TVar, class store, class front>
TPZParFrontMatrix<TVar, store, front>::TPZParFrontMatrix(int64_t globalsize) :
TPZRegisterClassId(&TPZParFrontMatrix::ClassId), TPZFrontMatrix<TVar, store, front>(globalsize),
fFinish(0)
{
	fEqnStack.Resize(0);
	pthread_mutex_t mlocal = PTHREAD_MUTEX_INITIALIZER;
	fwritelock = mlocal;
	pthread_cond_t clocal = PTHREAD_COND_INITIALIZER;
	fwritecond = clocal;
}

template<class TVar, class store, class front>
TPZParFrontMatrix<TVar, store, front>::~TPZParFrontMatrix(){
}

template<class TVar, class store, class front>
void TPZParFrontMatrix<TVar, store, front>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec < int64_t > & destinationindex)
{
	
	// message #1.3 to fFront:TPZFront
	this->fFront.AddKel(elmat, destinationindex);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Frondwidth after AddKel "<< this->fFront.FrontSize();
		LOGPZ_INFO(loggerfw,sout.str())
	}
#endif

	int64_t mineq, maxeq;
	this->EquationsToDecompose(destinationindex, mineq, maxeq);
	if(maxeq >= mineq) {

		TPZEqnArray<TVar> *AuxEqn = new TPZEqnArray<TVar>;
		
		this->fFront.DecomposeEquations(mineq,maxeq,*AuxEqn);
		this->CheckCompress();
		PZ_PTHREAD_MUTEX_LOCK(&fwritelock,"TPZParFrontMatrix<...>::AddKel()");
		fEqnStack.Push(AuxEqn);
		if(maxeq == this->Rows()-1){
			cout << "Decomposition finished" << endl;
			cout.flush();
			FinishWriting();
			//fStorage.ReOpen();
		}
		PZ_PTHREAD_MUTEX_UNLOCK(&fwritelock,"TPZParFrontMatrix<...>::AddKel()");
		PZ_PTHREAD_COND_SIGNAL(&fwritecond,"TPZParFrontMatrix<...>::AddKel()");
	}
	this->fDecomposed = this->fFront.GetDecomposeType();
} 
template<class TVar, class store, class front>
void TPZParFrontMatrix<TVar, store, front>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec < int64_t > & sourceindex, TPZVec < int64_t > & destinationindex)
{
	this->fFront.AddKel(elmat, sourceindex, destinationindex);
#ifdef LOG4CXX
    if (loggerfw->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Frondwidth after AddKel "<< this->fFront.FrontSize();
		LOGPZ_DEBUG(loggerfw,sout.str())
	}
#endif
	int64_t mineq, maxeq;
	this->EquationsToDecompose(destinationindex, mineq, maxeq);

	if(maxeq >= mineq) {
		TPZEqnArray<TVar> *AuxEqn = new TPZEqnArray<TVar>;

		this->fFront.DecomposeEquations(mineq,maxeq,*AuxEqn);
		this->CheckCompress();
		PZ_PTHREAD_MUTEX_LOCK(&fwritelock,"TPZParFrontMatrix<...>::AddKel()");
		fEqnStack.Push(AuxEqn);
		if(maxeq == this->Rows()-1){
            //check if writeing is over and closes file
			cout << endl << "Decomposition finished" << endl;
			cout.flush();
			FinishWriting();
			this->fFront.Reset(0);
			//fStorage.ReOpen();
		}
		PZ_PTHREAD_MUTEX_UNLOCK(&fwritelock,"TPZParFrontMatrix<...>::AddKel()");
		PZ_PTHREAD_COND_SIGNAL(&fwritecond,"TPZParFrontMatrix<...>::AddKel()");
	}
	this->fDecomposed = this->fFront.GetDecomposeType();
}

template<class TVar, class store, class front>
void TPZParFrontMatrix<TVar, store, front>::FinishWriting(){
	// FinishWriting already has a lock
	//PZ_PTHREAD_MUTEX_LOCK(&fwritelock);
	cout << endl << "FinishWriting" << endl;
	cout.flush();     
	fFinish = 1;
	// FinishWriting already has a lock
	//pthread_mutex_unlock(&fwritelock);
	PZ_PTHREAD_COND_SIGNAL(&fwritecond,"TPZParFrontMatrix<...>::FinishWriting()");
} 

template<class TVar, class store, class front>
void * TPZParFrontMatrix<TVar, store, front>::WriteFile(void *t){
	TPZParFrontMatrix<TVar, store, front> *parfront = (TPZParFrontMatrix<TVar, store, front>*) t;    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Entering WriteFile thread execution";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	cout << endl << "Entering Decomposition" << endl;
	cout.flush();
	while(1){
		TPZStack<TPZEqnArray<TVar> *> local;
		PZ_PTHREAD_MUTEX_LOCK(&parfront->fwritelock,"TPZParFrontMatrix<...>::WriteFile()");
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
		{
			std::stringstream sout;
			sout << "Acquired writelock";
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		if(parfront->fEqnStack.NElements() == 0){
			if(parfront->fFinish == 1) {
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
				{
					std::stringstream sout;
					sout << "Terminating WriteFile thread execution";
					LOGPZ_DEBUG(logger,sout.str())
				}
#endif
				cout << "Leaving WHILE" << endl;
				cout.flush();
				break;
			}
#ifdef LOG4CXX
            if (logger->isDebugEnabled())
			{
				std::stringstream sout;
				sout << "Entering cond_wait on fwritecond variable";
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			PZ_PTHREAD_COND_WAIT(&parfront->fwritecond, &parfront->fwritelock,"TPZParFrontMatrix<...>::WriteFile()");
		}
		
		local = parfront->fEqnStack;
		parfront->fEqnStack.Resize(0);
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
		{
			std::stringstream sout;
			sout << "Copied the equation stack releasing the writelock";
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		
		PZ_PTHREAD_MUTEX_UNLOCK(&parfront->fwritelock,"TPZParFrontMatrix<...>::WriteFile()");
		int64_t neqn = local.NElements();

		int64_t eq;
		for(eq=0; eq<neqn; eq++) {
			parfront->fStorage.AddEqnArray(local[eq]);
			delete local[eq];
		}
	}
	parfront->fStorage.FinishWriting();
	parfront->fStorage.ReOpen();
	parfront->fFinish = 0;
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Releasing writelock";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	PZ_PTHREAD_MUTEX_UNLOCK(&parfront->fwritelock,"TPZParFrontMatrix<...>::WriteFile()");
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Falling through on the write thread";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	std::cout << "Terminating write thread\n";
	return (0);
} 

template<class TVar>
class TPZStackEqnStorage;
template<class TVar>
class TPZFileEqnStorage;
template<class TVar>
class TPZFrontSym;
template<class TVar>
class TPZFrontNonSym;


template class TPZParFrontMatrix<float, TPZStackEqnStorage<float>, TPZFrontSym<float> >;
template class TPZParFrontMatrix<float, TPZFileEqnStorage<float>, TPZFrontSym<float> >;
template class TPZParFrontMatrix<float, TPZStackEqnStorage<float>, TPZFrontNonSym<float> >;
template class TPZParFrontMatrix<float, TPZFileEqnStorage<float>, TPZFrontNonSym<float> >;

template class TPZParFrontMatrix<double, TPZStackEqnStorage<double>, TPZFrontSym<double> >;
template class TPZParFrontMatrix<double, TPZFileEqnStorage<double>, TPZFrontSym<double> >;
template class TPZParFrontMatrix<double, TPZStackEqnStorage<double>, TPZFrontNonSym<double> >;
template class TPZParFrontMatrix<double, TPZFileEqnStorage<double>, TPZFrontNonSym<double> >;

template class TPZParFrontMatrix<long double, TPZStackEqnStorage<long double>, TPZFrontSym<long double> >;
template class TPZParFrontMatrix<long double, TPZFileEqnStorage<long double>, TPZFrontSym<long double> >;
template class TPZParFrontMatrix<long double, TPZStackEqnStorage<long double>, TPZFrontNonSym<long double> >;
template class TPZParFrontMatrix<long double, TPZFileEqnStorage<long double>, TPZFrontNonSym<long double> >;

template class TPZParFrontMatrix<std::complex<float>, TPZStackEqnStorage<std::complex<float> >, TPZFrontSym<std::complex<float> > >;
template class TPZParFrontMatrix<std::complex<float>, TPZFileEqnStorage<std::complex<float> >, TPZFrontSym<std::complex<float> > >;
template class TPZParFrontMatrix<std::complex<float>, TPZStackEqnStorage<std::complex<float> >, TPZFrontNonSym<std::complex<float> > >;
template class TPZParFrontMatrix<std::complex<float>, TPZFileEqnStorage<std::complex<float> >, TPZFrontNonSym<std::complex<float> > >;

template class TPZParFrontMatrix<std::complex<double>, TPZStackEqnStorage<std::complex<double> >, TPZFrontSym<std::complex<double> > >;
template class TPZParFrontMatrix<std::complex<double>, TPZFileEqnStorage<std::complex<double> >, TPZFrontSym<std::complex<double> > >;
template class TPZParFrontMatrix<std::complex<double>, TPZStackEqnStorage<std::complex<double> >, TPZFrontNonSym<std::complex<double> > >;
template class TPZParFrontMatrix<std::complex<double>, TPZFileEqnStorage<std::complex<double> >, TPZFrontNonSym<std::complex<double> > >;

template class TPZParFrontMatrix<std::complex<long double>, TPZStackEqnStorage<std::complex<long double> >, TPZFrontSym<std::complex<long double> > >;
template class TPZParFrontMatrix<std::complex<long double>, TPZFileEqnStorage<std::complex<long double> >, TPZFrontSym<std::complex<long double> > >;
template class TPZParFrontMatrix<std::complex<long double>, TPZStackEqnStorage<std::complex<long double> >, TPZFrontNonSym<std::complex<long double> > >;
template class TPZParFrontMatrix<std::complex<long double>, TPZFileEqnStorage<std::complex<long double> >, TPZFrontNonSym<std::complex<long double> > >;


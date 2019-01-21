/**
 * @file
 * @brief Contains the TPZDohrPrecond class which implements a matrix which computes the preconditioner developed by Dohrmann. \n
 * Also contains TPZDohrPrecondThreadV1Data, TPZDohrPrecondV2SubData and TPZDohrPrecondV2SubDataList structure.
 */

#ifndef TPZDOHRPRECOND_H
#define TPZDOHRPRECOND_H

#include "pzmatrix.h"
#include "tpzautopointer.h"
#include "tpzdohrsubstruct.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrassembly.h"
#include "tpzdohrassemblelist.h"

#include <list>
#include <semaphore.h>

#include "pz_pthread.h"

/**
 * \addtogroup substructure
 * @{
 */
/**
 * @brief Implements a matrix which computes the preconditioner developed by Dohrmann. \ref substructure "Sub Structure"
 * @author Philippe Devloo
 * @since 2006
 */
template <class TVar, class TSubStruct> 
class TPZDohrPrecond : public TPZMatrix<TVar>
{
	/** @brief The matrix class is a placeholder for a list of substructures */
	std::list<TPZAutoPointer<TSubStruct> > fGlobal;
	/** @brief The global matrix associated with the coarse degrees of freedom */
	TPZStepSolver<TVar> * fCoarse; //K(c)

	/** @brief Size of the coarse system */
	int64_t fNumCoarse; //n(c)
	
	/** @brief Number of threads used during preconditioning */
	int fNumThreads;
	
	TPZAutoPointer<TPZDohrAssembly<TVar> > fAssemble;
	
    //
    
public:
    /** @brief Constructor with matrix */
    TPZDohrPrecond(TPZDohrMatrix<TVar, TSubStruct> &origin, TPZAutoPointer<TPZDohrAssembly<TVar> > assemble);
	/** @brief Copy constructor */
	TPZDohrPrecond(const TPZDohrPrecond &copy);
    /** @brief Empty constructor for restoring */
    TPZDohrPrecond();
	
    ~TPZDohrPrecond();
    
	// CLONEDEF(TPZDohrPrecond)
	virtual TPZMatrix<TVar>*Clone() const { return new TPZDohrPrecond(*this); }
    
    /** @brief The matrix class is a placeholder for a list of substructures */
	std::list<TPZAutoPointer<TSubStruct> > &Global()
    {
        return fGlobal;
    }
	
	/** @brief Initialize the necessary datastructures */
	/** It will compute the coarse matrix, coarse residual and any other necessary data structures */
	void Initialize();
    
        /**
         * Specify the number of threads for preconditioning
         */
        void SetNumThreads(int numthreads)
        {
	  fNumThreads = numthreads;
        }

	/**
	 * @brief The only method any matrix class needs to implement
	 * @param x Is x on the above operation. It must be a vector!
	 * @param y Is y on the above operation
	 * @param z Is z on the above operation
	 * @param alpha Is alpha on the above operation
	 * @param beta Is beta on the above operation
	 * @param opt Indicates if is Transpose or not
	 */
	/**
	 * In this case the variable x represents the rhs and z the result of the preconditioning \n
	 * When used as a preconditioner y will be zero
	 * In fact, it will compute \f$v1+v2+v3\f$ \n
	 * It computes \f$ z = beta * y + alpha * opt(this)*x\f$ but z and x can not overlap in memory.
	 */
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z, const TVar alpha,const TVar beta,const int opt) const;
	
    /** Copy of the MultAdd using TBB */
    virtual void MultAddTBB(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z, const TVar alpha,const TVar beta,const int opt) const;
    
	/** @brief Specify the solution process for the coarse matrix */
	void SetSolver(TPZSolver<TVar> &solver);
	
	/** @brief Compute the contribution of the coarse matrix */
	void ComputeV1(const TPZFMatrix<TVar> &x, TPZFMatrix<TVar> &v1) const;
	/** @brief Compute the contribution of the sub domains */
	void ComputeV2(const TPZFMatrix<TVar> &x, TPZFMatrix<TVar> &v2) const;
    
	/** @brief Routines to send and receive messages */
	public:
virtual int ClassId() const;
	
    /**
	 * @brief Unpacks the object structure from a stream of bytes
	 * @param buf The buffer containing the object in a packed form
	 * @param context 
	 */
	virtual void  Read(TPZStream &buf, void *context );
	/**
	 * @brief Packs the object structure in a stream of bytes
	 * @param buf Buffer which will receive the bytes
	 * @param withclassid
	 */
	virtual void Write( TPZStream &buf, int withclassid );

};

template <class TVar, class TSubStruct> 
int TPZDohrPrecond<TVar, TSubStruct>::ClassId() const{
    return Hash("TPZDohrPrecond") ^ TPZMatrix<TVar>::ClassId() << 1 ^ TSubStruct().ClassId() << 2;
}

#include <pthread.h>

/**
 * @brief Auxiliar structure with thread to compute the preconditioner developed by Dohrmann. \ref substructure "Sub Structure"
 */
template <class TVar, class TSubStruct> 
struct TPZDohrPrecondThreadV1Data {
	TPZDohrPrecondThreadV1Data() : fDohrMatrix(0), fInput(0), fOutput(0)
	{
	}
	TPZDohrPrecondThreadV1Data(const TPZDohrPrecond<TVar, TSubStruct> *ptr, const TPZFMatrix<TVar> &input, TPZFMatrix<TVar> &output) : fDohrMatrix(ptr),
	fInput(&input), fOutput(&output)
	{
	}
	/** @brief Pointer to the dohr matrix */
	const TPZDohrPrecond<TVar, TSubStruct> * fDohrMatrix;
	/** @brief Input matrix */
	const TPZFMatrix<TVar> * fInput;
	/** @brief Matrix where the coarse solution will be contributed */
	TPZFMatrix<TVar> *fOutput;
	/** @brief Compute the contribution of the coarse matrix */
	static void *ComputeV1(void *dataptr)
	{
		TPZDohrPrecondThreadV1Data<TVar, TSubStruct> *ptr = (TPZDohrPrecondThreadV1Data<TVar, TSubStruct> *) dataptr;
		ptr->fDohrMatrix->ComputeV1(*(ptr->fInput),*(ptr->fOutput));
		return dataptr;
	}
};

/**
 * @brief Auxiliar structure for v2 vector to compute the preconditioner developed by Dohrmann. \ref substructure "Sub Structure"
 */
template <class TVar, class TSubStruct> 
struct TPZDohrPrecondV2SubData {
	
	TPZDohrPrecondV2SubData() : fSubStructure(0), fInput_local(0), fv2_local(0)
	{
	}
	
	TPZDohrPrecondV2SubData(int subindex, const TPZAutoPointer<TSubStruct> &substruct, TPZAutoPointer<TPZFMatrix<TVar> > res_local) : fSubStructure(substruct),
	fInput_local(res_local)
	{
		fv2_local = new TPZDohrAssembleItem<TVar>(subindex, res_local->Rows());
	}
	/** @note Protect ourselves from default copy constructors */
	TPZDohrPrecondV2SubData(const TPZDohrPrecondV2SubData<TVar, TSubStruct> &copy) : fSubStructure(copy.fSubStructure), fInput_local(copy.fInput_local),
	fv2_local(copy.fv2_local)
	{
	}
	
	TPZDohrPrecondV2SubData &operator=(const TPZDohrPrecondV2SubData &copy)
	{
		fSubStructure = copy.fSubStructure;
		fInput_local = copy.fInput_local;
		fv2_local = copy.fv2_local;
		return *this;
	}
	
	bool IsValid()
	{
		return fSubStructure;
	}
	
	~TPZDohrPrecondV2SubData()
	{
	}
	
	/** @brief Pointer to the dohr matrix */
	TPZAutoPointer<TSubStruct> fSubStructure;
	/** @brief Input matrix */
	TPZAutoPointer<TPZFMatrix<TVar> > fInput_local;
	/** @brief The local contribution to the v2 vector */
	TPZAutoPointer<TPZDohrAssembleItem<TVar> > fv2_local;
};

template<class TVar>
struct TPZDohrAssembleList;

/**
 * @brief Auxiliar structure with list for v2 vector data. \ref substructure "Sub Structure"
 */
template <class TVar, class TSubStruct> 
struct TPZDohrPrecondV2SubDataList {
	TPZDohrPrecondV2SubDataList(TPZAutoPointer<TPZDohrAssembleList<TVar> > &assemble) : fAssemblyStructure(assemble)
	{
	  PZ_PTHREAD_MUTEX_INIT(&fAccessLock, 0, "TPZDohrPrecondV2SubDataList::TPZDohrPrecondV2SubDataList()");
	}
	~TPZDohrPrecondV2SubDataList()
	{
	  PZ_PTHREAD_MUTEX_DESTROY(&fAccessLock, "TPZDohrPrecondV2SubDataList::~TPZDohrPrecondV2SubDataList()");
	}
	
    /** @brief Mutex which will enable the access protection of the list */
	pthread_mutex_t fAccessLock;
	
    /** @brief The list of structures which need to be computed */
	std::list<TPZDohrPrecondV2SubData<TVar, TSubStruct> > fWork;
    
	/** @brief Interface to add items in a thread safe way */
	void AddItem(TPZDohrPrecondV2SubData<TVar, TSubStruct> &data)
	{
	  PZ_PTHREAD_MUTEX_LOCK(&fAccessLock, "TPZDohrPrecondV2SubDataList::AddItem()");
	  fWork.push_back(data);
	  PZ_PTHREAD_MUTEX_UNLOCK(&fAccessLock, "TPZDohrPrecondV2SubDataList::AddItem()");
	}
	/** @brief Interface to pop an item in a thread safe way */
	TPZDohrPrecondV2SubData<TVar, TSubStruct> PopItem()
	{
		TPZDohrPrecondV2SubData<TVar, TSubStruct> result;
		PZ_PTHREAD_MUTEX_LOCK(&fAccessLock, "TPZDohrPrecondV2SubDataList::PopItem()");
		if (fWork.size()) {
			result = *fWork.begin();
			fWork.pop_front();
		}
		PZ_PTHREAD_MUTEX_UNLOCK(&fAccessLock, "TPZDohrPrecondV2SubDataList::PopItem()");
		return result;
	}
	
	/** @brief The local contribution to the v2 vector */
	TPZAutoPointer<TPZDohrAssembleList<TVar> > fAssemblyStructure;
	
	/** @brief The procedure which executes the lengthy process */
	static void *ThreadWork(void *voidptr);
	
};

/** @} */

#endif

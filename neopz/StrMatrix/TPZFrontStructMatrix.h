/**
 * @file
 * @brief Contains the TPZFrontStructMatrix class which responsible for a interface among Finite Element Package and Matrices package to frontal method.
 */

#ifndef TPZFRONTSTRUCTMATRIX_H
#define TPZFRONTSTRUCTMATRIX_H

#include "pzstrmatrix.h"
#include "pzcmesh.h" 

#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"

#include "pzelmat.h"

#include "pzmatrix.h"
#include "pzfmatrix.h"

/**
 * @brief Responsible for a interface among Finite Element Package and Matrices package to frontal method. \ref structural "Structural Matrix" \ref frontal "Frontal"
 * @ingroup structural frontal
 * @note Type parameter for TPZFrontStructMatrix frontal matrix. \n It can assume values TPZFrontSym and TPZFrontNonSym for symmetric and non symmetric matrices
 */
/**
 * Prevents users from all the necessary information to work with all matrices classes \n
 * It facilitates considerably the use of TPZAnalysis
 */
template<class front> 
class TPZFrontStructMatrix : public TPZStructMatrix {
	
protected:
	/** @brief This vector contains an ordered list */
	/** The elements must be asssembled in that order so the frontal works on its best performance */
	TPZVec<int> fElementOrder;
    int f_quiet;
	
	/**
	 * @brief Returns a vector containing all elements connected to a degree of freedom.
	 * @param numelconnected Vector containing the number of connections for every ith dof
	 */
	void GetNumElConnected(TPZVec <int> &numelconnected);

	/** @brief It is applied over fElementOrder putting it in the correct order. */
	void OrderElement();//TPZVec <int> &elorder);
	
	/** @brief Resequence the connects according to the element order*/
	void AdjustSequenceNumbering();
	
    /** @brief Used Decomposition method */
    DecomposeType fDecomposeType;

public:

    /**
     * @brief Class constructor
     * < href="http://www.fec.unicamp.br/~longhin">link text</a>
     * < href="http://www.fec.unicamp.br/~phil">link text</a>
     */ 
	TPZFrontStructMatrix(TPZCompMesh *);
	
	TPZFrontStructMatrix(const TPZFrontStructMatrix &copy) : TPZStructMatrix(copy), fElementOrder(copy.fElementOrder),f_quiet(copy.f_quiet), fDecomposeType(copy.fDecomposeType)
	{
	}
	
    /// Set the decomposition type
    virtual void SetDecomposeType(DecomposeType dectype)
    {
        fDecomposeType = dectype;
    }
    

	static int main();
	
    /** @brief Class destructor */ 
	virtual ~TPZFrontStructMatrix();
    
    
	
	/** @brief Returns a pointer to TPZMatrix<STATE> */
	TPZMatrix<STATE> * Create();
	
	/** @brief Clones a TPZFrontStructMatrix */
	TPZStructMatrix * Clone();
	
	/**
	 * @brief Assemble a stiffness matrix according to rhs
	 * @param stiffness Stiffness matrix to assembled
	 * @param rhs Vector containing loads
	 * @param guiInterface pointer to user interface
	 */ 	
	void AssembleNew(TPZMatrix<STATE> & stiffness
					 , TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	/**
	 * @brief Assemble a stiffness matrix.
	 * @param stiffness Stiffness matrix to assembled
	 * @param rhs Vector containing loads
	 * @param guiInterface pointer to user interface
	 */ 	
	void Assemble(TPZMatrix<STATE> & stiffness
				  , TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface);
	
	/**
	 * @brief Computes element matrices.
	 * @param el Actual element being computed
	 * @param ek Formed element matrix
	 * @param ef Global element load matrix
	 * @param stiffness Global stiffness matrix
	 * @param rhs Global load matrix
	 */
	/** 
	 * Each computed element matrices would then be added to Stiffness matrix
	 */
	void AssembleElement(TPZCompEl *el, TPZElementMatrix & ek
						 , TPZElementMatrix & ef, TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs); 
	
    using TPZStructMatrix::CreateAssemble;
	/**
	 * @brief Returns a pointer to TPZMatrix.
	 * @param rhs Load matrix
	 * @param guiInterface pointer to user interface
	 */
	/** 
	 * This is a mandatory function, it is neded by all StructMatrix. \n
	 * Except in frontal matrices, the returned matrix is not in its decomposed form.
	 */
	TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    void SetQuiet(int quiet);
    
private:
    TPZFrontStructMatrix();
	
    friend TPZPersistenceManager;
	
};

#endif //TPZFRONTSTRUCTMATRIX_H

/**
 * @file
 * @brief Contains the TPZVTKGraphMesh class which implements the graphical mesh to VTK environment.
 */

#ifndef PZVTKMESH
#define PZVTKMESH

#include "pzgraphmesh.h"
#include "pzvec.h"

template<class TVar>
class TPZBlock;

/**
 * @ingroup post
 * @brief To export a graphical mesh to VTK environment. \ref post "Post processing"
 */
class TPZVTKGraphMesh : public TPZGraphMesh {
	
public:
	
	/** @brief Constructor for graphical mesh using VTK format */
	TPZVTKGraphMesh(TPZCompMesh *cmesh, int dimension, TPZMaterial * mat, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames);
    /** @brief Constructor for graphical mesh using VTK format with tensor variables */
    TPZVTKGraphMesh(TPZCompMesh *cmesh, int dimension, TPZMaterial * mat, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensnames);
	/** @brief Copy constructor for graphical mesh using VTK format */
	TPZVTKGraphMesh(TPZCompMesh *cmesh,int dim,TPZVTKGraphMesh *graph,TPZMaterial * mat);
	
	virtual void DrawMesh(int numcases);
	virtual void DrawNodes();
	virtual void DrawConnectivity(MElementType type);
	virtual void DrawSolution(int step, REAL time);
	virtual void DrawSolution(TPZBlock<REAL> &Sol);
	virtual void DrawSolution(char *var = 0);
	
protected:
	virtual void SequenceNodes();
	int fNumCases;
	int fNumSteps;
	
};

#endif

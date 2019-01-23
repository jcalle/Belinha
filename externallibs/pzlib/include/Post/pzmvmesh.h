/**
 * @file
 * @brief Contains the TPZMVGraphMesh class which implements graphical mesh to MVGraph package.
 */

#ifndef MVGRIDH
#define MVGRIDH

#include "pzgraphmesh.h"
#include "pzvec.h"

template<class TVar>
class TPZBlock;

/**
 * @ingroup post
 * @brief Implements graphical mesh to MVGraph package. \ref post "Post processing"
 */
/** 
 * MVGraph: Multivariate Interactive Visualization can to be obtained from: \n
 * <a href="http://cran.r-project.org/web/packages/mvgraph/index.html">Lattes</a>
 */
class TPZMVGraphMesh : public TPZGraphMesh {
	
public:
	
	/** @brief Constructor for graphical mesh using MVGraph format */
    TPZMVGraphMesh(TPZCompMesh *cmesh, int dimension, TPZMaterial * mat, const TPZVec<std::string> &scalarnames, const TPZVec<std::string> &vecnames);
	/** @brief Copy constructor for graphical mesh using MVGraph format */
	TPZMVGraphMesh(TPZCompMesh *cmesh,int dim,TPZMVGraphMesh *graph,TPZMaterial * mat);
	
	/** @brief Draw graphical mesh */
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


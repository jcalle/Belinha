/**
 * @file
 * @brief Contains the implementation of the TPZFStructMatrix methods. 
 */

#include "pzfstrmatrix.h"
#include "pzfmatrix.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include <sstream>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.tpzfstructmatrix"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
#endif


using namespace std;

TPZMatrix<STATE> * TPZFStructMatrix::Create(){
	int64_t neq = fEquationFilter.NActiveEquations();
    
	return new TPZFMatrix<STATE>(neq,neq,0.);
}

TPZFStructMatrix::TPZFStructMatrix() : TPZStructMatrix()
{
}

TPZFStructMatrix::TPZFStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}

TPZFStructMatrix::TPZFStructMatrix(TPZAutoPointer<TPZCompMesh> mesh) : TPZStructMatrix(mesh)
{
}

TPZStructMatrix * TPZFStructMatrix::Clone(){
    return new TPZFStructMatrix(*this);
}

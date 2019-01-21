/**
 * @file
 * @brief Contains implementations of the TPZSubMeshAnalysis methods: implementation of the TPZSubMeshAnalysis class.
 */

#include "pzsmanal.h"
#include "pzsubcmesh.h"
#include "pzfmatrix.h"
#include "pzstrmatrix.h"
#include "pzsolve.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.analysis.pzsmanalysis"));

#endif

using namespace std;

// Construction/Destruction

TPZSubMeshAnalysis::TPZSubMeshAnalysis(TPZSubCompMesh *mesh) : TPZRegisterClassId(&TPZSubMeshAnalysis::ClassId),
TPZAnalysis(mesh,true), fReducableStiff(0){
	fMesh = mesh;
    if (fMesh)
    {
        fReferenceSolution.Redim(fCompMesh->NEquations(),1);
    }
}

TPZSubMeshAnalysis::~TPZSubMeshAnalysis()
{
	
}

/** @brief Set the computational mesh of the analysis. */
void TPZSubMeshAnalysis::SetCompMesh(TPZCompMesh * mesh, bool mustOptimizeBandwidth)
{
    TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(mesh);
    if (submesh) {
        fMesh = submesh;
    }
    else
    {
        DebugStop();
    }
    TPZAnalysis::SetCompMesh(mesh, mustOptimizeBandwidth);
    if (fCompMesh) {
        fReferenceSolution.Redim(fCompMesh->NEquations(), 1);
    }
}



void TPZSubMeshAnalysis::Assemble(){
	

    TPZCompMesh *mesh = Mesh();
    TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(mesh);
    if (!submesh) {
        DebugStop();
    }
    TPZCompMesh *fathermesh = submesh->Mesh();
    if (fathermesh->NElements() < 50)
    {
        std::cout << "Assembling the SubCompMesh index " << fMesh->Index() << std::endl;
    }
	int numeq = fCompMesh->NEquations();
	int numinternal = fMesh->NumInternalEquations();
	fReferenceSolution.Redim(numeq,1);
	fRhs.Redim(numeq,1);
    if(!fReducableStiff) 
    {
        fReducableStiff = new TPZMatRed<STATE, TPZFMatrix<STATE> > ();
    }
	fReducableStiff->Redim(numeq,numinternal);
	TPZMatRed<STATE, TPZFMatrix<STATE> > *matred = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (fReducableStiff.operator->());
    if(!fSolver->Matrix())
    {
        if (fStructMatrix->HasRange()) {
            DebugStop();
        }
        fStructMatrix->SetEquationRange(0, numinternal);
        fSolver->SetMatrix(fStructMatrix->Create());	
        fStructMatrix->EquationFilter().Reset();
    }
    matred->SetMaxNumberRigidBodyModes(fMesh->NumberRigidBodyModes());
	//	fReducableStiff.SetK00(fSolver->Matrix());
	// this will initialize fK00 too
	matred->SetSolver(dynamic_cast<TPZMatrixSolver<STATE> *>(fSolver->Clone()));
	//	TPZStructMatrix::Assemble(fReducableStiff,fRhs, *fMesh);
//	time_t before = time (NULL);
	fStructMatrix->Assemble(fReducableStiff,fRhs,fGuiInterface);
//	time_t after = time(NULL);
//	double diff = difftime(after, before);
//	std::cout << __PRETTY_FUNCTION__ << " tempo " << diff << std::endl;
}

void TPZSubMeshAnalysis::Run(std::ostream &out){
	
	//fReducableStiff.Print("Reducable stiff before assembled");
	fReferenceSolution = fSolution;
//	time_t tempo = time(NULL);
	Assemble();
//	time_t tempodepois = time(NULL);
//	double elapsedtime = difftime(tempodepois, tempo);
	
//	std::cout << "Tempo para assemblagem " << elapsedtime << std::endl;
	if (!fReducableStiff) {
        DebugStop();
    }
	TPZMatRed<STATE, TPZFMatrix<STATE> > *matred = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (fReducableStiff.operator->());
	matred->SetF(fRhs);
}
void TPZSubMeshAnalysis::CondensedSolution(TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
//	time_t tempo = time(NULL);
	if (!fReducableStiff) {
        DebugStop();
    }
	TPZMatRed<STATE, TPZFMatrix<STATE> > *matred = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (fReducableStiff.operator->());
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        matred->Print("Before = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    matred->K11Reduced(ek, ef);
	
//	time_t tempodepois = time(NULL);
//	double elapsedtime = difftime(tempodepois, tempo);
	
//	std::cout << "Tempo para inversao " << elapsedtime << std::endl;
	
}

/** @brief compute the reduced right hand side using the current stiffness. Abort if there is no stiffness computed */
void TPZSubMeshAnalysis::ReducedRightHandSide(TPZFMatrix<STATE> &rhs)
{
    if (!fReducableStiff) {
        DebugStop();
    }
    TPZMatRed<STATE, TPZFMatrix<STATE> > *matred = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (fReducableStiff.operator->());
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        matred->Print("Before = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    matred->SetF(fRhs);
    matred->F1Red(rhs);
    
}


void TPZSubMeshAnalysis::LoadSolution(const TPZFMatrix<STATE> &sol)
{
	
	//	sol.Print("sol");
	int numinter = fMesh->NumInternalEquations();
	int numeq = fMesh->TPZCompMesh::NEquations();
	TPZFMatrix<STATE> soltemp(numeq-numinter,1,0.);
	int i;
	for(i=0; i<numeq-numinter; i++) {
		soltemp(i,0) = sol.GetVal(numinter+i,0)-fReferenceSolution(numinter+i,0);
	}
	TPZFMatrix<STATE> uglobal(numeq,1,0.);
    if(fReducableStiff)
    {
        TPZMatRed<STATE, TPZFMatrix<STATE> > *matred = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (fReducableStiff.operator->());
        matred->UGlobal(soltemp,uglobal);        
        fSolution = fReferenceSolution + uglobal;
    }
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        soltemp.Print("External DOF Solution",sout);
        uglobal.Print("Expanded solution",sout);
        fSolution.Print("fSolution",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	TPZAnalysis::LoadSolution();
}

int TPZSubMeshAnalysis::ClassId() const{
    return Hash("TPZSubMeshAnalysis") ^ TPZAnalysis::ClassId() << 1;
}

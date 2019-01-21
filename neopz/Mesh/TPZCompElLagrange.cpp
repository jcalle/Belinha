//
//  TPZCompElLagrange.cpp
//  PZ
//
//  Created by Philippe Devloo on 11/2/13.
//
//

#include "TPZCompElLagrange.h"
#include "TPZMaterial.h"
#include "pzelmat.h"

TPZCompElLagrange::~TPZCompElLagrange()
{
    
}

/**
 * @brief Method for creating a copy of the element in a patch mesh
 * @param mesh Patch clone mesh
 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
 * @param gl2lcElMap map the computational elements
 */
/**
 * Otherwise of the previous clone function, this method don't
 * copy entire mesh. Therefore it needs to map the connect index
 * from the both meshes - original and patch
 */
TPZCompEl *TPZCompElLagrange::ClonePatchEl(TPZCompMesh &mesh,
                                std::map<int64_t,int64_t> & gl2lcConMap,
                                std::map<int64_t,int64_t> & gl2lcElMap) const
{
    TPZCompElLagrange *newel = new TPZCompElLagrange(mesh,*this,gl2lcElMap);
    for (int64_t l=0; l<fDef.size(); l++) {
        for (int i=0; i<2; i++) {
            newel->fDef[l].fIdf[i] = fDef[l].fIdf[i];
            std::map<int64_t,int64_t>::iterator it = gl2lcConMap.find(fDef[l].fConnect[i]);
            if (it != gl2lcConMap.end()) {
                newel->fDef[l].fConnect[i] = it->second;
            }
            else
            {
                DebugStop();
            }
        }
    }
    return newel;
}

/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZCompElLagrange::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    InitializeElementMatrix(ek, ef);
#ifdef PZDEBUG
    if (ef.fMat.Cols() != 1) {
        DebugStop();
    }
#endif
    int64_t nlagrange = fDef.size();
    int64_t count = 0;
    for (int64_t l=0; l<nlagrange; l++)
    {
        TPZConnect &c0 = Connect(2*l);
        int blsize0 = c0.NShape()*c0.NState();
        TPZConnect &c1 = Connect(2*l+1);
        int blsize1 = c1.NShape()*c1.NState();
        ek.fMat(count+fDef[l].fIdf[0],count+fDef[l].fIdf[0]) = 1.;
        ek.fMat(count+fDef[l].fIdf[0],count+blsize0+fDef[l].fIdf[1]) = -1.;
        ek.fMat(count+blsize0+fDef[l].fIdf[1],count+fDef[l].fIdf[0]) = -1.;
        ek.fMat(count+blsize0+fDef[l].fIdf[1],count+blsize0+fDef[l].fIdf[1]) = 1.;
        const TPZBlock<STATE> &bl = Mesh()->Block();
        STATE diff = bl(c0.SequenceNumber(),0,fDef[l].fIdf[0],0)-bl(c1.SequenceNumber(),0,fDef[l].fIdf[1],0);
        ef.fMat(count+fDef[l].fIdf[0],0) = -diff;
        ef.fMat(count+blsize0+fDef[l].fIdf[1],0) = diff;
        count += blsize0+blsize1;
    }
}

/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
//void TPZCompElLagrange::CalcResidual(TPZElementMatrix &ef)
//{
//    
//}

void TPZCompElLagrange::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef){
    int numloadcases = 1;
	int numdof = 1;
    TPZMaterial *mat = this->Material();
    if (mat)
    {
        numloadcases = mat->NumLoadCases();
    }
	const int ncon = this->NConnects();
    
    ek.fMesh = Mesh();
    ek.fType = TPZElementMatrix::EK;
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
    
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fNumStateVars = numdof;
	ef.fNumStateVars = numdof;
	int i;
    int numeq=0;
	for(i=0; i<ncon; i++){
        TPZConnect &c = Connect(i);
        int nshape = c.NShape();
        int nstate = c.NState();
        
		ek.fBlock.Set(i,nshape*nstate);
		ef.fBlock.Set(i,nshape*nstate);
        numeq += nshape*nstate;
	}
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,numloadcases);
	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ef.fConnect)[i] = ConnectIndex(i);
		(ek.fConnect)[i] = ConnectIndex(i);
	}
}//void



int TPZCompElLagrange::ClassId() const{
    return Hash("TPZCompElLagrange") ^ TPZCompEl::ClassId() << 1;
}

#ifndef REFLECTEDSHOCKNONALIGNEDHPP
#define REFLECTEDSHOCKNONALIGNEDHPP

#include "pzreal.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "pzgmesh.h"
#include "pzflowcmesh.h"
#include "pzeulerconslaw.h"

void RSNAMeshPoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int64_t> > &elms);
TPZGeoMesh * CreateRSNAGeoMesh(TPZGeoMesh *gmesh, TPZVec< TPZVec< REAL > > & nodes,
                               TPZVec< TPZVec< int64_t > > & elms,
                               MElementType ElType, int matId,
                               TPZVec<TPZGeoEl *> & gEls,
                               int nSubdiv);
// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh * RSNACompMesh(TPZFlowCompMesh *cmesh, REAL CFL, REAL delta,
                               int degree, int nSubdiv,
                               TPZArtDiffType DiffType,
                               TPZTimeDiscr Diff_TD,
                               TPZTimeDiscr ConvVol_TD,
                               TPZTimeDiscr ConvFace_TD);


#endif
//
//  TPZCompMeshTools.cpp
//  PZ
//
//  Created by Philippe Devloo on 6/4/15.
//
//

#include "TPZCompMeshTools.h"
#include "pzelchdiv.h"
#include "pzshapepiram.h"
#include "TPZOneShapeRestraint.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzelementgroup.h"
#include "pzsubcmesh.h"
#include "pzcondensedcompel.h"
#include "pzmultiphysicselement.h"
#include "TPZMeshSolution.h"

#include <algorithm>

#include "pzmetis.h"
#include "pzsloan.h"
#include "TPZSloanRenumbering.h"
#include "TPZCutHillMcKee.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.analysis"));
#endif

#ifdef USING_BOOST
#include "TPZBoostGraph.h"
/**
 * @brief To renumbering will use boost library.
 * @ingroup analysis
 */
#define RENUMBER TPZSloanRenumbering()
#else
/**
 * @brief To renumbering will use sloan library.
 * @ingroup analysis
 */
#define RENUMBER TPZSloanRenumbering()
//#define RENUMBER TPZCutHillMcKee()
#endif

static TPZOneShapeRestraint SetupPyramidRestraint(TPZCompEl *cel, int side);

using namespace pzshape;

void TPZCompMeshTools::AddHDivPyramidRestraints(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZCompElHDiv< pzshape::TPZShapePiram > *pyram = dynamic_cast<TPZCompElHDiv< pzshape::TPZShapePiram > *>(cel);
        if (!pyram) {
            continue;
        }
        
        if (cel->Type() != EPiramide)
        {
            DebugStop();
        }
        TPZGeoEl *gel = cel->Reference();
        /// we need to identify the connect to restraint
        /// start with face 14
        int is;
        int numneigh = 0;
        for (is=14; is<18; is++) {
            TPZGeoElSide gelside(gel,is);
            TPZCompElSide large = gelside.LowerLevelCompElementList2(1);
            // if the side is constrained, it wont be used as a reference
            if (large) {
                numneigh++;
                continue;
            }
            // if the side has smaller elements connected to it, discard the side
            TPZStack<TPZCompElSide> celstack;
            gelside.HigherLevelCompElementList2(celstack, 1, 0);
            if (celstack.size()) {
                numneigh++;
                continue;
            }
            // see if there is a pyramid element  neighbour of the pyramid
            gelside.EqualLevelCompElementList(celstack, 1, 0);
            int localel;
            for (localel=0; localel<celstack.size(); localel++) {
                if (celstack[localel].Reference().Element()->Type() == ETetraedro || celstack[localel].Reference().Element()->Type() == ETriangle) break;
            }
            if (celstack.size()) {
                numneigh++;
            }
            TPZCompElHDiv<pzshape::TPZShapeTetra> *pir;
            /// the pyramid has a neighbour of pyramid type
            if (localel < celstack.size()) {
                // find out if the face is the restrained face

                pir = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeTetra> *>(celstack[localel].Element());
                if (pir == NULL) {
                    break;
                }
                const int restrainedface = pir->RestrainedFace();
                const int neighbourside = celstack[localel].Side();
                if (restrainedface == neighbourside) {
                    continue;
                }
                break;
            }
        }
        if (numneigh == 0) {
            is = 17;
        }
        if (is == 18 || is == 13) {
            cmesh->Print();
            DebugStop();
        }
        TPZOneShapeRestraint rest = SetupPyramidRestraint(pyram, is);
        pyram->AddShapeRestraint(rest);
        /// add the restraint datastructure to the neighbour
        TPZGeoElSide gelside(cel->Reference(),is);
        TPZGeoElSide neighbour = gelside.Neighbour();
        TPZCompElSide celneigh = neighbour.Reference();
        TPZInterpolatedElement *interp = dynamic_cast<TPZInterpolatedElement *>(celneigh.Element());
        interp->AddShapeRestraint(rest);
    }
    
        
    ExpandHDivPyramidRestraints(cmesh);
}

static int idf[6] = {2,1,1,2,0,0};

static int nextface(int side, int inc)
{
    return ((side-14)+inc)%4+14;
}

static void GetGeoNodeIds(TPZGeoEl *gel, int side, TPZVec<int64_t> &ids)
{
    int nsidenodes = gel->NSideNodes(side);
    if (nsidenodes != 3) {
        DebugStop();
    }
    for (int i=0; i<3; i++) {
        ids[i] = gel->SideNodePtr(side, i)->Id();
    }
}

static int SideIdf(TPZCompElHDiv<TPZShapePiram> *cel, int side)
{
    TPZManVector<int64_t,3> ids(3);
    GetGeoNodeIds(cel->Reference(),side,ids);
    int transformid = TPZShapeTriang::GetTransformId2dT(ids);
    return idf[transformid];
}
static TPZOneShapeRestraint SetupPyramidRestraint(TPZCompEl *cel, int side)
{
    TPZCompElHDiv<TPZShapePiram> *pir = dynamic_cast<TPZCompElHDiv<TPZShapePiram> *>(cel);
    if (!pir) {
        DebugStop();
    }
    TPZOneShapeRestraint result;
    int mult = 1;
    for (int fc=0; fc<4; fc++)
    {
        int face = nextface(side, fc);
        int64_t nodeindex = pir->SideConnectIndex(0, face);
        result.fFaces[fc] = std::make_pair(nodeindex, SideIdf(pir, face));
        result.fOrient[fc] = mult*pir->SideOrient(face-13);
        mult *= -1;
    }
    return result;
}

void TPZCompMeshTools::ExpandHDivPyramidRestraints(TPZCompMesh *cmesh)
{
    /// Add the shapeonerestraints of all restrained connects build a map indexed on the restrained connect
    std::map<int64_t, TPZOneShapeRestraint> AllRestraints;
    typedef std::list<TPZOneShapeRestraint>::iterator itlist;
    if (!cmesh) {
        DebugStop();
    }
    int64_t nelem = cmesh->NElements();
    for (int64_t iel=0; iel<nelem; iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        std::list<TPZOneShapeRestraint> ellist;
        ellist = cel->GetShapeRestraints();
        for (itlist it = ellist.begin(); it != ellist.end(); it++) {
            int64_t elindex = it->fFaces[0].first;
            AllRestraints[elindex] = *it;
        }
    }
    /// For each element that has a oneshape restraint, add the oneshape restraints of all connects
    for (int64_t iel=0; iel<nelem; iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        std::list<TPZOneShapeRestraint> ellist;
        std::map<int64_t, TPZOneShapeRestraint> restraintmap;
        std::set<int64_t> connectset;
        ellist = cel->GetShapeRestraints();
        // we will insert the restraints afterwards
        cel->ResetShapeRestraints();
        for (itlist it = ellist.begin(); it != ellist.end(); it++) {
//            std::cout << "1 including connect index " << it->fFaces[0].first << " to restraint map\n";
//            it->Print(std::cout);
            restraintmap[it->fFaces[0].first] = *it;
            for (int i=1; i<4; i++) {
                connectset.insert(it->fFaces[i].first);
            }
        }
        while (connectset.size()) {
            int64_t cindex = *connectset.begin();
            connectset.erase(cindex);
            if (AllRestraints.find(cindex) != AllRestraints.end()) {
                TPZOneShapeRestraint restloc = AllRestraints[cindex];
//                std::cout << "2 including connect index " << cindex << " to restraint map\n";
//                restloc.Print(std::cout);
                restraintmap[cindex] = restloc;
                for (int i=1; i<4; i++) {
                    int64_t locindex = restloc.fFaces[i].first;
                    std::cout << "3 verifying connect index " << locindex << "\n";
                    if (restraintmap.find(locindex) != restraintmap.end()) {
                        std::cout << "********************* CIRCULAR RESTRAINT\n";
                    }
                    else
                    {
                        connectset.insert(locindex);
                    }
                }
            }
        }
        for (std::map<int64_t,TPZOneShapeRestraint>::iterator it = restraintmap.begin(); it != restraintmap.end(); it++) {
            cel->AddShapeRestraint(it->second);
        }
    }
    
}

void TPZCompMeshTools::LoadSolution(TPZCompMesh *cpressure, TPZFunction<STATE> &Forcing)
{
    int64_t nel = cpressure->NElements();
    for (int64_t iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cpressure->Element(iel);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != 3) {
            continue;
        }
        if (gel->Type() == EPiramide) {
            TPZGeoNode *top = gel->NodePtr(4);
            //TPZManVector<REAL,3> topco(3),valvec(1);AQUIFRAN
            TPZManVector<REAL,3>topco(3);
            TPZManVector<STATE,3>valvec(1);
            top->GetCoordinates(topco);
            Forcing.Execute(topco, valvec);
            STATE topval = valvec[0];
            TPZConnect &c = cel->Connect(0);
            int64_t seqnum = c.SequenceNumber();
            for (int i=0; i<4; i++) {
                cpressure->Block()(seqnum,0,i,0) = topval;
            }
            for (int i=0; i<4; i++) {
                TPZConnect &c = cel->Connect(i+1);
                TPZGeoNode *no = gel->NodePtr(i);
                no->GetCoordinates(topco);
                Forcing.Execute(topco, valvec);
                STATE nodeval = valvec[0];
                seqnum = c.SequenceNumber();
                cpressure->Block()(seqnum,0,0,0) = nodeval-topval;
            }
        }
        else
        {
            int ncorner = gel->NCornerNodes();
            for (int i=0; i<ncorner; i++) {
                TPZConnect &c = cel->Connect(i);
                TPZGeoNode *no = gel->NodePtr(i);
                TPZManVector<REAL,3> topco(3);
                no->GetCoordinates(topco);
                TPZManVector<STATE,3> valvec(1);
                Forcing.Execute(topco, valvec);
                STATE nodeval = valvec[0];
                int64_t seqnum = c.SequenceNumber();
                cpressure->Block()(seqnum,0,0,0) = nodeval;
            }
        }
    }
}

void TPZCompMeshTools::GroupElements(TPZCompMesh *cmesh, std::set<int64_t> elbasis, std::set<int64_t> &grouped)
{
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    int dim = cmesh->Dimension();
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        if (elbasis.find(el) == elbasis.end()) {
            continue;
        }
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        std::set<int64_t> elgroup;
        
        std::set<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        std::cout << "cel " << cel->Index() << " connects ";
        for (std::set<int64_t>::iterator it=connectlist.begin(); it != connectlist.end(); it++) {
            std::cout << *it << " ";
        }
        std::cout << std::endl;
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            std::cout << "side " << is << std::endl;
            TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> celstack;
            gelside.ConnectedCompElementList(celstack, 0, 0);
            int64_t nelstack = celstack.size();
            for (int64_t elst=0; elst<nelstack; elst++) {
                TPZCompElSide celst=celstack[elst];
//                TPZGeoElSide gelst =celst.Reference();
                TPZCompEl *celsidelement = celst.Element();
                if (elbasis.find(celsidelement->Index()) == elbasis.end()) {
                    continue;
                }
                std::set<int64_t> smallset;
                celsidelement->BuildConnectList(smallset);
                std::cout << "neigh " << celsidelement->Index() << " connects ";
                for (std::set<int64_t>::iterator it=smallset.begin(); it != smallset.end(); it++) {
                    std::cout << *it << " ";
                }
                std::cout << std::endl;
                if (std::includes(connectlist.begin(), connectlist.end(), smallset.begin(), smallset.end()))
                {
                    std::cout << "Is included\n";
                    elgroup.insert(celsidelement->Index());
                    elbasis.erase(celsidelement->Index());
                }
            }
        }
        if (elgroup.size()) {
            elgroup.insert(el);
            int64_t grindex;
            TPZElementGroup *elgr = new TPZElementGroup(*cmesh,grindex);
            for (std::set<int64_t>::iterator it = elgroup.begin(); it != elgroup.end(); it++) {
                elgr->AddElement(cmesh->Element(*it));
            }
            grouped.insert(grindex);
        }
        else
        {
            grouped.insert(el);
        }
    }
}

void TPZCompMeshTools::GroupElements(TPZCompMesh *cmesh)
{
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    int dim = cmesh->Dimension();
    int64_t nel = cmesh->NElements();
    // all elements that have been put in a group
    std::set<int64_t> grouped;
    for (int64_t el=0; el<nel; el++) {
        if (grouped.find(el) != grouped.end()) {
            continue;
        }
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        std::set<int64_t> elgroup;
        
        std::set<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> celstack;
            gelside.ConnectedCompElementList(celstack, 0, 0);
            int64_t nelstack = celstack.size();
            for (int64_t elst=0; elst<nelstack; elst++) {
                TPZCompElSide celst=celstack[elst];
                TPZCompEl *celsidelement = celst.Element();
                if (grouped.find(celsidelement->Index()) != grouped.end()) {
                    continue;
                }
                std::set<int64_t> smallset;
                celsidelement->BuildConnectList(smallset);
                if (std::includes(connectlist.begin(), connectlist.end(), smallset.begin(), smallset.end()))
                {
                    elgroup.insert(celsidelement->Index());
                }
            }
        }
        if (elgroup.size()) {
            elgroup.insert(el);
            int64_t grindex;
            TPZElementGroup *elgr = new TPZElementGroup(*cmesh,grindex);
            for (std::set<int64_t>::iterator it = elgroup.begin(); it != elgroup.end(); it++) {
                elgr->AddElement(cmesh->Element(*it));
                grouped.insert(*it);
            }
            grouped.insert(grindex);
        }
    }
    cmesh->ComputeNodElCon();
}

/// ungroup all embedded elements of the computational mesh
void TPZCompMeshTools::UnGroupElements(TPZCompMesh *cmesh){
    
    int64_t nelem = cmesh->NElements();

    //unwrapping element groups
    for(int64_t i=0; i<nelem; i++){
        TPZCompEl *el = cmesh->ElementVec()[i];
        TPZElementGroup * group = dynamic_cast<TPZElementGroup*>(el);
        if(group){
            group->Unwrap();
        }
    }
    
}

/// uncondensed elements for the elements that have internal nodes
void TPZCompMeshTools::UnCondensedElements(TPZCompMesh *cmesh){
    
    int64_t nelem = cmesh->NElements();
    
    //unwrapping condensed compel
    for(int64_t i=0; i<nelem; i++){
        TPZCompEl *el = cmesh->ElementVec()[i];
        TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(el);
        if(cond){
            cond->Unwrap();
        }
    }
    
}

/// Put the element set into a subcompmesh and make the connects internal
void TPZCompMeshTools::PutinSubmeshes(TPZCompMesh *cmesh, std::map<int64_t,std::set<int64_t> >&elindices, std::map<int64_t,int64_t> &indices, bool KeepOneLagrangian)
{
    for (std::map<int64_t,std::set<int64_t> >::iterator it = elindices.begin(); it != elindices.end(); it++) {
        int64_t index;
        TPZSubCompMesh *subcmesh = new TPZSubCompMesh(*cmesh,index);
        indices[it->first] = index;
        for (std::set<int64_t>::iterator itloc = it->second.begin(); itloc != it->second.end(); itloc++) {
            subcmesh->TransferElement(cmesh, *itloc);
        }
    }
    cmesh->ComputeNodElCon();
    for (std::map<int64_t,int64_t>::iterator it = indices.begin(); it != indices.end(); it++) {
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cmesh->Element(it->second));
        if (!subcmesh) {
            DebugStop();
        }
        if (KeepOneLagrangian)
        {
            int64_t nconnects = subcmesh->NConnects();
            for (int64_t ic=0; ic<nconnects; ic++) {
                TPZConnect &c = subcmesh->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    break;
                }
            }
        }
        subcmesh->MakeAllInternal();
    }

    
    
}



/// Put the element set into a subcompmesh and make the connects internal
void TPZCompMeshTools::PutinSubmeshes(TPZCompMesh *cmesh, std::set<int64_t> &elindices, int64_t &index, bool KeepOneLagrangian)
{
    TPZSubCompMesh *subcmesh = new TPZSubCompMesh(*cmesh,index);
    for (std::set<int64_t>::iterator it = elindices.begin(); it != elindices.end(); it++) {
        subcmesh->TransferElement(cmesh, *it);
    }
    cmesh->ComputeNodElCon();
    if (KeepOneLagrangian)
    {
        int64_t nconnects = subcmesh->NConnects();
        for (int64_t ic=0; ic<nconnects; ic++) {
            TPZConnect &c = subcmesh->Connect(ic);
            if (c.LagrangeMultiplier() > 0) {
                c.IncrementElConnected();
                break;
            }
        }
    }
    subcmesh->MakeAllInternal();
//    int numthreads = 0;
//    subcmesh->SetAnalysisSkyline(numthreads, 0, 0);
}


/// created condensed elements for the elements that have internal nodes
void TPZCompMeshTools::CreatedCondensedElements(TPZCompMesh *cmesh, bool KeepOneLagrangian, bool keepmatrix)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        int nc = cel->NConnects();
        int ic;
        if (KeepOneLagrangian) {
            for (ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    break;
                }
            }
        }
        for (ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.HasDependency() || c.NElConnected() > 1) {
                continue;
            }
            break;
        }
        bool cancondense = (ic != nc);
        if(cancondense)
        {
            TPZCondensedCompEl *cond = new TPZCondensedCompEl(cel, keepmatrix);
        }
        
    }

    cmesh->CleanUpUnconnectedNodes();
    
}

static void ComputeError(TPZCompEl *cel, TPZFunction<STATE> &func, TPZCompMesh *mesh2, TPZVec<STATE> &square_errors);

static void ComputeError(TPZCondensedCompEl *cel, TPZFunction<STATE> &func, TPZCompMesh *mesh2, TPZVec<STATE> &square_errors)
{
    TPZCompEl *ref = cel->ReferenceCompEl();
    ComputeError(ref, func, mesh2, square_errors);
}

static void ComputeError(TPZMultiphysicsElement *cel, TPZFunction<STATE> &func, TPZCompMesh *mesh2, TPZVec<STATE> &square_errors)
{
    TPZManVector<STATE,3> errors(3,0.);
    bool store_error = false;
    cel->EvaluateError(func, errors, store_error);
    int64_t index = cel->Index();
    TPZCompMesh *mesh = cel->Mesh();
    for (int i=0; i<3; i++) {
        mesh->ElementSolution()(index,i) = errors[i];
        square_errors[i] += errors[i]*errors[i];
    }
}

static void ComputeError(TPZElementGroup *cel, TPZFunction<STATE> &func, TPZCompMesh *mesh2, TPZVec<STATE> &square_errors)
{
    TPZStack<TPZCompEl *,5> celstack;
    celstack = cel->GetElGroup();
    int64_t nel = celstack.size();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *subcel = celstack[el];
        ComputeError(subcel, func, mesh2, square_errors);
    }
}

static void ComputeError(TPZInterpolationSpace *cel, TPZFunction<STATE> &func, TPZCompMesh *mesh2, TPZVec<STATE> &square_errors)
{
    DebugStop();
}


static void ComputeError(TPZCompEl *cel, TPZFunction<STATE> &func, TPZCompMesh *mesh2, TPZVec<STATE> &square_errors)
{
    TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
    // acumulate the errors of the submeshes
    if (sub) {
        TPZCompMeshTools::ComputeDifferenceNorm(sub, mesh2, square_errors);
        return;
    }
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
    if(elgr)
    {
        ComputeError(elgr, func, mesh2, square_errors);
        return;
    }
    TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
    if (cond) {
        ComputeError(cond, func, mesh2, square_errors);
        return;
    }
    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
    if (intel) {
        ComputeError(intel, func, mesh2, square_errors);
        return;
    }
    TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if(mphys)
    {
        ComputeError(mphys, func, mesh2, square_errors);
        return;
    }
    
}
/// compute the norm of the difference between two meshes
/// square of the errors are computed for each element of mesh1
void TPZCompMeshTools::ComputeDifferenceNorm(TPZCompMesh *mesh1, TPZCompMesh *mesh2, TPZVec<STATE> &square_errors)
{
    int64_t nel = mesh1->NElements();
    int dim = mesh1->Dimension();
    if(square_errors.size() != 3)
    {
        DebugStop();
    }
    mesh1->ElementSolution().Redim(mesh1->NElements(), 3);
    
    int materialid = 1;
    TPZMeshSolution func(mesh2,materialid);
    
//    mesh2->Reference()->ResetReference();
//    mesh2->LoadReferences();
    if (nel >= 1000) {
        std::cout << "ComputeDifferenceNorm nelem " << nel << std::endl;
    }
    
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = mesh1->Element(el);
        if (!cel) {
            continue;
        }
        ComputeError(cel, func, mesh2, square_errors);

        if (nel >= 1000 && (el+1)%1000 == 0) {
            std::cout << "*";
        }
        if(el%(20*1000) == 0)
        {
            std::cout << square_errors <<  std::endl;
        }
    }
    if(nel >= 1000) std::cout << std::endl;
}

/// adjust the polynomial orders of the hdiv elements such that the internal order is higher than the sideorders
void TPZCompMeshTools::AdjustFluxPolynomialOrders(TPZCompMesh *fluxmesh, int hdivplusplus)
{
    int dim = fluxmesh->Dimension();
    /// loop over all the elements
    int64_t nel = fluxmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != dim) {
            continue;
        }
        // compute the maxorder
        int maxorder = -1;
        int ncon = intel->NConnects();
        for (int i=0; i<ncon-1; i++) {
            int conorder = intel->Connect(i).Order();
            maxorder = maxorder < conorder ? conorder : maxorder;
        }
        int nsides = gel->NSides();
        int nconside = intel->NSideConnects(nsides-1);
        // tive que tirar para rodar H1
        //        if (nconside != 1 || maxorder == -1) {
        //            DebugStop();
        //        }
        int64_t cindex = intel->SideConnectIndex(nconside-1, nsides-1);
        TPZConnect &c = fluxmesh->ConnectVec()[cindex];
        if (c.NElConnected() != 1) {
            DebugStop();
        }
        if (c.Order() != maxorder+hdivplusplus) {
            //            std::cout << "Changing the order of the central connect " << cindex << " from " << c.Order() << " to " << maxorder+hdivplusplus << std::endl;
            // change the internal connect order to be equal do maxorder
            intel->SetSideOrder(nsides-1, maxorder+hdivplusplus);
        }
    }
    fluxmesh->ExpandSolution();
}

/// set the pressure order acording to the order of internal connect of the elements of the fluxmesh
void TPZCompMeshTools::SetPressureOrders(TPZCompMesh *fluxmesh, TPZCompMesh *pressuremesh)
{
    // build a vector with the required order of each element in the pressuremesh
    // if an element of the mesh dimension of the fluxmesh does not have a corresponding element in the pressuremesh DebugStop is called
    int meshdim = fluxmesh->Dimension();
    pressuremesh->Reference()->ResetReference();
    pressuremesh->LoadReferences();
    TPZManVector<int64_t> pressorder(pressuremesh->NElements(),-1);
    int64_t nel = fluxmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != meshdim) {
            continue;
        }
        int nsides = gel->NSides();
        int64_t cindex = intel->SideConnectIndex(0, nsides-1);
        TPZConnect &c = fluxmesh->ConnectVec()[cindex];
        int order = c.Order();
        TPZCompEl *pressureel = gel->Reference();
        TPZInterpolatedElement *pintel = dynamic_cast<TPZInterpolatedElement *>(pressureel);
        if (!pintel) {
            DebugStop();
        }
        pressorder[pintel->Index()] = order;
    }
    pressuremesh->Reference()->ResetReference();
    nel = pressorder.size();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        TPZInterpolatedElement *pintel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!pintel) {
            continue;
        }
        if (pressorder[el] == -1) {
            continue;
        }
        pintel->PRefine(pressorder[el]);
    }
    
    pressuremesh->ExpandSolution();
}


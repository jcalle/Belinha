/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticPyramid methods. 
 */
#include "tpzquadraticpyramid.h"
#include "tpzgeoblend.h"
#include "tpzgeoelmapped.h"

#include "pzgeoelrefless.h.h"
#include "tpzgeoelrefpattern.h.h"
#include "pznoderep.h.h"
#include "pzshapepiram.h"
#include "tpzgeomid.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.specialmaps.quadraticpyramid"));
#endif

#ifdef _AUTODIFF
#include "fad.h"
#endif

using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

template<class T>
void TPZQuadraticPyramid::TShape(TPZVec<T> &par,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
    
    T qsi = par[0], eta = par[1], zeta = par[2];
    T check = zeta-1.;
    if(fabs(check) < 1.E-6)
    {
        phi(0,0)  = 0.;
        phi(1,0)  = 0.;
        phi(2,0)  = 0.;
        phi(3,0)  = 0.;
        phi(4,0)  = 1.;
        phi(5,0)  = 0.;
        phi(6,0)  = 0.;
        phi(7,0)  = 0.;
        phi(8,0)  = 0.;
        phi(9,0)  = 0.;
        phi(10,0) = 0.;
        phi(11,0) = 0.;
        phi(12,0) = 0.;
        
        dphi(0,0) = 0.;
        dphi(1,0) = 0.;
        dphi(2,0) = -0.25;
        
        dphi(0,1) = 0.;
        dphi(1,1) = 0.;
        dphi(2,1) = -0.25;
        
        dphi(0,2) = 0.;
        dphi(1,2) = 0.;
        dphi(2,2) = -0.25;
        
        dphi(0,3) = 0.;
        dphi(1,3) = 0.;
        dphi(2,3) = -0.25;
        
        dphi(0,4) = 0.;
        dphi(1,4) = 0.;
        dphi(2,4) = 3.;
        
        dphi(0,5) = 0.;
        dphi(1,5) = 0.5;
        dphi(2,5) = 0.5;
        
        dphi(0,6) = -0.5;
        dphi(1,6) = 0.;
        dphi(2,6) = 0.5;
        
        dphi(0,7) = 0.;
        dphi(1,7) = -0.5;
        dphi(2,7) = 0.5;
        
        dphi(0,8) = 0.5;
        dphi(1,8) = 0.;
        dphi(2,8) = 0.5;
        
        dphi(0,9) = -1.;
        dphi(1,9) = -1.;
        dphi(2,9) = -1.;
        
        dphi(0,10) = 1.;
        dphi(1,10) = -1.;
        dphi(2,10) = -1.;
        
        dphi(0,11) = 1.;
        dphi(1,11) = 1.;
        dphi(2,11) = -1.;
        
        dphi(0,12) = -1;
        dphi(1,12) = 1.;
        dphi(2,12) = -1.;
        
        return;
    }
    
    phi(0,0) = (1. - 2.*zeta)*(1. - zeta)*(-0.25*(1. - eta*eta/((1. - zeta)*(1. - zeta)))*(1. - qsi*qsi/((1. - zeta)*(1. - zeta))) + (eta*qsi*(-1. + eta + zeta)*(-1. + qsi + zeta))/(4.*((1. - zeta)*(1. - zeta)*(1. - zeta)*(1. - zeta))));
    phi(1,0) = (1. - 2.*zeta)*(1. - zeta)*(-0.25*(1. - eta*eta/((1. - zeta)*(1. - zeta)))*(1. - qsi*qsi/((1. - zeta)*(1. - zeta))) + (eta*qsi*(1. + qsi - zeta)*(-1. + eta + zeta))/(4.*((1. - zeta)*(1. - zeta)*(1. - zeta)*(1. - zeta))));
    phi(2,0) = (-0.25*(1. - eta*eta/((1. - zeta)*(1. - zeta)))*(1. - qsi*qsi/((1. - zeta)*(1. - zeta))) + (eta*qsi*(1. + eta - zeta)*(1. + qsi - zeta))/(4.*((1. - zeta)*(1. - zeta)*(1. - zeta)*(1. - zeta))))*(1. - 2.*zeta)*(1. - zeta);
    phi(3,0) = (1. - 2.*zeta)*(1. - zeta)*(-0.25*(1. - eta*eta/((1. - zeta)*(1. - zeta)))*(1. - qsi*qsi/((1. - zeta)*(1. - zeta))) + (eta*qsi*(1. + eta - zeta)*(-1. + qsi + zeta))/(4.*((1. - zeta)*(1. - zeta)*(1. - zeta)*(1. - zeta))));
    phi(4,0) = zeta*(-1. + 2.*zeta);
    phi(5,0) = (1. - 2.*zeta)*(1. - zeta)*(0.5*(1. - eta*eta/((1. - zeta)*(1. - zeta)))*(1. - qsi*qsi/((1. - zeta)*(1. - zeta))) + (eta*(1. - qsi*qsi/((1. - zeta)*(1. - zeta)))*(-1. + eta + zeta))/(2.*((1. - zeta)*(1. - zeta))));
    phi(6,0) = (0.5*(1. - eta*eta/((1. - zeta)*(1. - zeta)))*(1. - qsi*qsi/((1. - zeta)*(1. - zeta))) + (qsi*(1. - eta*eta/((1. - zeta)*(1. - zeta)))*(1. + qsi - zeta))/(2.*((1. - zeta)*(1. - zeta))))*(1. - 2.*zeta)*(1. - zeta);
    phi(7,0) = (0.5*(1. - eta*eta/((1. - zeta)*(1. - zeta)))*(1. - qsi*qsi/((1. - zeta)*(1. - zeta))) + (eta*(1. - qsi*qsi/((1. - zeta)*(1. - zeta)))*(1. + eta - zeta))/(2.*((1. - zeta)*(1. - zeta))))*(1. - 2.*zeta)*(1. - zeta);
    phi(8,0) = (1. - 2.*zeta)*(1. - zeta)*(0.5*(1. - eta*eta/((1. - zeta)*(1. - zeta)))*(1. - qsi*qsi/((1. - zeta)*(1. - zeta))) + (qsi*(1. - eta*eta/((1. - zeta)*(1. - zeta)))*(-1. + qsi + zeta))/(2.*((1. - zeta)*(1. - zeta))));
    phi(9,0) = (zeta*(-1. + eta + zeta)*(-1. + qsi + zeta))/(1. - zeta);
    phi(10,0) = -(((1. + qsi - zeta)*zeta*(-1. + eta + zeta))/(1. - zeta));
    phi(11,0) = ((1. + eta - zeta)*(1. + qsi - zeta)*zeta)/(1. - zeta);
    phi(12,0) = -(((1. + eta - zeta)*zeta*(-1. + qsi + zeta))/(1. - zeta));
    
    
    dphi(0,0) = (0.5*(-0.5 + 1.*zeta)*(1.*(eta*eta) + eta*(-1. + 2.*qsi + 1.*zeta) + qsi*(-2. + 2.*zeta)))/((-1. + zeta)*(-1. + zeta));
    dphi(1,0) = (1.*(-0.5 + 1.*zeta)*(qsi*(-0.5 + 0.5*qsi + 0.5*zeta) + eta*(-1. + 1.*qsi + 1.*zeta)))/((-1. + zeta)*(-1. + zeta));
    dphi(2,0) = ((qsi*qsi)*(0.25 - 0.25*zeta) - 1.*((-1. + 1.*zeta)*(-1. + 1.*zeta)*(-1. + 1.*zeta))*(-0.75 + 1.*zeta) +
                 (eta*eta)*(0.25 + (-0.25 - 0.5*qsi)*zeta) + eta*qsi*(0.25 + (-0.25 - 0.5*qsi)*zeta))/((-1. + 1.*zeta)*(-1. + 1.*zeta)*(-1. + 1.*zeta));
    
    dphi(0,1) = (-0.5*(-0.5 + 1.*zeta)*(1.*(eta*eta) + qsi*(2. - 2.*zeta) + eta*(-1. - 2.*qsi + 1.*zeta)))/((-1. + zeta)*(-1. + zeta));
    dphi(1,1) = (1.*(-0.5 + 1.*zeta)*(qsi*(0.5 + 0.5*qsi - 0.5*zeta) + eta*(-1. - 1.*qsi + 1.*zeta)))/((-1. + zeta)*(-1. + zeta));
    dphi(2,1) = ((qsi*qsi)*(0.25 - 0.25*zeta) - 1.*((-1. + 1.*zeta)*(-1. + 1.*zeta)*(-1. + 1.*zeta))*(-0.75 + 1.*zeta) +
                 eta*qsi*(-0.25 + (0.25 - 0.5*qsi)*zeta) + (eta*eta)*(0.25 + (-0.25 + 0.5*qsi)*zeta))/((-1. + 1.*zeta)*(-1. + 1.*zeta)*(-1. + 1.*zeta));
    
    dphi(0,2) = (0.5*(-0.5 + 1.*zeta)*(-1.*(eta*eta) + eta*(-1. - 2.*qsi + 1.*zeta) + qsi*(-2. + 2.*zeta)))/((-1. + zeta)*(-1. + zeta));
    dphi(1,2) = (1.*(-0.5 + 1.*zeta)*(qsi*(-0.5 - 0.5*qsi + 0.5*zeta) + eta*(-1. - 1.*qsi + 1.*zeta)))/((-1. + zeta)*(-1. + zeta));
    dphi(2,2) = ((qsi*qsi)*(0.25 - 0.25*zeta) - 1.*((-1. + 1.*zeta)*(-1. + 1.*zeta)*(-1. + 1.*zeta))*(-0.75 + 1.*zeta) + (eta*eta)*(0.25 + (-0.25 + 0.5*qsi)*zeta) + eta*qsi*(0.25 + (-0.25 + 0.5*qsi)*zeta))/((-1. + zeta)*(-1. + zeta)*(-1. + zeta));
    
    dphi(0,3) = (-0.5*(-0.5 + 1.*zeta)*(-1.*(eta*eta) + qsi*(2. - 2.*zeta) + eta*(-1. + 2.*qsi + 1.*zeta)))/((-1. + zeta)*(-1. + zeta));
    dphi(1,3) = (1.*(-0.5 + 1.*zeta)*(qsi*(0.5 - 0.5*qsi - 0.5*zeta) + eta*(-1. + 1.*qsi + 1.*zeta)))/((-1. + zeta)*(-1. + zeta));
    dphi(2,3) = ((qsi*qsi)*(0.25 - 0.25*zeta) - 1.*((-1. + 1.*zeta)*(-1. + 1.*zeta)*(-1. + 1.*zeta))*(-0.75 + 1.*zeta) +
                 (eta*eta)*(0.25 + (-0.25 - 0.5*qsi)*zeta) + eta*qsi*(-0.25 + (0.25 + 0.5*qsi)*zeta))/((-1. + 1.*zeta)*(-1. + 1.*zeta)*(-1. + 1.*zeta));
    
    dphi(0,4) = 0.0;
    dphi(1,4) = 0.0;
    dphi(2,4) = -1.0 + 4.0*zeta;
    
    
    dphi(0,5) = (-2.*qsi*(-0.5 + 1.*zeta)*(1.*((1. - 1.*zeta)*(1. - 1.*zeta)) + eta*(-1. + 1.*zeta)))/((-1. + zeta)*(-1. + zeta)*(-1. + zeta));
    dphi(1,5) = (1.*(-0.5 + 1.*zeta)*((qsi*qsi)*(1. - 1.*zeta) + 1.*((-1. + zeta)*(-1. + zeta)*(-1. + zeta))))/((-1. + zeta)*(-1. + zeta)*(-1. + zeta));
    dphi(2,5) = (0.5*(qsi*qsi)*((1. - 1.*zeta)*(1. - 1.*zeta)) +
                 2.0*((1.0 - 1.*zeta)*(1.0 - 1.*zeta)*(1.0 - 1.*zeta)*(1.0 - 1.*zeta))*(-0.75 + 1.*zeta) +
                 eta*(1. + zeta*(-4. + (qsi*qsi)*(-1. + 1.*zeta) + zeta*(6. + zeta*(-4.0 + 1.*zeta)))))/((-1. + zeta)*(-1. + zeta)*(-1. + zeta)*(-1. + zeta));
    
    
    dphi(0,6) = (1.*(-0.5 + 1.*zeta)*((eta*eta)*(-1. + 1.*zeta) - 1.*((-1. + 1.*zeta)*(-1. + 1.*zeta)*(-1. + 1.*zeta))))/((-1. + zeta)*(-1. + zeta)*(-1. + zeta));
    dphi(1,6) = (2.*eta*(-0.5 + 1.*zeta)*(-1.*((1. - 1.*zeta)*(1. - 1.*zeta)) + qsi*(-1. + 1.*zeta)))/((-1. + zeta)*(-1. + zeta)*(-1. + zeta));
    dphi(2,6) = (2.0*((1.0 - 1.*zeta)*(1.0 - 1.*zeta)*(1.0 - 1.*zeta)*(1.0 - 1.*zeta))*(-0.75 + 1.*zeta) +
                 (eta*eta)*(0.5 + (-1. + qsi*(1. - 1.*zeta) + 0.5*zeta)*zeta) +
                 qsi*(-1. + zeta*(4. + zeta*(-6. + (4.0 - 1.*zeta)*zeta))))/((-1. + zeta)*(-1. + zeta)*(-1. + zeta)*(-1. + zeta));
    
    dphi(0,7) = (2.*qsi*(-0.5 + 1.*zeta)*(-1.*((1. - 1.*zeta)*(1. - 1.*zeta)) + eta*(-1. + 1.*zeta)))/((-1. + zeta)*(-1. + zeta)*(-1. + zeta));
    dphi(1,7) = (1.*(-0.5 + 1.*zeta)*((qsi*qsi)*(-1. + 1.*zeta) - 1.*((-1. + 1.*zeta)*(-1. + 1.*zeta)*(-1. + 1.*zeta))))/((-1. + zeta)*(-1. + zeta)*(-1. + zeta));
    dphi(2,7) = (0.5*(qsi*qsi)*((1. - 1.*zeta)*(1. - 1.*zeta)) + 2.0*((1.0 - 1.*zeta)*(1.0 - 1.*zeta)*(1.0 - 1.*zeta)*(1.0 - 1.*zeta))*
                 (-0.75 + 1.*zeta) + eta*(-1. + zeta*(4. + (qsi*qsi)*(1. - 1.*zeta) + zeta*(-6. + (4.0 - 1.*zeta)*zeta))))/((-1. + zeta)*(-1. + zeta)*(-1. + zeta)*(-1. + zeta));
    
    dphi(0,8) = (1.*(-0.5 + 1.*zeta)*((eta*eta)*(1. - 1.*zeta) + 1.*((-1. + 1.*zeta)*(-1. + 1.*zeta)*(-1. + 1.*zeta))))/((-1. + zeta)*(-1. + zeta)*(-1. + zeta));
    dphi(1,8) = (-2.*eta*(-0.5 + 1.*zeta)*(-1. + 1.*qsi + 1.*zeta))/((1. - 1.*zeta)*(1. - 1.*zeta));
    dphi(2,8) = (2.0*((1.0 - 1.*zeta)*(1.0 - 1.*zeta)*(1.0 - 1.*zeta)*(1.0 - 1.*zeta))*(-0.75 + 1.*zeta) +
                 (eta*eta)*(0.5 + zeta*(-1. + 0.5*zeta + qsi*(-1. + 1.*zeta))) +
                 qsi*(1. + zeta*(-4. + zeta*(6. + zeta*(-4.0 + 1.*zeta)))))/((1.0 - 1.*zeta)*(1.0 - 1.*zeta)*(1.0 - 1.*zeta)*(1.0 - 1.*zeta));
    
    dphi(0,9) = (-1.*zeta*(-1. + eta + zeta))/(-1. + zeta);
    dphi(1,9) = (-1.*zeta*(-1. + qsi + zeta))/(-1. + zeta);
    dphi(2,9) = 1.0 - 1.*eta - 1.*qsi + (1.*eta*qsi)/((1. - 1.*zeta)*(1. - 1.*zeta)) - 2.*zeta;
    
    dphi(0,10) = (1.*zeta*(-1. + eta + zeta))/(-1. + zeta);
    dphi(1,10) = (-1.*zeta*(-1. - 1.*qsi + zeta))/(-1. + zeta);
    dphi(2,10) = 1.0 - 1.*eta + 1.*qsi - (1.*eta*qsi)/((1. - 1.*zeta)*(1. - 1.*zeta)) - 2.*zeta;
    
    dphi(0,11) = (1.*zeta*(-1. - 1.*eta + zeta))/(-1. + zeta);
    dphi(1,11) = (1.*zeta*(-1. - 1.*qsi + zeta))/(-1. + zeta);
    dphi(2,11) = 1.0 + 1.*eta + 1.*qsi + (1.*eta*qsi)/((1. - 1.*zeta)*(1. - 1.*zeta)) - 2.*zeta;
    
    dphi(0,12) = (-1.*zeta*(-1. - 1.*eta + zeta))/(-1. + zeta);
    dphi(1,12) = (1.*zeta*(-1. + qsi + zeta))/(-1. + zeta);
    dphi(2,12) = 1.0 + 1.*eta - 1.*qsi - (1.*eta*qsi)/((1. - 1.*zeta)*(1. - 1.*zeta)) - 2.*zeta;
    
}


template<class T>
void TPZQuadraticPyramid::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
    
    TPZFNMatrix<13,T> phi(NNodes,1);
    TPZFNMatrix<39,T> dphi(3,NNodes);
    TShape(loc,phi,dphi);
    int space = nodes.Rows();
    
    for(int i = 0; i < space; i++) {
        x[i] = 0.0;
        for(int j = 0; j < NNodes; j++) {
            x[i] += phi(j,0)*nodes.GetVal(i,j);
        }
    }
    
}

template<class T>
void TPZQuadraticPyramid::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
    
    gradx.Resize(3,3);
    gradx.Zero();
    int nrow = nodes.Rows();
    int ncol = nodes.Cols();
#ifdef PZDEBUG
    if(nrow != 3 || ncol  != 13){
        std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
        std::cout << "nodes matrix must be 3x13." << std::endl;
        DebugStop();
    }
    
#endif
    TPZFNMatrix<13,T> phi(NNodes,1);
    TPZFNMatrix<39,T> dphi(3,NNodes);
    TShape(loc,phi,dphi);
    for(int i = 0; i < NNodes; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            gradx(j,0) += nodes.GetVal(j,i)*dphi(0,i);
            gradx(j,1) += nodes.GetVal(j,i)*dphi(1,i);
            gradx(j,2) += nodes.GetVal(j,i)*dphi(2,i);
            
        }
    }
    
}


/**
 * Creates a geometric element according to the type of the father element
 */

TPZGeoEl *TPZQuadraticPyramid::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                                TPZVec<int64_t>& nodeindexes,
                                                int matid,
                                                int64_t& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

TPZGeoEl *TPZQuadraticPyramid::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) 
{
	int ns = orig->NSideNodes(side);
	TPZManVector<int64_t> nodeindices(ns);
	int in;
	for(in=0; in<ns; in++)
	{
		nodeindices[in] = orig->SideNodeIndex(side,in);
	}
	int64_t index;
	
	TPZGeoMesh *mesh = orig->Mesh();
	MElementType type = orig->Type(side);
	
	TPZGeoEl *newel = mesh->CreateGeoBlendElement(type, nodeindices, bc, index);
	TPZGeoElSide me(orig,side);
	TPZGeoElSide newelside(newel,newel->NSides()-1);
	
	newelside.InsertConnectivity(me);
	newel->Initialize();
	
	return newel;
}

/// create an example element based on the topology
/* @param gmesh mesh in which the element should be inserted
 @param matid material id of the element
 @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
 @param size (in) size of space where the element should be created
 */
#include "tpzchangeel.h"

void TPZQuadraticPyramid::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
{
    TPZManVector<REAL,3> co(3),shift(3),scale(3);
    TPZManVector<int64_t,4> nodeindexes(NCornerNodes);
    for (int i=0; i<3; i++) {
        scale[i] = size[i]/3.;
        shift[i] = size[i]/2.+lowercorner[i];
    }
    
    for (int i=0; i<NCornerNodes; i++) {
        ParametricDomainNodeCoord(i, co);
        co.Resize(3,0.);
        for (int j=0; j<3; j++) {
            co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
        }
        nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
    }
    int64_t index;
    CreateGeoElement(gmesh, EPiramide, nodeindexes, matid, index);
    TPZGeoEl *gel = gmesh.Element(index);
    int nsides = gel->NSides();
    for (int is=0; is<nsides; is++) {
        gel->SetSideDefined(is);
    }
    gel = TPZChangeEl::ChangeToQuadratic(&gmesh, index);
    for (int node = gel->NCornerNodes(); node < gel->NNodes(); node++) {
        TPZManVector<REAL,3> co(3);
        gel->NodePtr(node)->GetCoordinates(co);
        for (int i=0; i<3; i++) {
            co[i] += (0.2*rand())/RAND_MAX - 0.1;
        }
        gel->NodePtr(node)->SetCoord(co);
    }
}

//void TPZQuadraticPyramid::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
//{
//    if(node > this->NNodes)
//    {
//        DebugStop();
//    }
//    nodeCoord.Resize(Dimension, 0.);
//    switch (node) {
//        case (0):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (1):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (2):
//        {
//            nodeCoord[0] = 1.;
//            nodeCoord[1] = 1.;
//            nodeCoord[2] = 0.;
//            break;
//        }
//        case (3):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] =  1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (4):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 0.;
//            nodeCoord[2] = 1.;
//            break;
//        }
//        case (5):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (6):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (7):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] =  1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (8):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (9):
//        {
//            nodeCoord[0] = -0.5;
//            nodeCoord[1] = -0.5;
//            nodeCoord[2] =  0.5;
//            break;
//        }
//        case (10):
//        {
//            nodeCoord[0] =  0.5;
//            nodeCoord[1] = -0.5;
//            nodeCoord[2] =  0.5;
//            break;
//        }
//        case (11):
//        {
//            nodeCoord[0] =  0.5;
//            nodeCoord[1] =  0.5;
//            nodeCoord[2] =  0.5;
//            break;
//        }
//        case (12):
//        {
//            nodeCoord[0] = -0.5;
//            nodeCoord[1] =  0.5;
//            nodeCoord[2] =  0.5;
//            break;
//        }
//        default:
//        {
//            DebugStop();
//            break;
//        }
//    }
//}

///CreateGeoElement -> TPZQuadraticPyramid

int TPZQuadraticPyramid::ClassId() const{
    return Hash("TPZQuadraticPyramid") ^ TPZNodeRep<13,pztopology::TPZPyramid>::ClassId() << 1;
}

template class TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticPyramid>>;
template class pzgeom::TPZNodeRep<13,TPZQuadraticPyramid>;

namespace pzgeom {
    template void TPZQuadraticPyramid::X(const TPZFMatrix<REAL>&, TPZVec<REAL>&, TPZVec<REAL>&);
    template void TPZQuadraticPyramid::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc, TPZFMatrix<REAL> &gradx);

#ifdef _AUTODIFF
    template void TPZQuadraticPyramid::X(const TPZFMatrix<REAL>&, TPZVec<Fad<REAL> >&, TPZVec<Fad<REAL> >&);
    template void TPZQuadraticPyramid::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<Fad<REAL> > &loc, TPZFMatrix<Fad<REAL> > &gradx);
#endif

}

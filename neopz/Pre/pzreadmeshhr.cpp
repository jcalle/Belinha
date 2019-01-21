/**
 * @file
 * @brief Contains the implementation of the TPZReadMeshHR methods. 
 */

#include "pzreadmeshhr.h"

#include "pzgeoel.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzelasmat.h"
#include "pzmat2dlin.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "pzcompel.h"
#include "pzlog.h"

#include <iostream>
#include <sstream>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.pre.tpzmeshreaderhr"));
#endif

TPZReadMeshHR::TPZReadMeshHR(const char* infInputFile): TPZReadMesh(infInputFile)
{}

TPZReadMeshHR::~TPZReadMeshHR()
{}

TPZCompMesh* TPZReadMeshHR::ReadMesh()
{
	std::string numberOf;
	TPZGeoMesh *gmesh = readGeoMesh();
	TPZCompMesh *cmesh = new TPZCompMesh (gmesh);
	removeComents (numberOf);
	int nmat = atoi (numberOf.c_str());
	ReadMaterials(nmat,*cmesh);
	removeComents (numberOf);
	int nbcs = atoi (numberOf.c_str());
	ReadBCs(nbcs,*cmesh);
	cmesh->AutoBuild();
	return cmesh;
}

void TPZReadMeshHR::removeComents (std::string &NumberOf)
{
	if (!fInputFile) return;
	while (fInputFile)
	{
		char buf[512];
		fInputFile.getline(buf,512);
		std::string aux (buf);
		int pos = aux.find (":");
		if ((pos > -1 && pos < (int)aux.size()) || !aux.size()) continue; //encontrou um coment�io na linha
		NumberOf = aux;
		break;
	}
}

void TPZReadMeshHR::ReadNodes (int64_t NNos, TPZGeoMesh & GMesh)
{
	int64_t i, id;
	int c;
	TPZVec <REAL> coord(3,0.);
	for (i=0;i<NNos;i++)
	{
		fInputFile >> id;
		for (c=0;c<3;c++) fInputFile >> coord[c];
		int64_t index = GMesh.NodeVec().AllocateNewElement();
		GMesh.NodeVec()[index] = TPZGeoNode(id,coord,GMesh);
	}
}

void TPZReadMeshHR::ReadElements (int64_t NElem, TPZGeoMesh & GMesh)
{
	int64_t i, c, id;
	int type,matid;
	for (i=0;i<NElem;i++)
	{
		fInputFile >> id >> type >> matid;
		int nnos = -1;
		MElementType etype;
		switch (type)
		{
			case (0):
				nnos = 1;
				etype = EPoint;
				break;
			case (1):
				nnos = 2;
				etype = EOned;
				break;
			case (2):
				nnos = 3;
				etype = ETriangle;
				break;
			case (3):
				nnos = 4;
				etype = EQuadrilateral;
				break;
			case (4):
				nnos = 4;
				etype = ETetraedro;
				break;
			case (5):
				nnos = 5;
				etype = EPiramide;
				break;
			case (6):
				nnos = 6;
				etype = EPrisma;
				break;
			case (7):
				nnos = 8;
				etype = ECube;
				break;
			default:
				std::stringstream sout;
#ifndef WINDOWS
				sout << __PRETTY_FUNCTION__;
#endif
				sout << "Nao sei que elemento " << type << " eh esse indicado para o elemento " << id;
#ifdef LOG4CXX
				LOGPZ_WARN (logger,sout.str().c_str());
#else
				std::cout << sout.str().c_str() << std::endl;
#endif
				continue;
		}
		TPZVec<int64_t> corneridx(nnos,-1);
		int64_t id, idx;
		for (c=0;c<nnos;c++)
		{
			fInputFile >> id;
			idx = GetNodeIndex(&GMesh,id);
			corneridx [c] = idx;
		}
		GMesh.CreateGeoElement(etype,corneridx,matid,id,1);
	}
}

void TPZReadMeshHR::ReadMaterials(int NMat, TPZCompMesh & CMesh) {
    int i;
    int id, classId;
    for (i = 0; i < NMat; i++) {
        fInputFile >> id >> classId;
        if (classId == ClassIdOrHash<TPZElasticityMaterial>()) {
            double e, nu, px, py;
            fInputFile >> e >> nu >> px >> py;
            TPZMaterial * mat = new TPZElasticityMaterial(id, e, nu, px, py, 0);
            CMesh.InsertMaterialObject(mat);
            break;
        } else if (classId == ClassIdOrHash<TPZMat2dLin>()) {
            TPZMat2dLin *mat2d = new TPZMat2dLin(id);
            int nstate;
            fInputFile >> nstate;

            int ist, jst;
            TPZFMatrix<STATE> xk(nstate, nstate, 1.), xc(nstate, nstate, 0.), xf(nstate, 1, 0.);
            //xk
            for (ist = 0; ist < nstate; ist++) {
                fInputFile >> xf(ist, 0);
            }
            //xc
            for (ist = 0; ist < nstate; ist++) {
                for (jst = 0; jst < nstate; jst++) {
                    fInputFile >> xc(ist, jst);
                }
            }
            //xf
            for (ist = 0; ist < nstate; ist++) {
                for (jst = 0; jst < nstate; jst++) {
                    fInputFile >> xk(ist, jst);
                }
            }
            mat2d->SetMaterial(xk, xc, xf);
            TPZMaterial * mat = mat2d;
            CMesh.InsertMaterialObject(mat);
            break;
        } else if (classId == ClassIdOrHash<TPZMatPoisson3d>()) {
            TPZMatPoisson3d *mat3d = new TPZMatPoisson3d(id, 3);
            int nstate;
            fInputFile >> nstate;
            int ist;
            REAL diff, conv;
            TPZVec<REAL> dir(3, 0.);
            fInputFile >> diff >> conv;
            for (ist = 0; ist < 3; ist++) {
                fInputFile >> dir[ist];
            }
            mat3d->SetParameters(diff, conv, dir);
            TPZMaterial * mat = mat3d;
            CMesh.InsertMaterialObject(mat);
            break;
        } else {
            std::stringstream sout;
            sout << "Could not identify material of type: " << classId
                    << " check material identifier for material " << i;
#ifdef LOG4CXX
            LOGPZ_FATAL(logger, sout.str().c_str());
#endif
            std::cout << sout.str().c_str() << std::endl;
            DebugStop();
            break;
        }
    }
}

void TPZReadMeshHR::ReadBCs (int NMat, TPZCompMesh & CMesh)
{
	int i;
	int id, type;
	if(!CMesh.MaterialVec().size())
	{
		std::stringstream sout;
#ifndef WINDOWS
		sout << __PRETTY_FUNCTION__ << " no materials " << std::endl;
#endif
		sout << "\tN� encontrei material na malha!";
#ifdef LOG4CXX
		LOGPZ_ERROR (logger, sout.str().c_str());
#else
		std::cout << sout.str().c_str() << std::endl;
#endif
		return;
	}
	std::map<int, TPZMaterial * >::iterator matit = CMesh.MaterialVec().begin();
	TPZMaterial * mat = matit->second;
	if(!mat)
	{
		std::cout << " empty material " << std::endl;
		return;
	}
	
	
	for (i=0;i<NMat;i++)
	{
		fInputFile >> id >> type;
		TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
		fInputFile >> val1 (0,0) >> val1(0,1) >> val1(0,2)
        >> val1 (1,0) >> val1(1,1) >> val1(1,2)
        >> val1 (2,0) >> val1(2,1) >> val1(2,2);
		fInputFile >> val2(0,0) >> val2(1,0) >> val2(2,0);
		TPZMaterial * bnd;
		bnd = mat->CreateBC (mat,id,type,val1,val2);
		CMesh.InsertMaterialObject(bnd);
	}
}

int64_t TPZReadMeshHR::GetNodeIndex(TPZGeoMesh *GMesh,int64_t Id)
{
	TPZAdmChunkVector<TPZGeoNode> nodeVec = GMesh->NodeVec();
	int64_t vid,index,size = nodeVec.NElements();
	if (Id < size) index = Id;
	else index = size - 1;
	vid = nodeVec[index].Id();
	if (vid == Id) return index;
	
	int64_t i;
	for (i=index-1;i>-1;i--)
	{
		vid = nodeVec[i].Id();
		if (vid == Id) return i;
	}
	std::stringstream sout;
#ifndef WINDOWS
	sout << __PRETTY_FUNCTION__;
#endif
	sout << " N�" << Id << " n� encontrado!";
#ifdef LOG4CXX
	LOGPZ_WARN (logger, sout.str().c_str());
#else
	std::cout << sout.str().c_str() << std::endl;
#endif
	
	return -1;
}

TPZGeoMesh * TPZReadMeshHR::readGeoMesh()
{
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	std::string numberOf;
	removeComents (numberOf);
	int64_t nnos = atoi (numberOf.c_str());
	ReadNodes (nnos, *gmesh);
	removeComents (numberOf);
	int64_t nelem = atoi (numberOf.c_str());
	ReadElements(nelem, *gmesh);
	gmesh->BuildConnectivity();
	return gmesh;
}

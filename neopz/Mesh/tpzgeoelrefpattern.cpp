/**
 * @file
 * @brief Contains the implementation of the TPZGeoElRefPattern methods.
 */

#include "tpzgeoelrefpattern.h"
#include "tpzgeoelrefpattern.h.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "TPZGeoElement.h"
#include "TPZRefPattern.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "TPZRefPatternDataBase.h"

using namespace pzgeom;
using namespace pzshape;

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzgeoelrefpattern"));
#endif
TPZGeoEl *CreateGeoElementPattern(TPZGeoMesh &mesh, MElementType type,
                                  TPZVec<int64_t>& nodeindexes,
                                  int matid,
                                  int64_t& index)

{
	switch( type ){
		case 0://point
		{
			TPZGeoEl * gel =
			new TPZGeoElRefPattern<TPZGeoPoint>(nodeindexes, matid, mesh, index);
			return gel;
		}
		case 1://line
		{
			TPZGeoEl *gel =
			new TPZGeoElRefPattern< TPZGeoLinear >
			(nodeindexes, matid, mesh, index);
			return gel;
		}
		case 2://triangle
		{
			TPZGeoEl *gel =
			new TPZGeoElRefPattern< TPZGeoTriangle >
			(nodeindexes, matid, mesh, index);
			return gel;
		}
		case 3://quadrilatera
		{
			TPZGeoEl* gel =
			new TPZGeoElRefPattern< TPZGeoQuad >
			(nodeindexes, matid, mesh, index);
			return gel;
		}
		case 4://tetraedra
		{
			TPZGeoEl*gel =
			new TPZGeoElRefPattern< TPZGeoTetrahedra >
			(nodeindexes, matid, mesh, index);
			return gel;
		}
		case 5://pyramid
		{
			TPZGeoEl *gel =
			new TPZGeoElRefPattern< TPZGeoPyramid >
			(nodeindexes, matid, mesh, index);
			return gel;
		}
		case 6://prism
		{
			TPZGeoEl*gel =
			new TPZGeoElRefPattern< TPZGeoPrism >
			(nodeindexes, matid, mesh, index);
			return gel;
		}
		case 7://cube
		{
			TPZGeoEl*gel =
			new TPZGeoElRefPattern< TPZGeoCube >
			(nodeindexes, matid, mesh, index);
			return gel;
		}
		default:
		{
			PZError << "TPZGeoMesh::CreateGeoElementRefPattern type element not exists:"
			<< " type = " << type << std::endl;
			return NULL;
		}
	}
}



/** registration of the class in the TPZRestoreClass */

#ifndef BORLAND
template class TPZRestoreClass< TPZGeoElRefPattern<TPZGeoCube>>;
template class TPZRestoreClass< TPZGeoElRefPattern<TPZGeoLinear>>;
template class TPZRestoreClass<TPZGeoElRefPattern<TPZGeoQuad>>;
template class TPZRestoreClass<TPZGeoElRefPattern<TPZGeoTriangle>>;
template class TPZRestoreClass< TPZGeoElRefPattern<TPZGeoPrism>>;
template class TPZRestoreClass< TPZGeoElRefPattern<TPZGeoTetrahedra>>;
template class TPZRestoreClass< TPZGeoElRefPattern<TPZGeoPyramid>>;
template class TPZRestoreClass< TPZGeoElRefPattern<TPZGeoPoint>>;
#endif


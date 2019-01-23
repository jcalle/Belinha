/**
 * @file
 * @brief Contains the TPZReadMesh class which implements the interface for build a computational mesh from a file.
 */
/*****************************************************************************
 * O conte�do desse arquivo � de propriedade do LabMeC-DES-FEC-UNICAMP e do
 * CENPES-Petrobras. 
 * O uso de qualquer parte ou do todo est� condicionado � expressa autoriza��o
 * dos propriet�rios.
 *****************************************************************************/

#ifndef PZREADMESH_H
#define PZREADMESH_H

#include <fstream>
class TPZCompMesh;

/**
 * @ingroup pre
 * @brief Virtual class that implements the interface for build a computational mesh from a file. \ref pre "Getting Data"
 * @author Edimar Cesar Rylo
 * @since September, 2006
 */
class TPZReadMesh
{
public:
	/**
	 * @brief Default constructor
	 * @param inFile [in] contains a full path to the input file
	 */
	TPZReadMesh(const char * inFile);
	
	/** @brief Default destructor */
	virtual ~TPZReadMesh() {
	}
	
	virtual TPZCompMesh *ReadMesh() = 0;
	
protected:
	/** @brief Input file */
	std::ifstream fInputFile;
};

inline TPZReadMesh::TPZReadMesh(const char * FileFullPath) : fInputFile(FileFullPath)
{
}

#endif

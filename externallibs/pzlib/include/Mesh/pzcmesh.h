/**
 * @file
 * @brief Contains declaration of TPZCompMesh class which is a repository for computational elements, nodes and material objects.
 */

#ifndef PZCMESHHPP
#define PZCMESHHPP


#include <stddef.h>               // for NULL
#include <iostream>               // for operator<<, string, cout, ostream
#include <map>                    // for map
#include <set>                    // for set
#include "pzadmchunk.h"           // for TPZAdmChunkVector
#include "pzblock.h"              // for TPZBlock
#include "TPZChunkVector.h"              // for TPZChunkVector
#include "pzconnect.h"            // for TPZConnect
#include "pzcreateapproxspace.h"  // for TPZCreateApproximationSpace
#include "pzgmesh.h"              // for TPZGeoMesh
#include "pzmatrix.h"             // for TPZFMatrix, TPZMatrix
#include "pzreal.h"               // for STATE, REAL
#include "TPZSavable.h"          // for TPZSavable
#include "pzstack.h"              // for TPZStack
#include "pzvec.h"                // for TPZVec
#include "tpzautopointer.h"       // for TPZAutoPointer
#include "pzcheckgeom.h"		  // for TPZCheckGeom
#include <functional>
#include "pzcompel.h"

class TPZGeoEl;
class TPZMaterial;
class TPZStream;
template <class TVar> class TPZTransfer;


/**
 * @brief Implements computational mesh. \ref CompMesh "Computational Mesh"
 * @ingroup CompMesh
*/ 
/**
 * The computational mesh is a repository for computational elements, nodes and
 * material objects \n
 * The computational mesh also contains the current solution of the mesh and an
 * elementwise solution vector \n
 * The data structure of this object is rather simple
 */
class TPZCompMesh : public virtual TPZSavable {
	
protected:
	/** @brief Geometric grid to which this grid refers */
	TPZGeoMesh	*fReference;
    
    /** @brief Autopointer to the geometric mesh used in case the user has passed an autopointer */
    TPZAutoPointer<TPZGeoMesh> fGMesh;
	
	/** @brief Grid name for model identification */
	std::string fName;
	
	/** @brief List of pointers to elements */
	TPZAdmChunkVector<TPZCompEl *>		fElementVec;	
	/** @brief List of pointers to nodes */
	TPZAdmChunkVector<TPZConnect>			fConnectVec;
	
	/** @brief Map of pointers to materials */
	std::map<int, TPZMaterial * >	fMaterialVec;
	/** @brief Block structure of the solution vector ???? */
	//TPZBlock<REAL>		fSolutionBlock;
	TPZBlock<STATE>		fSolutionBlock;
	
	/** @brief Solution vector */
	//TPZFMatrix<REAL>	fSolution;
	TPZFMatrix<STATE>	fSolution;
	
	/** @brief Block structure to right construction of the stiffness matrix and load vector */
	//TPZBlock<REAL>		fBlock;
	TPZBlock<STATE>		fBlock;
	
	/** @brief Solution vectors organized by element */
	TPZFMatrix<STATE> fElementSolution;
	
	/* @brief set the dimension of the simulation or the model */
	int fDimModel;
	
	/** @brief Default order for all elements of this mesh */
	int fDefaultOrder;
    
    /** @brief The object which defines the type of space being created */
    TPZCreateApproximationSpace fCreate;
    
    /// Modify the permute vector swapping the lagrangeq with maxeq and shifting the intermediate equations
    void ModifyPermute(TPZVec<int64_t> &permute, int64_t lagrangeq, int64_t maxeq);
    
    //quantity of meshes associated with a mesh multiphysics
	int64_t fNmeshes;

public:
	
	/**
	 * @brief Constructor from geometrical mesh
	 * @param gr pointer to geometrical reference mesh
	 */
	TPZCompMesh(TPZGeoMesh* gr=0);
    /** @brief Constructor based on an autopointer to a geometric mesh */
    TPZCompMesh(TPZAutoPointer<TPZGeoMesh> &gmesh);
	/** @brief Copy constructor */
	TPZCompMesh(const TPZCompMesh &copy);
    /** @brief copy the content of the mesh */
    TPZCompMesh &operator=(const TPZCompMesh &copy);
	
	/** @brief Simple Destructor */
	virtual ~TPZCompMesh();
	
	/**
	 * @brief This method will initiate the comparison between the current computational
	 * mesh and the mesh which is referenced by the geometric mesh
	 * @param var state variable number
	 * @param matname pointer to material name
	 */
	virtual REAL CompareMesh(int var, char *matname);
	/**
	 * @brief Gives all patches of the mesh
	 * @param grpatch - stack where will be inserted the geometric elements
	 */
	void GetRefPatches(std::set<TPZGeoEl *> &grpatch);
	/** @brief Gives the conects graphs */
	void GetNodeToElGraph(TPZVec<int64_t> &nodtoelgraph, TPZVec<int64_t> &nodtoelgraphinde,TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindexx);
	
    /** @brief Gives the element patch */
	void GetElementPatch(TPZVec<int64_t> nodtoelgraph, TPZVec<int64_t> nodtoelgraphindex, TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex,int64_t elind ,TPZStack<int64_t> &patch);
	
	/** @brief Set the mesh name */
	void SetName(const std::string &nm);
	
	/** @brief Set de dimension of the domain of the problem*/
	void SetDimModel(int dim) { fDimModel = dim; }
	
	/** @brief Sets the geometric reference mesh */
	void SetReference(TPZGeoMesh * gmesh);
	
    /** @brief Sets the geometric reference mesh */
    void SetReference(TPZAutoPointer<TPZGeoMesh> &gmesh);
    
	/** @brief Returns the dimension of the simulation*/
	int Dimension() const {return fDimModel;}

	/** @brief Returns the mesh name */
	std::string &Name() {return fName;}
	
	/** @brief Delete all the dynamically allocated data structures */
	void CleanUp();
	
	/**
	 * @name Access_Private_Data
	 * @brief Routines for access to private data
	 * @{
	 */
	
	/** @brief Number of connects allocated including free nodes */
	int64_t NConnects() const {return fConnectVec.NElements();}
	
	/** @brief Number of independent connect objects */
	int64_t NIndependentConnects();
	
	/** @brief Number of computational elements allocated */
	int64_t NElements() const {return fElementVec.NElements();}
	
	/** @brief Number of materials */
	size_t NMaterials() const {return fMaterialVec.size();}
	
	//@}
	
	/**
	 * @name Access_Internal_Data_Structures
	 * @brief Methods for accessing the internal data structures
	 * @{
	 */
	
	/** @brief Returns a reference to the element pointers vector */
	TPZAdmChunkVector<TPZCompEl *>   &ElementVec() { return fElementVec; }
    
    TPZCompEl * Element(int64_t iel)
    {
        return fElementVec[iel];
    }

    const TPZCompEl * Element(int64_t iel) const {
        return fElementVec[iel];
    }
	
	/** @brief Returns a reference to the element pointers vector */
	const TPZAdmChunkVector<TPZCompEl *>   &ElementVec() const { return fElementVec; }
	
	/** @brief Return a reference to the connect pointers vector */
	TPZAdmChunkVector<TPZConnect>    &ConnectVec() { return fConnectVec; }
    
	const TPZAdmChunkVector<TPZConnect> &ConnectVec() const { return fConnectVec; }
	
	/** @brief Returns a reference to the material pointers vector */
	std::map<int ,TPZMaterial * >	&MaterialVec() { return fMaterialVec; }

	/** @brief Returns a pointer to the geometrical mesh associated */
	TPZGeoMesh *Reference() const { return fReference; }
	
	/** @brief Access the block structure of the solution vector */
	//const TPZBlock<REAL> &Block() const { return fBlock;}
	const TPZBlock<STATE> &Block() const { return fBlock;}
	
	/** @brief Access the block structure of the solution vector */
	TPZBlock<STATE> &Block() { return fBlock;}
	
	/** @brief Access the solution vector */
	TPZFMatrix<STATE> &Solution(){ return fSolution;}
	
	/** @brief Access method for the element solution vectors */
	TPZFMatrix<STATE> &ElementSolution() { return fElementSolution;}
	
	/** @} */
	
	/**
	 * @name Modify_Mesh_Structure
	 * @brief Methods for modify mesh structure, materials, elements, node etc
	 * @{
	 */
	
	/**
	 * @brief Returns an index to a new connect
	 * @param nshape Number of function shapes
	 * @param nstate Number of variable state, depending of the material.
	 * @param order Order of approximation
	 */
	virtual  int64_t AllocateNewConnect(int nshape, int nstate, int order);
	
	/**
	 * @brief Returns an index to a new connect
	 * @param connect Connect from which attributes should be copied
	 */
	virtual  int64_t AllocateNewConnect(const TPZConnect &connect);
	
	/**
	 * @brief Insert a material object in the datastructure
	 * @param mat pointer to the material
	 */
	int InsertMaterialObject(TPZMaterial * mat);
	
	/** @brief Resequence the block object, remove unconnected connect objects and reset the dimension of the solution vector */
	void InitializeBlock();
	
	/** @brief Adapt the solution vector to new block dimensions */
	virtual void ExpandSolution();
	
	/** @} */
	
	/**
	 * @name Access_Solution
	 * @brief Methods for access and manipulates solution
	 * @{
	 */
	
	/** @brief Set a ith element solution, expanding the element-solution matrix if necessary */
	void SetElementSolution(int64_t i, TPZVec<STATE> &sol);
	
	/** @} */
	
	/**
	 * @name Print
	 * @brief Print methods
	 * @{
	 */
	
	/**
	 * @brief Prints mesh data
	 * @param out indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream & out = std::cout) const;
	
	/**
	 * @brief Print the solution by connect index
	 * @param out indicates the device where the data will be printed
	 */
	void ConnectSolution(std::ostream & out);

	/** @} */
	
	/**
	 * @brief Find the material with identity id
	 * @param id material id to be found
	 */
	TPZMaterial * FindMaterial(int id);
	
	/** @brief Map this grid in the geometric grid */
	void LoadReferences();
	
	/** @brief Compute the number of elements connected to each connect object */
	virtual void ComputeNodElCon();
	
	/** @brief Compute the number of elements connected to each connect object */
	virtual void ComputeNodElCon(TPZVec<int> &nelconnected) const;
	
	/** @brief Delete the nodes which have no elements connected to them */
	virtual void CleanUpUnconnectedNodes();

	/**
	 * @brief Computes the connectivity graph of the elements, as appropriate for the TPZRenumbering class
	 * @param elgraph stack of elements to create the grapho????
	 * @param elgraphindex graphos indexes vector
	 */
	void ComputeElGraph(TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex);

    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const;
    

	/**
	 * @name Submesh_methods
	 * @brief Methods required by submeshes
	 * @{
	 */
	
	/** @brief Get the father meshes stack */
	virtual TPZCompMesh *FatherMesh() const { return NULL; }
	
	/**
	 * @brief Makes a specified connection a internal mesh connection.
	 * @param local connection local number to be processed
	 */
	virtual void MakeInternal(int64_t local) {;}
	
	/** @brief Make all mesh connections internal mesh connections. Connects to an internal connection */
	virtual void MakeAllInternal() {;}
	
	/** @brief Returns the rootmesh who have the specified connection. */
	virtual TPZCompMesh* RootMesh (int64_t local) { return this; }
	
	/**
	 * @brief Transfer one element from a specified mesh to the current submesh.
	 * @param mesh pointer to the mesh whose the element from
	 * @param elindex element index to transfer
	 */
	virtual int64_t TransferElementFrom(TPZCompMesh *mesh, int64_t elindex) { return elindex; }
	
	/**
	 * @brief Transfer one element from a submesh to another mesh.
	 * @param mesh mesh pointer who will receive the element
	 * @param elindex element index to transfer
	 */
	virtual int64_t TransferElementTo(TPZCompMesh *mesh, int64_t elindex) { return elindex; }
	
	/**
	 * @brief Transfer one element form a submesh to another mesh.
	 * @param mesh mesh pointer who will receive the element
	 * @param elindex element index to transfer
	 */
	virtual int64_t TransferElement(TPZCompMesh *mesh, int64_t elindex) { return elindex; }
	
	/**
	 * @brief Put an local connection in the supermesh - Supermesh is one
	 * mesh who contains the analised submesh.
	 * @param local local index of the element to be trasnfered
	 * @param super pointer to the destination mesh
	 */
	virtual int64_t PutinSuperMesh (int64_t local, TPZCompMesh *super);
	
	/**
	 * @brief Get an external connection from the supermesh - Supermesh is one
	 * mesh who contains the analised submesh.
	 * @param superind index of the element to be trasnfered
	 * @param super pointer to the destination mesh
	 */
	virtual int64_t GetFromSuperMesh (int64_t superind, TPZCompMesh *super);
    
    /**
	 * @brief Gives the commom father mesh of the specified mesh and the current submesh.
	 * @param mesh pointer to other mesh whosw want to know the commom father mesh
	 */
	virtual TPZCompMesh * CommonMesh (TPZCompMesh *mesh);
	


	/** @} */

	int GetDefaultOrder()
	{
		return fDefaultOrder;
	}
	
	void SetDefaultOrder( int order )
	{
		fDefaultOrder = order;
	}

	/**
	 * @name SCIENTIFIC_ROUTINES
	 * @brief Scientific manipulating routines
	 * @{
	 */
	
	/** @brief This computes the number of equations associated with non-restrained nodes */
	int64_t NEquations();
	
	/** @brief This method computes the bandwidth of the system of equations */
	int BandWidth();
	
	/**
	 * @brief This method computes the skyline of the system of equations
	 * @param skyline vector where the skyline will be computed
	 */
	virtual void Skyline(TPZVec<int64_t> &skyline);
	
	/**
	 * @brief Assemble the vector with errors estimators
	 * @param estimator vector where will be assembled the errors
	 * @param errorid index for dual or wheeler estimator
	 */
	void AssembleError(TPZFMatrix<REAL> &estimator, int errorid);
	
	/**
	 * @brief Builds the transfer matrix from the current grid to the coarse grid
	 * @param coarsemesh grid for where the matrix will be transfered
	 * @param transfer transfer matrix between the current mesh and the coarse mesh
	 */
	void BuildTransferMatrix(TPZCompMesh &coarsemesh, TPZTransfer<STATE> &transfer);
	
	/** @brief To discontinuous elements */
	void BuildTransferMatrixDesc(TPZCompMesh &transfermesh,TPZTransfer<STATE> &transfer);
	void ProjectSolution(TPZFMatrix<STATE> &projectsol);
    
    //set nummber of meshs
    void SetNMeshes(int64_t nmeshes){
        fNmeshes = nmeshes;
    }
    
    int64_t GetNMeshes(){
        return fNmeshes;
    }
	
private:
	
	/** @brief Creates the computational elements, and the degree of freedom nodes */ 
	/** If MaterialIDs is passed, only element of material id in the set<int> will be created */
	virtual void AutoBuild(const std::set<int> *MaterialIDs);
	
public:
	
    /** @brief Create a computational element based on the geometric element */
    TPZCompEl *CreateCompEl(TPZGeoEl *gel, int64_t &index)
    {
        return fCreate.CreateCompEl(gel, *this, index);
    }
    
	/** @brief Creates the computational elements, and the degree of freedom nodes */ 
	/** Only element of material id in the set<int> will be created */
	virtual void AutoBuild(const std::set<int> &MaterialIDs){
		fCreate.BuildMesh(*this,MaterialIDs);
	}
	
	/** @brief Creates the computational elements, and the degree of freedom nodes */
	virtual void AutoBuild(){
#ifdef PZDEBUG
        {
            TPZGeoMesh *gmesh = Reference();
            TPZCheckGeom check(gmesh);
            check.CheckUniqueId();
        }
#endif
		fCreate.BuildMesh(*this);
	}
		
	/** @brief Creates the computational elements, and the degree of freedom nodes */
	/**
	 * Elements created may be TPZInterpolatedElement or TPZCompElDisc. \n
	 * indices contains the type of the element. Element type are given by the enumerate MCreationType.
	 */
	virtual void AutoBuildContDisc(const TPZVec<TPZGeoEl*> &continuous, const TPZVec<TPZGeoEl*> &discontinuous);
    
    TPZCreateApproximationSpace &ApproxSpace()
    {
        return fCreate;
    }
	
    void SetAllCreateFunctionsDiscontinuous()
    {
        fCreate.SetAllCreateFunctionsDiscontinuous();
    }
    
    void SetAllCreateFunctionsContinuous()
    {
        fCreate.SetAllCreateFunctionsContinuous();
    }
    
    void SetAllCreateFunctionsDiscontinuousReferred()
    {
        fCreate.SetAllCreateFunctionsDiscontinuousReferred();
    }
    
    void SetAllCreateFunctionsContinuousReferred()
    {
        fCreate.SetAllCreateFunctionsContinuousReferred();
    }
    
    void SetAllCreateFunctionsHDiv()
    {
        fCreate.SetAllCreateFunctionsHDiv(Dimension());
    }
	
	//space for full basis for quadrilateral element
	void SetAllCreateFunctionsHDivFull()
    {
        fCreate.SetAllCreateFunctionsHDivFull(Dimension());
    }
	

#ifndef STATE_COMPLEX
    void SetAllCreateFunctionsHDivPressure()
    {
        fCreate.SetAllCreateFunctionsHDivPressure(Dimension());
    }
#endif
    
    void SetAllCreateFunctions(TPZCompEl &cel)
    {
        fCreate.SetAllCreateFunctions(cel, this);
    }

    void SetAllCreateFunctionsMultiphysicElem()
    {
        fCreate.SetAllCreateFunctionsMultiphysicElem();
    }

    void SetAllCreateFunctionsMultiphysicElemWithMem()
    {
        fCreate.SetAllCreateFunctionsMultiphysicElemWithMem();
    }
	
    void SetAllCreateFunctionsContinuousWithMem()
    {
        fCreate.SetAllCreateFunctionsContinuousWithMem();
    }
		
	/** @brief Will build the list of element boundary conditions build the list of connect boundary conditions. */
	/** Put material pointers into the elements. Check on the number of dof of the connects */
	int Consolidate();
	
	/**
	 * ??
	 */
	int Check();

	/**
	 * @brief Given the solution of the global system of equations, computes and stores the solution for the restricted nodes
	 * @param sol given solution matrix
	 */
	void LoadSolution(const TPZFMatrix<STATE> &sol);
    
    /**
     * @brief Transfer multiphysics mesh solution
     */
    void TransferMultiphysicsSolution();
	
	/**
	 * @brief Divide the element corresponding to index
	 * @param index element to divide index
	 * @param subindex divided elements indexes vector
	 * @param interpolate ????
	 */
	void Divide(int64_t index,TPZVec<int64_t> &subindex,int interpolate = 0);
	
	/**
	 * @brief Create a computational element father for the comp. elem. into elements
	 * @param elements indices of subelements to be agrouped
	 * @param index of new coarse element
	 * @param CreateDiscontinuous indicates a TPZInterpolatedElement must be created. True indicates a TPZCompElDisc
	 */
	void Coarsen(TPZVec<int64_t> &elements, int64_t &index, bool CreateDiscontinuous = false);
	
	/** @brief Deletes all interfaces and rebuild them all */
	void RemakeAllInterfaceElements();
	
	/**
	 * @brief Will refine the elements associated with a boundary condition till there are
	 * no elements constrained by boundary condition elements
	 */
	void AdjustBoundaryElements();
	
	/**
	 * @brief Permute the sequence number of the connect objects
	 * It is a permute gather operation
	 * @param permute vector of elements to permute
	 */
	void Permute(TPZVec<int64_t> &permute);
    
    /** @brief Put the sequence number of the pressure connects after the seq number of the flux connects */
    void SaddlePermute();
    void SaddlePermute2();

    /// Integrate the variable name over the mesh
    /*
     * @param varname name of the variable that will be integrated
     * @param matids ids of the materials that will contribute to the integral
     */
    TPZVec<STATE> Integrate(const std::string &varname, const std::set<int> &matids);

	
	/**
	 * @brief Evaluates the error given the two vectors of the analised parameters
	 * @param fp pointer for the function with following arguments:
	 * @note Parameter loc - local vector of the analised parameter
	 * @note Parameter val - given vector to compare
	 * @note Parameter deriv - ????
	 * @param errorSum - return the L1 error
	 */
    void EvaluateError(std::function<void (const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv)> fp, bool store_error,
					   TPZVec<REAL> &errorSum);
	
	/** @brief This method compute the jump solution of interface and convert discontinuous elements with jump less than eps in continuous elements. */
	/**
	 * It may be compared the following values to eps: \n
	 * If \f$ opt = 0 \f$ then  \f$ eps \approx \sqrt{ \int { (leftsol - rightsol)^2 } } \f$ \n
	 * If \f$ opt = 1 \f$ then \f$ eps \approx Max{ \| leftsol - rightsol } \| \f$
	 *
	 * @param eps Tolerancy of the jump to cancel and to convert discontinuous element in continuous
	 * @param opt Option by type of norm (\f$ L_2 \f$ norm or \f$ L \infty \f$ norm).
	 * @param dim Dimension of the working discontinuous elements
	 * @param celJumps Vector to store the diference between the values from right and left elements connected on the interface
	 */
	void ConvertDiscontinuous2Continuous(REAL eps, int opt, int dim, TPZVec<STATE> &celJumps);
	
	/**
	 * @brief This method convert a discontinuous element with index disc_index in continuous element
	 * @param disc_index Index of the discontinuous element to be converted
	 * @param new_index Returns the index of the new continuous element created
	 */
	void Discontinuous2Continuous(int64_t disc_index, int64_t &new_index);
	
	/** @} */
	
	
	/** @brief Clone this mesh */
	TPZCompMesh * Clone() const;
	
	/** @brief Copies the materials of this mesh to the given mesh */
	void CopyMaterials(TPZCompMesh &mesh) const ;
	
	REAL DeltaX();
	
	REAL MaximumRadiusOfMesh();
	
	REAL LesserEdgeOfMesh();
	
	/** @brief Creates a mesh inserting into accumlist the element list to more refined mesh */
	TPZCompMesh *ComputeMesh(TPZVec<int64_t> &accumlist,int64_t numaggl);
	
	/**
	 * @brief This method will fill the matrix passed as parameter with a representation of the fillin of the global stiffness matrix, based on the sequence number of the connects
	 * @param resolution Number of rows and columns of the matrix
	 * @param fillin Matrix which is mapped onto the global system of equations and represents the fillin be assigning a value between 0. and 1. in each element 
	 */
	void ComputeFillIn(int64_t resolution, TPZFMatrix<REAL> &fillin);
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	public:
virtual int ClassId() const;

	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);

};

inline int64_t TPZCompMesh::AllocateNewConnect(int nshape, int nstate, int order) {
	int64_t connectindex = fConnectVec.AllocateNewElement();
    TPZConnect &c = fConnectVec[connectindex];
    c.SetNShape(nshape);
    c.SetNState(nstate);
    c.SetOrder(order,connectindex);
    c.SetLagrangeMultiplier(0);
	int64_t blocknum = fBlock.NBlocks();
	fBlock.SetNBlocks(blocknum+1);
	fBlock.Set(blocknum,nshape*nstate);
	c.SetSequenceNumber(blocknum);
	return connectindex;
}

/**
 * @brief Returns an index to a new connect
 * @param connect Connect from which attributes should be copied
 */
inline int64_t TPZCompMesh::AllocateNewConnect(const TPZConnect &connect)
{
   	int64_t connectindex = fConnectVec.AllocateNewElement();
    TPZConnect &c = fConnectVec[connectindex];
    c = connect;
    c.RemoveDepend();
    int nshape = c.NShape();
    int nstate = c.NState();
	int64_t blocknum = fBlock.NBlocks();
	fBlock.SetNBlocks(blocknum+1);
	fBlock.Set(blocknum,nshape*nstate);
	c.SetSequenceNumber(blocknum);
	return connectindex;
}


inline void TPZCompMesh::SetReference(TPZGeoMesh * gmesh){
	this->fReference = gmesh;
    fGMesh = TPZAutoPointer<TPZGeoMesh>(0);
}

inline void TPZCompMesh::SetReference(TPZAutoPointer<TPZGeoMesh> & gmesh){
    this->fReference = gmesh.operator->();
    fGMesh = gmesh;
}

#endif

/**
 * @file
 * @brief Contains declaration of TPZCompElHDivPressure class which implements a generic computational element (HDiv scope).
 */

#ifndef PZHDIVPRESSUREHTT
#define PZHDIVPRESSUREHTT

#include "pzelchdiv.h"


/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivPressure : public TPZCompElHDiv<TSHAPE> {
	
	/** @brief Defines the interpolation order for pressure variable*/
	int fPressureOrder;
	
	/** @brief To append vectors */
	void Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12);
public:
	
	TPZCompElHDivPressure(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
	TPZCompElHDivPressure(TPZCompMesh &mesh, const TPZCompElHDivPressure<TSHAPE> &copy);
	
	/**
	 * @brief Constructor used to generate patch mesh... generates a map of connect index from
	 * global mesh to clone mesh
	 */
	TPZCompElHDivPressure(TPZCompMesh &mesh,
						  const TPZCompElHDivPressure<TSHAPE> &copy,
						  std::map<int64_t,int64_t> & gl2lcConMap,
						  std::map<int64_t,int64_t> & gl2lcElMap);
	
	TPZCompElHDivPressure();
	
	virtual ~TPZCompElHDivPressure();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		return new TPZCompElHDivPressure<TSHAPE> (mesh, *this);
	}
	
	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes
	 * mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap) const
	{
		return new TPZCompElHDivPressure<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
    /** @brief Set create function in TPZCompMesh to create elements of this type */
#ifndef STATE_COMPLEX
  virtual void SetCreateFunctions(TPZCompMesh *mesh){
		mesh->SetAllCreateFunctionsHDivPressure();
	}
#endif
	
	virtual MElementType Type();
	
	virtual int NConnects() const;
	
	virtual void SetConnectIndex(int i, int64_t connectindex);
	
	virtual int NConnectShapeF(int connect, int order) const;
	
	virtual int Dimension() const {
		return TSHAPE::Dimension;
	}
	
	virtual int NCornerConnects() const {
		return 0;
	}
	
	virtual int64_t ConnectIndex(int node) const;
    
    /** @brief returns the index of the pressure connect
     * returns -1 if their is no pressure connect
     */
    virtual int PressureConnectIndex() const
    {
        return NConnects()-1;
    }
	
	/** @brief Identifies the interpolation order for pressure variable*/
	virtual void SetPressureOrder(int ord);
	/** @brief Returns the interpolation order to dual variable */
	int DualOrder();
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord);
	
	/** @brief Sets the preferred interpolation order along a side
	 
	 This method only updates the datastructure of the element
	 In order to change the interpolation order of an element, use the method PRefine*/
	virtual void SetPreferredOrder(int order);
	
	virtual int ConnectOrder(int connect) const;		
	
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions
	 */
	virtual void InitMaterialData(TPZMaterialData &data);
	
	/** @brief Computes the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	
	/** 
	 * @brief Compute the shape functions corresponding to the dual space
	 */
	virtual void ShapeDual(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
	
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
    
    ///Compute the solution for a given variable
	virtual void Solution( TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol);
    
private:
	virtual	void ComputeSolution(TPZVec<REAL> &qsi, TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes);
    
public:
    
	/** 
	 * @brief Compute shape functions based on master element in the classical FEM manne. 
	 * @param[in] qsi point in master element coordinates 
	 * @param[in] data stores all input data
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data);
	
	void ComputeSolutionPressureHDiv(TPZVec<REAL> &qsi, TPZMaterialData &data);
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
                                 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);	
	
    /** @brief Compute the solution using Hdiv structure */
	void ComputeSolutionPressureHDiv(TPZMaterialData &data);
	
	
	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);
	
	/** @brief Returns the transformation which transform a point from the side to the interior of the element */
	TPZTransform<> TransformSideToElement(int side);
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */

virtual int ClassId() const;

	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
};

template<class TSHAPE>
int TPZCompElHDivPressure<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHDivPressure") ^ TPZCompElHDiv<TSHAPE>::ClassId() << 1;
}

/** @brief Creates computational point element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressurePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational linear element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressureLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressureQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressureTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational cube element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressureCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational prismal element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressurePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational pyramidal element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressurePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational tetrahedral element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressureTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @} */

#endif

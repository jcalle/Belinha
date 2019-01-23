/**
 * @file
 * @brief Contains the declaration of the Reduced Space class.
 * @author Philippe Devloo
 * @since 7/30/12.
 */
#ifndef PZ_pzreducedspace_h
#define PZ_pzreducedspace_h

#include "pzinterpolationspace.h"

/**
 * @ingroup CompElement
 */
class TPZReducedSpace : public TPZInterpolationSpace
{
public:
    /** @brief Default constructor */
	TPZReducedSpace();
	
	/** @brief Default destructor */
	virtual ~TPZReducedSpace();
	
	/** @brief Puts a copy of the element in the referred mesh */
	TPZReducedSpace(TPZCompMesh &mesh, const TPZReducedSpace &copy);
	
	/** @brief Puts a copy of the element in the patch mesh */
	TPZReducedSpace(TPZCompMesh &mesh, const TPZReducedSpace &copy, std::map<int64_t,int64_t> &gl2lcElMap);
	
	/** @brief Copy of the element in the new mesh whit alocated index */
	TPZReducedSpace(TPZCompMesh &mesh, const TPZReducedSpace &copy, int64_t &index);
	
	/**
	 * @brief Create a computational element within mesh
	 * @param mesh mesh where will be created the element
	 * @param gel geometrical element to insert
	 * @param index new elemen index
	 */
	/** Inserts the element within the data structure of the mesh */
	TPZReducedSpace(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
    static void SetAllCreateFunctionsReducedSpace(TPZCompMesh *cmesh);

    /** @brief Returns the number of nodes of the element */
	virtual int NConnects() const 
    {
        return 1;
    }
    
    /** @brief Returns the number of dof nodes along side iside*/
    virtual int NSideConnects(int iside) const
    {
        return NConnects();
    }
    
    /**
     * @brief Returns the local node number of icon along is
     * @param icon connect number along side is
     * @param is side which is being queried
     */
    virtual int SideConnectLocId(int icon,int is) const
    {
#ifdef PZDEBUG
        if (icon != 0) {
            DebugStop();
        }
#endif
        return 0;
    }
    

	
	/** @brief It returns the shapes number of the element */
	virtual int NShapeF() const;
	
	/** @brief Returns the number of shapefunctions associated with a connect*/
	virtual int NConnectShapeF(int inod, int order) const;
	
	/** @brief Returns the max order of interpolation. */
	virtual int MaxOrder();
	
	/** 
	 * @brief Computes the shape function set at the point x. 
	 * @param qsi point in master element coordinates
	 * @param phi vector of values of shapefunctions, dimension (numshape,1)
	 * @param dphi matrix of derivatives of shapefunctions, dimension (dim,numshape)
	 */
	/**
	 * This method uses the order of interpolation
	 * of the element along the sides to compute the number of shapefunctions
	 */
	virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	
	/** 
	 * @brief Computes the shape function set at the point x. 
	 * @param qsi point in master element coordinates
	 * @param phi vector of values of shapefunctions, dimension (numshape,1)
	 * @param dphix matrix of derivatives of shapefunctions, dimension (dim,numshape)
     * @param axes axes indicating the direction of the derivatives
	 */
	/**
	 * This method uses the order of interpolation
	 * of the element along the sides to compute the number of shapefunctions
	 */
	virtual void ShapeX(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphix, TPZFMatrix<REAL> &axes);
    
    //by Agnaldo
    virtual void ShapeX(TPZVec<REAL> &qsi,TPZMaterialData &data);
    
    virtual void ComputeShape(TPZVec<REAL> &qsi,TPZMaterialData &data);
    virtual void ComputeSolution(TPZVec<REAL> &qsi,TPZMaterialData &data);

	/** 
	 * @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions
	 */
	virtual void InitMaterialData(TPZMaterialData &data);
	
	/** @brief Compute and fill data with requested attributes */
	virtual void ComputeRequiredData(TPZMaterialData &data,
									 TPZVec<REAL> &qsi);

    /**
     * @brief Computes solution and its derivatives in local coordinate qsi
     * @param qsi master element coordinate
     * @param phi matrix containing shape functions compute in qsi point
     * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
     * @param axes axes indicating the direction of the derivatives
     * @param sol finite element solution
     * @param dsol solution derivatives
     */
    void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
                         const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);
    
	/** @brief Initialize element matrix in which is computed CalcStiff */
	void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);
	
	/** @brief Initialize element matrix in which is computed in CalcResidual */
	void InitializeElementMatrix(TPZElementMatrix &ef);
	
	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
    
    
    virtual void PRefine ( int order ){
        DebugStop();
    }
    
    virtual void SetConnectIndex(int inode, int64_t index) {
        DebugStop();
    }
    
    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const;
    
    virtual const TPZIntPoints &GetIntegrationRule() const
    {
        TPZInterpolationSpace *intel = ReferredIntel();
        return intel->GetIntegrationRule();
    }
    
    virtual TPZIntPoints &GetIntegrationRule()
    {
        TPZInterpolationSpace *intel = ReferredIntel();
        return intel->GetIntegrationRule();
    }
    
    virtual int Dimension() const {
        TPZInterpolationSpace *intel = ReferredIntel();
        return intel->Dimension();
    }
	
    virtual TPZCompEl * ClonePatchEl (TPZCompMesh &mesh, std::map< int64_t, int64_t > &gl2lcConMap, std::map< int64_t, int64_t > &gl2lcElMap) const;
    
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const{
        
    }
    
    virtual int64_t ConnectIndex(int i) const {
        if (i != 0) {
            DebugStop();
        }
        return 0;
    }
    
    virtual void SetPreferredOrder ( int order ){
        PZError <<"This method was not implemented";
        DebugStop();
    }
    
    void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);
    public:
virtual int ClassId() const;

private:

    TPZInterpolationSpace *ReferredIntel() const;
};


#endif

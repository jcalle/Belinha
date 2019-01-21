/**
 * @file
 * @brief Contains declaration of TPZCompelDisc class which implements a computational element for discontinuous interpolation space.
 */

/** Discontinous Elements */

#ifndef ELCOMPDISCHPP
#define ELCOMPDISCHPP

#include <iostream>
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzcompel.h"
#include "pzinterpolationspace.h"
#include "pzgeoel.h"
#include "pzreal.h"
#include "TPZShapeDisc.h"
#include "tpzautopointer.h"
#include "pzquad.h"
#include "pzfunction.h"
#include "pzfmatrix.h"

struct TPZElementMatrix;
class TPZBndCond;
class TPZConnect;
class TPZMaterial;
class TPZGeoEl;
class TPZCompMesh;

/**
 * @brief This class implements a discontinuous element (for use with discontinuous Galerkin). \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
class TPZCompElDisc : public TPZInterpolationSpace {
	
private:
	
	TPZAutoPointer<TPZIntPoints> CreateIntegrationRule() const;
	
protected:
	
	TPZAutoPointer<TPZIntPoints> fIntRule;
	
	/** @brief Shape function type used by the element */
	pzshape::TPZShapeDisc::MShapeType fShapefunctionType;
	
public:
	
	/** @brief Sets tensorial shape functions for all Discontinuous elements in cmesh. */
	static void SetTensorialShape(TPZCompMesh * cmesh);
	
	/** @brief Set total order shape functions for all Discontinuous elements in cmesh. */
	static void SetTotalOrderShape(TPZCompMesh * cmesh);
	
	/** @brief Set tensorial shape functions. */
	void SetTensorialShape();
	
	/** @brief Set total order shape functions. */
	void SetTotalOrderShape();
	
	/** @brief Set tensorial shape functions with many derivatives. */
	/** Available only for 2D shape functions. */
	void SetTensorialShapeFull();
	
	/** @brief Set total order shape functions. */
	/** Available only for 2D shape functions. */
	void SetTotalOrderShapeFull();
	
protected:
	
	/** @brief It preserves index of connect associated to the element */
	int64_t fConnectIndex;
	
	/** @brief Normalizing constant for shape functions */
	REAL fConstC;
    
    /**  @brief Variable to choose the  qsi point or the X point  in the calculus of the phis and dphis */
    bool fUseQsiEta;
	
	/** @brief A pz function to allow the inclusion of extra shape functions which are defined externally. */
	TPZAutoPointer<TPZFunction<STATE> > fExternalShape;
	
protected:
	
	/** @brief It keeps the interior point coordinations of the element */
	TPZManVector<REAL,3> fCenterPoint;
	
	/**
	 * @brief It creates new conect that it associates the degrees of freedom of the
	 * element and returns its index
	 */
	virtual int64_t CreateMidSideConnect();
	
public:
	
	/** @brief Define external shape functions which are stored in class attribute fExternalShape */
	void SetExternalShapeFunction(TPZAutoPointer<TPZFunction<STATE> > externalShapes);
	
	/** @brief Return whether element has external shape functions set to */
	bool HasExternalShapeFunction();
	
	int GetMaterial( const TPZGeoElSide& gside );
	
	/** @brief Creates discontinuous computational element */
	static TPZCompEl *CreateDisc(TPZGeoEl *geo, TPZCompMesh &mesh, int64_t &index);
	
	/**
	 * @brief Sets the orthogonal function which will be used throughout the program.
	 * @param orthogonal pointer to a function which will be used to generate the shape functions
	 */
	static void SetOrthogonalFunction(void (*orthogonal)(REAL C, REAL x0, REAL x,int degree, TPZFMatrix<REAL> & phi, TPZFMatrix<REAL> & dphi, int n)){
		pzshape::TPZShapeDisc::fOrthogonal = orthogonal;
	}
	
	/** @brief Sets the inner radius value. */
	virtual void SetInnerRadius(REAL InnerRadius) {PZError << "TPZCompElDisc::SetInnerRadius - This method should never be called because the inner" << std::endl
		<< "radius is not stored in TPZCompElDisc. It is stored in TPZAgglomerateElement." << std::endl;}
	
	/** @brief Sets element's number of interfaces. */
	virtual void SetNInterfaces(int nfaces) {PZError << "TPZCompElDisc::SetNFaces - This method should never be called because the number of interfaces" << std::endl
		<< "is not stored in TPZCompElDisc. It is only stored by TPZAgglomerateElement." << std::endl;}
	
	/** @brief Returns the number of interfaces. */
	virtual int NInterfaces();
	
	/** @brief Default constructor */
	TPZCompElDisc();
	/** @brief Constructor of the discontinuous element associated with geometric element */
	TPZCompElDisc(TPZCompMesh &mesh,TPZGeoEl *ref,int64_t &index);//original
	/** @brief Constructor */
	TPZCompElDisc(TPZCompMesh &mesh,int64_t &index);//construtor do aglomerado
	/** @brief Copy constructor */
	TPZCompElDisc(TPZCompMesh &mesh, const TPZCompElDisc &copy);

	/**
	 * @brief Creates a clone of the given element in a pathc mesh
	 * @param mesh patch mesh
	 * @param copy element to be copied
	 * @param gl2lcConMap map between the connect indexes in original and patch mesh
	 * @param gl2lcElMap map between the element indexes in original an patch mesh
	 */
	TPZCompElDisc(TPZCompMesh &mesh,
				  const TPZCompElDisc &copy,
				  std::map<int64_t,int64_t> &gl2lcConMap,
				  std::map<int64_t,int64_t> &gl2lcElMap);
	
	TPZCompElDisc(TPZCompMesh &mesh, const TPZCompElDisc &copy,int64_t &index);
	
	/** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh);
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		return new TPZCompElDisc(mesh,*this);
	}
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh,int64_t &index) const {
		return new TPZCompElDisc(mesh,*this,index);
	}
	
	/**
	 * @see class TPZCompEl
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
									std::map<int64_t,int64_t> & gl2lcConMap,
									std::map<int64_t,int64_t> & gl2lcElMap) const {
		return new TPZCompElDisc(mesh,*this,gl2lcConMap,gl2lcElMap);
	}
	/** @brief Default destructor */
	virtual ~TPZCompElDisc();
	
	/** @brief Divide the computational element */
	void Divide(int64_t index, TPZVec<int64_t> &subindex, int interpolate = 0);
	
	/**
	 * @brief Computes the shape function set at the point x. This method uses the order of interpolation
	 * of the element along the sides to compute the number of shapefunctions
	 * @param qsi point in master element coordinates
	 * @param phi vector of values of shapefunctions, dimension (numshape,1)
	 * @param dphi matrix of derivatives of shapefunctions, dimension (dim,numshape)
	 */
	virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphidxi);
	
	
protected:
	virtual void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidx);
    
    public:
    
    /** 
	 * @brief Compute shape functions based on master element in the classical FEM manne. 
	 * @param[in] intpoint point in master element coordinates 
	 * @param[in] data stores all input data
	 */
    virtual void ComputeShape(TPZVec<REAL> &intpoint,TPZMaterialData &data);
	
	/**
	 */
	virtual void ShapeX(TPZVec<REAL> &X, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
	
	/** @brief Add extenal shape function into already computed phi and dphi discontinuous functions. */
	void AppendExternalShapeFunctions(TPZVec<REAL> &X, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
	
	/** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
	virtual const TPZIntPoints &GetIntegrationRule() const;
	
	/** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
	virtual TPZIntPoints &GetIntegrationRule();
	
	/** @brief Type of the element */
	virtual MElementType Type() {return EDiscontinuous;}
	
	/** @brief It returns the constant that normalizes the bases of the element */
	REAL ConstC() const {return fConstC;}
	
	void SetConstC(REAL c){fConstC = c;}
    
    void SetTrueUseQsiEta(){
        fUseQsiEta = true;
        fCenterPoint.Fill(0.);
        fConstC = 1.;

        TPZGeoEl *gel = Reference();
        if (gel)
        {
            int ns = Reference()->NSides();
            this->Reference()->CenterPoint(ns-1, fCenterPoint);
        }
    }

    void SetFalseUseQsiEta(){
        fUseQsiEta = false;
        Reference()->CenterPoint(Reference()->NSides()-1,fCenterPoint);
        TPZVec<REAL> csi(fCenterPoint);
        Reference()->X(csi,fCenterPoint);
        fConstC = NormalizeConst();
    }

	/** @brief Returns the center point */
	void InternalPoint(TPZVec<REAL> &point);
	
	/** @brief Prints the features of the element */
	virtual void Print(std::ostream & out = std::cout) const;
	
	/** @brief Returns the degree of interpolation of the element */
	virtual int Degree() const;
	
	/** @brief Assigns the degree of the element */
	virtual void SetDegree(int degree);
	/** @brief Returns the number of connects */
	virtual int NConnects() const;
	
	/** @brief Amount of vertices of the element */
	int NCornerConnects() const { return Reference()->NNodes();}
    
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const;
	
	/** @brief Returns dimension from the element */
	int Dimension() const { return Reference()->Dimension();}
	
	/** @brief Calculates the normalizing constant of the bases of the element */
	virtual REAL NormalizeConst();
	
	/** @brief Returns the connect index from the element */
	int64_t ConnectIndex(int side = 0) const;
	void  SetConnectIndex(int /*inode*/, int64_t index) {fConnectIndex = index;}
    
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

	
	/** @brief Returns the shapes number of the element */
	virtual int NShapeF() const;
	
	/** @brief Returns the max order of interpolation. */
	virtual int MaxOrder();
	
	/** @brief Returns the max order of interpolation excluding external shape order */
	int MaxOrderExceptExternalShapes();
	
	/** @brief Returns the number of shapefunctions associated with a connect*/
	virtual int NConnectShapeF(int inod, int order) const;
	
	REAL CenterPoint(int index) const {return fCenterPoint[index];}
	
	virtual void CenterPoint(TPZVec<REAL> &center);
	
	void SetCenterPoint(int i,REAL x){fCenterPoint[i] = x;}
	
	REAL SizeOfElement();
	
	/** @brief Returns the volume of the geometric element associated. */
	virtual  REAL VolumeOfEl() { return Reference()->Volume(); }
	
	/**
	 * @brief Creates corresponding graphical element(s) if the dimension matches
	 * graphical elements are used to generate output files
	 * @param grmesh graphical mesh where the element will be created
	 * @param dimension target dimension of the graphical element
	 */
	void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension);
	
	/** @deprecated
	 * @brief Computes the solution in function of a point in cartesian space
	 */
	/**
	 * Deprecated shape function method. It is kept because of TPZAgglomerateElement. \n
	 * It does not include singular shape functions if they exist.
	 */
	void SolutionX(TPZVec<REAL> &x,TPZVec<STATE> &uh);
	
    public:
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi master element coordinate
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 * @param axes axes associated with the derivative of the solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> & axes);
    
public:
    
    /** 
	 * @brief Compute shape functions based on master element in the classical FEM manne. 
	 * @param[in] qsi point in master element coordinates 
	 * @param[in] data stores all input data
	 */
    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data);
	
	/**
	 * @brief Computes solution and its derivatives in local coordinate qsi
	 * @param qsi master element coordinate
	 * @param phi matrix containing shape functions compute in qsi point
	 * @param dphix matrix containing the derivatives of shape functions with respect of global coordinates: D[phi,x], D[phi,y], D[phi,z]
	 * @param axes axes associated with the derivative of the solution
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
								 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);
	
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi. \n
	 * This method will function for both volumetric and interface elements
	 * @param qsi master element coordinate of the interface element
	 * @param normal unitary normal vector
	 * @param leftsol finite element solution
	 * @param dleftsol solution derivatives
	 * @param leftaxes axes associated with the left solution
	 * @param rightsol finite element solution
	 * @param drightsol solution derivatives
	 * @param rightaxes axes associated with the right solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZVec<REAL> &normal,
								 TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
								 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes);
	
	virtual void AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight);
	
	/** @brief Accumulates the vertices of the agglomerated elements */
	virtual void AccumulateVertices(TPZStack<TPZGeoNode *> &nodes);
	
	int NSides();
	
	void BuildTransferMatrix(TPZCompElDisc &coarsel, TPZTransfer<STATE> &transfer);
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	public:
virtual int ClassId() const;

	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
	/** @brief Define the desired order for entire element. */
	virtual void SetPreferredOrder ( int order ) { SetDegree( order ); }
	
	/**
	 * @brief Change the preferred order for the element and proceed the
	 * adjust of the aproximation space taking in acount the type \n
	 * of formulation and the neighbours of the element
	 */
	virtual void PRefine ( int order ) { SetDegree( order ); }
	
	/** @brief Compute the integral of the square residual over the element domain. */
	/** 
	 * For instance, let us take the Poisson's equation: -Laplac(u) = f.
	 * Then this method compute the integral of ( -Laplac(u) - f )^2.
	 * For the given example it is observed that approximation orders < 2 will
	 * not work properly since Laplac(u) is always = 0.
	 * It works only for elements that lay in the same space dimension of the master element.
	 */
	static REAL EvaluateSquareResidual2D(TPZInterpolationSpace *cel);
	
	/** @brief Evaluates the square residual for every element in mesh. */
	static void EvaluateSquareResidual2D(TPZCompMesh &cmesh, TPZVec<REAL> &error, bool verbose = false);
    
	
};

inline TPZCompEl *TPZCompElDisc::CreateDisc(TPZGeoEl *geo, TPZCompMesh &mesh, int64_t &index) {
	if(!geo->Reference() && geo->NumInterfaces() == 0)
		return new TPZCompElDisc(mesh,geo,index);
	return NULL;
}
//Exemplo do quadrilatero:
//acessar com -> TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);

#endif

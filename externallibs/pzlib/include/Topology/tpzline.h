/**
 * @file
 * @brief Contains the TPZLine class which defines the topology of a line element.
 */

#ifndef PZTOPOLOGYTPZLINE_H
#define PZTOPOLOGYTPZLINE_H

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzquad.h"
#include "pzeltype.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

class TPZIntPoints;
class TPZInt1d;
class TPZGraphEl1dd;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

/// Groups all classes defining the structure of the master element
namespace pztopology {
	
	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of a line element. \ref topology "Topology"
	 * Sides 0 and 1 are vertices, side 2 is the line. 
	 */
	class TPZLine : public TPZSavable {
	public:
		
		/** @brief Enumerate for topological characteristics */
		enum {NCornerNodes = 2, NSides = 3, Dimension = 1, NFaces = 2};
		
                public:
                virtual int ClassId() const;
                void Read(TPZStream& buf, void* context);
                void Write(TPZStream& buf, int withclassid) const;

		/** @brief Default constructor */
        TPZLine() : TPZRegisterClassId(&TPZLine::ClassId){
		}
		
		/** @brief Default destructor */
		virtual ~TPZLine() {
		}

		/** @name About sides of the topological element
		 * @{ */

		/** @brief Returns the dimension of the side */
		static int SideDimension(int side);

		/** @brief Get all sides with lower dimension on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides);
		/** @brief Get all sides with lower dimension but equal to DimTarget on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget);
		
		/**
		 * @brief Returns all sides whose closure contains side
		 * @param side smaller dimension side
		 * @param high vector which will contain all sides whose closure contain sidefrom
		 */
		static void HigherDimensionSides(int side, TPZStack<int> &high);

		/** @brief Returns the number of nodes (not connectivities) associated with a side */
		static int NSideNodes(int side);

		/** @brief Returns the local node number of the node "node" along side "side" */
		static int SideNodeLocId(int side, int node);

		/** @brief Returns the number of connects for a set dimension */
		static int NumSides();
		/** @brief Returns the number of connects for a set dimension */
		static int NumSides(int dimension);
		
		/** @brief Returns the number of nodes (not connectivities) associated with a side */
		static int NContainedSides(int side);

		/** @brief Returns the local connect number of the connect "c" along side "side" */
		static int ContainedSideLocId(int side, int c);
		
		/** @} */
		
		/** @name About points at the parametric spaces
		 * @{ */

		/** @brief Returns the barycentric coordinates in the master element space of the original element */
		static void CenterPoint(int side, TPZVec<REAL> &center);
		
		/** @brief Verifies if the parametric point pt is in the element parametric domain */
		static bool IsInParametricDomain(TPZVec<REAL> &pt, REAL tol = 1e-6);
        
        /** @brief Generates a random point in the master domain */
        static void RandomPoint(TPZVec<REAL> &pt);
        
        template<class T>
        static bool MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide);
        
        static void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);

		/** @} */

		/** @name About type of the topological element
		 * @{ */
		
		/** @brief Returns the type of the element as specified in file pzeltype.h */
		static MElementType Type();
		/** @brief Returns the type of the element side as specified in file pzeltype.h */
		static MElementType Type(int side);

		/** @} */
		
		/** @name About Transformations
		 * @{ */

		/**
		 * @brief Returns the transformation which takes a point from the side sidefrom to the side sideto
		 * @param sidefrom Side where the point resides
		 * @param sideto Side whose closure contains sidefrom
		 */
		static TPZTransform<> SideToSideTransform(int sidefrom, int sideto);

		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side Side from which the point will be tranformed (0<=side<=2)
		 * @return TPZTransform<> object
		 */
		static TPZTransform<> TransformSideToElement(int side);

		/**
		 * @brief Returns the transformation which transform a point from the interior of the element to the side
		 * @param side Side to which the point will be tranformed (0<=side<=2)
		 * @return TPZTransform<> object
		 */
		static TPZTransform<> TransformElementToSide(int side);
		
		/**
		 * @brief Method which identifies the transformation based on the IDs of the corner nodes
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */
		static int GetTransformId(TPZVec<int64_t> &id);
		
		/**
		 * @brief Method which identifies the transformation of a side based on the IDs
		 * of the corner nodes
		 * @param side Index of side
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */	
		static int GetTransformId(int side, TPZVec<int64_t> &id);
		
		/** @} */
		
		/** @name Methods related over numeric integration
		 * @{ */

		/**
		 * @brief Create an integration rule over side
		 * @param side Side to create integration rule
		 * @param order Order of the integration rule to be created
		 */
		static TPZIntPoints *CreateSideIntegrationRule(int side, int order);
		
		/** @brief Typedef to numerical integration rule */
		typedef TPZInt1d IntruleType;
		/** @brief Typedef to graphical element type */
		typedef TPZGraphEl1dd GraphElType;
		
		/** @} */
		
		/**
		 * @brief Identifies the permutation of the nodes needed to make neighbouring elements compatible 
		 * in terms of order of shape functions
		 * @param side Side for which the permutation is needed
		 * @param id Ids of the corner nodes of the elements
		 * @param permgather Permutation vector in a gather order
		 */
		static void GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather);

		/** @brief Volume of the master element (measure of the element) */
		static REAL RefElVolume() { return 2.0; }
        
        /* Given side and gradx the method returns directions needed for Hdiv space */
        static void ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors);
        static void GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao);
        static void GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao, TPZVec<int> &sidevectors);
        
        /// Compute the directions of the HDiv vectors
        static void ComputeDirections(TPZFMatrix<REAL> &gradx, REAL detjac, TPZFMatrix<REAL> &directions)
        {
            for (int i=0; i<3; i++) {
                directions(i,0) = -gradx(i,0)/detjac;
                directions(i,1) = gradx(i,0)/detjac;
                directions(i,2) = gradx(i,0)/detjac;
            }
        }
        

        
        /**
         * Returns the number of bilinear sides to this shape. Needed to compute the number shapefunctions( NConnectShapeF )
         */
        static int NBilinearSides();
		
	};
	
}

#endif

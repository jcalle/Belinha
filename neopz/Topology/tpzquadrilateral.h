/**
 * @file
 * @brief Contains the TPZQuadrilateral class which defines the topology of a quadrilateral element. 
 */

#ifndef PZTOPOLOGYTPZQUADRILATERAL_H
#define PZTOPOLOGYTPZQUADRILATERAL_H

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzquad.h"
#include "pzeltype.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

class TPZIntPoints;
class TPZIntQuad;
class TPZGraphElQ2dd;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

/// Groups all classes defining the structure of the master element
namespace pztopology {
	
	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of a quadrilateral element. \ref topology "Topology"
	 * Sides 0 to 3 are vertices, sides 4 to 7 are lines, side 8 is the quadrilateral. 
	 */
	class TPZQuadrilateral : public TPZSavable {
	public:

		/** @brief Enumerate for topological characteristics */
		enum {NSides = 9, NCornerNodes = 4, Dimension = 2, NFaces = 4};

            virtual int ClassId() const;
            void Read(TPZStream& buf, void* context);
            void Write(TPZStream& buf, int withclassid) const;

                
		/** @brief Default constructor */
        TPZQuadrilateral() : TPZRegisterClassId(&TPZQuadrilateral::ClassId){
		}
		
		/** @brief Default destructor */
		virtual ~TPZQuadrilateral() {
		}
		
		/** @name About sides of the topological element
		 * @{ */

		/** @brief returns the dimension of the side */
		static int SideDimension(int side);
		
		/** @brief Get all sides with lower dimension on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides);
		/** @brief Get all sides with lower dimension but equal to DimTarget on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget);
		
		/**
		 * @brief returns all sides whose closure contains side
		 * @param side smaller dimension side
		 * @param high vector which will contain all sides whose closure contain sidefrom
		 */
		static void HigherDimensionSides(int side, TPZStack<int> &high);
		
		/** @brief return the number of vertices (not connectivities) associated with a side */
		static int NSideNodes(int side);
		/** @brief returns the local node number of the node "node" along side "side" */
		static int SideNodeLocId(int side, int node);

		/** @brief Returns number of sides of the element (9) */
		static int NumSides() {
			return NSides;
		}
		/** @brief return the number of sides of a given dimension */
		static int NumSides(int dimension);

		/** @brief return the number of sides contained in the closure of a side */
		static int NContainedSides(int side);
		
		/** @brief returns the local side number of the side "c" contained in the closure of side "side" */
		static int ContainedSideLocId(int side, int c);
		
		/** @} */
		
		/** @name About points at the parametric spaces
		 * @{ */

		/** @brief returns the barycentric coordinates in the master element space of the original element */
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
		 * @brief returns the transformation which takes a point from the side sidefrom to the side sideto
		 * @param sidefrom side where the point resides
		 * @param sideto side whose closure contains sidefrom
		 * @see the class TPZTransform
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
		 * @param side Side to which the point will be tranformed (0<=side<=8)
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
        
        /**
         * @brief return the vector which permutes the connects according to the transformation id
         */
        static void GetGatherPermute(int transformid, TPZVec<int> &permute);
		
		/** @} */
		
		/** @name Methods related over numeric integration
		 * @{ */
		
		/**
		 * @brief Create an integration rule over side
		 * @param side Side to create integration rule
		 * @param order Order of the integration rule to be created
		 */		
		static TPZIntPoints * CreateSideIntegrationRule(int side, int order);
		
		/** @brief Typedef to numerical integration rule */
		typedef TPZIntQuad IntruleType;
		/** @brief Typedef to graphical element type */
		typedef TPZGraphElQ2dd GraphElType;
		
		/** @} */

		/**
		 * @brief Identifies the permutation of the nodes needed to make neighbouring elements compatible 
		 * in terms of order of shape functions
		 * @param side Side for which the permutation is needed
		 * @param id Ids of the corner nodes of the elements
		 * @param permgather Permutation vector in a gather order
		 */
		static void GetSideHDivPermutation(int transformid, TPZVec<int> &permgather);
		
		/** @brief Volume of the master element (measure) */
		static REAL RefElVolume() { return 4.0L; }
        
        /* Given gradx the method returns directions needed for Hdiv space */
        static void ComputeDirections(int is, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors);
        static void GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao);
        static void GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao, TPZVec<int> &sidevectors);
        
        /// Compute the directions of the HDiv vectors
        static void ComputeDirections(TPZFMatrix<REAL> &gradx, REAL detjac, TPZFMatrix<REAL> &directions);
        
        /**
         * Returns the number of bilinear sides to this shape. Needed to compute the number shapefunctions( NConnectShapeF )
         */
        static int NBilinearSides();

	};
	
}

#endif

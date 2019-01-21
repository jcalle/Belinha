//
//  pzmatmixedpoisson3d.h
//  PZ
//
//  Created by Agnaldo Farias on 03/11/14.
//
//

#ifndef __PZ__pzmatmixedpoisson3d__
#define __PZ__pzmatmixedpoisson3d__

#include <stdio.h>
#include <iostream>
#include "pzdiscgal.h"
#include "TPZMaterial.h"
#include "TPZLagrangeMultiplier.h"

/**
 * @ingroup material
 * @author Douglas
 * @since 11/03/2014
 * @brief Material to solve a mixed poisson problem 3D by multiphysics simulation
 * @brief Pressure(p): uses L2 space.  Velocity (Q): uses Hdiv space
 */
/**
 * \f$ Q = -K*Grad(p)  ==> Int{K^-1*Q.v}dx_e + Int{Grad(p)*v}dx_e - Int{sign(vn)*p*vn}ds_e = - Int{pD*vn}ds (Eq. 1)  \f$
 *
 * \f$ div(Q) = f  ==> Int{v.Grad(w)}dx_e - Int{sign(Qn)Qn*w}ds_e = - Int{f*w}dx_e  (Eq. 2) \f$
 *
 * \f$ p = pD in Dirichlet boundary and Q.n = qN in Neumann boundary\f$
 */

class TPZMatMixedPoisson3D : public TPZMaterial {
    
protected:
    /** Material Id */
    int fMatId;
    
    /** Valor da funcao de carga */
    REAL fF; //fF
    
    /** Dimensao do dominio */
    int fDim;
    
    /** Coeficiente que multiplica o gradiente */
    REAL fK;
    
     /** Coeficiente que multiplica a pressão: termo de reacao */
    REAL falpha;
    
    /** @brief fluid viscosity*/
    REAL fvisc;
    
    /** @brief permeability tensor. Coeficient which multiplies the gradient operator*/
    TPZFMatrix<REAL> fTensorK;
    
    /** @brief inverse of the permeability tensor.*/
    TPZFMatrix<REAL> fInvK;
    
    /** @brief Pointer to forcing function, it is the Permeability and its inverse */
    TPZAutoPointer<TPZFunction<STATE> > fPermeabilityFunction;
    
    /** @brief Pointer to forcing function, it is the reaction term */
    TPZAutoPointer<TPZFunction<STATE> > fReactionTermFunction;
    
    //object to material lagrange multiplier
    TPZLagrangeMultiplier *fmatLagr;
    
    /** @brief Parameter to choose the second integration by parts in the variational formulation*/
    bool fSecondIntegration;
    
    bool fReactionTerm;
    
public:
    
    TPZMatMixedPoisson3D();
    
    TPZMatMixedPoisson3D(int matid, int dim);
    
    virtual ~TPZMatMixedPoisson3D();
    
    /** @brief copy constructor */
    TPZMatMixedPoisson3D(const TPZMatMixedPoisson3D &copy);
    
    TPZMatMixedPoisson3D &operator=(const TPZMatMixedPoisson3D &copy);
    
    virtual std::string Name() { return "TPZMatMixedPoisson3D"; }
    
    int Dimension() const {return fDim;}
    
    int MatId()
    {
        return fMatId;
    }
    
    virtual int NStateVariables();
    
    void SetPermeability(REAL perm) {
        fK = perm;
    }
    
    //Set the permeability tensor and inverser tensor
    void SetPermeabilityTensor(TPZFMatrix<REAL> K, TPZFMatrix<REAL> invK){
        
        //        if(K.Rows() != fDim || K.Cols() != fDim) DebugStop();
        //        if(K.Rows()!=invK.Rows() || K.Cols()!=invK.Cols()) DebugStop();
        
        fTensorK = K;
        fInvK = invK;
    }
    
    void SetReactionTerm(REAL alpha)
    {
        fReactionTerm = true;
        falpha = alpha;
    }
    
    void SetViscosity(REAL visc) {
        fvisc = visc;
    }
    
    void GetPermeability(REAL &perm) {
        perm = fK;
    }
    
    void SetInternalFlux(REAL flux) {
        fF = flux;
    }
    
    void SetPermeabilityFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
        fPermeabilityFunction = fp;
    }
    
    TPZAutoPointer<TPZFunction<STATE> > PermeabilityFunction()
    {
        return fPermeabilityFunction;
    }
    
    void SetfReactionTermFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
        fReactionTermFunction = fp;
    }
    
    TPZAutoPointer<TPZFunction<STATE> > ReactionTermFunction()
    {
        return fReactionTermFunction;
    }
    
    
    
    /** @brief Gets the order of the integration rule necessary to integrate an element with polinomial order p */
    virtual int IntegrationRuleOrder(int elPMaxOrder) const {
        
        return 20;
    }
    
    /** @brief Gets the order of the integration rule necessary to integrate an element multiphysic */
    virtual int IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const
    {
        int polorder = elPMaxOrder[0]*2;
        int forceorder = 0;
        if (fForcingFunction) {
            forceorder = fForcingFunction->PolynomialOrder();
        }
        if (forceorder > elPMaxOrder[0]) {
            polorder = forceorder+elPMaxOrder[0];
        }
        return polorder;
    }
    
    void Print(std::ostream &out);
    
    /** @name Contribute methods
     * @{
     */
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * This method use normalized piola contravariant mapping for nonlinear mappings. With second integration by parts
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since June 2, 2014
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    ///This method use piola contravariant mapping for nonlinear mappings
    void ContributeWithoutSecondIntegration(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since June 2, 2014
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
        DebugStop();
    }
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since June 2, 2014
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since June 2, 2014
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    //    /**
    //     * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
    //     * @param datavec [in]
    //     * @param dataleft [in]
    //     * @param weight [in]
    //     * @param ek [out] is the stiffness matrix
    //     * @param ef [out] is the load vector
    //     * @param bc [in] is the boundary condition object
    //     * @since June 2, 2014
    //     */
    //    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    //
    //    /**
    //     * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
    //     * @param data [in]
    //     * @param dataleft [in]
    //     * @param weight [in]
    //     * @param ek [out] is the stiffness matrix
    //     * @param ef [out] is the load vector
    //     * @param bc [in] is the boundary condition object
    //     * @since June 2, 2014
    //     */
    //    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    //    void         ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    //
    //    /**
    //     * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
    //     * @param datavec [in]
    //     * @param dataleft [in]
    //     * @param dataright [in]
    //     * @param weight [in]
    //     * @param ek [out] is the stiffness matrix
    //     * @param ef [out] is the load vector
    //     * @since June 2, 2014
    //     */
    //    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    //
    //    /**
    //     * @brief It computes a contribution to stiffness matrix and load vector at one integration point
    //     * @param data [in]
    //     * @param dataleft [in]
    //     * @param dataright [in]
    //     * @param weight [in]
    //     * @param ek [out] is the stiffness matrix
    //     * @param ef [out] is the load vector
    //     * @since June 2, 2014
    //     */
    //    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {DebugStop();}
    
    /**
     * @brief It return a solution to multiphysics simulation.
     * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */
    //virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    /**
     * @brief .
     * @param datavec [in] Data material vector
     * @param DivergenceofPhi Divergence of Phi over the master element
     * @param DivergenceofU Divergence of flux over the master element
     */
    void ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU);
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    // metodo para gerar vtk
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    // metodo para computar erros Pressao
    void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
    // metodo para computar erros Hdiv
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    // metodo para computar erros Pressao
    void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
                TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &/*flux*/,
                TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
    // metodo para computar erros Hdiv
    void ErrorsHdiv(TPZMaterialData &data,TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
    
    void UseSecondIntegrationByParts(){
        fSecondIntegration=true;
    }
    
    bool IsUsedSecondIntegration(){
        return fSecondIntegration;
    }
    public:
virtual int ClassId() const;

};

#endif /* defined(__PZ__pzmatmixedpoisson3d__) */

//
//  TPZSBFemElementGroup.hpp
//  PZ
//
//  Created by Philippe Devloo on 4/4/16.
//
//

#ifndef TPZSBFemElementGroup_hpp
#define TPZSBFemElementGroup_hpp

#include <stdio.h>

#include "pzelementgroup.h"
#include "TPZSBFemVolume.h"


class TPZSBFemElementGroup : public TPZElementGroup
{
    
public:
    enum EComputationMode {EStiff, EOnlyMass, EMass};
    
private:
    
    /// Matrix of eigenvectors which compose the stiffness matrix
    TPZFMatrix<std::complex<double> > fPhi;
    
    /// Inverse of the eigenvector matrix (transfers eigenvector coeficients to side shape coeficients)
    TPZFMatrix<std::complex<double> > fPhiInverse;
    
    /// Vector of eigenvalues of the SBFem analyis
    TPZManVector<std::complex<double> > fEigenvalues;
    
    /// Multiplying coefficients of each eigenvector
    TPZFMatrix<std::complex<double> > fCoef;
    
    TPZFMatrix<STATE> fMassMatrix;
    
    EComputationMode fComputationMode = EStiff;
    
    /// multiplier to multiply the mass matrix
    REAL fMassDensity = 1.;
    
    /// timestep coeficient
    REAL fDelt = 1.;
    
    /// Compute the mass matrix based on the value of M0 and the eigenvectors
    void ComputeMassMatrix(TPZElementMatrix &M0);
    
public:
    
    TPZSBFemElementGroup() : TPZElementGroup()
    {
        
    }
    
    /// constructor
    TPZSBFemElementGroup(TPZCompMesh &mesh, int64_t &index) : TPZElementGroup(mesh,index)
    {
        
    }
    
    /** @brief add an element to the element group
     */
    virtual void AddElement(TPZCompEl *cel);
    

    /// Compute the SBFem matrices
    /// method to assemble E0, E1, E2
    void ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2, TPZElementMatrix &M0);
    
    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
    
    /// set the density or specific heat of the material
    void SetDensity(REAL density)
    {
        fMassDensity = density;
    }
    /// Set the element to compute the mass matrix
    void SetComputeOnlyMassMatrix()
    {
        if(fMassMatrix.Rows() == 0)
        {
            DebugStop();
        }
        fComputationMode = EOnlyMass;
    }
    
    /// Set the element to compute stiffness plus mass
    void SetComputeTimeDependent(REAL delt)
    {
        fDelt = delt;
        fComputationMode = EMass;
    }
    
    void SetComputeStiff()
    {
        fComputationMode = EStiff;
    }
    /**
     * @brief Prints element data
     * @param out Indicates the device where the data will be printed
     */
    virtual void Print(std::ostream &out = std::cout) const
    {
        out << __PRETTY_FUNCTION__ << std::endl;
        TPZElementGroup::Print(out);
        int nel = fElGroup.size();
        out << "Element indexes of the volume elements ";
        for (int el=0; el<nel; el++) {
            out << fElGroup[el]->Index() << " ";
        }
        out << std::endl;
        out << "Indices of the associated computational skeleton elements\n";
        for (int el=0; el<nel; el++) {
            TPZCompEl *cel = fElGroup[el];
            TPZSBFemVolume *vol = dynamic_cast<TPZSBFemVolume *>(cel);
            if(!vol) DebugStop();
            out << vol->SkeletonIndex() << " ";
        }
        out << std::endl;
        out << "Connect indexes of the contained elements\n";
        for (int el=0; el<nel; el++) {
            TPZCompEl *cel = fElGroup[el];
            int nc = cel->NConnects();
            for (int ic=0; ic<nc; ic++) {
                out << cel->ConnectIndex(ic) << " ";
            }
            out << std::endl;
        }

/*
        out << "EigenVectors for displacement\n";
        fPhi.Print("Phi = ",out,EMathematicaInput);
        out << "Inverse EigenVectors\n";
        fPhiInverse.Print("PhiInv = ",out,EMathematicaInput);
        out << "EigenValues " << fEigenvalues << std::endl;
        out << "Mass Matrix\n";
        fMassMatrix.Print("Mass = ",out);
        out << "Solution Coeficients\n";
        fCoef.Print("Coef ",out);
        for (int el=0; el<nel; el++) {
            fElGroup[el]->Print(out);
        }
 */
        out << "End of " << __PRETTY_FUNCTION__ << std::endl;
    }
    

    
    /**
     * @brief Computes the element right hand side
     * @param ef element load vector(s)
     */
    virtual void CalcResidual(TPZElementMatrix &ef)
    {
        TPZElementMatrix ek(Mesh(),TPZElementMatrix::EK);
        CalcStiff(ek,ef);
    }
    

    /** @brief Loads the solution within the internal data structure of the element */
    /**
     * Is used to initialize the solution of connect objects with dependency. \n
     * Is also used to load the solution within SuperElements
     */
    virtual void LoadSolution();

    /** @brief Loads the geometric element referece */
    virtual void LoadElementReference()
    {
        for (int64_t i = 0; i < fElGroup.size(); i++) {
            fElGroup[i]->LoadElementReference();
        }
    }
    


    int64_t NumEigenValues()
    {
        return fEigenvalues.size();
    }
    
    /// Load the coeficients such that we visualize an eigenvector
    void LoadEigenVector(int64_t eig);
    /// method to compute the stiffness
    /// method to compute the solution
    
    TPZVec<STATE> MultiplyingCoeficients()
    {
        int nel = fCoef.Rows();
        TPZVec<STATE> result(nel,0.);
        for (int ir=0; ir<nel; ir++) {
            result[ir] = fCoef(ir,0).real();
        }
        return result;
    }
    
    TPZVec<std::complex<double> > &EigenValues()
    {
        return fEigenvalues;
    }
    
    TPZFMatrix<std::complex<double> > &Phi()
    {
        return fPhi;
    }
    
    TPZFMatrix<std::complex<double> > &PhiInverse()
    {
        return fPhiInverse;
    }
    
    TPZFMatrix<STATE> &MassMatrix()
    {
        return fMassMatrix;
    }
    
    TPZFMatrix<std::complex<double> > Coeficients()
    {
        return fCoef;
    }
    
    TPZFMatrix<double> PhiReal()
    {
        int64_t rows = fPhi.Rows(),cols = fPhi.Cols();
        TPZFMatrix<double> phireal(rows,cols);
        for(int64_t i=0; i<rows; i++)
        {
            for(int64_t j=0; j<cols; j++)
            {
                phireal(i,j) = fPhi(i,j).real();
            }
        }
        return phireal;
    }
    
    TPZManVector<double> EigenvaluesReal()
    {
        int64_t nel = fEigenvalues.NElements();
        TPZManVector<double> eig(nel);
        for(int64_t el=0; el<nel; el++)
        {
            eig[el] = fEigenvalues[el].real();
        }
        return eig;
    }
    
    TPZFMatrix<double> CoeficientsReal()
    {
        int64_t rows = fCoef.Rows(),cols = fCoef.Cols();
        TPZFMatrix<double> coefreal(rows,cols);
        for(int64_t i=0; i<rows; i++)
        {
            for(int64_t j=0; j<cols; j++)
            {
                coefreal(i,j) = fCoef(i,j).real();
            }
        }
        return coefreal;
    }

};

#endif /* TPZSBFemElementGroup_hpp */

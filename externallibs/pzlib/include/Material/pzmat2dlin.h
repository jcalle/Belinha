/**
 * @file
 * @brief Contains the TPZMat2dLin class which implements a bi-dimensional linear problem.
 */

#ifndef MAT2DLINHPP
#define MAT2DLINHPP

#include <iostream>
#include "TPZMaterial.h"
#include "pzfmatrix.h"

struct TPZElementMatrix;
class TPZBndCond;
template<class T>
class TPZVec;

/**
 * @ingroup material
 * @brief Implements a bi-dimensional linear problem.
 */
class TPZMat2dLin : public TPZMaterial {
	
	TPZFMatrix<STATE>    fKxx, fKxy, fKyx, fKyy, fKx0, fK0x, fKy0, fK0y, fK00, fXf;
	public :
	
    TPZMat2dLin(int num = 1) : TPZRegisterClassId(&TPZMat2dLin::ClassId),
    TPZMaterial(num), fKxx(), fKxy(),
	fKyx() , fKyy(), fKx0(), fK0x(), fKy0(), fK0y(), fK00(), fXf() {
    }
	
	TPZMat2dLin(TPZMat2dLin &copy) : TPZRegisterClassId(&TPZMat2dLin::ClassId),
    TPZMaterial(copy),
	fKxx(copy.fKxx), fKxy(copy.fKxy), fKyx(copy.fKyx), fKyy(copy.fKyy),
	fKx0(copy.fKx0), fK0x(copy.fK0x), fKy0(copy.fKy0),
	fK0y(copy.fK0y), fK00(copy.fK00), fXf(copy.fXf){}
	
	virtual int NStateVariables() { return fKxx.Rows(); }
	
	virtual int Dimension() const { return 2; }
	
	void Print(std::ostream & out = std::cout);
	
	void SetMaterial(TPZFMatrix<STATE> &xkin,TPZFMatrix<STATE> &xcin,TPZFMatrix<STATE> &xfin){
		int r = xkin.Rows();
		fKxx = xkin;
		fKyy = xkin;
		fK00 = xcin;
		fXf = xfin;
		fKxy.Redim(r,r);
		fKyx.Redim(r,r);
		fKx0.Redim(r,r);
		fK0x.Redim(r,r);
		fKy0.Redim(r,r);
		fK0y.Redim(r,r);
	}
	
	void ConvectionDiffusion(REAL angle,REAL diff);
	
	TPZFMatrix<STATE> &Xk() {return fKxx;}
	TPZFMatrix<STATE> &Ck() {return fK00;}
	TPZFMatrix<STATE> &Xf() {return fXf;}
	
	virtual std::string Name() { return "TPZMat2dLin"; }
	
	virtual void Contribute(TPZMaterialData &data,REAL weight,
							TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
	
	virtual void Contribute(TPZMaterialData &data,REAL weight,
							TPZFMatrix<STATE> &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC(TPZMaterialData &data, REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
	virtual void ContributeBC(TPZMaterialData &data, REAL weight,
							  TPZFMatrix<STATE> &ef,TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int NFluxes();
	
	virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &fl);
	
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &dudx,TPZFMatrix<REAL> &axes,TPZVec<STATE> &flux,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int index);
	
protected:
	void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes, int var,TPZVec<STATE> &Solout);
public:
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	/** @brief Creates a copy of the material object */
	virtual TPZMaterial * NewMaterial();
	
	TPZBndCond *OutflowFlux(TPZMaterial * &reference, int bc);
	
    /** @{
     * @name Save and Load methods
     */
    
	/** @brief returns the unique identifier for reading/writing objects to streams */
	public:
virtual int ClassId() const;

    
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
    /**
     * @}
     */
};

#endif

/**
 * @file
 * @brief Contains the declaration of the TPZCompElWithMem class, it is as TPZCompEl with enable material memory feature.
 */

#ifndef PZCOMPELWITHMEM_H
#define PZCOMPELWITHMEM_H

class TPZMaterialData;

#include "pzinterpolationspace.h"
#include "pzstack.h"
#include "pzcmesh.h"
#include "pzquad.h"
#include "TPZMaterial.h"
#include "pzelctemp.h"
#include "pzmultiphysicscompel.h"

//#include "tpzpoint.h"

#include "pzlog.h"



#ifdef LOG4CXX
static LoggerPtr CompElWMemlogger(Logger::getLogger("mesh.TPZCompElWithMem"));
#endif

/**
 * @brief This class implements the TPZCompEl structure to enable material memory feature. \n
 * It should be instantiated using one of TPZCompEl bottom classes as parent in the template parameter.
 * @ingroup CompElement
 * @since July, 23 2008
 */

template <class TBASE>
class TPZCompElWithMem : public TBASE
{
public:
  
  TPZCompElWithMem();
  
  virtual ~TPZCompElWithMem();
  
  TPZCompElWithMem(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
  
  TPZCompElWithMem(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index, TPZCompElSide left, TPZCompElSide right);
  
  TPZCompElWithMem(TPZCompMesh &mesh, const TPZCompElWithMem<TBASE> &copy);
  
  /** @brief used to generate patch mesh... generates a map of connect index from global mesh to clone mesh */
  TPZCompElWithMem(TPZCompMesh &mesh,
                   const TPZCompElWithMem<TBASE> &copy,
                   std::map<int64_t,int64_t> & gl2lcConMap,
                   std::map<int64_t,int64_t> & gl2lcElMap);
  
  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override{
    return new TPZCompElWithMem<TBASE> (mesh, *this);
  }
  
  /** @brief Set create function in TPZCompMesh to create elements of this type */
  virtual void SetCreateFunctions(TPZCompMesh *mesh) override{
    mesh->SetAllCreateFunctionsContinuousWithMem();
  }
  
  /**
   * @brief Create a copy of the given element. The clone copy have the connect indexes mapped to the local clone connects by the given map
   * @param mesh Patch clone mesh
   * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
   * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
   */
  virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap) const override
  {
    return new TPZCompElWithMem<TBASE> (mesh, *this, gl2lcConMap, gl2lcElMap);
  }
  
 
  virtual void ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &qsi) override;
  
  virtual void ComputeRequiredData(TPZVec<REAL> &intpointtemp, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialData> &datavec) override;
 
  int64_t GetGlobalIntegrationPointIndex(TPZMaterialData &data);
  
protected:
  
  /** @brief PrepareIntPtIndices initializes the material damage varibles memory in the proper material class. */
  virtual void PrepareIntPtIndices() override;
  
  /** @brief PrepareIntPtIndices initializes the material damage varibles memory in the proper material class. */
  virtual void ForcePrepareIntPtIndices() override;
  
  /** @brief PrepareIntPtIndices initializes the material damage varibles memory in the proper material class. */
  virtual void SetMemoryIndices(TPZVec<int64_t> &indices) override;
  
  /** @brief Frees the material damage varibles memory in the proper material class. */
  virtual void SetFreeIntPtIndices() override;
  
  void CopyIntPtIndicesFrom(const TPZCompElWithMem<TBASE> & copy);
  
public:
  
  /** @brief Get the indices of the vector of element memory associated with the integration points */
  /**
   * Will return an empty vector if no memory is associated with the integration point
   * Is implemented in TPZCompElWithMem
   */
  void GetMemoryIndices(TPZVec<int64_t> &indices) const override;
  
  /// Modify the maximum order an integration rule can integrate
  virtual void SetIntegrationRule(int ord) override;
  
  /** @brief Saves the element data to a stream */
  virtual void Write(TPZStream &buf, int withclassid) const override;
  
  /** @brief Reads the element data from a stream */
  virtual void Read(TPZStream &buf, void *context) override;
  
  /** @brief ClassId of the class. Is implemented for each type of compel in this .h */
  public:
virtual int ClassId() const override;

  /**
   * @name Print
   * @brief Methods for print data structure
   * @{
   */
  
  /**
   * @brief Prints element data
   * @param out indicates the device where the data will be printed
   */
  virtual void Print(std::ostream & out = std::cout) const override
  {
    TBASE::Print(out);
    out << "Integration point indexes " << fIntPtIndices << std::endl;
  }
  
  /** @} */
  
private:
  
  TPZStack<int64_t,128> fIntPtIndices;
  
};

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem() : TPZRegisterClassId(&TPZCompElWithMem::ClassId),
TBASE() {
  //PrepareIntPtIndices();
}

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElWithMem::ClassId),
TBASE(mesh, gel, index){
  PrepareIntPtIndices();
}

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index, TPZCompElSide left, TPZCompElSide right) :
TPZRegisterClassId(&TPZCompElWithMem::ClassId),
TBASE(mesh, ref, index, left, right){
  PrepareIntPtIndices();
}

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem(TPZCompMesh &mesh, const TPZCompElWithMem<TBASE> &copy) :
TPZRegisterClassId(&TPZCompElWithMem::ClassId),
TBASE(mesh, copy) {
  CopyIntPtIndicesFrom(copy);
}


template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem(TPZCompMesh &mesh,
                                          const TPZCompElWithMem<TBASE> &copy,
                                          std::map<int64_t,int64_t> & gl2lcConMap,
                                          std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElWithMem::ClassId),
TBASE(mesh,copy,gl2lcConMap,gl2lcElMap)
{
  CopyIntPtIndicesFrom(copy);
}



template <class TBASE>
inline void TPZCompElWithMem<TBASE>::ForcePrepareIntPtIndices() {
  
  TPZMaterial * material = TBASE::Material();
  if(!material){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
    return;
  }
  
  if (this->NumberOfCompElementsInsideThisCompEl() == 0) {
    // This is suposed to happen if in the constructor of a multiphysics element. The CompEl vector is only initialized after the autobuild
    return;
  }
  
  const TPZIntPoints &intrule = TBASE::GetIntegrationRule();
  
  int intrulepoints = intrule.NPoints();
  
  fIntPtIndices.Resize(intrulepoints);
  
  for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){
    fIntPtIndices[int_ind] = this->Material()->PushMemItem();
    // Pushing a new entry in the material memory
  } //Loop over integratin points generating a reference vector of memory
  //entries in the related pzmatwithmem for further use.
  
}

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::SetMemoryIndices(TPZVec<int64_t> &indices) {
    int n = indices.size();
    fIntPtIndices.Resize(n);
    for (int i = 0; i < n; i++) {
        fIntPtIndices[i] = indices[i];
    }
}

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::SetFreeIntPtIndices() {

    TPZMaterial * material = TBASE::Material();

    if (material) {
        int64_t n = fIntPtIndices.NElements();

        for (int64_t i = 0; i < n; i++) {
            this->Material()->FreeMemItem(fIntPtIndices[i]);
        }
    }
    fIntPtIndices.Resize(0);
}

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::SetIntegrationRule(int ord)
{
  TBASE::SetIntegrationRule(ord);
  // verify if the number of integration points changed
  const TPZIntPoints &intrule = TBASE::GetIntegrationRule();
  int intrulepoints = intrule.NPoints();
  if (intrulepoints != fIntPtIndices.size()) {
    SetFreeIntPtIndices();
    PrepareIntPtIndices();
  }
  
}


/** @brief Get the indices of the vector of element memory associated with the integration points */
/**
 * Will return an empty vector if no memory is associated with the integration point
 * Is implemented in TPZCompElWithMem
 */
template<class TBASE>
inline void TPZCompElWithMem<TBASE>::GetMemoryIndices(TPZVec<int64_t> &indices) const
{
  indices = fIntPtIndices;
}


template <class TBASE>
void TPZCompElWithMem<TBASE>::CopyIntPtIndicesFrom(const TPZCompElWithMem<TBASE> & copy)
{
  
  TPZMaterial * material = TBASE::Material();
  if(!material){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
    return;
  }
  
  int64_t i, n = copy.fIntPtIndices.NElements();
  fIntPtIndices.Resize(n);
  
  for(i = 0; i < n; i++)
  {
    fIntPtIndices[i] = copy.fIntPtIndices[i];
  }
}

/** Save the element data to a stream */
template <class TBASE>
inline void TPZCompElWithMem<TBASE>::ComputeRequiredData(TPZMaterialData &data,
                                                         TPZVec<REAL> &qsi){
  TBASE::ComputeRequiredData(data, qsi);
  data.intGlobPtIndex = GetGlobalIntegrationPointIndex(data);
  //material index for the n-th CompEl integration point
}

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::ComputeRequiredData(TPZVec<REAL> &intpointtemp, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialData> &datavec)
{
  TBASE::ComputeRequiredData(intpointtemp,trvec,datavec);
  
  int nelofthismphysics = this->NumberOfCompElementsInsideThisCompEl();
  for (int icel = 0; icel < nelofthismphysics; icel++) {
    datavec[icel].intGlobPtIndex = GetGlobalIntegrationPointIndex(datavec[icel]);
  }
}

template <class TBASE>
inline int64_t TPZCompElWithMem<TBASE>::GetGlobalIntegrationPointIndex(TPZMaterialData &data)
{
  int64_t glIntegralPt = -1;
  if (data.intLocPtIndex >= 0) {
    glIntegralPt = fIntPtIndices[ data.intLocPtIndex ]; // returning the
  }
  return glIntegralPt;
}

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::Write(TPZStream &buf, int withclassid) const
{
  TBASE::Write(buf,withclassid);
  buf.Write( fIntPtIndices);
  int classid = ClassId();
  buf.Write(&classid);
}

/** Read the element data from a stream */
template <class TBASE>
inline void TPZCompElWithMem<TBASE>::Read(TPZStream &buf, void *context)
{
  TBASE::Read(buf,context);
  buf.Read( fIntPtIndices);
  int classid;
  buf.Read(&classid);
  if (classid != ClassId()) {
    DebugStop();
  }
}

template <class TBASE>
int TPZCompElWithMem<TBASE>::ClassId() const{
    return Hash("TPZCompElWithMem") ^ TBASE::ClassId() << 1;
}

#endif

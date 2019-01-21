#ifndef PZEQUATIONFILTERHPP
#define PZEQUATIONFILTERHPP

#include "pzvec.h"
#include "pzfmatrix.h"

#include <set>

class TPZEquationFilter : public TPZSavable {
public:

    TPZEquationFilter() {} 
    
    TPZEquationFilter(int64_t numeq) : fNumEq(numeq), fIsActive(false), fActiveEqs(), fDestIndices()
    {

    }

    TPZEquationFilter(const TPZEquationFilter&cp):fNumEq(cp.fNumEq),fIsActive(cp.fIsActive),
           fActiveEqs(cp.fActiveEqs),fDestIndices(cp.fDestIndices)
    {
      ///nothing here
    }

    ~TPZEquationFilter()
    {
      ///nothing here
    }

    TPZEquationFilter & operator=(const TPZEquationFilter&cp)
    {
        this->fNumEq = cp.fNumEq;
        this->fIsActive = cp.fIsActive;
        this->fActiveEqs = cp.fActiveEqs;
        this->fDestIndices = cp.fDestIndices;
        return *this;
    }

    int ClassId() const{
        return Hash("TPZEquationFilter");
    }
    
    void Read(TPZStream& buf, void* context){
        buf.Read(&fNumEq);
        buf.Read(fIsActive);
        buf.Read(fActiveEqs);
        buf.Read(fDestIndices);
    }
    
    void Write(TPZStream& buf, int withclassid) const{
        buf.Write(&fNumEq);
        buf.Write(fIsActive);
        buf.Write(fActiveEqs);
        buf.Write(fDestIndices);
    }

    ///Define as equacoes ativas de [mineq, maxeq)
    void SetMinMaxEq(int64_t mineq, int64_t maxeq)
    {
      if (mineq < 0 || mineq > fNumEq ||
          maxeq < 0 || maxeq > fNumEq ||
          mineq > maxeq){
          DebugStop();
      }

      const int64_t n = maxeq-mineq;
      TPZVec<int64_t> activeEquations(n);
      for(int64_t i = 0; i < n; i++){
        activeEquations[i] = i + mineq;
      }
      this->SetActiveEquations( activeEquations );
    }

    ///Define as equacoes ativas
    void SetActiveEquations(TPZVec<int64_t> &active)
    {
        if(fActiveEqs.NElements()) DebugStop();///oops, call reset first

        fIsActive = true;
        ///removendo duplicados e reordenando
        std::set<int64_t> activeset;
        int64_t neq = active.size();
        if (neq) {
            activeset.insert(&active[0], &active[neq-1]+1);
        }

        fDestIndices.Resize(fNumEq);
        fDestIndices.Fill(-1);
        int64_t count = 0;
        fActiveEqs.Resize(activeset.size());
        for (std::set<int64_t>::iterator it=activeset.begin(); it != activeset.end(); it++) {
            fActiveEqs[count] = *it;
            fDestIndices[*it] = count++;
        }
    }

    /// Reset method
    void Reset()
    {
        fIsActive = false;
        fActiveEqs.Resize(0);
        fDestIndices.Resize(0);
    }

    /** Filtra as equações:
     * @param orig [in][out] - remove de orig equacoes nao ativas
     * @param dest [in][out] - remove de dest as equcoes nao ativas
     */
    void Filter(TPZVec<int64_t> &orig, TPZVec<int64_t> &dest) const
    {
        if (fDestIndices.size() == 0) {
            return;
        }
        else {
            int64_t count = 0;
            int64_t numeq = dest.size();
            for (int64_t i=0; i<numeq; i++) {
                if (fDestIndices[dest[i]] != -1) {
                    orig[count] = orig[i];
                    dest[count] = fDestIndices[dest[i]];
                    count++;
                }
            }
            orig.Resize(count);
            dest.Resize(count);
        }
    }

    /** Filtra as equações:
      * @param dest [in][out] - remove de dest as equacoes nao ativas
      */
    void Filter(TPZVec<int64_t> &dest) const
    {
        if (fDestIndices.size() == 0) {
            return;
        }
        else{
            int64_t count = 0;
            int64_t numeq = dest.size();
            for (int64_t i=0; i<numeq; i++) {
                if (fDestIndices[dest[i]] != -1) {
                    dest[count] = fDestIndices[dest[i]];
                    count++;
                }
            }
            dest.Resize(count);
        }
    }

    ///Retorna o numero de equacoes ativas do sistema
    int64_t NActiveEquations() const
    {
        if (IsActive()) {
            return fActiveEqs.size();
        }
        else
        {
            return fNumEq;
        }
    }

    ///Retorna o numero de equacoes do sistema original
    int64_t NEqExpand() const
    {
        return fNumEq;
    }

	/**
	 * Returns true if the filter is active
	 */
    bool IsActive() const
    {
        return fIsActive;

    }

	/**
	 * Expands the vector small to a original system, fill zeros into the no active equations.
	 */
	template<class TVar>
    void Scatter(const TPZFMatrix<TVar> &vsmall, TPZFMatrix<TVar> &vexpand) const
    {
        int64_t neqcondense = this->NActiveEquations();
        if(vsmall.Rows() != neqcondense || vexpand.Rows() != fNumEq)
        {
            DebugStop();
        }
        if(! IsActive())
        {
            vexpand = vsmall;
            return;
        }
        vexpand.Zero();

#ifdef PZDEBUG
        {
            for(int64_t i=0; i<neqcondense; i++)
            {
                if(fActiveEqs[i] >= fNumEq)
                {
                    DebugStop();
                }
            }
        }
#endif
        for(int64_t i=0; i<neqcondense; i++) vexpand(fActiveEqs[i],0) = vsmall.GetVal(i,0);
    }

    /**
     * @brief Reduce the vector to the number of active equations.
     */
    template<class T>
    void Gather(const TPZFMatrix<T> &large, TPZFMatrix<T> &gathered) const
    {
        int64_t neqcondense = this->NActiveEquations();
        if(gathered.Rows() != neqcondense || large.Rows() != fNumEq)
        {
            DebugStop();
        }
        if(! IsActive())
        {
            gathered = large;
            return;
        }
        gathered.Zero();
        for(int64_t i=0; i<neqcondense; i++) gathered(i,0) = large.GetVal(fActiveEqs[i],0);
    }

    /**
     * @brief Returns the number of active equations between [minindex,maxindex]
     */
    int64_t NumActive(int64_t minindex, int64_t maxindex) const
    {
        if (minindex < 0 || maxindex < 0 || minindex > fNumEq || maxindex > fNumEq ||
            maxindex < minindex) {
            DebugStop();
        }
        if (!IsActive()) {
            return maxindex-minindex;
        }
        int numactive = 0;
        for (int64_t i=minindex; i<maxindex; i++) {
            if (fDestIndices[i] != -1) {
                numactive++;
            }
        }
        return numactive;
    }
    
    void FilterSkyline(TPZVec<int64_t> &skyline) const
    {
        if (!IsActive()) {
            return;
        }

        for (int64_t ieq = 0; ieq<fActiveEqs.size(); ieq++)
        {
            int64_t skyl = skyline[fActiveEqs[ieq]];
            while (fDestIndices[skyl] == -1 && skyl < fNumEq) {
                skyl++;
            }
#ifdef PZDEBUG
            // all active equations should have a destination
            if (skyl > fActiveEqs[ieq] || fDestIndices[skyl] < 0) {
                DebugStop();
            }
#endif
            skyline[ieq] = fDestIndices[skyl];
        }
        skyline.Resize(fActiveEqs.size());

    }

    void SetNumEq(const int64_t numEq){
        fNumEq = numEq;
    }

private:

    /// Numero de equacoes do sistema original
    int64_t fNumEq;
    
    /// Flag indicating whether the filter is active
    bool fIsActive;
    
    /// Equacoes ativas
    TPZVec<int64_t> fActiveEqs;

    /// Posicao das equacoes originais no sistema reduzido
    TPZVec<int64_t> fDestIndices;
    
};

#endif

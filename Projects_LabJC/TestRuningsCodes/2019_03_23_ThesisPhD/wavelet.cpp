#include "pzcompel.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "wavelet.h"
#include "pzmanvector.h"
#include "ttimeanalysis.h"
#include "pzvec.h"
#include "pzcmesh.h"
#include "pzgraphmesh.h"
#include "pzdxmesh.h"
#include "myheader.h"
#include "hadaptive.h"


TWavelet::TWavelet(double thread,double maxcoef) : THREADHOLDER(thread),
    MAXWAVELETCOEFFICIENT(maxcoef) {
}

void TWavelet::PrintDecomposition(int var,TPZVec<REAL> &coefs,int nsubel, std::ostream &out) {
  int i, j;
  int nfathers = coefs.NElements()/nsubel;
  out <<std::endl << "WAVELET DECOMPOSITION" <<std::endl;
  out << "To variable : " << var <<std::endl <<std::endl;
  for(i=0;i<nfathers;i++) {
    for(j=0;j<nsubel;j++)
      out << coefs[i*nsubel+j] << "\t";
    out <<std::endl;
  }
}

THaar1D::THaar1D(double thread,double maxcoef) : TWavelet(thread,maxcoef) {
}

REAL THaar1D::OneLevelEven(TPZVec<REAL> &values) {
  return .5*(values[0]+values[1]);
}
REAL THaar1D::OneLevelDetails(TPZVec<REAL> &values) {
  return .5*(values[0]-values[1]);
}

int THaar1D::Decomposition(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
     TPZManVector<int> &elsons,int MinLevel) {

  int nelem = cmesh.NElements(), nfathers = elfathers.NElements(), nsons = elsons.NElements();
  /**Para marcar elementos ja trabalhados*/
  TPZVec<int> indexel(nelem,-1);
  int indexintel, indexbrother;
  /**Variaveis para calculo dos coeficientes wavelets*/
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  REAL djk;
  int i, indexfather = nfathers, indexson = nsons;
//  if(nfathers || nsons) PZError << "HaarDecomposition. Already exist fathers or sons.\n";
  TPZVec<REAL> mean(2,0.);
  TPZCompEl *el;
  TPZInterpolatedElement *intel, *brotherel;

  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el || el->IsInterface() || el->Dimension() != 1) continue;
    intel = (TPZInterpolatedElement *)el;
    indexintel = intel->Index();
    if(!indexel[indexintel]) continue;
    int mask = 1;
    TPZGeoEl *gel = intel->Reference()->Father();
    if(!gel) continue;
    if(!gel->HasSubElement()) {
      PZError << "HaarDecomposition. This father element has no sons.\n";
      continue;
    }
    brotherel = (TPZInterpolatedElement *)gel->SubElement(1)->Reference();
    if(brotherel==intel) {
      brotherel = (TPZInterpolatedElement *)gel->SubElement(0)->Reference();
      mask = 0;
    }
    mean[0] = intel->MeanSolution(var);

    if(!brotherel) {
      /**Entao o subelemento foi dividido, gel armazena o subelemento geometrico.
         Calcula-se entao a media celular de gel desde os subelementos computacionais.*/
      gel = gel->SubElement(mask);
      mean[1] = MeanSolutionFromSubElements(gel,var);
      indexbrother = -1;
    }
    else {
      /**Elemento irmao existe*/
      indexbrother = brotherel->Index();
      mean[1] = brotherel->MeanSolution(var);
      /** Para evitar o calculo novamente no sub-elemento irmao de um outro ja calculado*/
      indexel[indexbrother] = 0;
    }

    /** Determinando os coeficientes wavelet para realizar threadholder */
    djk = fabs(OneLevelDetails(mean));
    if(djk > THREADHOLDER) {
      if(djk > MAXWAVELETCOEFFICIENT) ReAssembling = 1;
//      if(intel->Reference()->Level() > MaxLevelRefinements-1 && ReAssembling) continue;
      if(nfathers-(2)<indexfather) {
        nfathers += 100;
        elfathers.Resize(nfathers);
      }
    /**Pode-se evitar o indice se no Divide do elemento computacional tira-se o parametro index*/
      elfathers[indexfather++] = indexintel;
      if(brotherel) elfathers[indexfather++] = indexbrother;
    }
    else if(IsZeroToCoarsing(djk) && indexbrother!=-1) {
      /**If level == 2 we can not to coarse element, because we can not
         to apply decomposition wavelet at next step*/
      if(intel->Reference()->Level()<=MinLevel) continue;
      if(nsons<indexson+2) {
        nsons += 100;
        elsons.Resize(nsons);
      }
      elsons[indexson++] = indexintel;
      elsons[indexson++] = indexbrother;
    }
  }
  elfathers.Resize(indexfather);
  elsons.Resize(indexson);
	return 0;
}

void THaar1D::DecompositionOnlyOrderZero(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
     TPZManVector<int> &elsons,int MinLevel) {

  int nelem = cmesh.NElements(), nfathers = elfathers.NElements(), nsons = elsons.NElements();
  /**Para marcar elementos ja trabalhados*/
  TPZVec<int> indexel(nelem,-1);
  int indexintel, indexbrother;
  /**Variaveis para calculo dos coeficientes wavelets*/
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  REAL djk;
  int i, indexfather = nfathers, indexson = nsons;
  TPZVec<REAL> mean(2,0.);
  TPZCompEl *el;
  TPZInterpolatedElement *intel, *brotherel;

  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el || el->IsInterface() || el->Dimension() != 1) continue;
    intel = (TPZInterpolatedElement *)el;
    indexintel = intel->Index();
    if(!indexel[indexintel]) continue;
    int mask = 1, nsides = intel->Reference()->NSides();
    /**Apenas para elementos com ordem de interpolacao zero.*/
    if(intel->PreferredSideOrder(nsides-1)) continue;
    TPZGeoEl *gel = intel->Reference()->Father();
    if(!gel) continue;
    if(!gel->HasSubElement()) {
      PZError << "HaarDecomposition. This father element has no sons.\n";
      continue;
    }
    brotherel = (TPZInterpolatedElement *)gel->SubElement(1)->Reference();
    if(brotherel==intel) {
      brotherel = (TPZInterpolatedElement *)gel->SubElement(0)->Reference();
      mask = 0;
    }
    mean[0] = intel->MeanSolution(var);

    if(!brotherel) {
      /**Entao o subelemento foi dividido, gel armazena o subelemento geometrico.
         Calcula-se entao a media celular de gel desde os subelementos computacionais.*/
      gel = gel->SubElement(mask);

      mean[1] = MeanSolutionFromSubElements(gel,var);
      indexbrother = -1;
    }
    else {
      /**Elemento irmao existe*/
      indexbrother = brotherel->Index();
      mean[1] = brotherel->MeanSolution(var);
      /** Para evitar o calculo novamente no sub-elemento irmao de um outro ja calculado*/
      indexel[indexbrother] = 0;
    }

    /** Determinando os coeficientes wavelet para realizar threadholder */
    djk = fabs(OneLevelDetails(mean));
    if(djk > THREADHOLDER) {
      if(djk > MAXWAVELETCOEFFICIENT) ReAssembling = 1;
//      if(intel->Reference()->Level() > MaxLevelRefinements-1) continue;
      if(nfathers-(2)<indexfather) {
        nfathers += 100;
        elfathers.Resize(nfathers);
      }
    /**Pode-se evitar o indice se no Divide do elemento computacional tira-se o parametro index*/
      elfathers[indexfather++] = indexintel;
      if(brotherel) elfathers[indexfather++] = indexbrother;
    }
/*    if(djk > THREADHOLDER) {
      if(intel->Reference()->Level() > MaxLevelRefinements-1) continue;
      int ii, mcount = 0;
      TPZManVector<int> indexneigh(0);
      if(segurityzone) {
        /**Procura elementos vizinhos de maior ou igual nivel para refinar-os
        TPZStack<TPZCompElSide> elstack(0);
        for(ii=0;ii<nsides;ii++) {
          TPZCompElSide intelside(intel,ii);
          TPZCompElSide brotherside(brotherel,ii);
          TPZCompElSide large;
          intelside.EqualLevelElementList(elstack,1,0);
          brotherside.EqualLevelElementList(elstack,1,0);
          large = intelside.LowerLevelElementList(1);
          if(large.Exists()) elstack.Push(large);
          large = brotherside.LowerLevelElementList(1);
          if(large.Exists()) elstack.Push(large);
        }
        mcount = elstack.NElements();
        if(mcount) indexneigh.Resize(mcount,-1);
        for(ii=0;ii<mcount;ii++)
          indexneigh[ii] = elstack[ii].Element()->Index();
      }
      if(nfathers-(2+mcount)<indexfather) {
        nfathers += 100;
        elfathers.Resize(nfathers);
      }
    /**Pode-se evitar o indice se no Divide do elemento computacional tira-se o parametro index
      elfathers[indexfather++] = indexintel;
      if(indexbrother!=-1) elfathers[indexfather++] = indexbrother;
      for(ii=0;ii<mcount;ii++) elfathers[indexfather++] = indexneigh[ii];
    }*/
    else if(IsZeroToCoarsing(djk) && indexbrother!=-1) {
      /**If level == 2 we can not to coarse element, because we can not
         to apply decomposition wavelet at next step*/
      if(intel->Reference()->Level()<=MinLevel) continue;
      if(nsons<indexson+2) {
        nsons += 100;
        elsons.Resize(nsons);
      }
      elsons[indexson++] = indexintel;
      elsons[indexson++] = indexbrother;
    }
  }
  elfathers.Resize(indexfather);
  elsons.Resize(indexson);
}

void THaar1D::Decomposition(TPZCompMesh &cmesh,int var,TPZManVector<REAL> &coef,
			 TPZManVector<TPZGeoEl *> &idgfather,int level) {
  TPZManVector<TPZGeoEl *> idgeoel(idgfather);
  int indexfather = 0;
  int k = idgeoel.NElements();
  int nelem = cmesh.NElements();
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  int i, l, j, m, id;
  TPZCompEl *el;
  TPZGeoEl *gel;
  if(!k) {
    int ngelem = cmesh.Reference()->NElements();
    idgeoel.Resize(ngelem,0);
    TPZVec<int> maskgel(ngelem,-1);
    /** Determine all the appropriate elements for decomposition */
    for(i=0;i<nelem;i++) {
      el = elvec[i];
      if(!el || el->IsInterface() || el->Dimension()!=1) continue;
      gel = el->Reference();
      for(j=0;j<level;j++)
        gel = gel->Father();
      if(!gel) {
        PZError << "Haar1d. Bad parameter level, no father.\n";
        continue;
      }
      id = gel->Id();
      if(maskgel[id]!=-1) continue;
      maskgel[id] = id;
      idgeoel[k++] = gel;
    }
    /** Redimension of the coefficient vector */
    idgeoel.Resize(k);
  }
  /** Number of subelements of the geometrical father */
  id = idgeoel[0]->NSubElements();
  int power = (int)pow(id,level);
  int newk = k*power;
  coef.Resize(newk);
  idgfather.Resize(newk,0);
  REAL djk;
  int highid, lowerpower;
  TPZVec<TPZGeoEl *> gelson(power);
  /** Put the values over elements in first part of the vector
      put coefficient wavelets in second part of the vector */
  for(i=0;i<k;i++) {
    power /= id;
    highid = id;
    gel = idgeoel[i];
    for(j=0;j<id;j++) gelson[j*power] = gel->SubElement(j);
    for(j=1;j<level;j++) {
      lowerpower = power/id;
      for(m=0;m<highid;m++) {
        gel = gelson[m*power];
        for(l=0;l<id;l++) gelson[m*power+l*lowerpower] = gel->SubElement(l);
      }
      power /= id;
      highid *= id;
    }
    TPZVec<REAL> tempcoef(highid);
    for(j=0;j<highid;j++)
      tempcoef[j] = ((TPZInterpolatedElement *)gelson[j]->Reference())->MeanSolution(var);
    power = highid;
    int indextemp, index = (i+1)*power-1;
    for(j=0;j<level;j++) {
      indextemp = 0;
      int mask = 0;
      for(l=0;l<highid;l+=2) {
        djk = fabs((tempcoef[l]-tempcoef[l+1])*.5);
        if(djk > THREADHOLDER) coef[index--] = 0.;
        else {
          coef[index--] = djk;
          /** Apenas armazenamos pais quando os coef. wavelet do maior nivel sao grandes (uma vez soh)*/
          if(!j && !mask) {
            for(m=0;m<gelson.NElements();m++)
              idgfather[indexfather++] = gelson[m];
          }
          mask = 1;
        }
        tempcoef[indextemp++] = (tempcoef[l]+tempcoef[l+1])/2;
      }
      highid /= id;
    }
    coef[index] = tempcoef[0];
  }
  idgfather.Resize(indexfather);
}

void THaar1D::DrawWavelets(int var,TTimeAnalysis &an,int posdim,int &step,
		      REAL &time,TPZVec<std::string> &scal,TPZVec<std::string> &vec,int last) {
  TPZCompMesh *cmesh = an.Mesh();
  if(!an.GraphMesh(posdim) || !step) return;
  int i, j;
  int seqnum, nsub = cmesh->Reference()->ElementVec()[0]->NSubElements();
  int64_t index, newindex = 0, nelem = cmesh->NElements();
  TPZVec<int> NewIndex(nelem,-1);
  TPZVec<REAL> Wav(nelem,0.);
  TPZVec<REAL> Val(nelem,0.);
  TPZBlock<STATE> &block = cmesh->Block();
  TPZManVector<TPZCompEl *> ElVec(nelem,0);
  TPZAdmChunkVector<TPZCompEl *> &MeshVec = cmesh->ElementVec();
  TPZCompEl *cel;
  TPZGeoEl *gel;
  for(i=0;i<nelem;i++) {
    cel = MeshVec[i];
    if(!cel || cel->IsInterface() || cel->Dimension()!=1) continue;
    ElVec[i] = cel;
  }

  TPZManVector<int64_t> subel(nsub,-1);  // index of the sons computational elements
  TPZVec<REAL> wav(nsub);
  REAL value;
  TPZInterpolatedElement *intel, *intel2;
  TPZFMatrix<STATE> Solution(cmesh->Solution());
  /**Taking the indexes of sub-elements with same father and filling comp
     mesh with wavelets basis of the previous level*/
  for(i=0;i<nelem;i++) {
    intel = (TPZInterpolatedElement *)ElVec[i];
    if(!intel) continue;
    gel = intel->Reference()->Father();
    intel = (TPZInterpolatedElement *)gel->SubElement(0)->Reference();
    NewIndex[2*newindex] = intel->Index();
    ElVec[NewIndex[2*newindex]] = 0;
    wav[0] = intel->MeanSolution(var);
    intel2 = (TPZInterpolatedElement *)gel->SubElement(1)->Reference();
    wav[1] = intel2->MeanSolution(var);
    NewIndex[2*newindex+1] = intel2->Index();
    ElVec[NewIndex[2*newindex+1]] = 0;
    value = Wav[newindex] = OneLevelDetails(wav);
    Val[newindex] = OneLevelEven(wav);
    int nconnect = intel->Reference()->NSides();
    seqnum = intel->Connect(nconnect).SequenceNumber();
    block(seqnum,0,var,0) = value;
    seqnum = intel2->Connect(nconnect).SequenceNumber();
    block(seqnum,0,var,0) = -value;
    newindex++;
  }
	an.GraphMesh(posdim)->SetNames(scal,vec);
  an.GraphMesh(posdim)->DrawSolution(step++,time);
  time += an.DeltaT();
  cmesh->Solution() = Solution;
  for(i=0;i<newindex;i++) {
    for(j=0;j<nsub;j++) subel[j] = NewIndex[nsub*i+j];
    cmesh->Coarsen(subel,index);
    NewIndex[i] = index;
  }
  cmesh->InitializeBlock();
  for(i=0;i<newindex;i++) {
    intel = (TPZInterpolatedElement *)MeshVec[NewIndex[i]];
    TPZConnect &con = intel->Connect(intel->Reference()->NSides());
    seqnum = con.SequenceNumber();
    block(seqnum,0,var,0) = Val[i];
  }
  an.ReDefineGraphMesh(posdim); //,scal,vec,an.GraphMesh(posdim));
  an.GraphMesh(posdim)->DrawMesh(MAXITERATIONS);
	an.GraphMesh(posdim)->SetNames(scal,vec);
  an.GraphMesh(posdim)->DrawSolution(step++,time);
  time += an.DeltaT();
  Solution = cmesh->Solution();
  for(i=0;i<newindex;i++) {
    intel =(TPZInterpolatedElement *)MeshVec[NewIndex[i]];
    TPZConnect &con = intel->Connect(intel->Reference()->NSides());
    seqnum = con.SequenceNumber();
    block(seqnum,0,var,0) = Wav[i];
  }
  if(!last)
    ((TPZDXGraphMesh *)(an.GraphMesh(posdim)))->SetNumCases(step+1);
	an.GraphMesh(posdim)->SetNames(scal,vec);
  an.GraphMesh(posdim)->DrawSolution(step++,time);
  time += an.DeltaT();
  cmesh->Solution() = Solution;
}
/*
REAL HaarOneLevelMean(TPZVec<REAL> &values) {
  return (values[0]+values[1])/M_SQRT2;
}
REAL HaarOneLevelDetails(TPZVec<REAL> &values) {
  return (values[0]-values[1])/M_SQRT2;
}
*/

THaar2D::THaar2D(double thread,double maxcoef) : TWavelet(thread,maxcoef) {
}

int THaar2D::Decomposition(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
        TPZManVector<int> &elsons,int MinLevel) {
  /**Considerar especialmente os elementos fronteira para ser refinados quando o
     seu vizinho do dominio interior eh de maior nivel que o seu geometrico*/
  int i, nelem, nfathers = elfathers.NElements(), nsons = elsons.NElements();
  int dimmat;
  std::map<int, TPZMaterial * >::const_iterator it;
  for (it = cmesh.MaterialVec().begin(); it != cmesh.MaterialVec().end(); it++) {
	  if (it->second->Id() < 0) continue;
	  dimmat = it->second->Dimension();
	  break;
  }
  nelem = cmesh.NElements();
  /**Para marcar elementos ja trabalhados*/
  TPZVec<int> indexel(nelem,-1);
  int nsub;
  int indexintel, indexbrother;
  /**Variaveis para calculo dos coeficientes wavelets*/
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  REAL djk;
  int j, p, indexfather = nfathers, indexson = nsons, mask;
//  if(nfathers || nsons) PZError << "HaarDecomposition. Already exist fathers or sons.\n";
  TPZVec<REAL> mean(2,0.);
  TPZCompEl *el;
  TPZInterpolatedElement *intel, *brotherel[4];

  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el || el->IsInterface() || el->Dimension()<dimmat) continue;
    intel = (TPZInterpolatedElement *)el;
    indexintel = intel->Index();
    if(!indexel[indexintel]) continue;
    mask = 0;   // suposse number of intel as subelement of the gel(father)
    TPZGeoEl *gel, *father = intel->Reference()->Father();
    if(!father) continue;
    if(!father->HasSubElement()) {
      PZError << "HaarDecomposition2d. This father element has no sons ("<<i<<".\n";
      continue;
    }
    /**Achando os comp elements into the more fine level.
       mask is the subelement number of the intel*/
    nsub = father->NSubElements();
    int coarse = 0, mcount = 0, indexcompel;
    for(j=0;j<nsub;j++) {
      brotherel[j] = (TPZInterpolatedElement *)father->SubElement(j)->Reference();
      if(brotherel[j]==intel) mask = j;
    }
    if(mask==nsub-1) {
      mean[0] = intel->MeanSolution(var);
      indexbrother = indexintel;
    }
    else if(!brotherel[nsub-1]) {
      gel = father->SubElement(nsub-1);
      mean[0] = MeanSolutionFromSubElements(gel,var);
      indexbrother = -1;
      coarse--;   //porque nao pode acontecer coarsing, nao existe um subelemento
    }
    else {
      indexbrother = brotherel[nsub-1]->Index();
      mean[0] = brotherel[nsub-1]->MeanSolution(var);
    /** Para evitar o calculo novamente neste sub-elemento*/
      indexel[indexbrother] = 0;
    }

    /**Determinando os coeficientes wavelet para realizar threadholder */
    TPZManVector<int> celindexes(10);
    for(j=0;j<nsub-1;j++) {
      TPZInterpolatedElement *compel = brotherel[j];
      if(!compel) {   // compel nao pode ser intel, intel sempre existe!!!
        gel = father->SubElement(j);
        mean[1] = MeanSolutionFromSubElements(gel,var);
        coarse--;   //porque nao pode acontecer coarsing, nao existe um subelemento
      }
      else {
        mean[1] = compel->MeanSolution(var);
        indexcompel = compel->Index();
        indexel[indexcompel] = 0;
      }
      djk = fabs(OneLevelDetails(mean));
      if(djk > MAXWAVELETCOEFFICIENT) ReAssembling = 1;   // ???
      if(djk > THREADHOLDER) {
//        if(intel->Reference()->Level() > MaxLevelRefinements-1) continue;
        if(brotherel[j]) celindexes[mcount++] = brotherel[j]->Index();
      }
      else if(IsZeroToCoarsing(djk) && intel->Reference()->Level()> MinLevel) coarse++;
    }
    /**Se existe algum subelemento refinado deve-se refinar o subelemento pivo tambem*/
    if(nsub>2 && mcount && brotherel[nsub-1]) celindexes[mcount++] = brotherel[nsub-1]->Index();
    if(nfathers<mcount+indexfather) {
      nfathers += 100;
      elfathers.Resize(nfathers);
    }
  /**passando os indices dos elementos a subdividir para elfathers*/
    for(p=0;p<mcount;p++)
      elfathers[indexfather++] = celindexes[p];
  /**preenchendo no vetor de elementos a emgrossar caso coarse seja 3*/
    if(coarse==nsub-1) {
      if(nsons<indexson+nsub) {
        nsons += 100;
        elsons.Resize(nsons);
      }
      for(p=0;p<nsub;p++)
        elsons[indexson++] = brotherel[p]->Index();
    }
  }
  elsons.Resize(indexson);
  elfathers.Resize(indexfather);
  return 0;
}

void THaar2D::DecompositionOnlyOrderZero(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
		 TPZManVector<int> &elsons,int MinLevel) {
  /**Considerar especialmente os elementos fronteira para ser refinados quando o
     seu vizinho do dominio interior eh de maior nivel que o seu geometrico*/
  int i, nelem, nfathers = elfathers.NElements(), nsons = elsons.NElements();
  int dimmat;
  std::map<int, TPZMaterial * >::const_iterator it;
  for (it = cmesh.MaterialVec().begin(); it != cmesh.MaterialVec().end(); it++) {
	  if (it->second->Id() < 0) continue;
	  dimmat = it->second->Dimension();
	  break;
  }
  nelem = cmesh.NElements();
  /**Para marcar elementos ja trabalhados*/
  TPZVec<int> indexel(nelem,-1);
  int nsub;
  int indexintel, indexbrother;
  /**Variaveis para calculo dos coeficientes wavelets*/
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  REAL djk;
  int j, p, indexfather = nfathers, indexson = nsons, mask;
  TPZVec<REAL> mean(2,0.);
  TPZCompEl *el;
  TPZInterpolatedElement *intel, *brotherel[4];

  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el || el->IsInterface() || el->Dimension()<dimmat) continue;
    intel = (TPZInterpolatedElement *)el;
    int nsides = intel->Reference()->NSides();
    /** Apenas serao decompostos os elementos de ordem zero */
    if(intel->PreferredSideOrder(nsides-1)) continue;
    indexintel = intel->Index();
    if(!indexel[indexintel]) continue;
    mask = 0;   // suposse number of intel as subelement of the gel(father)
    TPZGeoEl *gel, *father = intel->Reference()->Father();
    if(!father) continue;
    if(!father->HasSubElement()) {
      PZError << "HaarDecomposition2d. This father element has no sons ("<<i<<".\n";
      continue;
    }
    /**Achando os comp elements into the more fine level.
       mask is the subelement number of the intel*/
    nsub = father->NSubElements();
    int coarse = 0, mcount = 0, indexcompel;
    for(j=0;j<nsub;j++) {
      brotherel[j] = (TPZInterpolatedElement *)father->SubElement(j)->Reference();
      if(brotherel[j]==intel) mask = j;
    }
    if(mask==nsub-1) {
      mean[0] = intel->MeanSolution(var);
      indexbrother = indexintel;
    }
    else if(!brotherel[nsub-1]) {
      gel = father->SubElement(nsub-1);
      mean[0] = MeanSolutionFromSubElements(gel,var);
      indexbrother = -1;
      coarse--;   //porque nao pode acontecer coarsing, nao existe um subelemento
    }
    else {
      indexbrother = brotherel[nsub-1]->Index();
      mean[0] = brotherel[nsub-1]->MeanSolution(var);
    /** Para evitar o calculo novamente neste sub-elemento*/
      indexel[indexbrother] = 0;
    }

    /**Determinando os coeficientes wavelet para realizar threadholder */
    TPZManVector<int> celindexes(10);
    for(j=0;j<nsub-1;j++) {
      TPZInterpolatedElement *compel = brotherel[j];
      if(!compel) {   // compel nao pode ser intel, intel sempre existe!!!
        gel = father->SubElement(j);
        mean[1] = MeanSolutionFromSubElements(gel,var);
        coarse--;   //porque nao pode acontecer coarsing, nao existe um subelemento
      }
      else {
        mean[1] = compel->MeanSolution(var);
        indexcompel = compel->Index();
        indexel[indexcompel] = 0;
      }
      djk = fabs(OneLevelDetails(mean));
      if(djk > MAXWAVELETCOEFFICIENT) ReAssembling = 1;   // ???
      if(djk > THREADHOLDER) {
//        if(intel->Reference()->Level() > MaxLevelRefinements-1) continue;
        if(brotherel[j]) celindexes[mcount++] = brotherel[j]->Index();
      }
      else if(IsZeroToCoarsing(djk) && intel->Reference()->Level()>MinLevel) coarse++;
    }
    /**Se existe algum subelemento refinado deve-se refinar o subelemento pivo tambem*/
    if(nsub>2 && mcount && brotherel[nsub-1]) celindexes[mcount++] = brotherel[nsub-1]->Index();
    if(nfathers<mcount+indexfather) {
      nfathers += 100;
      elfathers.Resize(nfathers);
    }
  /**passando os indices dos elementos a subdividir para elfathers*/
    for(p=0;p<mcount;p++)
      elfathers[indexfather++] = celindexes[p];
  /**preenchendo no vetor de elementos a emgrossar caso coarse seja 3*/
    if(coarse==nsub-1) {
      if(nsons<indexson+nsub) {
        nsons += 100;
        elsons.Resize(nsons);
      }
      for(p=0;p<nsub;p++)
        elsons[indexson++] = brotherel[p]->Index();
    }
  }
  elsons.Resize(indexson);
  elfathers.Resize(indexfather);
}

REAL THaar2D::OneLevelEven(TPZVec<REAL> &values) {
  return .25*(values[0]+values[1]+values[2]+values[3]);
}
REAL THaar2D::OneLevelDetails(TPZVec<REAL> &values) {
  return .375*(values[1]-values[0]);  // di = 3/8*(ai - a0)
}

TSchauder1D::TSchauder1D(double thread,double maxcoef) : TWavelet(thread,maxcoef) {
}

REAL TSchauder1D::OneLevelEven(REAL value) {
  return value;
}
REAL TSchauder1D::OneLevelDetails(TPZVec<REAL> &values) {
  return (values[1] - .5*(values[0]+values[2]));
}

int TSchauder1D::Decomposition(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
     TPZManVector<int> &elsons,int MinLevel) {
  int i, order, nelem, nfathers = elfathers.NElements(), nsons = elsons.NElements();
  int dimmat,existorderzero = 0;
  std::map<int, TPZMaterial * >::const_iterator it;
  for (it = cmesh.MaterialVec().begin(); it != cmesh.MaterialVec().end(); it++) {
	  if (it->second->Id() < 0) continue;
	  dimmat = it->second->Dimension();
	  order = it->second->NStateVariables();
	  break;
  }
  nelem = cmesh.NElements();
  /**Para marcar elementos ja trabalhados*/
  TPZVec<int> indexel(nelem,-1);
  int indexintel, indexbrother;
  /**Variaveis para calcular os coeficientes wavelets*/
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  REAL djk;
//  int nconnect;
  int seqnum, indexvar, indexfather = nfathers, indexson = nsons;
  TPZVec<REAL> wav(3,0.);
  TPZCompEl *el;
  TPZInterpolatedElement *intel, *brotherel, *brother;
  TPZBlock<STATE> &block = cmesh.Block();

  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el || el->IsInterface() || el->Dimension() != 1) continue;
    intel = (TPZInterpolatedElement *)el;
    indexintel = intel->Index();
    if(!indexel[indexintel]) continue;
    int mask = 1, nsides = intel->Reference()->NSides();
    /** Apenas serao decompostos os elementos de ordem nao zero */
    if(!intel->PreferredSideOrder(nsides-1)) {
      existorderzero = 1;   //Eh melhor trocar por um vetor passado para Haar soh
      continue;
    }
    TPZGeoEl *gel = intel->Reference()->Father();
    if(!gel) continue;
    if(!gel->HasSubElement()) {
      PZError << "SchauderDecomposition. Geometrical father without sons.\n";
      continue;
    }
    brotherel = (TPZInterpolatedElement *)gel->SubElement(1)->Reference();
    if(brotherel==intel) { //intel->ref tem lado 1 comum com o seu father (mask=0)
      brotherel = (TPZInterpolatedElement *)gel->SubElement(0)->Reference();
      mask = 0;
    }
    brother = brotherel;
    /**Determinando o valor esquina do subelemento irmao*/
    if(!brotherel) {
      gel = gel->SubElement(mask);
      do {
        gel = gel->SubElement(mask);
      } while(!gel->Reference());   // if(!gel) error
      brotherel = (TPZInterpolatedElement *)gel->Reference();
    }
    else indexbrother = brother->Index();
    /** Determinando os coeficientes wavelet para realizar threadholder */
//    if(intel->IsConnectContinuous(0)) {
//      seqnum = intel->Connect(0).SequenceNumber();
//      indexvar = 0;
//    }
//    else {
      seqnum = intel->Connect(3).SequenceNumber();
      indexvar = order;
 //   }
    wav[(mask+1)%2] = block(seqnum,0,var,0);
//    if(intel->IsConnectContinuous(1)) {
 //     seqnum = intel->Connect(1).SequenceNumber();
   //   wav[(mask+1)%2+1] = block(seqnum,0,var,0);
     // indexvar = 0;
//    }
  //  else {
      seqnum = intel->Connect(3).SequenceNumber();
      wav[(mask+1)%2+1] = block(seqnum,0,indexvar+var,0);
   // }
    indexvar = 0;
  //  if(brotherel->IsConnectContinuous(mask))
    //  seqnum = brotherel->Connect(mask).SequenceNumber();
//    else {
//      if(mask && !brotherel->IsConnectContinuous(0))
	if(mask) indexvar = order;
      seqnum = brotherel->Connect(3).SequenceNumber();
 //   }
    wav[2*mask] = block(seqnum,0,indexvar+var,0);
    djk = fabs(OneLevelDetails(wav));
    if(djk > MAXWAVELETCOEFFICIENT) ReAssembling = 1;   // ???
    if(djk > THREADHOLDER) {
/**Tarefa deixada para CheckVectors*/
//      if(intel->Reference()->Level() > MaxLevelRefinements-1) continue;
      if(nfathers-(2)<indexfather) {
        nfathers += 100;
        elfathers.Resize(nfathers);
      }
    /**Pode-se evitar o indice se no Divide do elemento computacional tira-se o parametro index*/
      elfathers[indexfather++] = indexintel;
      if(brother) elfathers[indexfather++] = indexbrother;
    }
    else if(IsZeroToCoarsing(djk) && brother) {
      if(intel->Reference()->Level() <= MinLevel) continue;
      if(nsons<indexson+2) {
        nsons += 100;
        elsons.Resize(nsons);
      }
      elsons[indexson++] = indexintel;
      elsons[indexson++] = indexbrother;
    }
  }
  elsons.Resize(indexson);
  elfathers.Resize(indexfather);
  return existorderzero;
}


void TSchauder1D::DrawWavelets(int var,TTimeAnalysis &an,int posdim,int &step,
		      REAL &time,TPZVec<std::string> &scal,TPZVec<std::string> &vec,int last) {
  if(an.Mesh()->GetDefaultOrder()<1) {
    PZError << "DrawSchauderWavelets. TPZCompEl::gOrder is lower than 1.\n";
    return;
  }
  TPZCompMesh* cmesh = an.Mesh();
  if(!an.GraphMesh(posdim) || !step) return;
  int i, j;
  int seqnum, order = an.OrderLaw(), nsub = cmesh->Reference()->ElementVec()[0]->NSubElements();
  int64_t index, newindex = 0, nelem = cmesh->NElements();
  TPZVec<REAL> Wav(nelem,0.);
  TPZVec<REAL> Val(2*nelem,0.);
  TPZBlock<STATE> & block = cmesh->Block();
  TPZManVector<TPZCompEl *> ElVec(nelem,0);
  TPZAdmChunkVector<TPZCompEl *> &MeshVec = cmesh->ElementVec();
  TPZCompEl *cel;
  TPZGeoEl *gel;
  for(i=0;i<nelem;i++) {
    cel = MeshVec[i];
    if(!cel || cel->IsInterface() || cel->Dimension()!=1) continue;
    ElVec[i] = cel;
  }

  TPZManVector<int64_t> subel(nsub,-1);  // index of the sons computational elements
  TPZVec<REAL> wav(nsub+1);
  int nconnect, seqnumwav, varwav;
  TPZManVector<int> IntelIndexes(nelem,-1);
  TPZVec<int> seqnumzero(cmesh->NConnects(),-1);
  int seqnumindex = 0;
  TPZInterpolatedElement *intel, *intel2;
  TPZFMatrix<STATE> Solution(cmesh->Solution());
  /**Taking the indexes of sub-elements with same father and filling comp
     mesh with wavelets basis of the previous level*/
  for(i=0;i<nelem;i++) {
    intel = (TPZInterpolatedElement *)ElVec[i];
    if(!intel) continue;
    gel = intel->Reference()->Father();
    intel = (TPZInterpolatedElement *)gel->SubElement(0)->Reference();
    IntelIndexes[2*newindex] = intel->Index();
    nconnect = intel->Reference()->NSides();
    int indexvar = 0, seqnumdisc = -1;
 //   if(intel->IsConnectContinuous(0)) {
 //     seqnum = intel->Connect(0).SequenceNumber();
 //     seqnumzero[seqnumindex++] = seqnum;
 //     wav[0] = block(seqnum,0,var,0);
  //  }
    //else {
      nconnect = intel->Reference()->NSides();
      seqnum = intel->Connect(nconnect).SequenceNumber();
      indexvar = order;
      wav[0] = block(seqnum,0,var,0);
      block(seqnum,0,var,0) = 0.;
   // }
    Val[2*newindex] = OneLevelEven(wav[0]);
//    if(intel->IsConnectContinuous(1)) {
  //    seqnumwav = intel->Connect(1).SequenceNumber();
 //     varwav = var;
   // }
   // else {
      seqnumwav = intel->Connect(nconnect).SequenceNumber();
      varwav = indexvar + var;
   // }
    wav[1] = block(seqnumwav,0,varwav,0);
    intel2 = (TPZInterpolatedElement *)gel->SubElement(1)->Reference();
    IntelIndexes[2*newindex+1] = intel2->Index();
//    if(intel2->IsConnectContinuous(0)) indexvar = 0;
  //  else {
      indexvar = order;
      seqnumdisc = intel2->Connect(nconnect).SequenceNumber();
    //}
//    if(intel2->IsConnectContinuous(1)) {
 //     seqnum = intel2->Connect(1).SequenceNumber();
 //     seqnumzero[seqnumindex++] = seqnum;
 //     wav[2] = block(seqnum,0,var,0);
 //   }
   // else {
      seqnum = intel2->Connect(nconnect).SequenceNumber();
      wav[2] = block(seqnum,0,indexvar+var,0);
      block(seqnum,0,indexvar+var,0) = 0.;
    //}
    Val[2*newindex+1] = OneLevelEven(wav[2]);
    Wav[newindex] = block(seqnumwav,0,varwav,0) = OneLevelDetails(wav);
    if(seqnumdisc!=-1) block(seqnumdisc,0,var,0) = Wav[newindex];
    ElVec[IntelIndexes[2*newindex]] = 0;
    ElVec[IntelIndexes[2*newindex+1]] = 0;
    newindex++;
  }
  for(i=0;i<seqnumindex;i++)
    block(seqnumzero[i],0,var,0) = 0.;
	an.GraphMesh(posdim)->SetNames(scal,vec);
  an.GraphMesh(posdim)->DrawSolution(step++,time);
  time += an.DeltaT();
  cmesh->Solution() = Solution;
  for(i=0;i<newindex;i++) {
    for(j=0;j<nsub;j++) subel[j] = IntelIndexes[nsub*i+j];
    cmesh->Coarsen(subel,index);
    IntelIndexes[i] = index;
  }
  cmesh->InitializeBlock();
  for(i=0;i<newindex;i++) {
    intel = (TPZInterpolatedElement *)MeshVec[IntelIndexes[i]];
    int indexvar = 0;
//    if(intel->IsConnectContinuous(0))
  //    seqnum = intel->Connect(0).SequenceNumber();
    //else {
      seqnum = intel->Connect(nconnect).SequenceNumber();
      indexvar = order;
    //}
    block(seqnum,0,var,0) = Val[2*i];
//    if(intel->IsConnectContinuous(1))
  //    seqnum = intel->Connect(1).SequenceNumber();
    //else
      seqnum = intel->Connect(nconnect).SequenceNumber();
    block(seqnum,0,indexvar+var,0) = Val[2*i+1];
  }
  an.ReDefineGraphMesh(posdim); //,scal,vec,an.GraphMesh(posdim));
  an.GraphMesh(posdim)->DrawMesh(MAXITERATIONS);
	an.GraphMesh(posdim)->SetNames(scal,vec);
  an.GraphMesh(posdim)->DrawSolution(step++,time);
  time += an.DeltaT();
  Solution = cmesh->Solution();
  for(i=0;i<newindex;i++) {
    intel =(TPZInterpolatedElement *)MeshVec[IntelIndexes[i]];
    int indexvar = 0;
 //   if(intel->IsConnectContinuous(0))
   //   seqnum = intel->Connect(0).SequenceNumber();
   // else {
      seqnum = intel->Connect(nconnect).SequenceNumber();
      indexvar = order;
   // }
    block(seqnum,0,var,0) = Wav[i];
 //   if(intel->IsConnectContinuous(1))
   //   seqnum = intel->Connect(1).SequenceNumber();
   // else
      seqnum = intel->Connect(nconnect).SequenceNumber();
    block(seqnum,0,indexvar+var,0) = Wav[i];
  }
  if(!last)
    ((TPZDXGraphMesh *)(an.GraphMesh(posdim)))->SetNumCases(step+1);
	an.GraphMesh(posdim)->SetNames(scal,vec);
  an.GraphMesh(posdim)->DrawSolution(step++,time);
  time += an.DeltaT();
  cmesh->Solution() = Solution;
}

void TSchauder1D::Preconditioning(int /*var*/,int /*dim*/,TPZFMatrix<STATE> &/*pre*/,int /*nlevels*/) {

}

/*
REAL SchauderOneLevelEven(REAL value) {
  return M_SQRT2*value;
}
REAL SchauderOneLevelDetails(TPZVec<REAL> &values) {
  return M_SQRT2*(values[1] - .5*(values[0]+values[2]));
}
*/

TSchauder2D::TSchauder2D(double thread,double maxcoef) : TWavelet(thread,maxcoef) {
}

REAL TSchauder2D::OneLevelEven(REAL value) {
  return value;
}
REAL TSchauder2D::OneLevelDetails(TPZVec<REAL> &values) {
  return (values[1] - .5*(values[0]+values[2]));
}

int TSchauder2D::Decomposition(TPZCompMesh &cmesh,int var,TPZManVector<int> &elfathers,
        TPZManVector<int> &elsons,int MinLevel) {
  int i, j, order, nelem = cmesh.NMaterials(), nfathers = elfathers.NElements(), nsons = elsons.NElements();
  int existorderzero = 0, dimmat;
  std::map<int, TPZMaterial * >::const_iterator it;
  for (it = cmesh.MaterialVec().begin(); it != cmesh.MaterialVec().end(); it++) {
	  if (it->second->Id() < 0) continue;
	  dimmat = it->second->Dimension();
	  order = it->second->NStateVariables();
	  break;
  }
  nelem = cmesh.NElements();
  /**Para marcar elementos ja trabalhados*/
  TPZVec<int> indexel(nelem,-1);
  int nsub, numberconnect;
  int indexintel;
  /**Variaveis para calcular os coeficientes wavelets*/
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  REAL djk;
  int p, seqnum, indexvar, indexfather = nfathers, indexson = nsons;
  TPZVec<REAL> wav(3,0.);
  TPZCompEl *el;
  TPZInterpolatedElement *intel, *brotherel[4];
  TPZBlock<STATE> &block = cmesh.Block();

  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el || el->IsInterface() || el->Dimension()<dimmat) continue;
    intel = (TPZInterpolatedElement *)el;
    indexintel = intel->Index();
    if(!indexel[indexintel]) continue;
//    int mask = 0;
    int nsides = intel->Reference()->NSides();
    /** Apenas serao decompostos os elementos de ordem nao zero */
    if(!intel->PreferredSideOrder(nsides-1)) {   // Deve ser feito para cada sub-elemento
      existorderzero = 1;   //Eh melhor trocar por um vetor passado para Haar soh
      continue;
    }
    TPZGeoEl *gel, *father = intel->Reference()->Father();
    if(!father) continue;
    if(!father->HasSubElement()) {
      PZError << "SchauderDecomposition2d. It has no different levels.\n";
      continue;
    }
    nsub = father->NSubElements();
    if(nsub==2) numberconnect = 3;
    else numberconnect = 7;
    int coarse = 0, mcount = 0, indexcompel;
    for(j=0;j<nsub;j++) {
      brotherel[j] = (TPZInterpolatedElement *)father->SubElement(j)->Reference();
//      if(brotherel[j]==intel) mask = j;
    }
    /** Determinando os coeficientes wavelet para realizar threadholder */
    TPZManVector<int> celindexes(10);
    for(j=0;j<nsub-1;j++) {
      TPZInterpolatedElement *compel = brotherel[j];
      if(!compel) {   //compel nao pode ser intel
        gel = father->SubElement(j);
        do {
          gel = gel->SubElement(j);
        }while(!gel->Reference());
        compel = (TPZInterpolatedElement *)gel->Reference();
        indexcompel = -1;
        coarse--;
      }
      else {
        indexcompel = compel->Index();    // ???
        indexel[indexcompel] = 0;
      }
//      if(compel->IsConnectContinuous(j)) {
  //      seqnum = compel->Connect(j).SequenceNumber();
    //    wav[0] = block(seqnum,0,var,0);
      //}
 //     else {
        indexvar = 0;
        for(p=0;p<j;p++)
      //    if(!compel->IsConnectContinuous(p))
			indexvar += order;
        seqnum = compel->Connect(numberconnect).SequenceNumber();
        wav[0] = block(seqnum,0,indexvar+var,0);
   //   }
      int num;
      if(nsub==2) {
        compel = brotherel[1];
        num = 0;   // precisa-se os valores nos pontos medios dos segmentos
        if(!compel) { coarse--;
          gel = father->SubElement(1);
          do { gel = gel->SubElement(num);
          }while(!gel->Reference());
          compel = (TPZInterpolatedElement *)gel->Reference();
        }
				else
				  indexel[compel->Index()] = 0;
      }
      else {
        compel = brotherel[3];
        num = (j+2)%3;
        if(!compel) {
          coarse--;
          gel = father->SubElement(3);
          do {
            gel = gel->SubElement(num);
          }while(!gel->Reference());
          compel = (TPZInterpolatedElement *)gel->Reference();
        }
				else indexel[compel->Index()] = 0;
      }
//      if(compel->IsConnectContinuous(num)) {
  //      seqnum = compel->Connect(num).SequenceNumber();
    //    wav[1] = block(seqnum,0,var,0);
    //  }
    //  else {
        indexvar = 0;
        for(p=0;p<num;p++)
       //   if(!compel->IsConnectContinuous(p))
			indexvar += order;
        seqnum = compel->Connect(numberconnect).SequenceNumber();
        wav[1] = block(seqnum,0,indexvar+var,0);
      //}
      if(nsub==2) num = 1;
      else num = (j+1)%3;
      compel = brotherel[num];
      if(!compel) {
        coarse--;
        gel = father->SubElement(num);
        do {
          gel = gel->SubElement(num);
        }while(!gel->Reference());
        compel = (TPZInterpolatedElement *)gel->Reference();
      }
			else indexel[compel->Index()] = 0;
//      if(compel->IsConnectContinuous(num)) {
  //      seqnum = compel->Connect(num).SequenceNumber();
    //    wav[2] = block(seqnum,0,var,0);
      //}
    //  else {
        indexvar = 0;
        for(p=0;p<num;p++)
         // if(!compel->IsConnectContinuous(p)) 
			  indexvar += order;
        seqnum = compel->Connect(numberconnect).SequenceNumber();
        wav[2] = block(seqnum,0,indexvar+var,0);
    //  }
      djk = fabs(OneLevelDetails(wav));
      if(djk > MAXWAVELETCOEFFICIENT) ReAssembling = 1;   // ???
      if(djk > THREADHOLDER) {
//        if(intel->Reference()->Level() > MaxLevelRefinements-1) continue;
        if(brotherel[j]) celindexes[mcount++] = brotherel[j]->Index();
        if(brotherel[(j+1)%3]) celindexes[mcount++] = brotherel[(j+1)%3]->Index();
      }
      else if(IsZeroToCoarsing(djk) && intel->Reference()->Level()>MinLevel) coarse++;
    }
    /**Se existe algum subelemento refinado deve-se refinar o subelemento pivo tambem*/
    if(nsub>2 && mcount && brotherel[nsub-1]) celindexes[mcount++] = brotherel[nsub-1]->Index();
    if(nfathers<mcount+indexfather) {
      nfathers += 100;
      elfathers.Resize(nfathers);
    }
  /**passando os indices dos elementos a subdividir para elfathers*/
    for(p=0;p<mcount;p++)
      elfathers[indexfather++] = celindexes[p];
  /**preenchendo no vetor de elementos a emgrossar caso coarse seja 3*/
    if(coarse==nsub-1) {
      if(nsons<nsub+indexson) {
        nsons += 100;
        elsons.Resize(nsons);
      }
      for(p=0;p<nsub;p++)
        elsons[indexson++] = brotherel[p]->Index();
    }
  }
  elsons.Resize(indexson);
  elfathers.Resize(indexfather);
  return existorderzero;
}

void TSchauder2D::DrawWavelets(int var,TTimeAnalysis &an,int posdim,int &step,
			  REAL &time,TPZVec<std::string> &scal,TPZVec<std::string> &vec,int last) {
  TPZCompMesh *cmesh = an.Mesh();
  if (!cmesh || cmesh->GetDefaultOrder() < 1) {
	  PZError << "DrawSchauder2dWavelets. TPZCompEl::gOrder is lower than 1.\n";
	  return;
  }
  if(!an.GraphMesh(posdim) || !step) return;
  int i, j;
  int seqnum, nsub = 4;
  int64_t index, newindex = 0, nelem = cmesh->NElements();
  TPZVec<REAL> Wav(nelem,0.);
  TPZVec<REAL> Val(nelem,0.);
  TPZBlock<STATE> &block = cmesh->Block();
  TPZBlock<STATE> bloco(cmesh->Block());
  TPZManVector<TPZCompEl *> ElVec(nelem,0);
  TPZInterpolatedElement *intel, *intel2;
  TPZAdmChunkVector<TPZCompEl *> &MeshVec = cmesh->ElementVec();
  TPZCompEl *cel;
  TPZGeoEl *gel;
  int64_t indice;
  for(i=0;i<nelem;i++) {
    cel = MeshVec[i];
    if(!cel || cel->IsInterface() || cel->Dimension()!=2) continue;
    ElVec[i] = cel;
    index = cel->Reference()->NSides();
	an.Mesh()->Discontinuous2Continuous(cel->Index(), indice);
//    for(j=0;j<index;j++) ((TPZInterpolatedElement *)cel)->MakeConnectContinuous(j);
  }

  TPZManVector<int64_t> subel(nsub,-1);  // index of the sons computational elements
  TPZVec<REAL> wav(nsub-1);
  TPZVec<REAL> value(nsub-1);
  TPZVec<REAL> wacopy(nsub-1);
  TPZManVector<int> IntelIndexes(nelem,-1);
  TPZFMatrix<STATE> Solution(cmesh->Solution());
  bloco.SetMatrix(&Solution);
  /**Taking the indexes of sub-elements with same father and filling comp
     mesh with wavelets basis of the previous level*/
  for(i=0;i<nelem;i++) {
    intel = (TPZInterpolatedElement *)ElVec[i];
    if(!intel) continue;
    gel = intel->Reference()->Father();
    intel2 = (TPZInterpolatedElement *)gel->SubElement(nsub-1)->Reference();
    for(j=0;j<nsub-1;j++) {
      intel = (TPZInterpolatedElement *)gel->SubElement(j)->Reference();
      IntelIndexes[nsub*newindex+j] = intel->Index();
      seqnum = intel->Connect(j).SequenceNumber();
      wav[j] = block(seqnum,0,var,0);
      bloco(seqnum,0,var,0) = 0.;
      Val[nsub*newindex+j] = OneLevelEven(wav[j]);
      seqnum = intel2->Connect(j).SequenceNumber();
      value[j] = block(seqnum,0,var,0);
    }
    IntelIndexes[nsub*newindex+j] = intel2->Index();
    for(j=0;j<nsub-1;j++) {
      wacopy[0] = wav[j];
      wacopy[1] = value[j];
      wacopy[2] = wav[(j+1)%3];
      seqnum = intel2->Connect(j).SequenceNumber();
      bloco(seqnum,0,var,0) = Wav[3*newindex+j] = OneLevelDetails(wacopy);
    }
    for(j=0;j<nsub;j++)
      ElVec[IntelIndexes[nsub*newindex+j]] = 0;
    newindex++;
  }
  TPZFMatrix<STATE> Sol(cmesh->Solution());
  cmesh->Solution() = Solution;
	an.GraphMesh(posdim)->SetNames(scal,vec);
  an.GraphMesh(posdim)->DrawSolution(step++,time);
  time += an.DeltaT();
  cmesh->Solution() = Sol;
  for(i=0;i<newindex;i++) {
    for(j=0;j<nsub;j++) subel[j] = IntelIndexes[nsub*i+j];
    cmesh->Coarsen(subel,index);
    IntelIndexes[i] = index;
  }
  cmesh->InitializeBlock();
  for(i=0;i<newindex;i++) {
    intel = (TPZInterpolatedElement *)MeshVec[IntelIndexes[i]];
    for(j=0;j<3;j++) {
      seqnum = intel->Connect(j).SequenceNumber();
      block(seqnum,0,var,0) = Val[nsub*i+j];
    }
  }
  an.ReDefineGraphMesh(posdim);
  an.GraphMesh(posdim)->DrawMesh(MAXITERATIONS);
/*  an.GraphMesh(posdim)->DrawSolution(step++,time,scal,vec);
  time += an.DeltaT();
  Solution = cmesh->Solution();
  for(i=0;i<newindex;i++) {
    intel =(TPZInterpolatedElement *)MeshVec[IntelIndexes[i]];
    int indexvar = 0;
    if(intel->IsConnectContinuous(0))
      seqnum = intel->Connect(0).SequenceNumber();
    else {
      seqnum = intel->Connect(nconnect).SequenceNumber();
      indexvar = order;
    }
    block(seqnum,0,var,0) = Wav[i];
    if(intel->IsConnectContinuous(1))
      seqnum = intel->Connect(1).SequenceNumber();
    else
      seqnum = intel->Connect(nconnect).SequenceNumber();
    block(seqnum,0,indexvar+var,0) = Wav[i];
  }*/
  if(!last)
    ((TPZDXGraphMesh *)(an.GraphMesh(posdim)))->SetNumCases(step+1);
	an.GraphMesh(posdim)->SetNames(scal,vec);
  an.GraphMesh(posdim)->DrawSolution(step++,time);
  time += an.DeltaT();
//  cmesh->Solution() = Solution;
}


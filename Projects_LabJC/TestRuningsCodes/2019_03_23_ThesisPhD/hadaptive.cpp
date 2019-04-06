#include <stdlib.h>
#include <math.h>

#include "hadaptive.h"
#include "myheader.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzmanvector.h"
#include "pzelg1d.h"
#include "pzelgt2d.h"
#include "pzelgpoint.h"
#include "elcpointgd.h"
#include "interface.h"
#include "elc1dgd.h"
#include "elct2dgd.h"
#include "pzadmchunk.h"
#include "pzadmchunk.h"
#include "tarungek.h"
#include "pzgraphmesh.h"
#include "pzdxmesh.h"
#include "conslaw.h"
#include "wavelet.h"

TAdaptive::TAdaptive(std::istream &input,int dim,int hpadap) : elfathers(0), elsons(0),
    hpadaptive(hpadap), fWaveletSec(NULL) {

  double thread, maxcoef;
  /** Reading PARAMETERS */
  /**To Refinement - Maxime level of refinements*/
  GetDataCommented(input,MaxLevelRefinements);
  /** Parameters to refinement wavelet techniques */
  GetDataCommented(input,thread);
  GetDataCommented(input,maxcoef);

	int index;
  /**The wavelet type is used to modified numerical flux over interfaces */
  GetDataCommented(input,index);   // (0) Haar wavelet, (1) Schauder wavelet
	if(index==1) {
	  if(dim==2) {
		  fWavelet = new TSchauder2D(thread,maxcoef);
			fWaveletSec = new THaar2D(thread,maxcoef);
		}
		else {
		  fWavelet = new TSchauder1D(thread,maxcoef);
  		fWaveletSec = new THaar1D(thread,maxcoef);
		}
	}
	else {
	  if(dim==2) fWavelet = new THaar2D(thread,maxcoef);
		else fWavelet = new THaar1D(thread,maxcoef);
	}

  GetDataCommented(input,segurityzone);   // segurityzone!=1 garante considerar vizinhos
  GetDataCommented(input,p_decrement);
  GetDataCommented(input,p_increment);
	fSteadyState = 1;
}

TAdaptive::~TAdaptive() {
  if(fWavelet) delete fWavelet;
}

/** idgfather has the index of the elements that must to be refined */
void TAdaptive::Refinements(TPZCompMesh &cmesh,TPZManVector<TPZGeoEl *> &idgfather,int ncycles) {
  int ngelem = idgfather.NElements();
  if(!ngelem || !ncycles) return;
  int i, mask, index = 0;
  elfathers.Resize(ngelem);
  TPZCompEl *cel;
  /** We can to have geo. elements duplicates, then we need to mask them */
  TPZManVector<int> cels(cmesh.NElements(),-1);
  for(i=0;i<ngelem;i++) {
    cel = idgfather[i]->Reference();
    if(!cel) continue;
    mask = cel->Index();
    if(cels[mask]==-1) {
      cels[mask] = mask;
      elfathers[index++] = mask;
    }
  }
  elfathers.Resize(index);
  Refinements(cmesh,ncycles);
}

/**Necessita testar nivel de refinamento para refinar. O nivel nao pode ser maior de
   MAXLEVELREFINEMENTS definido no arquivo hadaptive.h*/
void TAdaptive::Refinements(TPZCompMesh &cmesh,int ncycles) {
  if(!ncycles) return;
  int nelem = elfathers.NElements(), nbounds = 0;
  if(!nelem) return;
  int i, j, k, n = 0, index, nsub, order = cmesh.GetDefaultOrder();
#ifndef NOTDEBUG
  /**Limpando indexcel*/
  TPZManVector<int> indexcelcopy(elfathers);
  for(i=0;i<nelem;i++) {
    index = indexcelcopy[i];
    if(index == -1) continue;
    elfathers[n++] = index;
    /**Limpando os indeces repetidos a index*/
    for(j=i+1;j<nelem;j++)
      if(index==indexcelcopy[j]) indexcelcopy[j]=-1;
  }
  /**Entao existem n indeces de elem. comp. diferentes e todos validos*/
  if(!n) return;
  nelem = n;
  elfathers.Resize(nelem);
#endif

  /**Depending on the material dimension choose the maxime number of sub-elements*/
  n = cmesh.NMaterials();
  TPZMaterial *mat = 0;
  int dimmat;
  std::map<int, TPZMaterial * >::const_iterator mit;
  for (mit = cmesh.MaterialVec().begin(); mit != cmesh.MaterialVec().end(); mit++) {
	  TPZMaterial *mat = mit->second;
    if(mat->Id()<0) continue;
    dimmat = mat->Dimension();
    if(dimmat==1) nsub = 2;
    else nsub = 4;
    break;
  }
  if(!mat) PZError << "Refinements. Is not found materials.\n";
  /**Inicia refinamento*/
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  TPZCompEl *cel;
  int ncon, nnewcel = nsub*nelem;
  elfathers.Resize(nnewcel,-1);
  TPZVec<int64_t> indexes(0);
  TPZVec<int> celbounds(nnewcel,-1);
  TPZInterpolatedElement *intel;

  for(j=0;j<nelem;j++) {
    index = elfathers[j];
    cel = elvec[index];
    if(!cel || cel->IsInterface() || cel->Dimension()<dimmat) {
      elfathers[j] = -1;
      continue;
    }
    intel = (TPZInterpolatedElement *)cel;
    if(intel->Reference()->Level() > MaxLevelRefinements - 1) {
      elfathers[j] = -1;
      continue;
    }
    ncon = intel->Reference()->NSides()-1;
    for(i=0;i<ncon;i++) {
      TPZStack<TPZCompElSide> celneighs;
      TPZCompElSide celside(intel,i);
      if(celside.Reference().Dimension()<dimmat-1) continue;   //apenas para elementos de codimensao 1
      celside.EqualLevelElementList(celneighs,1,0);
      n = celneighs.NElements();
      if(!n) continue;
      for(k=0;k<n;k++) {
        TPZCompEl *compelstack = celneighs[k].Element();
        if(compelstack->Material()->Id()<0) {
          celbounds[nbounds++] = compelstack->Index();
	  if(nbounds==nnewcel) {
	    nnewcel *= 2;
	    PZError << "Refinements. Elementos fronteira excedem capacidade." << std::endl;
	    celbounds.Resize(nnewcel,-1);
	  }
	}
      }
    }
    nsub = intel->Reference()->NSubElements();
    cmesh.SetDefaultOrder(intel->PreferredSideOrder(ncon));
//		REAL coef = intel->DiffusionCoefficient();
    cmesh.Divide(index,indexes,1);
    for(k=0;k<nsub;k++) {
//      ((TPZInterpolatedElement *)elvec[indexes[k]])->SetDiffusionCoefficient(coef);
      elfathers[k*nelem+j] = indexes[k];
   }
  }
  /**Subdividindo os elementos fronteira achados */
  for(j=0;j<nbounds;j++) {
    index = celbounds[j];
    intel = (TPZInterpolatedElement *)(elvec[index]);
    if(!intel) continue;
    ncon = intel->Reference()->NSides()-1;
//		REAL coef = intel->DiffusionCoefficient();
	cmesh.SetDefaultOrder(intel->PreferredSideOrder(ncon));
	nsub = intel->Reference()->NSubElements();
    cmesh.Divide(index,indexes,1);
//    for(k=0;k<nsub;k++)
//      ((TPZInterpolatedElement *) elvec[indexes[k]])->SetDiffusionCoefficient(coef);
  }
  /**Subdividindo os elementos fronteira apenas se os elementos vizinhos no
     interior do dominio tem sido ja subdivididos
  for(j=0;j<nelem;j++) {
    index = celbounds[j];
    if(index==-1) continue;
    intel = (TPZInterpolatedElement *)elvec[index];
    TPZCompElSide celside(intel,intel->NSides()-1);
    TPZStack<TPZCompElSide> celneighs(0);
    celside.EqualLevelElementList(celneighs,1,0);
    if(celneighs.NElements()) continue;
    nsub = intel->Reference()->NSubElements();
    ncon = intel->NSides()-1;
    TPZCompEl::gOrder = intel->SideOrder(ncon);
    cmesh.Divide(index,indexes,1);
    for(k=0;k<nsub;k++) indexcel[k*nelem+j] = indexes[k];
  }*/
  for(i=1;i<ncycles;i++)
    Refinements(cmesh,1);

  elfathers.Resize(0);
  cmesh.SetDefaultOrder(order);
//  cmesh.CleanInterfaces();
  cmesh.InitializeBlock();
}

void TAdaptive::Refinements(TPZCompMesh &cmesh,int ncycles,int minelem,int maxelem) {
  if(!ncycles) return;
  int j,cycle;
  int nelem,ninterp;
  TPZManVector<int64_t> subindex;
  TPZCompEl *el;
  TPZVec<int> vec;
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  for(cycle=0;cycle<ncycles;cycle++) {
    if(cycle==ncycles-1) {
      TPZGeoEl1d::SetCreateFunction(&(TCompEl1dWI::CreateElDiscWI));
      TPZGeoElPoint::SetCreateFunction(&(TCompElPointWI::CreateElDiscWI));
      TPZGeoElT2d::SetCreateFunction(&(TCompElT2dWI::CreateElDiscWI));
    }
    ninterp = 0;
    nelem = cmesh.NElements();
    if(nelem<maxelem) PZError << "Refinement error.\n";
    vec.Resize(nelem);
    for(j=0;j<nelem;j++) {
      el = elvec[j];
      if(!el || el->IsInterface()) continue;
      if(el->Reference()->Level() >= MaxLevelRefinements) continue;
      vec[ninterp++] = j;
    }
    for(j=minelem;j<maxelem;j++)
      cmesh.Divide(vec[j],subindex,1);
  }
}

/** Necessita testar nivel de refinamento para fazer elemento coarse(grosso)*/
void TAdaptive::Coarsing(TPZCompMesh &cmesh) {
  int nelem = elsons.NElements(), dim;
  if(!nelem) return;
  int i, j, p, nsub;
  int64_t index, n = 0;
  /**Depending on the material dimension choose the maxime number of sub-elements*/
  std::map<int, TPZMaterial * >::const_iterator mit;
  TPZMaterial * mat = 0;
  for (mit = cmesh.MaterialVec().begin(); mit != cmesh.MaterialVec().end(); mit++) {
	  mat = mit->second;
	  if (mat->Id() < 0) continue;
	  dim = mat->Dimension();
	  if (dim == 1) nsub = 2;
	  else nsub = 4;
	  break;
  }
  if(!mat) { PZError << "Refinements. Not found materials.\n"; return; }
  /**Inicia coarsing*/
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  TPZCompEl *cel;
  TPZInterpolatedElement *intel, *el;
  TPZManVector<int64_t> indexes(nsub);
  int counter = 0, minorder, order, orderold= cmesh.GetDefaultOrder();
  int nvar = mat->NStateVariables();
  
  for(j=0;j<nelem;j+=nsub) {
    cel = elvec[elsons[j]];
    if(!cel || cel->IsInterface() || !cel->Dimension()) {
		nsub = 1;
		continue;
    }
    el = (TPZInterpolatedElement *)cel;
    TPZGeoEl *gel = el->Reference()->Father();
    nsub = el->Reference()->NSubElements();
    /**Se nao existe father ou ele eh o mais grosso nao faz coarse*/
    if(!gel || gel->Level()<MinLevel) continue;     //Jorge 17/08
    indexes.Resize(nsub,-1);
    int nsides = el->Reference()->NSides();
    if(el->PreferredSideOrder(nsides-1)) order = 1;   //Se algum elemento tem ordem zero order = 0
    else order = 1;
    minorder = 10;
    /**Verificando se os elementos sucessivos tem o mesmo elemento pai
       e averiguando se existe elemento com ordem de interpolacao zero*/
    for(i=1;i<nsub;i++) {
      int orderpartial;
      cel = elvec[elsons[j+i]];
      if(!cel || cel->IsInterface() ||
             gel != ((TPZInterpolatedElement* )cel)->Reference()->Father())
        goto MASK1;
      orderpartial = ((TPZInterpolatedElement *)cel)->PreferredSideOrder(nsides-1);
      if(!orderpartial) order = 0;
      if(minorder>orderpartial) minorder = orderpartial;
    }
    for(i=0;i<nsub;i++) indexes[i] = elsons[j+i];

    /** Quando a ordem de interpolacao eh zero tem tratamento diferenciado */
    if(!order) {
      TPZVec<REAL> values(nvar,0.);
      /** Tomando a media entre todos os sub-elementos, como valores do elemento grosso */
      for(i=0;i<nsub;i++) {
        intel = (TPZInterpolatedElement *)elvec[indexes[i]];
        for(p=0;p<nvar;p++)
				  values[p] += intel->MeanSolution(p);
      }
      for(i=0;i<nvar;i++) values[i] /= nsub;
      cmesh.SetDefaultOrder(0);
      cmesh.Coarsen(indexes,index);
      cmesh.ExpandSolution();
      intel = (TPZInterpolatedElement *)elvec[index];
      int seqnum = intel->Connect(nsides).SequenceNumber();
      /** Atribuindo as medias ao elemento grosso */
      for(p=0;p<nvar;p++)
        cmesh.Block()(seqnum,0,p,0) = values[p];
    }
    else {
    /** Desconsidera-se o connect sobre o element para ser continuo */
      int ncorners = el->NCornerConnects();
      TPZVec<REAL> values(ncorners*nvar,0.);
      int seqnum, nstart;
      /**Tomando os valores dos connects esquinas como valores do elem. grosso.*/
      for(i=0;i<ncorners;i++) {
        nstart = 0;
        intel = (TPZInterpolatedElement *)gel->SubElement(i)->Reference();
        /**Se connect continuous ele preenche valores do mesmo connect*/
//        if(intel->IsConnectContinuous(i)) continue;
        /*caso contrario resgatamos os valores discontinuous*/
        for(p=0;p<i;p++)   //anteriores connects esquina do subelement
 //         if(!intel->IsConnectContinuous(p))
            nstart += nvar;   // ajust the index into the discontinuous connect
        seqnum = intel->Connect(nsides).SequenceNumber();
        for(p=0;p<nvar;p++) values[i*nvar+p] = cmesh.Block()(seqnum,0,nstart+p,0);
        /**Fazemos todos os connects do subelemento discontinuous exeto o atual connect*/
        for(p=0;p<nsides;p++) {
          if(p==i) continue;
  //        intel->MakeConnectDiscontinuous(p);
        }   // Para que os provaveis vizinhos do mesmo nivel do subelemento
      }     // preservem os seus valores. Estes connects ficaram dependentes
      cmesh.SetDefaultOrder(minorder);
      cmesh.Coarsen(indexes,index);
//      ((TPZInterpolatedElement *) elvec[index])->SetDiffusionCoefficient(0.);
      cmesh.ExpandSolution();
      if(index==-1) PZError << "Coarsing. Coarse element is not make.\n";
      intel = (TPZInterpolatedElement *)elvec[index];
      nstart = 0;
      seqnum = intel->Connect(nsides).SequenceNumber();
      for(i=0;i<ncorners;i++) {
    //    if(intel->IsConnectContinuous(i)) continue;
        for(p=0;p<nvar;p++)
          cmesh.Block()(seqnum,0,nstart+p,0) = values[i*nvar+p];
        nstart += nvar;
      }
      for(;i<nsides;i++) {
     //   if(intel->IsConnectContinuous(i)) {
       //   seqnum = intel->Connect(i).SequenceNumber();
         // int size = cmesh.Block().Size(seqnum);
    //      for(p=0;p<size;p++) cmesh.Block()(seqnum,0,p,0) = 0.;
      //  }
        //else {
          seqnum = intel->Connect(nsides).SequenceNumber();
          //int size = nvar*(intel->NConnectShapeF(i));
     //     for(p=0;p<size;p++) cmesh.Block()(seqnum,0,nstart+p,0) = 0.;
     //     nstart += size;
    //    }
      }
    }
    elsons[counter++] = index;
    MASK1:
    intel = 0; el = 0;
  }
  elsons.Resize(0);
  cmesh.SetDefaultOrder(orderold);
//  cmesh.CleanInterfaces();
  cmesh.InitializeBlock();
}

/** Before to refine a group of geo elements, decrement its interpolation order*/
int TAdaptive::DecrementOrder(TPZCompMesh &cmesh,TPZManVector<int> &elfathersvec) {
  int i, j, nelem = elfathersvec.NElements();
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  int side, order, minorder=0;
  TPZInterpolatedElement *el, *newel;
  TPZBlock<STATE> &block = cmesh.Block();
  int oldgorder = cmesh.GetDefaultOrder();
  int64_t index;
	if(!oldgorder || !ReAssembling) {
	  elfathersvec.Resize(0);
		return 0;
	}
  for(i=0;i<nelem;i++) {
    el = (TPZInterpolatedElement *)elvec[elfathersvec[i]];
    side = el->Reference()->NSides();
    order = el->PreferredSideOrder(side-1);
    if(el->Reference()->Level()< MaxLevelRefinements-1 || !order) continue;
    /* Aqui a ordem tem que ser escolhida: se order>1 entao 1 se order =1 entao 0 */
/*    if(order>1) {
    /** Falta recuperar os valores dos connects respectivos e preencher apos PRefine
//      el->PRefine(1);
//      minorder = 1;
      el->PRefine(order-1);
      minorder = order-1;                            //TESTE!!!
    }
    else if(order==1) {
*/
    if(order==1) {
      for(j=0;j<side;j++) {
        /* Sera que esta-se limpando os blocos com os valores das variaveis??? */
     //   if(el->IsConnectContinuous(j)) {
       //   PZError << "DecrementOrder. Element has connect continuous.\n";
       //   break;
     //   }
      }
      if(j<side) continue;
      minorder = 0;
      int nvar = el->Material()->NStateVariables();
      TPZVec<REAL> values(nvar);
      for(j=0;j<nvar;j++) values[j] = el->MeanSolution(j);
      el->PRefine(0);
      int seqnum = el->Connect(side).SequenceNumber();
      if(nvar!=block.Size(seqnum))
        PZError << "DecrementOrder. Dimension of block internal is uncompatible.\n";
      for(j=0;j<nvar;j++) block(seqnum,0,j,0) = values[j];
    }
    else {
      cmesh.SetDefaultOrder(order-1);
      TPZGeoEl *gel = el->Reference();
      gel->ResetReference();
	  cmesh.CreateCompEl(gel, index);
      newel = (TPZInterpolatedElement *)(cmesh.ElementVec()[index]);
      newel->CheckConstraintConsistency();
      cmesh.ExpandSolution();
      newel->InterpolateSolution(*el);
      index = el->Index();
      delete el;
      cmesh.ElementVec()[index] = 0;
      cmesh.ElementVec().SetFree(index);
      gel->SetReference(newel);
    }
  }
  cmesh.SetDefaultOrder(oldgorder - 1);
  return minorder;
}

/** After to coarsing a group of geo elements, increment its interpolation order*/
int TAdaptive::IncrementOrder(TPZCompMesh &cmesh) {
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  int i, nelem = elsons.NElements();
  int side, order, maxorder = 0;
  TPZInterpolatedElement *el;
  for(i=0;i<nelem;i++) {
    el = (TPZInterpolatedElement *)elvec[elsons[i]];
    side = el->Reference()->NSides()-1;
    order = el->PreferredSideOrder(side)+1;
    if(order <= el->PreferredSideOrder(side)) {
      el->PRefine(order);
      if(maxorder<order) maxorder = order;
    }
  }
  return maxorder;
}

void TAdaptive::CleanVectors() {
  elfathers.Resize(0);
	elsons.Resize(0);
}

/** Here is called a wavelet decomposition to identify singularity regions*/
/** 1. Se level (da malha) = 1 entao refina toda a malha, (nao eh possivel aplicar
       wavelets em um nivel soh). Mas se onlyfluxadaptive!=0 retorna zero.
    2. se level > 1 faz uma decomposicao wavelet obtendo um vetor de elementos para
       refinar(elfathers) e um vetor de elementos para formar outros mais grossos(elsons).
    3. Estos vetores podem ser utilizados apenas para modificar o fluxo numerico a utilizar
       quando onlyfluxadaptive!=0. Aplica-se um fluxo de alta ordem nas interfaces dos
       elementos em elsons e um fluxo de baixa ordem nas interfaces dos elementos em elfathers.
    3. Quando segurityzone!=0 alem dos elementos detetados pelos wavelets para ser refinados
       sao incrementados no vetor elfathers os vizinhos deles.*/
int TAdaptive::IsNecessaryToRefine(TPZCompMesh &cmesh,int level,int onlyfluxadaptive) {

	int i,refine;
  int nvar;
	std::map<int, TPZMaterial * >::const_iterator mit;
	TPZMaterial * mat = 0;
	for (mit = cmesh.MaterialVec().begin(); mit != cmesh.MaterialVec().end(); mit++) {
		if (mit->second->Id() < 0) continue;
		mat = mit->second;
		break;
	}

  if(!mat) {
    PZError << "IsNecessaryToRefine. Computational mesh has no materials or bad level.\n";
    return 0;
  }

	nvar = mat->NStateVariables();
  /**First decomposition if mesh was refined. It is level no zero.*/
  if(level>1) {
    if(!nvar) {
      PZError << "IsNecessaryToRefine. Comp mesh has no material associated.\n";
      return nvar;
    }
    /**Vector to store the geometrical elements with large wavelet coefficients*/
    TPZManVector<int> idgeoel(0);
    /** It is not necessary the descomposition over all variables */
//    for(i=0;i<nvar;i++) {
//      WaveletDecomposition(cmesh,i,level,dim,segurityzone);
		if(nvar < 3) nvar = 0;
    WaveletDecomposition(cmesh,idgeoel,nvar,level);
    refine = idgeoel.NElements();
    if(refine) {
      int nelfathers = elfathers.NElements();
      refine += nelfathers;
      elfathers.Resize(refine);
      for(int j=nelfathers;j<refine;j++)
        elfathers[j] = idgeoel[j-nelfathers];
      idgeoel.Resize(0);
    }
    if(onlyfluxadaptive) {
      AdaptNumericalFlux(cmesh,((TConservationLaw *)mat)->fFluxTypeHigher,elsons);
      AdaptNumericalFlux(cmesh,((TConservationLaw *)mat)->fFluxTypeLower,elfathers);
      elfathers.Resize(0);
      elsons.Resize(0);
    }
  }
  else if(!onlyfluxadaptive) {
    Refinements(cmesh,1);
    IsNecessaryToRefine(cmesh,level+1,0);
  }
  return (elfathers.NElements() + elsons.NElements());
}

/** To adapt numerical flux for high or lower resolution.
    elfathers is a vector of computational elements over continuity region */
void TAdaptive::AdaptNumericalFlux(TPZCompMesh &cmesh,int fluxtype,TPZManVector<int> &vecelem) {
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  int nelemesh = cmesh.NElements();
  int i, j, k, nelem = vecelem.NElements();
	if(!nelem) return;
  TPZVec<int> indexel(nelemesh,0);
  int nside, nsideneigh, index, inter_face;
  TPZInterpolatedElement *el, *neigh;
  for(i=0;i<nelem;i++) {
    index = vecelem[i];
    if(index==-1) continue;
    el = (TPZInterpolatedElement *)elvec[index];
    if(!el) continue;
    nside = el->Reference()->NSides();
    for(j=0;j<nside;j++) {
		inter_face = el->Interface(j);
      if(inter_face ==-1) continue;
      if(inter_face ==-2) {
        /**Entao este elemento tem elementos de nivel maior conectados ao longo do lado side*/
        TPZStack<TPZCompElSide> elemvec;
        TPZCompElSide thisside(el,j);
        TPZGeoElSide gelside = thisside.Reference();
        gelside = gelside.Neighbour();
		TPZStack<TPZGeoElSide> gsmallvec;
		if (!gelside.Exists()) PZError << "TCompEl1dGD::MakeConnectContinuous.. Not exist neighboard.";
		gelside.GetSubElements2(gsmallvec);
		for (i = 0; i < gsmallvec.NElements(); i++)
			elemvec.Push(gsmallvec.Pop().Reference());

        int nneighs = elemvec.NElements();
        for(k=0;k<nneighs;k++) {
          neigh = (TPZInterpolatedElement *)elemvec[k].Element();
          nsideneigh = neigh->Reference()->NSides();
          for(int p=0;p<nsideneigh;p++) {
			  inter_face = neigh->Interface(p);
            if(inter_face == -1 || inter_face == -2) continue;
            indexel[inter_face] = 1;
          }
        }
      }
      else
        indexel[inter_face] = 1;
    }
  }
  for(i=0;i<nelemesh;i++) {
    if(indexel[i])
      ((TInterfaceElement *)elvec[i])->SetFluxType(fluxtype);
  }
}

/**Baseado no parametro WaveletIndex aplica sobre a malha(pelo menos uma vez refinada)
   uma decomposicao wavelet utilizando as bases de Haar(WaveletIndex==0) ou as bases de Schauder.
   No caso das bases de Schauder, se o elemento comp. tem ordem de interpolacao zero ele
   eh separado para aplicar nele as bases de Haar.
   1. Se algum elemento em elsons esta em elfathers (por efeito da segurityzone) ele eh
      apagado em elsons.
   2. Tiram-se de elfathers e elsons todo indice -1 e os elementos que nao podem
      ser mais refinados (em elfathers) ou que estao no nivel mais grosso(em elsons).*/
void TAdaptive::WaveletDecomposition(TPZCompMesh &cmesh,TPZManVector<int> &elfathersvec,
        int var,int level) {
  if(level>1) {
	  if(!fWavelet) {
		  std::cout << "TAdaptive::WaveletDecomposition. Wavelet pointer is Null.\n";
			return;
		}
    if(fWavelet->Decomposition(cmesh,var,elfathersvec,elsons,MinLevel))
      fWaveletSec->DecompositionOnlyOrderZero(cmesh,var,elfathersvec,elsons,MinLevel);

    /**Tirando todo elemento de elsons que tambem esta no elfathers e tirando -1
       e os elementos que ja estao no nivel mais grosso*/
    if(segurityzone) SegurityZone(cmesh,elfathersvec);
//    CheckVectorsToAdaptivity(cmesh,elfathersvec,elsons);
  }
}

void TAdaptive::AllIncrementInterpolation(TPZCompMesh &cmesh) {
	int nelem = cmesh.NElements(), nsons = 0, i;
	elsons.Resize(nelem);
  /**Variaveis para calculo dos coeficientes wavelets*/
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  TPZCompEl *el;

  for(i=0;i<nelem;i++) {
    el = elvec[i];
    if(!el || el->IsInterface() || el->Dimension() != 1) continue;
    elsons[nsons++] = el->Index();
  }
  elsons.Resize(nsons);
	i = IncrementOrder(cmesh);
	std::cout << "Nova ordem de interpolacao: " << i << std::endl;
}

/** Making adaptive refinement, based on high wavelet coefficients elfathers contains
    the element ids to refine, elsons contains the element ids to coarsing */
void TAdaptive::WaveletHPAdaptive(TPZCompMesh &cmesh) {

  /** Refining whether the wavelet coefficients are large (respec. elements) */
  if(!elfathers.NElements() && !elsons.NElements()) return;
	int i, n;
	std::map<int, TPZMaterial * >::const_iterator mit;
	TPZMaterial * mat = 0;
	for (mit = cmesh.MaterialVec().begin(); mit != cmesh.MaterialVec().end(); mit++) {
		if (mit->second->Id() < 0) continue;
		mat = mit->second;
	    break;
  }
  if(!mat) PZError << "WaveletHPAdaptive. Not found materials.\n";

  TPZManVector<int> smalls(0);
  CheckVectorsToAdaptivity(cmesh,smalls);

  int nrefinements = 1;
  ((TConservationLaw *)mat)->SetFluxTypeLower();

  TAdaptive::Refinements(cmesh,nrefinements);

  if(p_decrement) DecrementOrder(cmesh,smalls);
  smalls.Resize(0);
  /** Now, to coarse the element with small wavelet coefficients*/
  ((TConservationLaw *)mat)->SetFluxTypeHigher();
  if(p_increment) IncrementOrder(cmesh);

  if(!fSteadyState) Coarsing(cmesh);

  /** Clean the vectors of the coarsing or refined elements*/
  elfathers.Resize(0);
  elsons.Resize(0);
//  cmesh.CleanInterfaces();
  cmesh.InitializeBlock();
}

void TAdaptive::CheckVectorsToAdaptivity(TPZCompMesh &cmesh,TPZManVector<int> &smalls) {   //smalls has level MAXLEVELREFINEMENTS
  int indexfather = elfathers.NElements();
  if(!indexfather && !elsons.NElements()) return;
  int i, j, k, index, nsub, indexintel, indexson = elsons.NElements();
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();

  int dimmat;
  std::map<int, TPZMaterial * >::const_iterator mit;
  TPZMaterial * mat = 0;
  for (mit = cmesh.MaterialVec().begin(); mit != cmesh.MaterialVec().end(); mit++) {
	  if (mit->second->Id() < 0) continue;
	  mat = mit->second;
    dimmat = mat->Dimension();
    break;
  }
  if(dimmat==1) nsub = 2;
  else if(dimmat==2) nsub = 4;

  /**1. Tirando os -1 de elfathers e os elementos que ja nao podem ser refinados
     2. Tirando os elementos repetidos em elfathers*/
  TPZManVector<int> vetor(elfathers);
  TPZCompEl *cel;
  index = 0;
	int lenghtfather = indexfather-1;
  for(i=0;i<indexfather;i++) {
    indexintel = vetor[i];
    if(indexintel==-1) continue;
    cel = elvec[indexintel];
    if(!cel || cel->IsInterface()) continue;
    if(((TPZInterpolatedElement *)cel)->Reference()->Level() > MaxLevelRefinements-1)
      continue;

		if(index>=lenghtfather) {
		  lenghtfather += 100;
			elfathers.Resize(lenghtfather+1);
		}
    elfathers[index++] = indexintel;
    /**Limpando os indeces repetidos a indexintel*/
    for(j=i+1;j<indexfather;j++)
      if(indexintel==vetor[j]) vetor[j]=-1;
  }
  elfathers.Resize(index);
  indexfather = index;   // Para em elsons tirar os que estao em elfathers

  /**Procura por elementos com coef_wave grande que nao podem mais ser refinados
     entao serao interpolados com ordem atual menos 1, emquanto ReAssembling */
  if(ReAssembling && indexfather) {
    TPZInterpolatedElement *el;
    smalls.Resize(indexfather,-1);
    int indexsmalls = 0;
    for(i=0;i<indexfather;i++) {
      el = (TPZInterpolatedElement *)elvec[elfathers[i]];
			k = el->PreferredSideOrder(el->Reference()->NSides() - 1);
      if(el->Reference()->Level() > MaxLevelRefinements-1 && k)
        smalls[indexsmalls++] = el->Index();
    }
    smalls.Resize(indexsmalls);
  }

  TPZManVector<int> vectorson(elsons);
  TPZInterpolatedElement *intel;
  int lenghtson = indexson;
  index = 0;
  /**Checando o vetor de filhos. 1. Tira os indices -1.
     2. Verifica se os elementos consecutivos sao filhos do mesmo pai (depende de nsub)
        se nao existe pai esta no nivel mais grosso, nao pode ser filho.
     3. Tira os elementos em elsons que estejam em elfathers.*/
  int nsubcurrent = nsub;
  for(i=0;i<indexson;i+=nsubcurrent) {
    indexintel = vectorson[i];
    if(indexintel==-1) continue;
    cel = elvec[indexintel];
    if(!cel || cel->IsInterface()) continue;
    nsubcurrent = cel->Reference()->NSubElements();
    if(cel->Dimension()<dimmat) {
//      PZError << "CheckVectors. Error in vector of sons.\n";
      continue;
    }
    intel = (TPZInterpolatedElement *)cel;
    /**Verificando se os elementos sucessivos tem o mesmo elemento pai*/
    TPZGeoEl *gel = intel->Reference()->Father();
    if(!gel || gel->Level() < MinLevel) continue;   // Esta no nivel mais grosso Jorge 17/08
    for(j=1;j<nsubcurrent;j++) {
      cel = elvec[elsons[j+i]];
      if(!cel || cel->IsInterface() ||
                 gel != ((TPZInterpolatedElement* )cel)->Reference()->Father()) {
        j = 9;
        continue;
      }
    }
    if(j==10) continue;
    /**Tira de elsons aqueles que estao no elfathers*/
    for(k=0;k<nsubcurrent;) {
      indexintel = vectorson[i+k];
      for(j=0;j<indexfather;j++) {
        if(indexintel == elfathers[j]) {      // ??? Para o que ???
          k = nsubcurrent+1;
          break;
        }
      }
      if(k>nsubcurrent) break;
      k++;
//      if((i+k)<(indexson) && k<(nsub+1)) indexintel = vectorson[i+k];   // 16/04/2000
    }
    /**Preenchendo em elsons apenas aqueles que resultaram adequados*/
    if(k==nsub) {
      for(j=0;j<nsub;j++)
        elsons[index++] = vectorson[i+j];
      if(dimmat==2) {
      /**Se algum subelemento tem vizinho de dimensao menor e igual nivel, este
      tambem deve ser engrossado*/
        int ncorners = intel->NCornerConnects();
        for(k=0;k<ncorners;k++) {
          TPZGeoEl *neighson, *neigh = (gel->Neighbour(k+ncorners)).Element();
          if(!neigh || neigh->Dimension()>=dimmat) continue;
          int jj, nsub1 = neigh->NSubElements();
          TPZVec<TPZCompEl *> cels(nsub1,0);
          for(jj=0;jj<nsub1;jj++) {
            neighson = neigh->SubElement(jj);
            if(!neighson) break;
            cels[jj] = neighson->Reference();
            if(!cels[jj]) break;
          }
          if(jj!=nsub1) break;
          for(jj=0;jj<nsub1;jj++) {
            if(lenghtson<index+nsub1) {
              lenghtson+=100;
              elsons.Resize(lenghtson,-1);
            }
            elsons[index++] = cels[jj]->Index();
          }
        }
      }
/**          cel = elvec[elsons[index-4+k]];
          TPZCompElSide celside(cel,k+ncorners);
          TPZStack<TPZCompElSide> celneighs(0);
          celside.EqualLevelElementList(celneighs,1,0);
          if(cel->Dimension()<dimmat) continue;
          if(celneighs.NElements()==1 && celneighs[0].Element()->Dimension()==1) {
            TPZGeoElSide neighside = celneighs[0].Element()->Reference()->Neighbour(1);
            if(neighside.Reference().Exists()) {
              indexintel = celneighs[0].Element()->Index();
        /**Tirando elementos fronteiras de elfathers caso tiver
            for(j=0;j<indexfather;j++)
              if(indexintel == elfathers[j])
                { elfathers[j] = -1; break; }
            elsons[index++] = indexintel;
            indexintel = neighside.Reference().Element()->Index();
            for(j=0;j<indexfather;j++)
              if(indexintel == elfathers[j])
                { elfathers[j] = -1; break; }
            elsons[index++] = indexintel;
          }
        }
      }*/
    }
  }
  elsons.Resize(index);
}

int TAdaptive::SegurityZone(TPZCompMesh &cmesh,TPZManVector<int> &elfathersvec) {
  int i, ii, p, nsides, nelem = elfathersvec.NElements();
  int index = nelem, nfathers = nelem;
  TPZCompEl *cel;
  TPZInterpolatedElement *intel;
  /**Procura elementos vizinhos de maior ou igual nivel para refinar-os*/
  for(i=0;i<nfathers;i++) {
    cel = cmesh.ElementVec()[elfathersvec[i]];
    if(!cel || cel->IsInterface()) continue;
    intel = (TPZInterpolatedElement *)cel;
    nsides = intel->Reference()->NSides();
    TPZStack<TPZCompElSide> elstack;
    for(ii=0;ii<nsides;ii++) {
      TPZCompElSide intelside(intel,ii);
      TPZCompElSide large;
      intelside.EqualLevelElementList(elstack,1,0);
      large = intelside.LowerLevelElementList(1);
      if(large.Exists()) elstack.Push(large);
    }
    p = elstack.NElements();
    if(nelem<index+p) {
      nelem += 100;
      elfathersvec.Resize(nelem);
    }
    for(ii=0;ii<p;ii++)
      elfathersvec[index++] = elstack[ii].Element()->Index();
  }
  elfathersvec.Resize(index);
  return index-nfathers;
}

/*int IsThreadholder(REAL value) {
  return (( value > -.001 && value < .001 ) ? 0 : 1);
}
/** Making adaptive refinement, based on high wavelet coefficients elfathers contains
    the element ids to refine, elsons contains the element ids to coarsing 
void WaveletHPAdaptive(TPZCompMesh &cmesh,TPZManVector<int>
&elfathers,TPZManVector<int> &elsons,int decrement,int increment) {

  /** Refining whether the wavelet coefficients are large (respec. elements) 
  if(!elfathers.NElements()) return;

  /** nelem!=0 means the mesh will be refined 
  int nrefine = 0;
  while(nelem && maxitera) {
    int ncycles = 1;  // O numero de ciclos depende do tamanho dos coeficientes wavelets
    TPZManVector<int> sons(0);
    DecrementOrder(cmesh,elfathers);
    Refinements(cmesh,elfathers,ncycles);
    nrefine++;
    nelem = IsNecessaryToRefine(cmesh,elfathers,sons,out,1);
    maxitera--;
  } 
  if(decrement) DecrementOrder(cmesh,elfathers);
  Refinements(cmesh,elfathers,3);
  cmesh.InitializeBlock();
  /** Now, to coarse the element with small wavelet coefficients
  if(increment) IncrementOrder(cmesh,elsons);
  Coarsing(cmesh,elsons);
  cmesh.InitializeBlock();
  /** Clean the vectors of the coarsing or refined elements
  elfathers.Resize(0);
  elsons.Resize(0);
}
*/

/*    brotherel = (TPZInterpolatedElement *)gel->SubElement(1)->Reference();
    gel = gel->Father();
    if(!gel) PZError << "Haar1d. Bad parameter level.\n";
    id = gel->Id();
    if(maskgel[id]!=-1) {
      idgeoel[i] = gel->Id();
      maskgel[id] = id;
    }
    meaneven = intel->MeanSolution(var);
    meanodd = brotherel->MeanSolution(var);
    coef[i] = (meaneven + meanodd)/M_SQRT2;
    djk = (meaneven - meanodd)/M_SQRT2;
    if(IsThreadholder(djk)) djk = 0.;
    coef[iniwav++] = djk;
    }
    Pause(iniwav);
    k /= 2;
  }
  iniwav = k/2;
  for(i=0;i<k;i++) {
    gel = gelvec[idgeoel[i]];
    intel = (TPZInterpolatedElement *)gel->SubElement(0)->Reference();
    brotherel = (TPZInterpolatedElement *)gel->SubElement(1)->Reference();
    meaneven = intel->MeanSolution(var);
    meanodd = brotherel->MeanSolution(var);
    coef[i] = (meaneven + meanodd)/M_SQRT2;
    djk = (meaneven - meanodd)/M_SQRT2;
    if(IsThreadholder(djk)) djk = 0.;
    coef[iniwav++] = djk;
  }
}

/** To adapt numerical flux for high or lower resolution.
    elfathers is a vector of computational elements over continuity region 
void AdaptNumericalFlux(TPZCompMesh &cmesh,int fluxtype,TPZManVector<int> &elfathers) {
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  int nelemesh = cmesh.NElements();
  int i, j, k, nelem = elfathers.NElements();
  TPZVec<int> indexel(nelemesh,0);
  int nside, nsideneigh, index, interface;
  TPZInterpolatedElement *el, *neigh;
  for(i=0;i<nelem;i++) {
    index = elfathers[i];
    if(index==-1) continue;
    el = (TPZInterpolatedElement *)elvec[index];
    if(!el) continue;
    nside = el->NSides();
    for(j=0;j<nside;j++) {
      interface = el->Interface(j);
      if(interface==-1) continue;
      if(interface==-2) {
        /**Entao este elemento tem elementos de nivel maior conectados ao longo do lado side
        TPZStack<TPZCompElSide> elemvec(0);
        TPZCompElSide thisside(el,j);
        TPZGeoElSide gelside = thisside.Reference();
        int nneighs = 0;
        if(gelside.Dimension()) {
          gelside = gelside.Neighbour();
          TPZVec<TPZGeoElSide> subelements(10);
          TPZVec<TPZGeoElSide> subs;
          gelside.GetSubElements(subs);
          int r, nallsubs = 0, nsubs = subs.NElements();
          while(nsubs) {
            for(r=0;r<nsubs-1;r++) {
              if(subs[r].Element()->Reference()) elemvec.Push(subs[r].Reference());
              else subelements[nallsubs++] = subs[r];
            }
            if(subs[r].Element()->Reference()) {
              elemvec.Push(subs[r].Reference());
              nsubs = 0;
            }
            else {
              subs[r].GetSubElements(subs);
              nsubs = subs.NElements();
            }
          }
          for(r=0;r<nallsubs;r++) {
            subelements[r].GetSubElements(subs);
            nsubs = subs.NElements();
            if(!nsubs) PZError << "TPZInterpolatedElement::Divide. It has not sub-elements.\n";
            for(int p=0;p<nsubs;p++) {
              if(subs[p].Element()->Reference()) elemvec.Push(subs[p].Reference());
              else subelements[nallsubs++] = subs[p];
            }
          }
//          thisside.HigherLevelElementList(elvec,1,0);
          nneighs = elemvec.NElements();
        }
        else {
          gelside = gelside.Neighbour();
          TPZGeoEl *gel = gelside.Element();
          while(!nneighs) {
            gelside = gel->SideSubElement(gelside.Side(),0);
            gel = gelside.Element();
            if(!gel) {
              PZError << "PZIntel::Divide. Error found interfaces.\n";
              break;
            }
            if(gel->Reference()) {
              elemvec.Push(gelside.Reference());
              nneighs = elemvec.NElements();
            }
          }
        }
        for(k=0;k<nneighs;k++) {
          neigh = (TPZInterpolatedElement *)elemvec[k].Element();
          nsideneigh = neigh->NSides();
          for(int p=0;p<nsideneigh;p++) {
            interface = neigh->Interface(p);
            if(interface == -1 || interface == -2) continue;
            indexel[interface] = 1;
          }
        }
      }
      else
        indexel[interface] = 1;
    }
  }
  for(i=0;i<nelemesh;i++) {
    if(indexel[i])
      ((TInterfaceElement *)elvec[i])->SetFluxType(fluxtype);
  }
}
*/

/** 
1. Problema com MeanSolutionFromSubElements na decomposicao wavelet Haar
2. Quando aplica-se coarsing para elementos com ordem zero ele joga
valor zero ao elemento grosso.
3. Mudar HaarOnlyZeros Schauder e todos os bi-dimensional
*/

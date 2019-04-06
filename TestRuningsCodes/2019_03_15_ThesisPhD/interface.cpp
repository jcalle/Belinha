#include "interface.h"
#include "pzcmesh.h"
#include "pzerror.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzelgpoint.h"
#include "pzelgt2d.h"
#include "pzelg1d.h"
#include "pzelgc3d.h"
#include "pzelgpi3d.h"
#include "pzelgpr3d.h"
#include "pzelgq2d.h"
#include "pzelgt3d.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "conslaw.h"


/** Constructor to interface element. It element compute the contribution on the
calstiff matrix and load vector of the integral over the boundary elements. It permits
to compute conservative contribution over the common boundary between neighboard
elements. This interface element is constructed by elleft element (then elleft is discontinuous.*/
TInterfaceElement::TInterfaceElement(TPZCompMesh &mesh,TPZCompElSide &elleft,
        TPZCompElSide &elright,int fluxtype, int64_t &index) : TPZInterfaceElement(mesh,0,index) {
  fContribution = 1;
  /**Depending on the dimension of the elements it attributed the left element*/
  if(elleft.Element()->Dimension()<elright.Element()->Dimension()) {
    fElVec[0] = elright;
    fElVec[1] = elleft;
  }
  else {
    fElVec[0] = elleft; fElVec[1] = elright;
  }
	  
  /**Material is passed from left element*/
  fMaterial = fElVec[0].Element()->Material();
  /**Filling the dimension from first computational side element */
  fDimension = fElVec[0].Reference().Dimension();
  TPZVec<int> order;
  int ord = 100;
  fElVec[0].Element()->GetIntegrationRule().GetOrder(order);
  for (int ii; ii < order.NElements(); ii++)
	  ord = (ord < order[ii]) ? ord : order[ii];
  order.Fill(0);
  fElVec[1].Element()->GetIntegrationRule().GetOrder(order);
  for (int ii; ii < order.NElements(); ii++)
	  ord = (ord < order[ii]) ? ord : order[ii];
  SetIntegrationRule(ord);
 
  if(!fIntegrationRule) fContribution = 0;
  TPZMaterial *mat0 = ((TPZInterpolatedElement *)fElVec[0].Element())->Material();
  TPZMaterial *mat1 = ((TPZInterpolatedElement *)fElVec[1].Element())->Material();
	fFluxType = fluxtype;
  if(elleft.Element()->Material()==NULL || elright.Element()->Material()==NULL){
    PZError << "TInterfaceElement. No material for this element\n";
    return;
  }
  /**Verifying the materials on the neighboards. It is possible to have
     one boundary condition and the other no boundary condition*/
  if(mat0 != mat1) {
    if(mat0->Id()<0 && ((TPZBndCond *)mat0)->Material()==mat1) return;
    if(mat1->Id()<0 && ((TPZBndCond *)mat1)->Material()==mat0) return;
    PZError << "TInterfaceElement is created between two comp. elements with different materials.\n";
  }
}

/**When interface element is deleted it put -1 in the respective fInterface into
   two related computational elements*/
TInterfaceElement::~TInterfaceElement() {
}

/**Return the type of interface. If dimension is zero then it is point_interface,
if dim = 1 it is linear_interface and if dim = 2 it is surface_interface*/
MElementType TInterfaceElement::Type() {
  switch(fDimension) {
    case 0: return EInterfacePoint;
    case 1: return EInterfaceLinear;
    case 2: return EInterfaceSurface;
    default:
      PZError << "TInterfaceElement::Type. Type of interface is undefined.\n";
  }
  return EInterface;
}

/**Return number of side connects of the left element related plus number of
   side connects of the right element related*/
int TInterfaceElement::NConnects() {
  int ncon0 = ((TPZInterpolatedElement *)fElVec[0].Element())->NSideConnects(fElVec[0].Side());
  return ncon0+(((TPZInterpolatedElement *)fElVec[1].Element())->NSideConnects(fElVec[1].Side()));
}
/**Return number of connects of the one element related with this interface*/
int TInterfaceElement::NConnects(int elnumber) {
#ifndef NOTDEBUG
  if(elnumber && elnumber!=1)
    PZError << "TInterfaceElement::NConnects. Bad parameter elnumber.\n";
#endif NOTDEBUG
  int nn = fElVec[elnumber].Side();
  TPZInterpolatedElement *el = (TPZInterpolatedElement *)fElVec[elnumber].Element();
  return el->NSideConnects(nn);
}

/**The connects to the interface element are enumerated same as the side connects
on the left element and the interior connect of the left element and side connects
of the right element and for last the internal connect of right element, if it exists*/
TPZConnect &TInterfaceElement::Connect(int i) {
  int ncon0 = NConnects(0), n = NConnects(1);
#ifndef NOTDEBUG
  if(i>(ncon0+n)-1 || i<0) {
    PZError << "TInterfaceElement::Connect. Bad parameter i.\n";
  }
#endif NOTDEBUG
  if(i<ncon0) n = 0;
  else { n = 1; i -= ncon0; }
  TPZInterpolatedElement *cel = (TPZInterpolatedElement *)fElVec[n].Element();
  return (cel->SideConnect(i,fElVec[n].Side()));
}
/**Return icon-th connect in elnumber element of the fElVec*/
TPZConnect &TInterfaceElement::Connect(int i,int elnumber) {
#ifndef NOTDEBUG
  if(elnumber && elnumber!=1)
    PZError << "TInterfaceElement::Connect. Bad parameter elnumber.\n";
#endif NOTDEBUG
  return (((TPZInterpolatedElement *)fElVec[elnumber].Element())->SideConnect(i,fElVec[elnumber].Side()));
}

int64_t TInterfaceElement::ConnectIndex(int i) {
  int ncon0 = NConnects(0), ncon1 = NConnects(1);
#ifndef NOTDEBUG
  if(i>ncon0+ncon1-1 || i<0) {
    PZError << "TInterfaceElement::ConnectIndex. Bad parameter i.\n";
  }
#endif NOTDEBUG
  if(i<ncon0) return ((TPZInterpolatedElement *)fElVec[0].Element())->SideConnectIndex(i,fElVec[0].Side());
  return ((TPZInterpolatedElement *)fElVec[1].Element())->SideConnectIndex(i-ncon0,fElVec[1].Side());
}

/**Prints the information of the interface element*/
void TInterfaceElement::Print(std::ostream &out) {
  out << std::endl << "   Interface Element : ";
  switch(fDimension) {
    case 0:
      out << "Point Interface"; break;
    case 1:
      out << "Linear Interface"; break;
    case 2:
      out << "Surface Interface"; break;
    default:
      out << "Undefined Interface"; break;
  }
  out <<std::endl;
  if(!fContribution) out << "Inactive interface" <<std::endl;
  else out << "Active interface" <<std::endl;
  return;
	/**All the connects related with this interface*/
  int i, ncon = NConnects();
  out <<std::endl << "Connects related with the interface, indexes :  " <<std::endl;
  for(i=0;i<ncon;i++) out << "\t" << ConnectIndex(i);
  out <<std::endl << "Material Id = " << fMaterial->Id() <<std::endl;
}

void TInterfaceElement::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef) {
  if(!fContribution) {
    ek.fBlock = NULL;
    ek.fMat.Zero();
    ef.fBlock = NULL;
    ef.fMat.Zero();
    return;
  }
  if(fMaterial->Id()<0) {
		PZError << "Interface with boundary condition.\n";
		ek.fMat.Zero();
		ef.fMat.Zero();
    return;
  }

  TPZInterpolatedElement *cel0 = (TPZInterpolatedElement *)fElVec[0].Element();
  TPZInterpolatedElement *cel1 = (TPZInterpolatedElement *)fElVec[1].Element();
/**First compute the dimensions to matrices*/
  TPZBlock<STATE> &block = Mesh()->Block();
  int ndof = fMaterial->NStateVariables();
  int ncon0 = NConnects(0);
  int ncon1 = NConnects(1);
  int nshape0 = cel0->NSideShapeF(fElVec[0].Side());
  int nshape1 = cel1->NSideShapeF(fElVec[1].Side());
  int numeq = (nshape0+nshape1)*ndof;
  int type = 2;
	
  /**Dimensioning the element matrices*/
  ek.fMat.Redim(numeq,numeq);
  ef.fMat.Redim(numeq, 1);
  int i, n, ncon = ncon0 + ncon1;
	int dim0 = cel0->Dimension(), dim1 = cel1->Dimension();
  ek.fBlock.SetNBlocks(ncon);
  ef.fBlock.SetNBlocks(ncon);
  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);
  for(i=0;i<ncon0;i++) {
    n = cel0->NConnectShapeF(cel0->SideConnectLocId(i,fElVec[0].Side()), GetgOrder())*ndof;
    ek.fBlock.Set(i,n);
    ef.fBlock.Set(i,n);
  }
  for(;i<ncon;i++) {
    n = cel1->NConnectShapeF(cel1->SideConnectLocId(i-ncon0,fElVec[1].Side()), GetgOrder())*ndof;
    ek.fBlock.Set(i,n);
    ef.fBlock.Set(i,n);
  }
  for(i=0;i<ncon;++i){
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }

  /**Acerting the storage to lateral solutions, phis, points, jacobians, ...*/
  TPZFMatrix<STATE> phi0(nshape0,1), phi1(nshape1,1);
  TPZFMatrix<STATE> dphi(fDimension,(nshape0<nshape1)?nshape1:nshape0);
  TPZFMatrix<STATE> axes(3,3,0.);

  TPZFMatrix<STATE> jacobian(fDimension,fDimension);
  REAL detjac;
  TPZVec<REAL> x(3,0.), normal(3,0.);

  /**Remember: cel0 is the element with greater dimension, and it use to calcule the normal*/
  TPZVec<REAL> leftpar(3,0.), leftparel(3,0.);
  TPZVec<REAL> rightpar(3,0.), rightparel(3,0.);
  REAL weight = 0.;
  TPZVec<REAL> leftsolu(ndof,0.);
  TPZVec<REAL> rightsolu(ndof,0.);
  TPZConnect *df;

  /**Obtaining the transformation to compute the same point on the common side
     from left as from right element*/
  TPZGeoElSide gelside = fElVec[0].Reference();
  TPZGeoElSide thisgeoside(gelside);
  gelside = fElVec[1].Reference();
  TPZGeoElSide neighgeoside(gelside);
  int dim = thisgeoside.Dimension();

  TPZTransform<STATE> t(thisgeoside.Dimension());
	TPZTransform<STATE> tleft(t), tright(t);
	if (dim) thisgeoside.SideTransform3(neighgeoside, t);
  double area;
  if(dim0==dim1) {
    area = cel0->Reference()->Volume();
    double area1 = cel1->Reference()->Volume();
    area = (area < area1) ? area : area1;
  }
  else if(dim0<dim1) {
	  area = cel1->Reference()->Volume();
		if(cel0->Material()->Id() <0) {
      int tipo = ((TPZBndCond *)cel0->Material())->Type();
			type = tipo;
      if(tipo==8 || tipo==2) type *= -1;
    }
	}
  else {
	  area = cel0->Reference()->Volume();
		if(cel1->Material()->Id() <0) {
      int tipo = ((TPZBndCond *)cel1->Material())->Type();
			type = tipo;
      if(tipo==8 || tipo==2) type *= -1;
    }
	}

  int ipoint, npoints = fIntegrationRule->NPoints();
  int iv, dfseq, dfvar;

  /**Compute the integral applying numeric integration*/
  for(ipoint=0;ipoint<npoints;ipoint++) {

	  fIntegrationRule->Point(ipoint,leftpar,weight);  //The integration rule was obtained from left element
    /**Working from right element. Compute parameter over side of the right element*/
    t.Apply(leftpar,rightpar);
    /**Working from left element. First transform parameter over side to parameter over element*/
    tleft = cel0->TransformSideToElement(fElVec[0].Side());
		tleft.Apply(leftpar,leftparel);
    /**Exterior normal to left element*/
	switch (fDimension) {
	case 0:
		((TPZGeoEl1d *)(cel0->Reference()))->NormalVector(fElVec[0].Side(), leftparel, normal, axes, jacobian);
		detjac = 1.;
		break;
	case 1:
		if (cel0->Reference()->NCornerNodes() < 4)
			((TPZGeoElT2d *)(cel0->Reference()))->NormalVector(fElVec[0].Side(), leftparel, normal, axes, jacobian);
		else if (cel0->Reference()->NCornerNodes() < 5)
			((TPZGeoElQ2d *)(cel0->Reference()))->NormalVector(fElVec[0].Side(), leftparel, normal, axes, jacobian);
		else  // Cuidado --> Deve verificarse se pode ser cubo, tetrahedra, prisma ou piramide
			((TPZGeoElC3d *)(cel0->Reference()))->NormalVector(fElVec[0].Side(), leftparel, normal, axes, jacobian);

		detjac = jacobian(0, 0);
		break;
	case 2:
		// Verificar se é interface triangulo ou quadrado EEEEE?????
		//Cuidado --> Deve verificarse se pode ser cubo, tetrahedra, prisma ou piramide
		((TPZGeoElC3d *)(cel0->Reference()))->NormalVector(fElVec[0].Side(), leftparel, normal, axes, jacobian);
		detjac = jacobian(0, 0)*jacobian(1, 1) - jacobian(1, 0)*jacobian(0, 1);
		break;
    case 3:
    default:
      PZError << "TInterface::CalcStiff. Define detjac to dimension 3.\n";
      break;
    }
    weight *= fabs(detjac);
    cel0->Reference()->X(leftparel,x);
    /**Compute the matrix phi and solution from interior of left element*/
    cel0->SideShapeFunction(fElVec[0].Side(),leftpar,phi0,dphi);   //Jorge	23/10/99

    for(i=0;i<ndof;i++) {
      leftsolu[i] = 0.;
      rightsolu[i] = 0.;
    }
    iv=0;
    for(i=0;i<ncon0;i++) {
      df = &Connect(i,0);
      dfseq = df->SequenceNumber();
      dfvar = block.Size(dfseq);
      for(n=0;n<dfvar;n++) {
        leftsolu[iv%ndof] += phi0(iv/ndof,0)*block(dfseq,0,n,0);
        iv++;
      }
    }

    /**First transform parameter over side to parameter over element*/
    cel1->SideShapeFunction(fElVec[1].Side(),rightpar,phi1,dphi);   // Jorge 23/10/99
    iv=0;
    for(i=0;i<ncon1;i++) {
      df = &Connect(i,1);
      dfseq = df->SequenceNumber();
      dfvar = block.Size(dfseq);
      for(n=0;n<dfvar;n++) {
        rightsolu[iv%ndof] += phi1(iv/ndof,0)*block(dfseq,0,n,0);
        iv++;
      }
    }

    /**Quando os fluxos numericos precisam de maior informacao, incrementa-se
       nos vetores solucao*/
    if((fFluxType>9 && fFluxType<14) || fFluxType==26) {
			/** Working from right element */
			tright = cel1->TransformSideToElement(fElVec[1].Side());
			tright.Apply(rightpar,rightparel);
      cel0->Reference()->CenterPoint(cel0->Reference()->NSides()-1,leftparel);
      cel1->Reference()->CenterPoint(cel1->Reference()->NSides()-1,rightparel);
      int nshapeel0 = cel0->NShapeF(), nshapeel1 = cel1->NShapeF();
      TPZFMatrix<STATE> phiel0(nshapeel0,1), phiel1(nshapeel1,1);
      int nconel0 = cel0->NConnects(), nconel1 = cel1->NConnects();
      TPZFMatrix<STATE> dphi(fDimension+1,(nshapeel0<nshapeel1)?nshapeel1:nshapeel0);
      cel0->Shape(leftparel,phiel0,dphi);
      cel1->Shape(rightparel,phiel1,dphi);
      leftsolu.Resize(2*ndof);
      rightsolu.Resize(2*ndof);
      iv = 0;
      for(i=0;i<ndof;i++) leftsolu[ndof+i] = rightsolu[ndof+i] = 0.;
      for(i=0;i<nconel0;i++) {
        df = &(cel0->Connect(i));
        dfseq = df->SequenceNumber();
        dfvar = block.Size(dfseq);
        for(n=0;n<dfvar;n++) {
          leftsolu[ndof+iv%ndof] += phiel0(iv/ndof,0)*block(dfseq,0,n,0);
          iv++;
        }
      }
      for(i=0;i<nconel1;i++) {
        df = &(cel1->Connect(i));
        dfseq = df->SequenceNumber();
        dfvar = block.Size(dfseq);
        for(n=0;n<dfvar;n++) {
          rightsolu[ndof+iv%ndof] += phiel1(iv/ndof,0)*block(dfseq,0,n,0);
          iv++;
        }
      }
    }

    /**Compute the contribution of the connects of left element at element matrices*/
		/**If type<0 pode ser utilizado para fluxos especiais ou fluxos fronteira*/
    ((TConservationLaw *)fMaterial)->ContributeOverInterface(x,leftsolu,rightsolu,weight,area,type,axes,phi0,phi1,normal,ek.fMat,ef.fMat);
  }
}

void TInterfaceElement::CalcRhs(TPZElementMatrix &ef) {
  if(!fContribution) {
    ef.fMat = NULL;
    return;
  }
  TPZInterpolatedElement *cel0 = (TPZInterpolatedElement *)fElVec[0].Element();
  TPZInterpolatedElement *cel1 = (TPZInterpolatedElement *)fElVec[1].Element();
/**First compute the dimensions to matrices*/
  TPZBlock<STATE> &block = Mesh()->Block();
  /**Let to use a same diferential equation but different coefficient*/
  int ndof = fMaterial->NStateVariables();
  int ncon0 = NConnects(0);
  int ncon1 = NConnects(1);
  int nshape0 = cel0->NSideShapeF(fElVec[0].Side());
  int nshape1 = cel1->NSideShapeF(fElVec[1].Side());
  int numeq = (nshape0+nshape1)*ndof;
  int dim0 = cel0->Dimension(), dim1 = cel1->Dimension();
	
  /**Acerting the element matrices*/
  ef.fMat.Redim(numeq,1);
  int i, n, ncon = ncon0 + ncon1;
  ef.fBlock.SetNBlocks(ncon);
  ef.fConnect.Resize(ncon);
  for(i=0;i<ncon0;i++) {
    n = cel0->NConnectShapeF(cel0->SideConnectLocId(i,fElVec[0].Side()), GetgOrder())*ndof;
    ef.fBlock.Set(i,n);
  }
  for(;i<ncon;i++) {
    n = cel1->NConnectShapeF(cel1->SideConnectLocId(i-ncon0,fElVec[1].Side()), GetgOrder())*ndof;
    ef.fBlock.Set(i,n);
  }
  for(i=0;i<ncon;++i)
    (ef.fConnect)[i] = ConnectIndex(i);

  /**Acerting the storage to lateral solutions, phis, points, jacobians, ...*/
  TPZFMatrix<STATE> phi0(nshape0,1), phi1(nshape1,1);
  TPZFMatrix<STATE> dphi(fDimension,((nshape0<nshape1)?nshape1:nshape0));
  TPZFMatrix<REAL> axes(3,3,0.);

  TPZFMatrix<REAL> jacobian(fDimension,fDimension);
  REAL detjac;
  TPZVec<REAL> x(3,0.), normal(3,0.);
  /**Remember: cel0 is the element with greater dimension, and it use to calcule the normal*/
  TPZVec<REAL> leftpar(3,0.), leftparel(3,0.);
  TPZVec<REAL> rightpar(3,0.), rightparel(3,0.);
  REAL weight = 0.;
  TPZVec<REAL> leftsolu(ndof,0.),rightsolu(ndof,0.);
  TPZConnect *df;

  /**Obtaining the transformation to compute the same point from left as from right element*/
  TPZGeoElSide gelside = fElVec[0].Reference();
  TPZGeoElSide thisgeoside(gelside);
  gelside = fElVec[1].Reference();
  TPZGeoElSide neighgeoside(gelside);
  TPZTransform<STATE> t(thisgeoside.Dimension());
  thisgeoside.SideTransform3(neighgeoside,t);

  double area;
  if(dim0==dim1) {
    area = cel0->Reference()->Volume();
    double area1 = cel1->Reference()->Volume();
    area = (area < area1) ? area : area1;
  }
  else if(dim0<dim1) area = cel1->Reference()->Volume();
  else area = cel0->Reference()->Volume();
  int ipoint, npoints = fIntegrationRule->NPoints();
  int iv, dfseq, dfvar;

  /**Compute the integral applying numeric integration*/
  for(ipoint=0;ipoint<npoints;ipoint++) {

	  fIntegrationRule->Point(ipoint,leftpar,weight);  //The integration rule was obtained from left element
    /**Working from right element. Compute parameter over side of the right element*/
    t.Apply(leftpar,rightpar);
    /**Working from left element. First transform parameter over side to parameter over element*/
    t = cel0->TransformSideToElement(fElVec[0].Side());
		t.Apply(leftpar,leftparel);
    /**Exterior normal to left element*/
    cel0->Reference()->X(leftparel,x);
    switch(fDimension) {
    case 0:
		((TPZGeoEl1d *)(cel0->Reference()))->NormalVector(fElVec[0].Side(), leftparel, normal, axes, jacobian);
		detjac = 1.;
      break;
    case 1:
		if(cel0->Reference()->NCornerNodes() < 4)
			((TPZGeoElT2d *)(cel0->Reference()))->NormalVector(fElVec[0].Side(), leftparel, normal, axes, jacobian);
		else if(cel0->Reference()->NCornerNodes() < 5)
			((TPZGeoElQ2d *)(cel0->Reference()))->NormalVector(fElVec[0].Side(), leftparel, normal, axes, jacobian);
		else  // Cuidado --> Deve verificarse se pode ser cubo, tetrahedra, prisma ou piramide
			((TPZGeoElC3d *)(cel0->Reference()))->NormalVector(fElVec[0].Side(), leftparel, normal, axes, jacobian);

		detjac = jacobian(0,0);
      break;
    case 2:
		// Verificar se é interface triangulo ou quadrado EEEEE?????
		//Cuidado --> Deve verificarse se pode ser cubo, tetrahedra, prisma ou piramide
	  ((TPZGeoElC3d *)(cel0->Reference()))->NormalVector(fElVec[0].Side(), leftparel, normal, axes, jacobian);
	  detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
      break;
    case 3:
    default:
      PZError << "TInterface::CalcRhs. Define detjac to dimension 3.\n";
      break;
    }
    weight *= fabs(detjac);
    /**Compute the matrix phi and solution from interior of left element*/
    cel0->SideShapeFunction(fElVec[0].Side(),leftpar,phi0,dphi);
    iv=0;
    for(i=0;i<ndof;i++) {
      leftsolu[i] = 0.;
      rightsolu[i] = 0.;
    }
    for(i=0;i<ncon0;i++) {
      df = &Connect(i,0);
      dfseq = df->SequenceNumber();
      dfvar = block.Size(dfseq);
      for(n=0;n<dfvar;n++) {
        leftsolu[iv%ndof] += phi0(iv/ndof,0)*block(dfseq,0,n,0);
        iv++;
      }
    }

    /**First transform parameter over side to parameter over element*/
    cel1->SideShapeFunction(fElVec[1].Side(),rightpar,phi1,dphi);
    iv=0;
    for(i=0;i<ncon1;i++) {
      df = &Connect(i,1);
      dfseq = df->SequenceNumber();
      dfvar = block.Size(dfseq);
      for(n=0;n<dfvar;n++) {
        rightsolu[iv%ndof] += phi1(iv/ndof,0)*block(dfseq,0,n,0);
        iv++;
      }
    }

    /**Quando os fluxos numericos precisam de maior informacao, incrementa-se
       nos vetores solucao*/
    if((fFluxType>9 && fFluxType<14) || fFluxType==26) {
      cel0->Reference()->CenterPoint(cel0->Reference()->NSides()-1,leftparel);
      cel1->Reference()->CenterPoint(cel1->Reference()->NSides() - 1,rightparel);
      int nshapeel0 = cel0->NShapeF(), nshapeel1 = cel1->NShapeF();
      TPZFMatrix<STATE> phiel0(nshapeel0,1), phiel1(nshapeel1,1);
      int nconel0 = cel0->NConnects(), nconel1 = cel1->NConnects();
      TPZFMatrix<STATE> dphi(fDimension+1,(nshapeel0<nshapeel1)?nshapeel1:nshapeel0);
      cel0->Shape(leftparel,phiel0,dphi);
      cel1->Shape(rightparel,phiel1,dphi);
      leftsolu.Resize(2*ndof);
      rightsolu.Resize(2*ndof);
      iv = 0;
      for(i=0;i<ndof;i++) leftsolu[ndof+i] = rightsolu[ndof+i] = 0.;
      for(i=0;i<nconel0;i++) {
        df = &(cel0->Connect(i));
        dfseq = df->SequenceNumber();
        dfvar = block.Size(dfseq);
        for(n=0;n<dfvar;n++) {
          leftsolu[ndof+iv%ndof] += phiel0(iv/ndof,0)*block(dfseq,0,n,0);
          iv++;
        }
      }
      for(i=0;i<nconel1;i++) {
        df = &(cel1->Connect(i));
        dfseq = df->SequenceNumber();
        dfvar = block.Size(dfseq);
        for(n=0;n<dfvar;n++) {
          rightsolu[ndof+iv%ndof] += phiel1(iv/ndof,0)*block(dfseq,0,n,0);
          iv++;
        }
      }
    }

    /**Compute the contribution of the connects of left element at element matrices*/
    ((TConservationLaw *)fMaterial)->ContributeOverInterface(x,leftsolu,rightsolu,weight,area,
	      fFluxType,axes,phi0,phi1,normal,ef.fMat,ef.fMat);
  }
}
void TInterfaceElement::CalcResidual(TPZElementMatrix &ef){
	TPZElementMatrix ek;
	CalcStiff(ek,ef);
}

/** To fill connect indexes and related TPZCompElSide into the fConnecIndexes
    and fElVec vectors. The interface element can to consider only
    the internal connects of the discontinuous elements into elvec. At least
    must to exist equal number of computational element as internal connects
TPZConnect &TInterfaceElement::GetConnects(int icon,TPZStack<CompElSide> &elvec) {
#ifndef NOTDEBUG
  if(icon < 0 || icon > elside.Element()->NSideConnects(elside.Side())-1) {
    PZError << "TInterfaceElement::GetConnect. Parameter icon is bad.\n";
    //It must to be improved. Here we must to return a TPZConnect()
    return (TPZInterpolatedElement *)(elside.Element())->ConnectSide(0,elside.Side());
  }
#endif NOTDEBUG

  /**To detect the connects related with this interface element
  TPZInterpolatedElement *el;
  int i, ncon;
  for(i=0;i<nelsides;i++) {
    el = (TPZInterpolatedElement *)elvec[i].Element();
    ncon = el->NConnects();
#ifndef NOTDEBUG
    if(Material() != el->Material()) {
      PZError << "TInterfaceElement::GetConnectIndexes. Exist connected element with different material.\n";
      SetMaterial(el->Material());
    }
    if(ncon == el->Reference()->NSides()) continue;
#endif NOTDEBUG
    /** To fill internal connects from the all connected element in elvec
    fConnectIndexes.Push(el->ConnectIndex(ncon-1));
    /** Storing the element pointers and the respective side of the element
        connecting this interface
    fElVec.Push(elvec[i]);

  }

  if(onlyinternals) return;
  /** Filling connects over the common boundary of the neighboard elements into elvec
  int side = elvec[i].Side();
  ncon = el->NSideConnects(side)-1;  //Remember: The internal connect was already considered
  for(i=0;i<ncon;i++)
    fConnectIndexes.Push(el->SideConnectIndex(i,side));
}
*/


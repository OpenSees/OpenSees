/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.25 $
// $Date: 2003-05-15 21:30:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/element/NLBeamColumn3d.cpp,v $
                                                                        
                                                                        
// Written: Remo Magalhaes de Souza (rmsouza.ce.berkeley.edu) on 03/99 
// Revised: rms 06/99 (mass matrix)
//          rms 07/99 (using setDomain)
//          rms 08/99 (included P-Delta effect)
//	    fmk 10/99 setResponse() & getResponse()
//          rms 11/99 (included rigid joint offsets)
//          rms 04/00 (using transformation class w/ linear or corotational transf)
//          rms 04/00 (generalized to iterative/non-iterative algorithm)
//          mhs 06/00 (using new section class w/ variable dimensions)
//          rms 06/00 (torsional stiffness considered at the section level)
//          rms 06/00 (making copy of the sections)
//          rms 06/00 (storing section history variables at the element level)
//          rms 07/00 (new state determination procedure, no need to store fscommit) 
//
// Purpose: This file contains the implementation for the NLBeamColumn3d class.
//          NLBeamColumn3d is a materially nonlinear flexibility based frame element.

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <Information.h>
#include <NLBeamColumn3d.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

#define  NDM   3         // dimension of the problem (3d)
#define  NND   6         // number of nodal dof's
#define  NEGD 12         // number of element global dof's
#define  NEBD  6         // number of element dof's in the basic system

#define DefaultLoverGJ 1.0e-10

Matrix NLBeamColumn3d::theMatrix(12,12);
Vector NLBeamColumn3d::theVector(12);
double NLBeamColumn3d::workArea[200];
GaussLobattoQuadRule1d01 NLBeamColumn3d::quadRule;

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
NLBeamColumn3d::NLBeamColumn3d():
Element(0,ELE_TAG_NLBeamColumn3d), connectedExternalNodes(2), 
nSections(0), sections(0), crdTransf(0), 
rho(0), maxIters(0), tol(0), initialFlag(0), isTorsion(false),
load(NEGD),
kv(NEBD,NEBD), Se(NEBD), 
kvcommit(NEBD,NEBD), Secommit(NEBD),
fs(0), vs(0), Ssr(0), vscommit(0), sp(0), Ki(0)
{
  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;
  
  theNodes[0] = 0;
  theNodes[1] = 0;
}

// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points. 
// allocates the necessary space needed by each object
NLBeamColumn3d::NLBeamColumn3d (int tag, int nodeI, int nodeJ, 
                                int numSections, SectionForceDeformation *sectionPtrs[],
                                CrdTransf3d &coordTransf, double massDensPerUnitLength,
                                int maxNumIters, double tolerance):
Element(tag,ELE_TAG_NLBeamColumn3d), connectedExternalNodes(2), 
nSections(numSections), sections(0), crdTransf(0),
rho(massDensPerUnitLength), maxIters(maxNumIters), tol(tolerance),
initialFlag(0), isTorsion(false),
load(NEGD), 
kv(NEBD,NEBD), Se(NEBD),  
kvcommit(NEBD,NEBD), Secommit(NEBD),
fs(0), vs(0), Ssr(0), vscommit(0), sp(0), Ki(0)
{
   connectedExternalNodes(0) = nodeI;
   connectedExternalNodes(1) = nodeJ;    

   // get copy of the sections
   
   if (!sectionPtrs)
   {
       opserr << "Error: NLBeamColumn3d::NLBeamColumn3d:  invalid section pointer ";
       exit(-1);
   }	  
   
   sections = new SectionForceDeformation *[nSections];
   if (!sections)
   {
       opserr << "Error: NLBeamColumn3d::NLBeamColumn3d: could not alocate section pointer";
       exit(-1);
   }  

   for (int i = 0; i < nSections; i++)
   {
      if (!sectionPtrs[i])
      {
	  opserr << "Error: NLBeamColumn3d::NLBeamColumn3d: section pointer " << i << endln;
          exit(-1);
      }  
       
      sections[i] = sectionPtrs[i]->getCopy();
      if (!sections[i])
      {
	  opserr << "Error: NLBeamColumn3d::NLBeamColumn3d: could not create copy of section " << i << endln;
          exit(-1);
      }

	  int order = sections[i]->getOrder();
	  const ID &code = sections[i]->getType();
	  for (int j = 0; j < order; j++) {
		if (code(j) == SECTION_RESPONSE_T)
			isTorsion = true;
	  }
   }

   if (!isTorsion)
     opserr << "NLBeamColumn3d::NLBeamColumn3d -- no torsion detected in sections, " <<
       "continuing with element torsional stiffness GJ/L = " << 1.0/DefaultLoverGJ;
   
   // get copy of the transformation object   
   crdTransf = coordTransf.getCopy(); 
   if (!crdTransf)
   {
     opserr << "Error: NLBeamColumn3d::NLBeamColumn3d: could not create copy of coordinate transformation object" << endln;
     exit(-1);
   }

   // alocate section flexibility matrices and section deformation vectors
   fs  = new Matrix [nSections];
   if (!fs) {
     opserr << "NLBeamColumn3d::NLBeamColumn3d() -- failed to allocate fs array";
     exit(-1);
   }
   
   vs = new Vector [nSections];
   if (!vs) {
     opserr << "NLBeamColumn3d::NLBeamColumn3d() -- failed to allocate vs array";
     exit(-1);
   }

   Ssr = new Vector [nSections];
   if (!Ssr) {
     opserr << "NLBeamColumn3d::NLBeamColumn3d() -- failed to allocate Ssr array";
     exit(-1);
   }
 
   vscommit = new Vector [nSections];
   if (!vscommit)
     {
       opserr << "NLBeamColumn3d::NLBeamColumn3d() -- failed to allocate vscommit array";   
       exit(-1);
   }

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  theNodes[0] = 0;
  theNodes[1] = 0;
}


// ~NLBeamColumn3d():
// 	destructor
//      delete must be invoked on any objects created by the object
NLBeamColumn3d::~NLBeamColumn3d()
{
   if (sections)
   {
      for (int i=0; i < nSections; i++)
         if (sections[i])
            delete sections[i];
      delete [] sections;
   }  
   
   if (fs) 
     delete [] fs;

   if (vs) 
     delete [] vs;

   if (Ssr) 
     delete [] Ssr;

   if (vscommit) 
     delete [] vscommit;

   if (crdTransf)
     delete crdTransf;   

   if (sp != 0)
     delete sp;

   if (Ki != 0)
     delete Ki;
}



int
NLBeamColumn3d::getNumExternalNodes(void) const
{
  return 2;
}


const ID &
NLBeamColumn3d::getExternalNodes(void)
{
   return connectedExternalNodes;
}

Node **
NLBeamColumn3d::getNodePtrs(void) 
{
  return theNodes;
}

int
NLBeamColumn3d::getNumDOF(void) 
{
   return NEGD;
}

void
NLBeamColumn3d::setDomain(Domain *theDomain)
{
   // check Domain is not null - invoked when object removed from a domain
   if (theDomain == 0) {
      theNodes[0] = 0;
      theNodes[1] = 0;
      return;
   }

   // get pointers to the nodes
   
   int Nd1 = connectedExternalNodes(0);  
   int Nd2 = connectedExternalNodes(1);
   
   theNodes[0] = theDomain->getNode(Nd1);
   theNodes[1] = theDomain->getNode(Nd2);  

   if (theNodes[0] == 0)
   {
      opserr << "NLBeamColumn3d::setDomain: Nd1: ";
      opserr << Nd1 << "does not exist in model\n";
      exit(0);
   }

   if (theNodes[1] == 0) 
   {
      opserr << "NLBeamColumn3d::setDomain: Nd2: ";
      opserr << Nd2 << "does not exist in model\n";
      exit(0);
   }

   // call the DomainComponent class method 
   this->DomainComponent::setDomain(theDomain);
    
   // ensure connected nodes have correct number of dof's
   int dofNode1 = theNodes[0]->getNumberDOF();
   int dofNode2 = theNodes[1]->getNumberDOF();
   
   if ((dofNode1 != NND) || (dofNode2 != NND))
   {
      opserr << "NLBeamColumn3d::setDomain(): Nd2 or Nd1 incorrect dof ";
      exit(0);
   }
   

   // initialize the transformation
   if (crdTransf->initialize(theNodes[0], theNodes[1]))
   {
      opserr << "NLBeamColumn3d::setDomain(): Error initializing coordinate transformation";  
      exit(0);
   }
    
   // get element length
   double L = crdTransf->getInitialLength();
   if (L == 0.0)
   {
      opserr << "NLBeamColumn3d::setDomain(): Zero element length:" << this->getTag();  
      exit(0);
   }

   if (initialFlag != 2) 
     this->initializeSectionHistoryVariables();
}


int
NLBeamColumn3d::commitState()
{

   int err = 0;
   int i = 0;

   // call element commitState to do any base class stuff
   if ((err = this->Element::commitState()) != 0) {
     opserr << "NLBeamColumn3d::commitState () - failed in base class";
     return err;
   }    

   do {
      vscommit[i] = vs[i];
      err = sections[i++]->commitState();  
   } while (err == 0 && i < nSections);
   
   if (err)
      return err;
   
   // commit the transformation between coord. systems
   if ((err = crdTransf->commitState()) != 0)
      return err;
      
   // commit the element variables state
   kvcommit = kv;
   Secommit = Se;

   return err;
}


int 
NLBeamColumn3d::revertToLastCommit()
{
   int err;
   int i = 0;
      
   do {
     vs[i] = vscommit[i];
     err = sections[i]->revertToLastCommit();
     
     sections[i]->setTrialSectionDeformation(vs[i]);        
     Ssr[i] = sections[i]->getStressResultant();
     fs[i]  = sections[i]->getSectionFlexibility();
      
     i++;
   } while (err == 0 && i < nSections);
   
   if (err)
      return err;
   
   // revert the transformation to last commit
   if ((err = crdTransf->revertToLastCommit()) != 0)
      return err;
     
   // revert the element state to last commit
   Se   = Secommit;
   kv   = kvcommit;
   
   initialFlag = 0;
   // return this->update();

   return err;
}


int 
NLBeamColumn3d::revertToStart()
{
   // revert the sections state to start
   int err;
   int i = 0;
     
   do {
      fs[i].Zero();
      vs[i].Zero();
      Ssr[i].Zero();
      err = sections[i++]->revertToStart();
 
   } while (err == 0 && i < nSections);

   if (err)
      return err;
   
   // revert the transformation to start
   if ((err = crdTransf->revertToStart()) != 0)
      return err;
  
   // revert the element state to start
   Secommit.Zero();
   kvcommit.Zero();

   Se.Zero();
   kv.Zero();

   initialFlag = 0;
   // this->update(); 
   
   return err;
}



const Matrix &
NLBeamColumn3d::getTangentStiff(void)
{
  // Will remove once we clean up the corotational 3d transformation -- MHS
  crdTransf->update();

  return crdTransf->getGlobalStiffMatrix(kv, Se);
}


const Matrix &
NLBeamColumn3d::getInitialStiff(void)
{

  // check for quick return
  if (Ki != 0)
    return *Ki;

  // get integration point positions and weights
  const Matrix &xi_pt = quadRule.getIntegrPointCoords(nSections);
  const Vector &weight = quadRule.getIntegrPointWeights(nSections);
  
  static Matrix f(NEBD,NEBD);   // element flexibility matrix

  static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse
  int i;
  
  I.Zero();
  for (i=0; i<NEBD; i++)
    I(i,i) = 1.0;
  
  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;
  
  // initialize f and vr for integration
  f.Zero();

  for (i=0; i<nSections; i++) {
    int order      = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    
    Vector Ss(workArea, order);
    Vector dSs(&workArea[order], order);
    Vector dvs(&workArea[2*order], order);
    
    Matrix fb(&workArea[3*order], order, NEBD);
	
    double xL  = xi_pt(i,0);
    double xL1 = xL-1.0;
    
    // get section flexibility matrix
    const Matrix &fSec = sections[i]->getInitialFlexibility();
    
    // f = f + (b^ fs * b) * weight(i);
    //f.addMatrixTripleProduct(1.0, b[i], fs[i], weight(i));
    int jj, ii;
    fb.Zero();
    double tmp;
    for (ii = 0; ii < order; ii++) {
      switch(code(ii)) {
      case SECTION_RESPONSE_P:
	for (jj = 0; jj < order; jj++)
	  fb(jj,0) += fSec(jj,ii)*weight(i);
	break;
      case SECTION_RESPONSE_MZ:
	for (jj = 0; jj < order; jj++) {
	  tmp = fSec(jj,ii)*weight(i);
	  fb(jj,1) += xL1*tmp;
	  fb(jj,2) += xL*tmp;
	}
	break;
      case SECTION_RESPONSE_VY:
	for (jj = 0; jj < order; jj++) {
	  tmp = oneOverL*fSec(jj,ii)*weight(i);
	  fb(jj,1) += tmp;
	  fb(jj,2) += tmp;
	}
	break;
      case SECTION_RESPONSE_MY:
	for (jj = 0; jj < order; jj++) {
	  tmp = fSec(jj,ii)*weight(i);
	  fb(jj,3) += xL1*tmp;
	  fb(jj,4) += xL*tmp;
	}
	break;
      case SECTION_RESPONSE_VZ:
	for (jj = 0; jj < order; jj++) {
	  tmp = oneOverL*fSec(jj,ii)*weight(i);
	  fb(jj,3) += tmp;
	  fb(jj,4) += tmp;
	}
	break;
      case SECTION_RESPONSE_T:
	for (jj = 0; jj < order; jj++)
	  fb(jj,5) += fSec(jj,ii)*weight(i);
	break;
      default:
	break;
      }
    }
    for (ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
	for (jj = 0; jj < 6; jj++)
	  f(0,jj) += fb(ii,jj);
	break;
      case SECTION_RESPONSE_MZ:
	for (jj = 0; jj < 6; jj++) {
	  tmp = fb(ii,jj);
	  f(1,jj) += xL1*tmp;
	  f(2,jj) += xL*tmp;
	}
	break;
      case SECTION_RESPONSE_VY:
	for (jj = 0; jj < 6; jj++) {
	  tmp = oneOverL*fb(ii,jj);
	  f(1,jj) += tmp;
	  f(2,jj) += tmp;
	}
	break;
      case SECTION_RESPONSE_MY:
	for (jj = 0; jj < 6; jj++) {
	  tmp = fb(ii,jj);
	  f(3,jj) += xL1*tmp;
	  f(4,jj) += xL*tmp;
	}
	break;
      case SECTION_RESPONSE_VZ:
	for (jj = 0; jj < 6; jj++) {
	  tmp = oneOverL*fb(ii,jj);
	  f(3,jj) += tmp;
	  f(4,jj) += tmp;
	}
	break;
      case SECTION_RESPONSE_T:
	for (jj = 0; jj < 6; jj++)
	  f(5,jj) += fb(ii,jj);
	break;
      default:
	break;
      }
    }
  }

  f  *= L;
  
  if (!isTorsion)
    f(5,5) = DefaultLoverGJ;
  
  // calculate element stiffness matrix
  static Matrix kvInit(NEBD, NEBD);
  if (f.Solve(I,kvInit) < 0)
    opserr << "NLBeamColumn3d::updateElementState() - could not invert flexibility\n";


  // set Ki
  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kvInit));

  return *Ki;
}

const Vector &
NLBeamColumn3d::getResistingForce(void)
{
  // Will remove once we clean up the corotational 3d transformation -- MHS
  crdTransf->update();

  Vector p0Vec(p0, 5);

  return crdTransf->getGlobalResistingForce(Se, p0Vec);
}


void
NLBeamColumn3d::initializeSectionHistoryVariables (void)
{
    for (int i = 0; i < nSections; i++)
    {
	int order = sections[i]->getOrder();
	
	fs[i] = Matrix(order,order);
	vs[i] = Vector(order);
        Ssr[i] = Vector(order);
	
	vscommit[i] = Vector(order);
    }
}



int 
NLBeamColumn3d::update(void)
{

  // if have completed a recvSelf() - do a revertToLastCommit
  // to get Ssr, etc. set correctly
  if (initialFlag == 2)
    this->revertToLastCommit();

    // update the transformation
    crdTransf->update();
       
    // get basic displacements and increments
    static Vector v(NEBD);
    static Vector dv(NEBD);
     
    v = crdTransf->getBasicTrialDisp();    
    dv = crdTransf->getBasicIncrDeltaDisp();    

    // get integration point positions and weights
    const Matrix &xi_pt = quadRule.getIntegrPointCoords(nSections);
    const Vector &weight = quadRule.getIntegrPointWeights(nSections);
     
    static Vector vr(NEBD);       // element residual displacements
    static Matrix f(NEBD,NEBD);   // element flexibility matrix

    static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse
    double dW;                    // section strain energy (work) norm 
    int i;
    
    I.Zero();
    for (i=0; i<NEBD; i++)
      I(i,i) = 1.0;

    // calculate nodal force increments and update nodal forces
    static Vector dSe(NEBD);

    // dSe = kv * dv;
    dSe.addMatrixVector(0.0, kv, dv, 1.0);

    double L = crdTransf->getInitialLength();
    double oneOverL  = 1.0/L;

    for (int j=0; j < maxIters; j++) {
      Se += dSe;
  
      // initialize f and vr for integration
      f.Zero();
      vr.Zero();

      for (i=0; i<nSections; i++) {
	int order      = sections[i]->getOrder();
	const ID &code = sections[i]->getType();

	Vector Ss(workArea, order);
	Vector dSs(&workArea[order], order);
	Vector dvs(&workArea[2*order], order);
	
	Matrix fb(&workArea[3*order], order, NEBD);
	
	double xL  = xi_pt(i,0);
	double xL1 = xL-1.0;

	// calculate total section forces
	// Ss = b*Se + bp*currDistrLoad;
	//Ss.addMatrixVector(0.0, b[i], Se, 1.0);
	int ii;
	for (ii = 0; ii < order; ii++) {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_P:
	    Ss(ii) = Se(0);
	    break;
	  case SECTION_RESPONSE_MZ:
	    Ss(ii) = xL1*Se(1) + xL*Se(2);
	    break;
	  case SECTION_RESPONSE_VY:
	    Ss(ii) = oneOverL*(Se(1)+Se(2));
	    break;
	  case SECTION_RESPONSE_MY:
	    Ss(ii) = xL1*Se(3) + xL*Se(4);
	    break;
	  case SECTION_RESPONSE_VZ:
	    Ss(ii) = oneOverL*(Se(3)+Se(4));
	    break;
	  case SECTION_RESPONSE_T:
	    Ss(ii) = Se(5);
	    break;
	  default:
	    Ss(ii) = 0.0;
	    break;
	  }
	}

	// Add the effects of element loads, if present
	if (sp != 0) {
	  const Matrix &s_p = *sp;
	  for (ii = 0; ii < order; ii++) {
	    switch(code(ii)) {
	    case SECTION_RESPONSE_P:
	      Ss(ii) += s_p(0,i);
	      break;
	    case SECTION_RESPONSE_MZ:
	      Ss(ii) += s_p(1,i);
	      break;
	    case SECTION_RESPONSE_VY:
	      Ss(ii) += s_p(2,i);
	      break;
	    case SECTION_RESPONSE_MY:
	      Ss(ii) += s_p(3,i);
	      break;
	    case SECTION_RESPONSE_VZ:
	      Ss(ii) += s_p(4,i);
	      break;
	    default:
	      break;
	    }
	  }
	}

	  dSs = Ss;
	  dSs.addVector(1.0, Ssr[i], -1.0);  // dSs = Ss - Ssr[i];
 
	  // compute section deformation increments
	  //       vs += fs * dSs;     
	  dvs.addMatrixVector(0.0, fs[i], dSs, 1.0);


	  if (initialFlag != 0)
	    vs[i] += dvs;

	  sections[i]->setTrialSectionDeformation(vs[i]);
	  
	  // get section resisting forces
	  Ssr[i] = sections[i]->getStressResultant();

	  // get section flexibility matrix
          fs[i] = sections[i]->getSectionFlexibility();

	  // calculate section residual deformations
	  // dvs = fs * (Ss - Ssr);
	  dSs = Ss;
	  dSs.addVector(1.0, Ssr[i], -1.0);  // dSs = Ss - Ssr[i];

	  dvs.addMatrixVector(0.0, fs[i], dSs, 1.0);

          // integrate element flexibility matrix
	  // f = f + (b^ fs * b) * weight(i);
      	  //f.addMatrixTripleProduct(1.0, b[i], fs[i], weight(i));
	  int jj;
	  const Matrix &fSec = fs[i];
	  fb.Zero();
	  double tmp;
	  for (ii = 0; ii < order; ii++) {
	    switch(code(ii)) {
	    case SECTION_RESPONSE_P:
	      for (jj = 0; jj < order; jj++)
		fb(jj,0) += fSec(jj,ii)*weight(i);
	      break;
	    case SECTION_RESPONSE_MZ:
	      for (jj = 0; jj < order; jj++) {
		tmp = fSec(jj,ii)*weight(i);
		fb(jj,1) += xL1*tmp;
		fb(jj,2) += xL*tmp;
	      }
	      break;
	    case SECTION_RESPONSE_VY:
	      for (jj = 0; jj < order; jj++) {
		tmp = oneOverL*fSec(jj,ii)*weight(i);
		fb(jj,1) += tmp;
		fb(jj,2) += tmp;
	      }
	      break;
	    case SECTION_RESPONSE_MY:
	      for (jj = 0; jj < order; jj++) {
		tmp = fSec(jj,ii)*weight(i);
		fb(jj,3) += xL1*tmp;
		fb(jj,4) += xL*tmp;
	      }
	      break;
	    case SECTION_RESPONSE_VZ:
	      for (jj = 0; jj < order; jj++) {
		tmp = oneOverL*fSec(jj,ii)*weight(i);
		fb(jj,3) += tmp;
		fb(jj,4) += tmp;
	      }
	      break;
	    case SECTION_RESPONSE_T:
	      for (jj = 0; jj < order; jj++)
		fb(jj,5) += fSec(jj,ii)*weight(i);
	      break;
	    default:
	      break;
	    }
	  }
	  for (ii = 0; ii < order; ii++) {
	    switch (code(ii)) {
	    case SECTION_RESPONSE_P:
	      for (jj = 0; jj < 6; jj++)
		f(0,jj) += fb(ii,jj);
	      break;
	    case SECTION_RESPONSE_MZ:
	      for (jj = 0; jj < 6; jj++) {
		tmp = fb(ii,jj);
		f(1,jj) += xL1*tmp;
		f(2,jj) += xL*tmp;
	      }
	      break;
	    case SECTION_RESPONSE_VY:
	      for (jj = 0; jj < 6; jj++) {
		tmp = oneOverL*fb(ii,jj);
		f(1,jj) += tmp;
		f(2,jj) += tmp;
	      }
	      break;
	    case SECTION_RESPONSE_MY:
	      for (jj = 0; jj < 6; jj++) {
		tmp = fb(ii,jj);
		f(3,jj) += xL1*tmp;
		f(4,jj) += xL*tmp;
	      }
	      break;
	    case SECTION_RESPONSE_VZ:
	      for (jj = 0; jj < 6; jj++) {
		tmp = oneOverL*fb(ii,jj);
		f(3,jj) += tmp;
		f(4,jj) += tmp;
	      }
	      break;
	    case SECTION_RESPONSE_T:
	      for (jj = 0; jj < 6; jj++)
		f(5,jj) += fb(ii,jj);
	      break;
	    default:
	      break;
	    }
	  }

	  // integrate residual deformations
	  // vr += (b^ (vs + dvs)) * weight(i);
	  //vr.addMatrixTransposeVector(1.0, b[i], vs[i] + dvs, weight(i));
	  dvs.addVector(1.0, vs[i], 1.0);
	  double dei;
	  for (ii = 0; ii < order; ii++) {
	    dei = dvs(ii)*weight(i);
	    switch(code(ii)) {
	    case SECTION_RESPONSE_P:
	      vr(0) += dei; break;
	    case SECTION_RESPONSE_MZ:
	      vr(1) += xL1*dei; vr(2) += xL*dei; break;
	    case SECTION_RESPONSE_VY:
	      tmp = oneOverL*dei;
	      vr(1) += tmp; vr(2) += tmp; break;
	    case SECTION_RESPONSE_MY:
	      vr(3) += xL1*dei; vr(4) += xL*dei; break;
	    case SECTION_RESPONSE_VZ:
	      tmp = oneOverL*dei;
	      vr(3) += tmp; vr(4) += tmp; break;
	    case SECTION_RESPONSE_T:
	      vr(5) += dei; break;
	    default:
	      break;
	    }
	  }
      }  
      
      f  *= L;
      vr *= L;

      if (!isTorsion) {
	f(5,5) = DefaultLoverGJ;
	vr(5) = Se(5)*DefaultLoverGJ;
      }

      // calculate element stiffness matrix
      if (f.Invert(kv) < 0)
	opserr << "NLBeamColumn3d::updateElementState() - could not invert flexibility\n";

      // dv = v - vr;
      dv = v;
      dv.addVector(1.0, vr, -1.0);
      
      // dSe = kv * dv;
      dSe.addMatrixVector(0.0, kv, dv, 1.0);

      dW = dv^ dSe;
      if (fabs(dW) < tol)
        break;

      if (maxIters == 1) {
	opserr << "NLBeamColumn3d::updateElementState() - element: " << this->getTag() << " failed to converge but going on\n";
	break;
      }
      if (j == (maxIters-1)) {
	opserr << "NLBeamColumn3d::updateElementState() - element: " << this->getTag() << " failed to converge\n";
	opserr << "dW: " << dW  << "\n dv: " << dv << " dSe: " << dSe << endln;
	return -1;
      }
    }
      
    // determine resisting forces
    Se += dSe;

    initialFlag = 1;

    return 0;
}



void NLBeamColumn3d::getGlobalDispls(Vector &dg) const
{
   // determine global displacements
   const Vector &disp1 = theNodes[0]->getTrialDisp();
   const Vector &disp2 = theNodes[1]->getTrialDisp();

   for (int i = 0; i < NND; i++)
   {
      dg(i)     = disp1(i);
      dg(i+NND) = disp2(i);
   }
}



void NLBeamColumn3d::getGlobalAccels(Vector &ag) const
{
   // determine global displacements
   const Vector &accel1 = theNodes[0]->getTrialAccel();
   const Vector &accel2 = theNodes[1]->getTrialAccel();

   for (int i = 0; i < NND; i++)
   {
      ag(i)     = accel1(i);
      ag(i+NND) = accel2(i);
   }
}


void NLBeamColumn3d::getForceInterpolatMatrix(double xi, Matrix &b, const ID &code)
{
   b.Zero();
   double L = crdTransf->getInitialLength();
    
   for (int i = 0; i < code.Size(); i++)
   {
      switch (code(i))
      {
         case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
            b(i,1) = xi - 1.0;
	    b(i,2) = xi;
	    break;
	 case SECTION_RESPONSE_P:		// Axial, P, interpolation
	    b(i,0) = 1.0;
	    break;
	 case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
	    b(i,1) = b(i,2) = 1.0/L;
	    break;
	 case SECTION_RESPONSE_MY:		// Moment, My, interpolation
	    b(i,3) = xi - 1.0;
	    b(i,4) = xi;
	    break;
	 case SECTION_RESPONSE_VZ:		// Shear, Vz, interpolation
	    b(i,3) = b(i,4) = 1.0/L;
	    break;
	 case SECTION_RESPONSE_T:		// Torque, T, interpolation
	    b(i,5) = 1.0;
	    break;
	 default:
	    break;
      }	
   }
}


void NLBeamColumn3d::getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code)
{
   bp.Zero();
   double L = crdTransf->getInitialLength();

   for (int i = 0; i < code.Size(); i++)
   {
      switch (code(i))
      {
	 case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
	    bp(i,1) = xi*(xi-1)*L*L/2;
	    break;
	 case SECTION_RESPONSE_P:		// Axial, P, interpolation
	    bp(i,0) = (1-xi)*L;
	    break;
	 case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
	    bp(i,1) = (xi-0.5)*L;
	    break;
	 case SECTION_RESPONSE_MY:		// Moment, My, interpolation
	    bp(i,2) = xi*(1-xi)*L*L/2;
	    break;
	 case SECTION_RESPONSE_VZ:		// Shear, Vz, interpolation
	    bp(i,2) = (0.5-xi)*L;
	    break;
	 case SECTION_RESPONSE_T:		// Torsion, T, interpolation
	    break;
	 default:
	    break;
      }
   }
}

    
const Matrix &
NLBeamColumn3d::getMass(void)
{ 
  double L = crdTransf->getInitialLength();

  theMatrix(0,0) = theMatrix(1,1) = theMatrix(2,2) =
    theMatrix(6,6) = theMatrix(7,7) = theMatrix(8,8) = 0.5*rho*L;
  
  return theMatrix;
}


void 
NLBeamColumn3d::zeroLoad(void)
{
  if (sp != 0)
    sp->Zero();

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  load.Zero();
}

int
NLBeamColumn3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  
  if (sp == 0) {
    sp = new Matrix(5,nSections);
    if (sp == 0) {
      opserr << "NLBeamColumn3d::addLoad -- out of memory\n";
      exit(-1);
    }
  }

  const Matrix &xi_pt = quadRule.getIntegrPointCoords(nSections);
  double L = crdTransf->getInitialLength();

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0)*loadFactor;  // Transverse
    double wz = data(1)*loadFactor;  // Transverse
    double wx = data(2)*loadFactor;  // Axial

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < nSections; i++) {
      double x = xi_pt(i,0)*L;
      // Axial
      s_p(0,i) += wx*(L-x);
      // Moment
      s_p(1,i) += wy*0.5*x*(x-L);
      // Shear
      s_p(2,i) += wy*(x-0.5*L);
      // Moment
      s_p(3,i) += wz*0.5*x*(L-x);
      // Shear
      s_p(4,i) += wz*(x-0.5*L);
    }
  
    // Accumulate reactions in basic system
    p0[0] -= wx*L;
    double V;
    V = 0.5*wy*L;
    p0[1] -= V;
    p0[2] -= V;
    V = 0.5*wz*L;
    p0[3] -= V;
    p0[4] -= V;
  }
  else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py = data(0)*loadFactor;
    double Pz = data(1)*loadFactor;
    double N  = data(2)*loadFactor;
    double aOverL = data(3);
    double a = aOverL*L;

    double Vy2 = Py*aOverL;
    double Vy1 = Py-Vy2;

    double Vz2 = Pz*aOverL;
    double Vz1 = Pz-Vz2;

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < nSections; i++) {
      double x = xi_pt(i,0)*L;
      if (x <= a) {
	s_p(0,i) += N;
	s_p(1,i) -= x*Vy1;
	s_p(2,i) -= Vy1;
	s_p(3,i) -= x*Vz1;
	s_p(4,i) -= Vz1;
      }
      else {
	s_p(1,i) -= (L-x)*Vy2;
	s_p(2,i) += Vy2;
	s_p(3,i) -= (L-x)*Vz2;
	s_p(4,i) += Vz2;
      }
    }

    // Accumulate reactions in basic system
    p0[0] -= N;
    p0[1] -= Vy1;
    p0[2] -= Vy2;
    p0[3] -= Vz1;
    p0[4] -= Vz2;
  }

  else {
    opserr << "NLBeamColumn3d::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }

  return 0;
}


int 
NLBeamColumn3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);    

  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;

  load(0) -= m*Raccel1(0);
  load(1) -= m*Raccel1(1);
  load(2) -= m*Raccel1(2);
  load(6) -= m*Raccel2(0);
  load(7) -= m*Raccel2(1);
  load(8) -= m*Raccel2(2);

  return 0;
}


const Vector &
NLBeamColumn3d::getResistingForceIncInertia()
{	
  // Check for a quick return
  if (rho == 0.0)
    theVector = this->getResistingForce();

  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
    // Compute the current resisting force
    theVector = this->getResistingForce();
    
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;
    
    theVector(0) += m*accel1(0);
    theVector(1) += m*accel1(1);
    theVector(2) += m*accel1(2);
    theVector(6) += m*accel2(0);
    theVector(7) += m*accel2(1);
    theVector(8) += m*accel2(2);

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      theVector += this->getRayleighDampingForces();
    
  } else {

    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      theVector += this->getRayleighDampingForces();
  }

  return theVector;
}

int
NLBeamColumn3d::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j , k;
  int loc = 0;
  
  static ID idData(9);
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);
  idData(3) = nSections;
  idData(4) = maxIters;
  idData(5) = initialFlag;
  idData(6) = crdTransf->getClassTag();
  int crdTransfDbTag  = crdTransf->getDbTag();
  if (crdTransfDbTag  == 0) {
      crdTransfDbTag = theChannel.getDbTag();
      if (crdTransfDbTag  != 0) {
	   crdTransf->setDbTag(crdTransfDbTag);
       }
  }
  idData(7) = crdTransfDbTag;
  idData(8) = (isTorsion) ? 1 : 0;
  
  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
     opserr << "NLBeamColumn3d::sendSelf() - %s\n",
	     		     "failed to send ID data";
     return -1;
  }    
  
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "NLBeamColumn3d::sendSelf() - %s\n",
	     		     "failed to send crdTranf";
     return -1;
  }      
  
  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //

  ID idSections(2*nSections);
  loc = 0;
  for (i = 0; i<nSections; i++) {
    int sectClassTag = sections[i]->getClassTag();
    int sectDbTag = sections[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      sections[i]->setDbTag(sectDbTag);
    }

    idSections(loc) = sectClassTag;
    idSections(loc+1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0)  {
    opserr << "NLBeamColumn3d::sendSelf() - %s\n",
			    "failed to send ID data";
    return -1;
  }    

  //
  // send the sections
  //
  
  for (j = 0; j<nSections; j++) {
    if (sections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "NLBeamColumn3d::sendSelf() - section " << 
	j << "failed to send itself\n";
      return -1;
    }
  }
  
  // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
  int secDefSize = 0;
  for (i = 0; i < nSections; i++) {
     int size = sections[i]->getOrder();
     secDefSize   += size;
  }
    
  Vector dData(2+NEBD+NEBD*NEBD+secDefSize); 
  loc = 0;

  // place double variables into Vector
  dData(loc++) = rho;
  dData(loc++) = tol;
  
  // place kvcommit into vector
  for (i=0; i<NEBD; i++) 
    dData(loc++) = Secommit(i);

  // place kvcommit into vector
  for (i=0; i<NEBD; i++) 
     for (j=0; j<NEBD; j++)
        dData(loc++) = kvcommit(i,j);
  
  // place vscommit into vector
  for (k=0; k<nSections; k++)
     for (i=0; i<sections[k]->getOrder(); i++)
	dData(loc++) = (vscommit[k])(i);
  
  if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
     opserr << "NLBeamColumn3d::sendSelf() - %s\n",
	 		     "failed to send Vector data";
     return -1;
  }    

  return 0;
}


int
NLBeamColumn3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  //
  // into an ID of size 9 place the integer data
  //
  int dbTag = this->getDbTag();
  int i,j,k;

  static ID idData(9);

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    opserr << "NLBeamColumn3d::recvSelf() - %s\n",
			    "failed to recv ID data";
    return -1;
  }    

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);
  maxIters = idData(4);
  initialFlag = idData(5);
  isTorsion = (idData(8) == 1) ? true : false;
  
  int crdTransfClassTag = idData(6);
  int crdTransfDbTag = idData(7);

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;
      crdTransf = theBroker.getNewCrdTransf3d(crdTransfClassTag);
      if (crdTransf == 0) {
	opserr << "NLBeamColumn3d::recvSelf() - failed to obtain a CrdTrans object with classTag" <<
	  crdTransfClassTag << endln;
	return -2;	  
      }
  }
  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf obkject
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "NLBeamColumn3d::sendSelf() - failed to recv crdTranf\n";
	     		     
     return -3;
  }      
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "NLBeamColumn3d::recvSelf() - %s\n",
			    "failed to recv ID data";
    return -1;
  }    

  //
  // now receive the sections
  //
  if (nSections != idData(3)) {

    //
    // we do not have correct number of sections, must delete the old and create
    // new ones before can recvSelf on the sections
    //

    // delete the old
    if (nSections != 0) {
      for (int i=0; i<nSections; i++)
	delete sections[i];
      delete [] sections;
    }

    // create a new array to hold pointers
    sections = new SectionForceDeformation *[idData(3)];
    if (sections == 0) {
      opserr << "NLBeamColumn3d::recvSelf() - out of memory creating sections array of size" << idData(3) << endln;
      exit(-1);
    }    

    // create all those matrices and vectors
    nSections = idData(3);
    loc = 0;

    // Delete the old
    if (vscommit != 0)
      delete [] vscommit;
  
    // Allocate the right number
    vscommit = new Vector[nSections];
    if (vscommit == 0) {
      opserr << "%s -- failed to allocate vscommit array",
			      "NLBeamColumn3d::recvSelf";
      return -1;
    }
  
    // Delete the old
    if (fs != 0)
     delete [] fs;
    
    // Allocate the right number
    fs  = new Matrix[nSections];  
    if (fs == 0) {
     opserr << "%s -- failed to allocate fs array",
			     "NLBeamColumn3d::recvSelf";
     return -1;
    } 
    
    // Delete the old
    if (vs != 0)
      delete [] vs;
    
    // Allocate the right number
    vs = new Vector[nSections];  
    if (vs == 0) {
      opserr << "%s -- failed to allocate vs array",
			      "NLBeamColumn3d::recvSelf";
      return -1;
    }
   
    // Delete the old
    if (Ssr != 0)
      delete [] Ssr;
    
    // Allocate the right number
    Ssr = new Vector[nSections];  
    if (Ssr == 0) {
      opserr << "%s -- failed to allocate Ssr array",
			      "NLBeamColumn3d::recvSelf";
      return -1;
    }
    
    for (i=0; i<nSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      sections[i] = theBroker.getNewSection(sectClassTag);
      if (sections[i] == 0) {
	opserr << "NLBeamColumn3d::recvSelf() - " << 
	  "Broker could not create Section of class type " << sectClassTag << endln;
	return -1;
      }

      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "NLBeamColumn3d::recvSelf() - section " <<
	  i << "failed to recv itself\n";
	return -1;
      }     
    }

    // Set up section history variables 
    this->initializeSectionHistoryVariables();

  } else {

    // 
    // for each existing section, check it is of correct type
    // (if not delete old & create a new one) then recvSelf on it
    //
    
    loc = 0;
    for (i=0; i<nSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;

      // check of correct type
      if (sections[i]->getClassTag() !=  sectClassTag) {
	// delete the old section[i] and create a new one
	delete sections[i];
	sections[i] = theBroker.getNewSection(sectClassTag);
	if (sections[i] == 0) {
	  opserr << "NLBeamColumn3d::recvSelf() - " <<
	    "Broker could not create Section of class type" << sectClassTag << endln;
	  exit(-1);
	}
      }

      // recvvSelf on it
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "NLBeamColumn3d::recvSelf() - section " << 
	  i << "failed to recv itself\n";
	return -1;
      }     
    }
  }
  
  // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
  int secDefSize = 0;
  for (int ii = 0; ii < nSections; ii++)
  {
     int size = sections[ii]->getOrder();
     secDefSize   += size;
  }

  Vector dData(2+NEBD+NEBD*NEBD+secDefSize);   
  
  if (theChannel.recvVector(dbTag, commitTag, dData) < 0)  {
    opserr << "NLBeamColumn3d::sendSelf() - failed to recv Vector data\n";
    return -1;
  }    
  
  loc = 0;
  
  // place double variables into Vector
  rho = dData(loc++);
  tol = dData(loc++);
  
  // place kvcommit into vector
  for (i=0; i<NEBD; i++) 
    Secommit(i) = dData(loc++);
  
  // place kvcommit into vector
  for (i=0; i<NEBD; i++) 
    for (j=0; j<NEBD; j++)
        kvcommit(i,j) = dData(loc++);
  
  kv   = kvcommit;
  Se   = Secommit;

  for (k = 0; k < nSections; k++) {
    int order = sections[k]->getOrder();
    
    // place vscommit into vector
    vscommit[k] = Vector(order);
    for (i = 0; i < order; i++)
      (vscommit[k])(i) = dData(loc++);
  }

  initialFlag = 2;  

  return 0;
}


void 
NLBeamColumn3d::compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const
{
   // update the transformation
   crdTransf->update();
       
   // get basic displacements and increments
   static Vector ub(NEBD);
   ub = crdTransf->getBasicTrialDisp();    
  
   // get integration point positions and weights
   const Matrix &xi_pt  = quadRule.getIntegrPointCoords(nSections);

   // setup Vandermode and CBDI influence matrices
   int i;
   double xi;
 
   // get CBDI influence matrix
   Matrix ls(nSections, nSections);
   double L = crdTransf->getInitialLength();
   getCBDIinfluenceMatrix(nSections, xi_pt, L, ls);
     
   // get section curvatures
   Vector kappa_y(nSections);  // curvature
   Vector kappa_z(nSections);  // curvature
   static Vector vs;                // section deformations 
        
   for (i=0; i<nSections; i++) {
       // THIS IS VERY INEFFICIENT ... CAN CHANGE IF RUNS TOO SLOW
       int sectionKey1 = 0;
       int sectionKey2 = 0;
       const ID &code = sections[i]->getType();
       int j;
       for (j = 0; j < code.Size(); j++)
       {
	   if (code(j) == SECTION_RESPONSE_MZ)
	       sectionKey1 = j;
	   if (code(j) == SECTION_RESPONSE_MY)
	       sectionKey2 = j;
       }
       if (sectionKey1 == 0) {
	 opserr << "FATAL NLBeamColumn3d::compSectionResponse - section does not provide Mz response\n";
	 exit(-1);
       }
       if (sectionKey2 == 0) {
	 opserr << "FATAL NLBeamColumn3d::compSectionResponse - section does not provide My response\n";
	 exit(-1);
       }
       
       // get section deformations
       vs = sections[i]->getSectionDeformation();
       
       kappa_z(i) = vs(sectionKey1);
       kappa_y(i) = vs(sectionKey2); 
   }

   //cout << "kappa_y: " << kappa_y;   
   //cout << "kappa_z: " << kappa_z;   
   
   Vector v(nSections), w(nSections);
   static Vector xl(NDM), uxb(NDM);
   static Vector xg(NDM), uxg(NDM); 
   // double theta;                             // angle of twist of the sections

   // v = ls * kappa_z;  
   v.addMatrixVector (0.0, ls, kappa_z, 1.0);  
   // w = ls * kappa_y *  (-1);  
   w.addMatrixVector (0.0, ls, kappa_y, -1.0);
   
   for (i=0; i<nSections; i++)
   {
      xi = xi_pt(i,0);

      xl(0) = xi * L;
      xl(1) = 0;
      xl(2) = 0;

      // get section global coordinates
      sectionCoords[i] = crdTransf->getPointGlobalCoordFromLocal(xl);

      // compute section displacements
      //theta  = xi * ub(5); // consider linear variation for angle of twist. CHANGE LATER!!!!!!!!!!
      uxb(0) = xi * ub(0); // consider linear variation for axial displacement. CHANGE LATER!!!!!!!!!!
      uxb(1) = v(i);
      uxb(2) = w(i);
          
      // get section displacements in global system 
      sectionDispls[i] = crdTransf->getPointGlobalDisplFromBasic(xi, uxb);
   }	       
}


void
NLBeamColumn3d::Print(OPS_Stream &s, int flag)
{

  // flags with negative values are used by GSA
   if (flag == -1) { 
    int eleTag = this->getTag();
    s << "NL_BEAM\t" << eleTag << "\t";
    s << sections[0]->getTag() << "\t" << sections[nSections-1]->getTag(); 
    s  << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
    s << "\t0\t0.0000000\n";
   }  

   else if (flag < -1) {
      int eleTag = this->getTag();
      int counter = (flag +1) * -1;
      int i;
      const Vector &force = this->getResistingForce();
      s << "FORCE\t" << eleTag << "\t" << counter << "\t0";
      for (i=0; i<3; i++)
	s << "\t" << force(i);
      s << endln;
      s << "FORCE\t" << eleTag << "\t" << counter << "\t1";
      for (i=0; i<3; i++)
	s << "\t" << force(i+6);
      s << endln;
      s << "MOMENT\t" << eleTag << "\t" << counter << "\t0";
      for (i=3; i<6; i++)
	s << "\t" << force(i);
      s << endln;
      s << "MOMENT\t" << eleTag << "\t" << counter << "\t1";
      for (i=3; i<6; i++)
	s << "\t" << force(i+6);
      s << endln;
   }
   
   else if (flag == 1) {
     static Vector xAxis(3);
     static Vector yAxis(3);
     static Vector zAxis(3);
     
     crdTransf->getLocalAxes(xAxis, yAxis, zAxis);
                        
     s << "#NLBeamColumn3D\n";
     s << "#LocalAxis " << xAxis(0) << " " << xAxis(1) << " " << xAxis(2) 
       << " " << zAxis(0) << " " << zAxis(1) << " " << zAxis(2) << endln;

     const Vector &node1Crd = theNodes[0]->getCrds();
     const Vector &node2Crd = theNodes[1]->getCrds();	
     const Vector &node1Disp = theNodes[0]->getDisp();
     const Vector &node2Disp = theNodes[1]->getDisp();    
     
     s << "#NODE " << node1Crd(0) << " " << node1Crd(1) << " " << node1Crd(2)
       << " " << node1Disp(0) << " " << node1Disp(1) << " " << node1Disp(2)
       << " " << node1Disp(3) << " " << node1Disp(4) << " " << node1Disp(5) << endln;
     
     s << "#NODE " << node2Crd(0) << " " << node2Crd(1) << " " << node2Crd(2)
       << " " << node2Disp(0) << " " << node2Disp(1) << " " << node2Disp(2)
       << " " << node2Disp(3) << " " << node2Disp(4) << " " << node2Disp(5) << endln;

     // allocate array of vectors to store section coordinates and displacements
     static int maxNumSections = 0;
     static Vector *coords = 0;
     static Vector *displs = 0;
     if (maxNumSections < nSections) {
       if (coords != 0) 
	 delete [] coords;
       if (displs != 0)
	 delete [] displs;
       
       coords = new Vector [nSections];
       displs = new Vector [nSections];
       
       if (!coords) {
	 opserr << "NLBeamColumn3d::Print() -- failed to allocate coords array";   
	 exit(-1);
       }
       
       int i;
       for (i = 0; i < nSections; i++)
	 coords[i] = Vector(NDM);
       
       if (!displs) {
	 opserr << "NLBeamColumn3d::Print() -- failed to allocate coords array";   
	 exit(-1);
       }
       
       for (i = 0; i < nSections; i++)
	 displs[i] = Vector(NDM);
       
       maxNumSections = nSections;
     }
     
     // compute section location & displacements
     this->compSectionDisplacements(coords, displs);
     
     // spit out the section location & invoke print on the scetion
     for (int i=0; i<nSections; i++) {
       s << "#SECTION " << (coords[i])(0) << " " << (coords[i])(1) << " " << (coords[i])(2);       s << " " << (displs[i])(0) << " " << (displs[i])(1) << " " << (displs[i])(2) << endln;
       sections[i]->Print(s, flag); 
     }
   }
   
   else {
      s << "\nElement: " << this->getTag() << " Type: NLBeamColumn3d ";
      s << "\tConnected Nodes: " << connectedExternalNodes ;
      s << "\tNumber of Sections: " << nSections << endln;
      s << "\tElement End Forces (P MZ1 MZ2 MY1 MY2 T): " << Secommit;
      s << "\tResisting Force: " << this->getResistingForce();
   }
}


OPS_Stream &operator<<(OPS_Stream &s, NLBeamColumn3d &E)
{
    E.Print(s);
    return s;
}

int
NLBeamColumn3d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
   
  if (displayMode == 1) {

    // first determine the two end points of the element based on
    //  the display factor (a measure of the distorted image)
    static Vector v1(NDM), v2(NDM);
    
    const Vector &node1Crd = theNodes[0]->getCrds();
    const Vector &node2Crd = theNodes[1]->getCrds();	
    const Vector &node1Disp = theNodes[0]->getDisp();
    const Vector &node2Disp = theNodes[1]->getDisp();    
    
    int i;
    
    // allocate array of vectors to store section coordinates and displacements
    static int maxNumSections = 0;
    static Vector *coords = 0;
    static Vector *displs = 0;
    if (maxNumSections < nSections) {
      if (coords != 0) 
	delete [] coords;
      if (displs != 0)
	delete [] displs;
      
      coords = new Vector [nSections];
      displs = new Vector [nSections];
      
      if (!coords) {
	opserr << "NLBeamColumn3d::displaySelf() -- failed to allocate coords array";   
	exit(-1);
      }
      
      for (i = 0; i < nSections; i++)
	coords[i] = Vector(NDM);
      
      if (!displs) {
	opserr << "NLBeamColumn3d::displaySelf() -- failed to allocate coords array";   
	exit(-1);
      }
      
      for (i = 0; i < nSections; i++)
	displs[i] = Vector(NDM);
      
      maxNumSections = nSections;
    }
       
    // compute section location & displacements
    int error;
    this->compSectionDisplacements(coords, displs);

    // get global displacements and coordinates of each section          
    
    v1 = node1Crd + node1Disp*fact;
    
    // get global displacements and coordinates of each section          
    
    for (i = 0; i<nSections; i++) {
      v2 = coords[i] + displs[i]*fact;
      
	 error = theViewer.drawLine(v1, v2, 1.0, 1.0);
	 
	 if (error)
	   return error;
	 v1 = v2;
	 
    }  
    
    v2 = node2Crd + node2Disp*fact;
    
    error = theViewer.drawLine(v1, v2, 1.0, 1.0);
    
    if (error)
      return error;
  }
  return 0;
}

Response* 
NLBeamColumn3d::setResponse(const char **argv, int argc, Information &eleInformation)
{
    //
    // we compare argv[0] for known response types 
    //

    // global force - 
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
		return new ElementResponse(this, 1, theVector);

    // local force -
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
		return new ElementResponse(this, 2, theVector);

    else if ((strcmp(argv[0],"defoANDforce") == 0) ||
	     (strcmp(argv[0],"deformationANDforce") == 0) ||
	     (strcmp(argv[0],"deformationsANDforces") == 0))
      return new ElementResponse(this, 4, Vector(24));

    // section response -
    else if (strcmp(argv[0],"section") ==0) {
      if (argc <= 2)
	return 0;
	
      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= nSections)
	return sections[sectionNum-1]->setResponse(&argv[2], argc-2, eleInformation);
      else
	return 0;
    }
    
    else
      return 0;
}

int 
NLBeamColumn3d::getResponse(int responseID, Information &eleInfo)
{
  static Vector force(12);
  static Vector defoAndForce(24);
  double V, N, T, M1, M2;
  int i;
  double L = crdTransf->getInitialLength();
  
  switch (responseID) {      
  case 1:  // forces
    return eleInfo.setVector(this->getResistingForce());
    
    
  case 4:
    this->getGlobalDispls(force);
    this->getResistingForce();
    for (i = 0; i < 12; i++) {
      defoAndForce(i) = force(i);
      defoAndForce(i+12) = theVector(i);
    }
    return eleInfo.setVector(defoAndForce);
    
  case 2:
    // Axial
    N = Se(0);
    force(6) =  N;
    force(0) = -N+p0[0];
    
    // Torsion
    T = Se(5);
    force(9) =  T;
    force(3) = -T;
    
    // Moments about z and shears along y
    M1 = Se(1);
    M2 = Se(2);
    force(5)  = M1;
    force(11) = M2;
    V = (M1+M2)/L;
    force(1) =  V+p0[1];
    force(7) = -V+p0[2];
    
    // Moments about y and shears along z
    M1 = Se(3);
    M2 = Se(4);
    force(4)  = M1;
    force(10) = M2;
    V = -(M1+M2)/L;
    force(2) = -V+p0[3];
    force(8) =  V+p0[4];
      
    return eleInfo.setVector(force);
    
  default: 
    return -1;
  }
}

int
NLBeamColumn3d::setParameter (const char **argv, int argc, Information &info)
{
	int ok = -1;

	// If the parameter belongs to the element itself
	if (strcmp(argv[0],"rho") == 0) {
		info.theType = DoubleType;
		return 1;
	}

	// If the parameter is belonging to a section or lower
	if (strcmp(argv[0],"section") == 0) {

	  // For now, no parameters of the section itself:
	  if (argc<5) {
	    opserr << "For now: cannot handle parameters of the section itself." << endln;
	    return -1;
	  }

	  // Reveal section and material tag numbers
	  int paramSectionTag = atoi(argv[1]);
	  int paramMatTag     = atoi(argv[3]);

	  // Store section and material tag in theInfo
	  ID *theID = new ID(2);
	  (*theID)(0) = paramSectionTag;
	  (*theID)(1) = paramMatTag;
	  info.theID = theID;
	  
	  // Find the right section and call its setParameter method
	  for (int i=0; i<nSections; i++) {
	    if (paramSectionTag == sections[i]->getTag()) {
	      ok = sections[i]->setParameter(&argv[2], argc-2, info);
	    }
	  }
	  
	  if (ok < 0) {
	    opserr << "NLBeamColumn2d::setParameter() - could not set parameter. " << endln;
	    return -1;
	  }
	  else {
	    return ok + 100;
	  }
	}
	
	// otherwise parameter is unknown for the NLBeamColumn2d class
	else
	  return -1;
}

int
NLBeamColumn3d::updateParameter (int parameterID, Information &info)
{
	ID *paramIDPtr;
	int ok = -1;

	switch (parameterID) {
	case 1:
		this->rho = info.theDouble;
		return 0;
	default:
		if (parameterID >= 100) {
			paramIDPtr = info.theID;
			ID paramID = (*paramIDPtr);
			int paramSectionTag = paramID(0);
			for (int i=0; i<nSections; i++) {
				if (paramSectionTag == sections[i]->getTag()) {
					ok = sections[i]->updateParameter(parameterID-100, info);
				}
			}
			if (ok < 0) {
				opserr << "NLBeamColumn2d::updateParameter() - could not update parameter. " << endln;
				return ok;
			}
			else {
				return ok;
			}
		}
		else
			return -1;
	}       
}

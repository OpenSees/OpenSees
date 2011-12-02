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
                                                                        
// $Revision: 1.27 $
// $Date: 2007-06-09 17:16:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ForceBeamColumn3d.cpp,v $

/*
 * References
 *

State Determination Algorithm
---
Neuenhofer, A. and F. C. Filippou (1997). "Evaluation of Nonlinear Frame Finite
Element Models." Journal of Structural Engineering, 123(7):958-966.

Spacone, E., V. Ciampi, and F. C. Filippou (1996). "Mixed Formulation of
Nonlinear Beam Finite Element." Computers and Structures, 58(1):71-83.


Plastic Hinge Integration
---
Scott, M. H. and G. L. Fenves (2006). "Plastic Hinge Integration Methods for
Force-Based Beam-Column Elements." Journal of Structural Engineering,
132(2):244-252.


Analytical Response Sensitivity (DDM)
---
Scott, M. H., P. Franchin, G. L. Fenves, and F. C. Filippou (2004).
"Response Sensitivity for Nonlinear Beam-Column Elements."
Journal of Structural Engineering, 130(9):1281-1288.


Software Design
---
Scott, M. H., G. L. Fenves, F. T. McKenna, and F. C. Filippou (2007).
"Software Patterns for Nonlinear Beam-Column Models."
Journal of Structural Engineering, Approved for publication, February 2007.

 *
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <Information.h>
#include <Parameter.h>
#include <ForceBeamColumn3d.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>

#include <ElementResponse.h>
#include <ElementalLoad.h>

#define  NDM   3         // dimension of the problem (3d)
#define  NND   6         // number of nodal dof's
#define  NEGD  12        // number of element global dof's
#define  NEBD  6         // number of element dof's in the basic system

#define DefaultLoverGJ 1.0e-10

Matrix ForceBeamColumn3d::theMatrix(12,12);
Vector ForceBeamColumn3d::theVector(12);
double ForceBeamColumn3d::workArea[200];

Vector *ForceBeamColumn3d::vsSubdivide = 0;
Matrix *ForceBeamColumn3d::fsSubdivide = 0;
Vector *ForceBeamColumn3d::SsrSubdivide = 0;

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
ForceBeamColumn3d::ForceBeamColumn3d(): 
  Element(0,ELE_TAG_ForceBeamColumn3d), connectedExternalNodes(2), 
  beamIntegr(0), numSections(0), sections(0), crdTransf(0),
  rho(0.0), maxIters(0), tol(0.0),
  initialFlag(0),
  kv(NEBD,NEBD), Se(NEBD),
  kvcommit(NEBD,NEBD), Secommit(NEBD),
  fs(0), vs(0), Ssr(0), vscommit(0), sp(0), Ki(0), isTorsion(false)
{
  theNodes[0] = 0;  
  theNodes[1] = 0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;
  
  v0[0] = 0.0;
  v0[1] = 0.0;
  v0[2] = 0.0;
  v0[3] = 0.0;
  v0[4] = 0.0;

  if (vsSubdivide == 0)
    vsSubdivide  = new Vector [maxNumSections];
  if (fsSubdivide == 0) 
    fsSubdivide  = new Matrix [maxNumSections];
  if (SsrSubdivide == 0)
    SsrSubdivide  = new Vector [maxNumSections];
  if (!vsSubdivide || !fsSubdivide || !SsrSubdivide) {
    opserr << "ForceBeamColumn3d::ForceBeamColumn3d() -- failed to allocate Subdivide arrays";   
    exit(-1);
  }
}

// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points. 
// allocates the necessary space needed by each object
ForceBeamColumn3d::ForceBeamColumn3d (int tag, int nodeI, int nodeJ,
				      int numSec, SectionForceDeformation **sec,
				      BeamIntegration &bi,
				      CrdTransf3d &coordTransf, double massDensPerUnitLength,
				      int maxNumIters, double tolerance):
  Element(tag,ELE_TAG_ForceBeamColumn3d), connectedExternalNodes(2),
  beamIntegr(0), numSections(0), sections(0), crdTransf(0),
  rho(massDensPerUnitLength),maxIters(maxNumIters), tol(tolerance), 
  initialFlag(0),
  kv(NEBD,NEBD), Se(NEBD), 
  kvcommit(NEBD,NEBD), Secommit(NEBD),
  fs(0), vs(0),Ssr(0), vscommit(0), sp(0), Ki(0), isTorsion(false)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;    

  beamIntegr = bi.getCopy();
  if (beamIntegr == 0) {
    opserr << "Error: ForceBeamColumn3d::ForceBeamColumn3d: could not create copy of beam integration object" << endln;
    exit(-1);
  }

  // get copy of the transformation object   
  crdTransf = coordTransf.getCopy(); 
  if (crdTransf == 0) {
    opserr << "Error: ForceBeamColumn3d::ForceBeamColumn3d: could not create copy of coordinate transformation object" << endln;
    exit(-1);
  }

  
  this->setSectionPointers(numSec, sec);

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;
  
  v0[0] = 0.0;
  v0[1] = 0.0;
  v0[2] = 0.0;
  v0[3] = 0.0;
  v0[4] = 0.0;

  if (vsSubdivide == 0)
    vsSubdivide  = new Vector [maxNumSections];
  if (fsSubdivide == 0)
    fsSubdivide  = new Matrix [maxNumSections];
  if (SsrSubdivide == 0)
    SsrSubdivide  = new Vector [maxNumSections];
  if (!vsSubdivide || !fsSubdivide || !SsrSubdivide) {
    opserr << "ForceBeamColumn3d::ForceBeamColumn3d() -- failed to allocate Subdivide arrays";   
    exit(-1);
  }
}

// ~ForceBeamColumn3d():
// 	destructor
//      delete must be invoked on any objects created by the object
ForceBeamColumn3d::~ForceBeamColumn3d()
{
  if (sections != 0) {
    for (int i=0; i < numSections; i++)
      if (sections[i] != 0)
	delete sections[i];
    delete [] sections;
  }
  
  if (fs != 0) 
    delete [] fs;
  
  if (vs != 0) 
    delete [] vs;
  
  if (Ssr != 0) 
    delete [] Ssr;
  
  if (vscommit != 0) 
    delete [] vscommit;
  
  if (crdTransf != 0)
    delete crdTransf;

  if (beamIntegr != 0)
    delete beamIntegr;
  
  if (sp != 0)
    delete sp;

  if (Ki != 0)
    delete Ki;
}

int
ForceBeamColumn3d::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
ForceBeamColumn3d::getExternalNodes(void) 
{
  return connectedExternalNodes;
}

Node **
ForceBeamColumn3d::getNodePtrs()
{
  return theNodes;
}

int
ForceBeamColumn3d::getNumDOF(void) 
{
  return NEGD;
}

void
ForceBeamColumn3d::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    opserr << "ForceBeamColumn3d::setDomain:  theDomain = 0 ";
    exit(0); 
  }

  // get pointers to the nodes
  
  int Nd1 = connectedExternalNodes(0);  
  int Nd2 = connectedExternalNodes(1);
  
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);  
  
  if (theNodes[0] == 0) {
    opserr << "ForceBeamColumn3d::setDomain: Nd1: ";
    opserr << Nd1 << "does not exist in model\n";
    exit(0);
  }
  
  if (theNodes[1] == 0) {
    opserr << "ForceBeamColumn3d::setDomain: Nd2: ";
    opserr << Nd2 << "does not exist in model\n";
    exit(0);
  }
  
  // call the DomainComponent class method 
  this->DomainComponent::setDomain(theDomain);
  
  // ensure connected nodes have correct number of dof's
  int dofNode1 = theNodes[0]->getNumberDOF();
  int dofNode2 = theNodes[1]->getNumberDOF();
  
  if ((dofNode1 != NND) || (dofNode2 != NND)) {
    opserr << "ForceBeamColumn3d::setDomain(): Nd2 or Nd1 incorrect dof ";
    exit(0);
  }
   
  // initialize the transformation
  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    opserr << "ForceBeamColumn3d::setDomain(): Error initializing coordinate transformation";  
    exit(0);
  }
    
  // get element length
  double L = crdTransf->getInitialLength();
  if (L == 0.0) {
    opserr << "ForceBeamColumn3d::setDomain(): Zero element length:" << this->getTag();  
    exit(0);
  }

  if (initialFlag == 0) 
    this->initializeSectionHistoryVariables();
}

int
ForceBeamColumn3d::commitState()
{
  int err = 0;
  int i = 0;

  // call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0) {
    opserr << "ForceBeamColumn3d::commitState () - failed in base class";
  }    
  
  do {
    vscommit[i] = vs[i];
    err = sections[i++]->commitState();
    
  } while (err == 0 && i < numSections);
  
  if (err)
    return err;
  
  // commit the transformation between coord. systems
  if ((err = crdTransf->commitState()) != 0)
    return err;
  
  // commit the element variables state
  kvcommit = kv;
  Secommit = Se;
  
  //   initialFlag = 0;  fmk - commented out, see what happens to Example3.1.tcl if uncommented
  //                         - i have not a clue why, ask remo if he ever gets in contact with us again!
  
  return err;
}

int ForceBeamColumn3d::revertToLastCommit()
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
  } while (err == 0 && i < numSections);
  
  
  if (err)
    return err;
  
  // revert the transformation to last commit
  if ((err = crdTransf->revertToLastCommit()) != 0)
    return err;
  
  // revert the element state to last commit
  Se   = Secommit;
  kv   = kvcommit;
  
  initialFlag = 0;
  // this->update();
  
  return err;
}

int ForceBeamColumn3d::revertToStart()
{
  // revert the sections state to start
  int err;
  int i = 0;
  
  do {
    fs[i].Zero();
    vs[i].Zero();
    Ssr[i].Zero();
    err = sections[i++]->revertToStart();
    
  } while (err == 0 && i < numSections);
  
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
ForceBeamColumn3d::getInitialStiff(void)
{
  // check for quick return
  if (Ki != 0)
    return *Ki;

  static Matrix f(NEBD,NEBD);   // element flexibility matrix  
  this->getInitialFlexibility(f);
  
  static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse  
  I.Zero();
  for (int i=0; i<NEBD; i++)
    I(i,i) = 1.0;
  
  // calculate element stiffness matrix
  // invert3by3Matrix(f, kv);
  static Matrix kvInit(NEBD, NEBD);
  if (f.Solve(I, kvInit) < 0)
    opserr << "ForceBeamColumn3d::getInitialStiff() -- could not invert flexibility";

    Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kvInit));

    return *Ki;
  }

  const Matrix &
  ForceBeamColumn3d::getTangentStiff(void)
  {
    crdTransf->update();	// Will remove once we clean up the corotational 3d transformation -- MHS
    return crdTransf->getGlobalStiffMatrix(kv, Se);
  }

  const Vector &
  ForceBeamColumn3d::getResistingForce(void)
  {
    // Will remove once we clean up the corotational 3d transformation -- MHS
    crdTransf->update();

    Vector p0Vec(p0, 5);

    return crdTransf->getGlobalResistingForce(Se, p0Vec);
  }

  void
  ForceBeamColumn3d::initializeSectionHistoryVariables (void)
  {
    for (int i = 0; i < numSections; i++) {
      int order = sections[i]->getOrder();

      fs[i] = Matrix(order,order);
      vs[i] = Vector(order);
      Ssr[i] = Vector(order);

      vscommit[i] = Vector(order);
    }
  }

  /********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS ********************
   */
  int
  ForceBeamColumn3d::update()
  {
    // if have completed a recvSelf() - do a revertToLastCommit
    // to get Ssr, etc. set correctly
    if (initialFlag == 2)
      this->revertToLastCommit();

    // update the transformation
    crdTransf->update();

    // get basic displacements and increments
    const Vector &v = crdTransf->getBasicTrialDisp();    

    static Vector dv(NEBD);
    dv = crdTransf->getBasicIncrDeltaDisp();    

    if (initialFlag != 0 && dv.Norm() <= DBL_EPSILON && sp == 0)
      return 0;

    static Vector vin(NEBD);
    vin = v;
    vin -= dv;
    double L = crdTransf->getInitialLength();
    double oneOverL  = 1.0/L;  

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);

    double wt[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wt);

    static Vector vr(NEBD);       // element residual displacements
    static Matrix f(NEBD,NEBD);   // element flexibility matrix

    static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse
    double dW;                    // section strain energy (work) norm 
    int i, j;

    I.Zero();
    for (i=0; i<NEBD; i++)
      I(i,i) = 1.0;

    int numSubdivide = 1;
    bool converged = false;
    static Vector dSe(NEBD);
    static Vector dvToDo(NEBD);
    static Vector dvTrial(NEBD);
    static Vector SeTrial(NEBD);
    static Matrix kvTrial(NEBD, NEBD);

    dvToDo = dv;
    dvTrial = dvToDo;

    static double factor = 10;
    double dW0 = 0.0;

    maxSubdivisions = 10;

    // fmk - modification to get compatable ele forces and deformations 
    //   for a change in deformation dV we try first a newton iteration, if
    //   that fails we try an initial flexibility iteration on first iteration 
    //   and then regular newton, if that fails we use the initial flexiblity
    //   for all iterations.
    //
    //   if they both fail we subdivide dV & try to get compatable forces
    //   and deformations. if they work and we have subdivided we apply
    //   the remaining dV.

    while (converged == false && numSubdivide <= maxSubdivisions) {

      // try regular newton (if l==0), or
      // initial tangent iterations (if l==1), or
      // initial tangent on first iteration then regular newton (if l==2)

      for (int l=0; l<3; l++) {

	//      if (l == 1) l = 2;
	SeTrial = Se;
	kvTrial = kv;
	for (i=0; i<numSections; i++) {
	  vsSubdivide[i] = vs[i];
	  fsSubdivide[i] = fs[i];
	  SsrSubdivide[i] = Ssr[i];
	}

	// calculate nodal force increments and update nodal forces      
	// dSe = kv * dv;
	dSe.addMatrixVector(0.0, kvTrial, dvTrial, 1.0);
	SeTrial += dSe;

	if (initialFlag != 2) {

	  int numIters = maxIters;
	  if (l == 1) 
	    numIters = 10*maxIters; // allow 10 times more iterations for initial tangent

	  for (j=0; j <numIters; j++) {

	    // initialize f and vr for integration
	    f.Zero();
	    vr.Zero();

	    if (beamIntegr->addElasticFlexibility(L, f) < 0) {
	      vr(0) += f(0,0)*SeTrial(0);
	      vr(1) += f(1,1)*SeTrial(1) + f(1,2)*SeTrial(2);
	      vr(2) += f(2,1)*SeTrial(1) + f(2,2)*SeTrial(2);
	      vr(3) += f(3,3)*SeTrial(3) + f(3,4)*SeTrial(4);
	      vr(4) += f(4,3)*SeTrial(3) + f(4,4)*SeTrial(4);
	      vr(5) += f(5,5)*SeTrial(5);
	    }

	    // Add effects of element loads
	    vr(0) += v0[0];
	    vr(1) += v0[1];
	    vr(2) += v0[2];
	    vr(3) += v0[3];
	    vr(4) += v0[4];

	    for (i=0; i<numSections; i++) {

	      int order      = sections[i]->getOrder();
	      const ID &code = sections[i]->getType();

	      static Vector Ss;
	      static Vector dSs;
	      static Vector dvs;
	      static Matrix fb;

	      Ss.setData(workArea, order);
	      dSs.setData(&workArea[order], order);
	      dvs.setData(&workArea[2*order], order);
	      fb.setData(&workArea[3*order], order, NEBD);

	      double xL  = xi[i];
	      double xL1 = xL-1.0;
	      double wtL = wt[i]*L;

	      // calculate total section forces
	      // Ss = b*Se + bp*currDistrLoad;
	      // Ss.addMatrixVector(0.0, b[i], Se, 1.0);
	      int ii;
	      for (ii = 0; ii < order; ii++) {
		switch(code(ii)) {
		case SECTION_RESPONSE_P:
		  Ss(ii) = SeTrial(0);
		  break;
		case SECTION_RESPONSE_MZ:
		  Ss(ii) = xL1*SeTrial(1) + xL*SeTrial(2);
		  break;
		case SECTION_RESPONSE_VY:
		  Ss(ii) = oneOverL*(SeTrial(1)+SeTrial(2));
		  break;
		case SECTION_RESPONSE_MY:
		  Ss(ii) = xL1*SeTrial(3) + xL*SeTrial(4);
		  break;
		case SECTION_RESPONSE_VZ:
		  Ss(ii) = oneOverL*(SeTrial(3)+SeTrial(4));
		  break;
		case SECTION_RESPONSE_T:
		  Ss(ii) = SeTrial(5);
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

	      // dSs = Ss - Ssr[i];
	      dSs = Ss;
	      dSs.addVector(1.0, SsrSubdivide[i], -1.0);

	      // compute section deformation increments
	      if (l == 0) {

		//  regular newton 
		//    vs += fs * dSs;     

		dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);

	      } else if (l == 2) {

		//  newton with initial tangent if first iteration
		//    vs += fs0 * dSs;     
		//  otherwise regular newton 
		//    vs += fs * dSs;     

		if (j == 0) {
		  const Matrix &fs0 = sections[i]->getInitialFlexibility();

		  dvs.addMatrixVector(0.0, fs0, dSs, 1.0);
		} else
		  dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);

	      } else {

		//  newton with initial tangent
		//    vs += fs0 * dSs;     

		const Matrix &fs0 = sections[i]->getInitialFlexibility();
		dvs.addMatrixVector(0.0, fs0, dSs, 1.0);
	      }

	      // set section deformations
	      if (initialFlag != 0)
		vsSubdivide[i] += dvs;

	      sections[i]->setTrialSectionDeformation(vsSubdivide[i]);

	      // get section resisting forces
	      SsrSubdivide[i] = sections[i]->getStressResultant();

	      // get section flexibility matrix
	      // FRANK 
	      fsSubdivide[i] = sections[i]->getSectionFlexibility();

	      /*
	      const Matrix &sectionStiff = sections[i]->getSectionTangent();
	      int n = sectionStiff.noRows();
	      Matrix I(n,n); I.Zero(); for (int l=0; l<n; l++) I(l,l) = 1.0;
	      Matrix sectionFlex(n,n);
	      sectionStiff.SolveSVD(I, sectionFlex, 1.0e-6);
	      fsSubdivide[i] = sectionFlex;	    
	      */

	      // calculate section residual deformations
	      // dvs = fs * (Ss - Ssr);
	      dSs = Ss;
	      dSs.addVector(1.0, SsrSubdivide[i], -1.0);  // dSs = Ss - Ssr[i];

	      dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);

	      // integrate element flexibility matrix
	      // f = f + (b^ fs * b) * wtL;
	      //f.addMatrixTripleProduct(1.0, b[i], fs[i], wtL);
	      int jj;
	      const Matrix &fSec = fsSubdivide[i];
	      fb.Zero();
	      double tmp;
	      for (ii = 0; ii < order; ii++) {
		switch(code(ii)) {
		case SECTION_RESPONSE_P:
		  for (jj = 0; jj < order; jj++)
		    fb(jj,0) += fSec(jj,ii)*wtL;
		  break;
		case SECTION_RESPONSE_MZ:
		  for (jj = 0; jj < order; jj++) {
		    tmp = fSec(jj,ii)*wtL;
		    fb(jj,1) += xL1*tmp;
		    fb(jj,2) += xL*tmp;
		  }
		  break;
		case SECTION_RESPONSE_VY:
		  for (jj = 0; jj < order; jj++) {
		    tmp = oneOverL*fSec(jj,ii)*wtL;
		    fb(jj,1) += tmp;
		    fb(jj,2) += tmp;
		  }
		  break;
		case SECTION_RESPONSE_MY:
		  for (jj = 0; jj < order; jj++) {
		    tmp = fSec(jj,ii)*wtL;
		    fb(jj,3) += xL1*tmp;
		    fb(jj,4) += xL*tmp;
		  }
		  break;
		case SECTION_RESPONSE_VZ:
		  for (jj = 0; jj < order; jj++) {
		    tmp = oneOverL*fSec(jj,ii)*wtL;
		    fb(jj,3) += tmp;
		    fb(jj,4) += tmp;
		  }
		  break;
		case SECTION_RESPONSE_T:
		  for (jj = 0; jj < order; jj++)
		    fb(jj,5) += fSec(jj,ii)*wtL;
		  break;
		default:
		  break;
		}
	      }

	      for (ii = 0; ii < order; ii++) {
		switch (code(ii)) {
		case SECTION_RESPONSE_P:
		  for (jj = 0; jj < NEBD; jj++)
		    f(0,jj) += fb(ii,jj);
		  break;
		case SECTION_RESPONSE_MZ:
		  for (jj = 0; jj < NEBD; jj++) {
		    tmp = fb(ii,jj);
		    f(1,jj) += xL1*tmp;
		    f(2,jj) += xL*tmp;
		  }
		  break;
		case SECTION_RESPONSE_VY:
		  for (jj = 0; jj < NEBD; jj++) {
		    tmp = oneOverL*fb(ii,jj);
		    f(1,jj) += tmp;
		    f(2,jj) += tmp;
		  }
		  break;
		case SECTION_RESPONSE_MY:
		  for (jj = 0; jj < NEBD; jj++) {
		    tmp = fb(ii,jj);
		    f(3,jj) += xL1*tmp;
		    f(4,jj) += xL*tmp;
		  }
		  break;
		case SECTION_RESPONSE_VZ:
		  for (jj = 0; jj < NEBD; jj++) {
		    tmp = oneOverL*fb(ii,jj);
		    f(3,jj) += tmp;
		    f(4,jj) += tmp;
		  }
		  break;
		case SECTION_RESPONSE_T:
		  for (jj = 0; jj < NEBD; jj++)
		    f(5,jj) += fb(ii,jj);
		  break;
		default:
		  break;
		}
	      }

	      // integrate residual deformations
	      // vr += (b^ (vs + dvs)) * wtL;
	      //vr.addMatrixTransposeVector(1.0, b[i], vs[i] + dvs, wtL);
	      dvs.addVector(1.0, vsSubdivide[i], 1.0);
	      double dei;
	      for (ii = 0; ii < order; ii++) {
		dei = dvs(ii)*wtL;
		switch(code(ii)) {
		case SECTION_RESPONSE_P:
		  vr(0) += dei;
		  break;
		case SECTION_RESPONSE_MZ:
		  vr(1) += xL1*dei; vr(2) += xL*dei;
		  break;
		case SECTION_RESPONSE_VY:
		  tmp = oneOverL*dei;
		  vr(1) += tmp; vr(2) += tmp;
		  break;
		case SECTION_RESPONSE_MY:
		  vr(3) += xL1*dei; vr(4) += xL*dei;
		  break;
		case SECTION_RESPONSE_VZ:
		  tmp = oneOverL*dei;
		  vr(3) += tmp; vr(4) += tmp;
		  break;
		case SECTION_RESPONSE_T:
		  vr(5) += dei;
		  break;
		default:
		  break;
		}
	      }
	    }

	    if (!isTorsion) {
	      f(5,5) = DefaultLoverGJ;
	      vr(5) = SeTrial(5)*DefaultLoverGJ;
	    }

	    // calculate element stiffness matrix
	    // invert3by3Matrix(f, kv);	  
	    // FRANK
	    //	  if (f.SolveSVD(I, kvTrial, 1.0e-12) < 0)
	    if (f.Solve(I, kvTrial) < 0)
	      opserr << "ForceBeamColumn3d::update() -- could not invert flexibility\n";
	    

	    // dv = vin + dvTrial  - vr
	    dv = vin;
	    dv += dvTrial;
	    dv -= vr;

	    // dv.addVector(1.0, vr, -1.0);

	    // dSe = kv * dv;
	    dSe.addMatrixVector(0.0, kvTrial, dv, 1.0);

	    dW = dv ^ dSe; 
	    if (dW0 == 0.0) 
	      dW0 = dW;

	    SeTrial += dSe;

	    // check for convergence of this interval
	    if (fabs(dW) < tol) { 

	      // set the target displacement
	      dvToDo -= dvTrial;
	      vin += dvTrial;

	      // check if we have got to where we wanted
	      if (dvToDo.Norm() <= DBL_EPSILON) {
		converged = true;

	      } else {  // we convreged but we have more to do

		// reset variables for start of next subdivision
		dvTrial = dvToDo;
		numSubdivide = 1;  // NOTE setting subdivide to 1 again maybe too much
	      }

	      // set kv, vs and Se values
	      kv = kvTrial;
	      Se = SeTrial;

	      for (int k=0; k<numSections; k++) {
		vs[k] = vsSubdivide[k];
		fs[k] = fsSubdivide[k];
		Ssr[k] = SsrSubdivide[k];
	      }

	      // break out of j & l loops
	      j = numIters+1;
	      l = 4;

	    } else {   //  if (fabs(dW) < tol) { 

	      // if we have failed to convrege for all of our newton schemes
	      // - reduce step size by the factor specified
	      if (j == (numIters-1) && (l == 2)) {
		dvTrial /= factor;
		numSubdivide++;
	      }
	    }

	  } // for (j=0; j<numIters; j++)
	} // if (initialFlag != 2)
      } // for (int l=0; l<2; l++)
    } // while (converged == false)

    // if fail to converge we return an error flag & print an error message

    if (converged == false) {
      opserr << "WARNING - ForceBeamColumn3d::update - failed to get compatible ";
      opserr << "element forces & deformations for element: ";
      opserr << this->getTag() << "(dW: << " << dW << ", dW0: " << dW0 << ")\n";

      /*
      opserr << "Section Tangent Condition Numbers: ";
      for (int i=0; i<numSections; i++) {
	const Matrix &sectionStiff = sections[i]->getSectionTangent();
	double conditionNumber = sectionStiff.conditionNumber();
	opserr << conditionNumber << " ";
      }
      opserr << endln;
      */

      return -1;
    }

    initialFlag = 1;

    return 0;
  }

  void ForceBeamColumn3d::getForceInterpolatMatrix(double xi, Matrix &b, const ID &code)
  {
    b.Zero();

    double L = crdTransf->getInitialLength();
    for (int i = 0; i < code.Size(); i++) {
      switch (code(i)) {
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
      case SECTION_RESPONSE_MY:              // Moment, My, interpolation
	b(i,3) = xi - 1.0;
	b(i,4) = xi;
	break;
      case SECTION_RESPONSE_VZ:              // Shear, Vz, interpolation
	b(i,3) = b(i,4) = 1.0/L;
	break;
      case SECTION_RESPONSE_T:               // Torque, T, interpolation
	b(i,5) = 1.0;
	break;
      default:
	break;
      }
    }
  }

  void ForceBeamColumn3d::getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code)
  {
    bp.Zero();

    double L = crdTransf->getInitialLength();
    for (int i = 0; i < code.Size(); i++) {
      switch (code(i)) {
      case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
	bp(i,1) = xi*(xi-1)*L*L/2;
	break;
      case SECTION_RESPONSE_P:		// Axial, P, interpolation
	bp(i,0) = (1-xi)*L;
	break;
      case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
	bp(i,1) = (xi-0.5)*L;
	break;
      case SECTION_RESPONSE_MY:              // Moment, My, interpolation
	bp(i,2) = xi*(1-xi)*L*L/2;
	break;
      case SECTION_RESPONSE_VZ:              // Shear, Vz, interpolation
	bp(i,2) = (0.5-xi)*L;
	break;
      case SECTION_RESPONSE_T:               // Torsion, T, interpolation
	break;
      default:
	break;
      }
    }
  }

  const Matrix &
  ForceBeamColumn3d::getMass(void)
  { 
    theMatrix.Zero();

    double L = crdTransf->getInitialLength();
    if (rho != 0.0)
      theMatrix(0,0) = theMatrix(1,1) = theMatrix(2,2) =
	theMatrix(6,6) = theMatrix(7,7) = theMatrix(8,8) = 0.5*L*rho;

    return theMatrix;
  }

  void 
  ForceBeamColumn3d::zeroLoad(void)
  {
    if (sp != 0)
      sp->Zero();

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    p0[3] = 0.0;
    p0[4] = 0.0;

    v0[0] = 0.0;
    v0[1] = 0.0;
    v0[2] = 0.0;
    v0[3] = 0.0;
    v0[4] = 0.0;
  }

  int
  ForceBeamColumn3d::addLoad(ElementalLoad *theLoad, double loadFactor)
  {
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);

    if (sp == 0) {
      sp = new Matrix(5,numSections);
      if (sp == 0) {
	opserr << "ForceBeamColumn3d::addLoad -- out of memory\n";
	exit(-1);
      }
    }

    double L = crdTransf->getInitialLength();

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);

    // Accumulate elastic deformations in basic system
    beamIntegr->addElasticDeformations(theLoad, loadFactor, L, v0);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0)*loadFactor;  // Transverse
      double wz = data(1)*loadFactor;  // Transverse
      double wx = data(2)*loadFactor;  // Axial

      Matrix &s_p = *sp;

      // Accumulate applied section forces due to element loads
      for (int i = 0; i < numSections; i++) {
	double x = xi[i]*L;
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

      if (aOverL < 0.0 || aOverL > 1.0)
	return 0;

      double a = aOverL*L;

      double Vy2 = Py*aOverL;
      double Vy1 = Py-Vy2;

      double Vz2 = Pz*aOverL;
      double Vz1 = Pz-Vz2;

      Matrix &s_p = *sp;

      // Accumulate applied section forces due to element loads
      for (int i = 0; i < numSections; i++) {
	double x = xi[i]*L;
	if (x <= a) {
	  s_p(0,i) += N;
	  s_p(1,i) -= x*Vy1;
	  s_p(2,i) -= Vy1;
	  s_p(3,i) += x*Vz1;
	  s_p(4,i) -= Vz1;
	}
	else {
	  s_p(1,i) -= (L-x)*Vy2;
	  s_p(2,i) += Vy2;
	  s_p(3,i) += (L-x)*Vz2;
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
      opserr << "ForceBeamColumn3d::addLoad() -- load type unknown for element with tag: " <<
	this->getTag() << endln;

      return -1;
    }

    return 0;
  }

  int 
  ForceBeamColumn3d::addInertiaLoadToUnbalance(const Vector &accel)
  {
    // Check for a quick return
    if (rho == 0.0)
      return 0;

    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);    

    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;

    // Should be done through p0[0]
    /*
    load(0) -= m*Raccel1(0);
    load(1) -= m*Raccel1(1);
    load(2) -= m*Raccel1(2);
    load(6) -= m*Raccel2(0);
    load(7) -= m*Raccel2(1);
    load(8) -= m*Raccel2(2);
    */

    return 0;
  }

  const Vector &
  ForceBeamColumn3d::getResistingForceIncInertia()
  {	
    // Compute the current resisting force
    theVector = this->getResistingForce();

    if (rho != 0.0) {
      const Vector &accel1 = theNodes[0]->getTrialAccel();
      const Vector &accel2 = theNodes[1]->getTrialAccel();

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
  ForceBeamColumn3d::sendSelf(int commitTag, Channel &theChannel)
  {  
    // place the integer data into an ID
    int dbTag = this->getDbTag();
    int i, j , k;
    int loc = 0;

    static ID idData(11);  
    idData(0) = this->getTag();
    idData(1) = connectedExternalNodes(0);
    idData(2) = connectedExternalNodes(1);
    idData(3) = numSections;
    idData(4) = maxIters;
    idData(5) = initialFlag;
    idData(6) = (isTorsion) ? 1 : 0;

    idData(7) = crdTransf->getClassTag();
    int crdTransfDbTag  = crdTransf->getDbTag();
    if (crdTransfDbTag  == 0) {
      crdTransfDbTag = theChannel.getDbTag();
      if (crdTransfDbTag  != 0) 
	crdTransf->setDbTag(crdTransfDbTag);
    }
    idData(8) = crdTransfDbTag;


    idData(9) = beamIntegr->getClassTag();
    int beamIntegrDbTag  = beamIntegr->getDbTag();
    if (beamIntegrDbTag  == 0) {
      beamIntegrDbTag = theChannel.getDbTag();
      if (beamIntegrDbTag  != 0) 
	beamIntegr->setDbTag(crdTransfDbTag);
    }
    idData(10) = beamIntegrDbTag;

    if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
      opserr << "ForceBeamColumn3d::sendSelf() - failed to send ID data\n";
      return -1;
    }    

    // send the coordinate transformation

    if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
      opserr << "ForceBeamColumn3d::sendSelf() - failed to send crdTranf\n";
      return -1;
    }      

    if (beamIntegr->sendSelf(commitTag, theChannel) < 0) {
      opserr << "ForceBeamColumn3d::sendSelf() - failed to send beamIntegr\n";
      return -1;
    }      

    //
    // send an ID for the sections containing each sections dbTag and classTag
    // if section ha no dbTag get one and assign it
    //

    ID idSections(2*numSections);
    loc = 0;
    for (i = 0; i<numSections; i++) {
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
      opserr << "ForceBeamColumn3d::sendSelf() - failed to send ID data\n";
      return -1;
    }    

    //
    // send the sections
    //

    for (j = 0; j<numSections; j++) {
      if (sections[j]->sendSelf(commitTag, theChannel) < 0) {
	opserr << "ForceBeamColumn3d::sendSelf() - section " <<
	  j << "failed to send itself\n";
	return -1;
      }
    }

    // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
    int secDefSize = 0;
    for (i = 0; i < numSections; i++) {
       int size = sections[i]->getOrder();
       secDefSize   += size;
    }

    Vector dData(1+1+NEBD+NEBD*NEBD+secDefSize + 4); 
    loc = 0;

    // place double variables into Vector
    dData(loc++) = rho;
    dData(loc++) = tol;

    // put  distrLoadCommit into the Vector
    //  for (i=0; i<NL; i++) 
    //dData(loc++) = distrLoadcommit(i);

    // place kvcommit into vector
    for (i=0; i<NEBD; i++) 
      dData(loc++) = Secommit(i);

    // place kvcommit into vector
    for (i=0; i<NEBD; i++) 
       for (j=0; j<NEBD; j++)
	  dData(loc++) = kvcommit(i,j);

    // place vscommit into vector
    for (k=0; k<numSections; k++)
       for (i=0; i<sections[k]->getOrder(); i++)
	  dData(loc++) = (vscommit[k])(i);

    // send damping coefficients
    dData(loc++) = alphaM;
    dData(loc++) = betaK;
    dData(loc++) = betaK0;
    dData(loc++) = betaKc;
    
    if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
       opserr << "ForceBeamColumn3d::sendSelf() - failed to send Vector data\n";

       return -1;
    }    

    return 0;
  }    

  int
  ForceBeamColumn3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
  {
    //
    // receive the integer data containing tag, numSections and coord transformation info
    //
    int dbTag = this->getDbTag();
    int i,j,k;

    static ID idData(11); // one bigger than needed 

    if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
      opserr << "ForceBeamColumn3d::recvSelf() - failed to recv ID data\n";

      return -1;
    }    

    this->setTag(idData(0));
    connectedExternalNodes(0) = idData(1);
    connectedExternalNodes(1) = idData(2);
    maxIters = idData(4);
    initialFlag = idData(5);
    isTorsion = (idData(6) == 1) ? true : false;

    int crdTransfClassTag = idData(7);
    int crdTransfDbTag = idData(8);

    int beamIntegrClassTag = idData(9);
    int beamIntegrDbTag = idData(10);

    // create a new crdTransf object if one needed
    if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
	if (crdTransf != 0)
	    delete crdTransf;

	crdTransf = theBroker.getNewCrdTransf3d(crdTransfClassTag);

	if (crdTransf == 0) {
	    opserr << "ForceBeamColumn3d::recvSelf() - failed to obtain a CrdTrans object with classTag" <<
	      crdTransfClassTag << endln;
	    return -2;	  
	}
    }

    crdTransf->setDbTag(crdTransfDbTag);

    // invoke recvSelf on the crdTransf obkject
    if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0)  
    {
       opserr << "ForceBeamColumn3d::sendSelf() - failed to recv crdTranf\n";
       return -3;
    }      


    // create a new beamIntegr object if one needed
    if (beamIntegr == 0 || beamIntegr->getClassTag() != beamIntegrClassTag) {
	if (beamIntegr != 0)
	    delete beamIntegr;

	beamIntegr = theBroker.getNewBeamIntegration(beamIntegrClassTag);

	if (beamIntegr == 0) {opserr << "ForceBeamColumn3d::recvSelf() - failed to obtain the beam integration object with classTag" <<
	    beamIntegrClassTag << endln;
	  exit(-1);
	}
    }

    beamIntegr->setDbTag(beamIntegrDbTag);

    // invoke recvSelf on the beamIntegr object
    if (beamIntegr->recvSelf(commitTag, theChannel, theBroker) < 0)  
    {
       opserr << "ForceBeamColumn3d::sendSelf() - failed to recv beam integration\n";

       return -3;
    }      

    //
    // recv an ID for the sections containing each sections dbTag and classTag
    //

    ID idSections(2*idData(3));
    int loc = 0;

    if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
      opserr << "ForceBeamColumn3d::recvSelf() - failed to recv ID data\n";
      return -1;
    }    

    //
    // now receive the sections
    //
    if (numSections != idData(3)) {

      //
      // we do not have correct number of sections, must delete the old and create
      // new ones before can recvSelf on the sections
      //

      // delete the old
      if (numSections != 0) {
	for (int i=0; i<numSections; i++)
	  delete sections[i];
	delete [] sections;
      }

      // create a section and recvSelf on it
      numSections = idData(3);

      // Delete the old
      if (vscommit != 0)
	delete [] vscommit;

      // Allocate the right number
      vscommit = new Vector[numSections];
      if (vscommit == 0) {
	opserr << "ForceBeamColumn3d::recvSelf -- failed to allocate vscommit array\n";
	return -1;
      }

      // Delete the old
      if (fs != 0)
	delete [] fs;

      // Allocate the right number
      fs = new Matrix[numSections];  
      if (fs == 0) {
	opserr << "ForceBeamColumn3d::recvSelf -- failed to allocate fs array\n";
	return -1;
      }

      // Delete the old
      if (vs != 0)
	delete [] vs;

      // Allocate the right number
      vs = new Vector[numSections];  
      if (vs == 0) {
	opserr << "ForceBeamColumn3d::recvSelf -- failed to allocate vs array\n";
	return -1;
      }

      // Delete the old
      if (Ssr != 0)
	delete [] Ssr;

      // Allocate the right number
      Ssr = new Vector[numSections];  
      if (Ssr == 0) {
	opserr << "ForceBeamColumn3d::recvSelf -- failed to allocate Ssr array\n";

	return -1;
      }

      // create a new array to hold pointers
      sections = new SectionForceDeformation *[idData(3)];
      if (sections == 0) {
	opserr << "ForceBeamColumn3d::recvSelf() - out of memory creating sections array of size" <<
	  idData(3) << endln;
	exit(-1);
      }    

      loc = 0;

      for (i=0; i<numSections; i++) {
	int sectClassTag = idSections(loc);
	int sectDbTag = idSections(loc+1);
	loc += 2;
	sections[i] = theBroker.getNewSection(sectClassTag);
	if (sections[i] == 0) {
	  opserr << "ForceBeamColumn3d::recvSelf() - Broker could not create Section of class type" <<
	    sectClassTag << endln;
	  exit(-1);
	}
	sections[i]->setDbTag(sectDbTag);
	if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "ForceBeamColumn3d::recvSelf() - section " << 
	    i << " failed to recv itself\n";
	  return -1;
	}     
      }

      this->initializeSectionHistoryVariables();

    } else {

      // 
      // for each existing section, check it is of correct type
      // (if not delete old & create a new one) then recvSelf on it
      //

      loc = 0;
      for (i=0; i<numSections; i++) {
	int sectClassTag = idSections(loc);
	int sectDbTag = idSections(loc+1);
	loc += 2;

	// check of correct type
	if (sections[i]->getClassTag() !=  sectClassTag) {
	  // delete the old section[i] and create a new one
	  delete sections[i];
	  sections[i] = theBroker.getNewSection(sectClassTag);
	  if (sections[i] == 0) {
	    opserr << "ForceBeamColumn3d::recvSelf() - Broker could not create Section of class type " 
		   << sectClassTag << endln;;
	    exit(-1);
	  }
	}

	// recvvSelf on it
	sections[i]->setDbTag(sectDbTag);
	if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "ForceBeamColumn3d::recvSelf() - section " <<
	    i << " failed to recv itself\n";
	  return -1;
	}     
      }
    }

    // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
    int secDefSize = 0;
    for (int ii = 0; ii < numSections; ii++) {
       int size = sections[ii]->getOrder();
       secDefSize   += size;
    }

    Vector dData(1+1+NEBD+NEBD*NEBD+secDefSize+4);   

    if (theChannel.recvVector(dbTag, commitTag, dData) < 0)  {
      opserr << "ForceBeamColumn3d::sendSelf() - failed to send Vector data\n";
      return -1;
    }    

    loc = 0;

    // place double variables into Vector
    rho = dData(loc++);
    tol = dData(loc++);

    // put  distrLoadCommit into the Vector
    //for (i=0; i<NL; i++) 
    // distrLoad(i) = dData(loc++);

    // place kvcommit into vector
    for (i=0; i<NEBD; i++) 
      Secommit(i) = dData(loc++);

    // place kvcommit into vector
    for (i=0; i<NEBD; i++) 
       for (j=0; j<NEBD; j++)
	  kvcommit(i,j) = dData(loc++);

    kv   = kvcommit;
    Se   = Secommit;

    for (k = 0; k < numSections; k++) {
      int order = sections[k]->getOrder();

      // place vscommit into vector
      vscommit[k] = Vector(order);
      for (i = 0; i < order; i++)
	(vscommit[k])(i) = dData(loc++);
    }

    // set damping coefficients
    alphaM = dData(loc++);
    betaK = dData(loc++);
    betaK0 = dData(loc++);
    betaKc = dData(loc++);
    
    initialFlag = 2;  

    return 0;
  }

  int
  ForceBeamColumn3d::getInitialFlexibility(Matrix &fe)
  {
    fe.Zero();

    double L = crdTransf->getInitialLength();
    double oneOverL  = 1.0/L;  

    // Flexibility from elastic interior
    beamIntegr->addElasticFlexibility(L, fe);

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);

    double wt[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wt);

    for (int i = 0; i < numSections; i++) {

      int order      = sections[i]->getOrder();
      const ID &code = sections[i]->getType();

      Matrix fb(workArea, order, NEBD);

      double xL  = xi[i];
      double xL1 = xL-1.0;
      double wtL = wt[i]*L;

      const Matrix &fSec = sections[i]->getInitialFlexibility();
      fb.Zero();
      double tmp;
      int ii, jj;
      for (ii = 0; ii < order; ii++) {
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  for (jj = 0; jj < order; jj++)
	    fb(jj,0) += fSec(jj,ii)*wtL;
	  break;
	case SECTION_RESPONSE_MZ:
	  for (jj = 0; jj < order; jj++) {
	    tmp = fSec(jj,ii)*wtL;
	    fb(jj,1) += xL1*tmp;
	    fb(jj,2) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VY:
	  for (jj = 0; jj < order; jj++) {
	    tmp = oneOverL*fSec(jj,ii)*wtL;
	    fb(jj,1) += tmp;
	    fb(jj,2) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_MY:
	  for (jj = 0; jj < order; jj++) {
	    tmp = fSec(jj,ii)*wtL;
	    fb(jj,3) += xL1*tmp;
	    fb(jj,4) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VZ:
	  for (jj = 0; jj < order; jj++) {
	    tmp = oneOverL*fSec(jj,ii)*wtL;
	    fb(jj,3) += tmp;
	    fb(jj,4) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_T:
	  for (jj = 0; jj < order; jj++)
	    fb(jj,5) += fSec(jj,ii)*wtL;
	  break;
	default:
	  break;
	}
      }
      for (ii = 0; ii < order; ii++) {
	switch (code(ii)) {
	case SECTION_RESPONSE_P:
	  for (jj = 0; jj < NEBD; jj++)
	    fe(0,jj) += fb(ii,jj);
	  break;
	case SECTION_RESPONSE_MZ:
	  for (jj = 0; jj < NEBD; jj++) {
	    tmp = fb(ii,jj);
	    fe(1,jj) += xL1*tmp;
	    fe(2,jj) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VY:
	  for (jj = 0; jj < NEBD; jj++) {
	    tmp = oneOverL*fb(ii,jj);
	    fe(1,jj) += tmp;
	    fe(2,jj) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_MY:
	  for (jj = 0; jj < NEBD; jj++) {
	    tmp = fb(ii,jj);
	    fe(3,jj) += xL1*tmp;
	    fe(4,jj) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VZ:
	  for (jj = 0; jj < NEBD; jj++) {
	    tmp = oneOverL*fb(ii,jj);
	    fe(3,jj) += tmp;
	    fe(4,jj) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_T:
	  for (jj = 0; jj < NEBD; jj++)
	    fe(5,jj) += fb(ii,jj);
	  break;
	default:
	  break;
	}
      }
    }

    if (!isTorsion)
      fe(5,5) = DefaultLoverGJ;

    return 0;
  }

  void
  ForceBeamColumn3d::compSectionDisplacements(Vector sectionCoords[],
					      Vector sectionDispls[]) const
  {
     // get basic displacements and increments
     static Vector ub(NEBD);
     ub = crdTransf->getBasicTrialDisp();    

     double L = crdTransf->getInitialLength();

     // get integration point positions and weights
     static double pts[maxNumSections];
     beamIntegr->getSectionLocations(numSections, L, pts);

     // setup Vandermode and CBDI influence matrices
     int i;
     double xi;

     // get CBDI influence matrix
     Matrix ls(numSections, numSections);
     getCBDIinfluenceMatrix(numSections, pts, L, ls);

     // get section curvatures
     Vector kappa_y(numSections);  // curvature
     Vector kappa_z(numSections);  // curvature
     static Vector vs;                // section deformations 

     for (i=0; i<numSections; i++) {
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
	   opserr << "FATAL ForceBeamColumn3d::compSectionResponse - section does not provide Mz response\n";
	   exit(-1);
	 }
	 if (sectionKey2 == 0) {
	   opserr << "FATAL ForceBeamColumn3d::compSectionResponse - section does not provide My response\n";
	   exit(-1);
	 }

	 // get section deformations
	 vs = sections[i]->getSectionDeformation();

	 kappa_z(i) = vs(sectionKey1);
	 kappa_y(i) = vs(sectionKey2); 
     }

     //cout << "kappa_y: " << kappa_y;   
     //cout << "kappa_z: " << kappa_z;   

     Vector v(numSections), w(numSections);
     static Vector xl(NDM), uxb(NDM);
     static Vector xg(NDM), uxg(NDM); 
     // double theta;                             // angle of twist of the sections

     // v = ls * kappa_z;  
     v.addMatrixVector (0.0, ls, kappa_z, 1.0);  
     // w = ls * kappa_y *  (-1);  
     w.addMatrixVector (0.0, ls, kappa_y, -1.0);

     for (i=0; i<numSections; i++)
     {
	xi = pts[i];

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
    return;	       
  }

  void
  ForceBeamColumn3d::Print(OPS_Stream &s, int flag)
  {
    // flags with negative values are used by GSA
    if (flag == -1) { 
      int eleTag = this->getTag();
      s << "EL_BEAM\t" << eleTag << "\t";
      s << sections[0]->getTag() << "\t" << sections[numSections-1]->getTag(); 
      s  << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
      s << "\t0\t0.0000000\n";
    }  

    // flags with negative values are used by GSA  
    else if (flag < -1) {
      int eleTag = this->getTag();
      int counter = (flag +1) * -1;
      double P  = Secommit(0);
      double MZ1 = Secommit(1);
      double MZ2 = Secommit(2);
      double MY1 = Secommit(3);
      double MY2 = Secommit(4);
      double L = crdTransf->getInitialLength();
      double VY = (MZ1+MZ2)/L;
      theVector(1) =  VY;
      theVector(4) = -VY;
      double VZ = (MY1+MY2)/L;
      double T  = Secommit(5);

      s << "FORCE\t" << eleTag << "\t" << counter << "\t0";
      s << "\t" << -P+p0[0] << "\t"  <<  VY+p0[1] << "\t"  << -VZ+p0[3]  << endln;
      s << "FORCE\t" << eleTag << "\t" << counter << "\t1";
      s << "\t"  << P  << ' '  << -VY+p0[2] << ' ' << VZ+p0[4] << endln;
      s << "MOMENT\t" << eleTag << "\t" << counter << "\t0";
      s << "\t" << -T << "\t"  << MY1 << "\t" << MZ1 << endln;
      s << "MOMENT\t" << eleTag << "\t" << counter << "\t1";
      s << "\t" << T << ' ' << MY2 << ' '  <<  MZ2 << endln;
    }

    // flag set to 2 used to print everything .. used for viewing data for UCSD renderer  
     else if (flag == 2) {
       static Vector xAxis(3);
       static Vector yAxis(3);
       static Vector zAxis(3);

       crdTransf->getLocalAxes(xAxis, yAxis, zAxis);

       s << "#ForceBeamColumn3D\n";
       s << "#LocalAxis " << xAxis(0) << " " << xAxis(1) << " " << xAxis(2);
       s << " " << yAxis(0) << " " << yAxis(1) << " " << yAxis(2) << " ";
       s << zAxis(0) << " " << zAxis(1) << " " << zAxis(2) << endln;

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

       double P  = Secommit(0);
       double MZ1 = Secommit(1);
       double MZ2 = Secommit(2);
       double MY1 = Secommit(3);
       double MY2 = Secommit(4);
       double L = crdTransf->getInitialLength();
       double VY = (MZ1+MZ2)/L;
       theVector(1) =  VY;
       theVector(4) = -VY;
       double VZ = (MY1+MY2)/L;
       double T  = Secommit(5);
       s << "#END_FORCES " << -P+p0[0] << ' '  <<  VY+p0[1] << ' '  << -VZ+p0[3] << ' ' 
	 << -T << ' '  << MY1 << ' ' << MZ1 << endln;
       s << "#END_FORCES "  << P  << ' '  << -VY+p0[2] << ' ' << VZ+p0[4] << ' '  
	 << T << ' ' << MY2 << ' '  <<  MZ2 << endln;

       // plastic hinge rotation
       static Vector vp(6);
       static Matrix fe(6,6);
       this->getInitialFlexibility(fe);
       vp = crdTransf->getBasicTrialDisp();
       vp.addMatrixVector(1.0, fe, Se, -1.0);
       s << "#PLASTIC_HINGE_ROTATION " << vp[1] << " " << vp[2] << " " << vp[3] << " " << vp[4] 
	 << " " << 0.1*L << " " << 0.1*L << endln;

       // allocate array of vectors to store section coordinates and displacements
       static int maxNumSections = 0;
       static Vector *coords = 0;
       static Vector *displs = 0;
       if (maxNumSections < numSections) {
	 if (coords != 0) 
	   delete [] coords;
	 if (displs != 0)
	   delete [] displs;

	 coords = new Vector [numSections];
	 displs = new Vector [numSections];

	 if (!coords) {
	   opserr << "ForceBeamColumn3d::Print() -- failed to allocate coords array";   
	   exit(-1);
	 }

	 int i;
	 for (i = 0; i < numSections; i++)
	   coords[i] = Vector(NDM);

	 if (!displs) {
	   opserr << "ForceBeamColumn3d::Print() -- failed to allocate coords array";   
	   exit(-1);
	 }

	 for (i = 0; i < numSections; i++)
	   displs[i] = Vector(NDM);

	 maxNumSections = numSections;
       }

       // compute section location & displacements
       this->compSectionDisplacements(coords, displs);

       // spit out the section location & invoke print on the scetion
       for (int i=0; i<numSections; i++) {
	 s << "#SECTION " << (coords[i])(0) << " " << (coords[i])(1) << " " << (coords[i])(2);       
	 s << " " << (displs[i])(0) << " " << (displs[i])(1) << " " << (displs[i])(2) << endln;
	 sections[i]->Print(s, flag); 
       }
     }

     else {
       s << "\nElement: " << this->getTag() << " Type: ForceBeamColumn3d ";
       s << "\tConnected Nodes: " << connectedExternalNodes ;
       s << "\tNumber of Sections: " << numSections;
       s << "\tMass density: " << rho << endln;
       beamIntegr->Print(s, flag);
       double P  = Secommit(0);
       double MZ1 = Secommit(1);
       double MZ2 = Secommit(2);
       double MY1 = Secommit(3);
       double MY2 = Secommit(4);
       double L = crdTransf->getInitialLength();
       double VY = (MZ1+MZ2)/L;
       theVector(1) =  VY;
       theVector(4) = -VY;
       double VZ = (MY1+MY2)/L;
       double T  = Secommit(5);
       s << "\tEnd 1 Forces (P MZ VY MY VZ T): "
	 << -P+p0[0] << " " << MZ1 << " " <<  VY+p0[1] << " " 
	 << MY1 << " " << -VZ+p0[3] << " " << T << endln;
       s << "\tEnd 2 Forces (P MZ VY MY VZ T): "
	 << P        << " " << MZ2 << " " << -VY+p0[2] << " " 
	 << MY2 << " " <<  VZ+p0[4] << " " << -T << endln;

       if (flag == 1) { 
	 for (int i = 0; i < numSections; i++)
	   s << "\numSections "<<i<<" :" << *sections[i];
       }
     }
  }

  OPS_Stream &operator<<(OPS_Stream &s, ForceBeamColumn3d &E)
  {
    E.Print(s);
    return s;
  }

  int
  ForceBeamColumn3d::displaySelf(Renderer &theViewer, int displayMode, float fact)
  {
    // first determine the end points of the beam based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    static Vector v1(3);
    static Vector v2(3);

    if (displayMode >= 0) {
      const Vector &end1Disp = theNodes[0]->getDisp();
      const Vector &end2Disp = theNodes[1]->getDisp();

      for (int i = 0; i < 3; i++) {
	v1(i) = end1Crd(i) + end1Disp(i)*fact;
	v2(i) = end2Crd(i) + end2Disp(i)*fact;    
      }
    } else {
      int mode = displayMode * -1;
      const Matrix &eigen1 = theNodes[0]->getEigenvectors();
      const Matrix &eigen2 = theNodes[1]->getEigenvectors();
      if (eigen1.noCols() >= mode) {
	for (int i = 0; i < 3; i++) {
	  v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	  v2(i) = end2Crd(i) + eigen2(i,mode-1)*fact;    
	}    
      } else {
	for (int i = 0; i < 3; i++) {
	  v1(i) = end1Crd(i);
	  v2(i) = end2Crd(i);
	}    
      }
    }

    return theViewer.drawLine (v1, v2, 1.0, 1.0);
  }

  Response*
  ForceBeamColumn3d::setResponse(const char **argv, int argc, OPS_Stream &output)
  {

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","ForceBeamColumn2d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    //
    // we compare argv[0] for known response types 
    //

    // global force - 
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");


      theResponse = new ElementResponse(this, 1, theVector);

    // local force -
    }  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

      output.tag("ResponseType","N_ 1");
      output.tag("ResponseType","Vy_1");
      output.tag("ResponseType","Vz_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Tz_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","T_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");

      theResponse = new ElementResponse(this, 2, theVector);

    // chord rotation -
    }  else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0 
	      || strcmp(argv[0],"basicDeformation") == 0) {

      output.tag("ResponseType","eps");
      output.tag("ResponseType","thetaZ_1");
      output.tag("ResponseType","thetaZ_2");
      output.tag("ResponseType","thetaY_1");
      output.tag("ResponseType","thetaY_2");
      output.tag("ResponseType","thetaX");

      theResponse = new ElementResponse(this, 3, Vector(6));

    // plastic rotation -
    } else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {

    output.tag("ResponseType","epsP");
    output.tag("ResponseType","thetaZP_1");
    output.tag("ResponseType","thetaZP_2");
    output.tag("ResponseType","thetaYP_1");
    output.tag("ResponseType","thetaYP_2");
    output.tag("ResponseType","thetaXP");

    theResponse = new ElementResponse(this, 4, Vector(6));
  
  // point of inflection
  } else if (strcmp(argv[0],"inflectionPoint") == 0) {
    theResponse = new ElementResponse(this, 5, Vector(2));
  
  // tangent drift
  } else if (strcmp(argv[0],"tangentDrift") == 0) {
    theResponse = new ElementResponse(this, 6, Vector(4));

  // section response -
  } else if (strcmp(argv[0],"section") ==0) { 
    if (argc > 2) {
    
      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections) {
	double xi[maxNumSections];
	double L = crdTransf->getInitialLength();
	beamIntegr->getSectionLocations(numSections, L, xi);
	
	output.tag("GaussPointOutput");
	output.attr("number",sectionNum);
	output.attr("eta",2.0*xi[sectionNum-1]-1.0);
	theResponse =  sections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
	
	output.endTag();
      }
    }
  }
  
  output.endTag();
  return theResponse;
}

int 
ForceBeamColumn3d::getResponse(int responseID, Information &eleInfo)
{
  static Vector vp(6);
  static Matrix fe(6,6);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());
  
  else if (responseID == 2) {
    // Axial
    double N = Se(0);
    theVector(6) =  N;
    theVector(0) = -N+p0[0];
    
    // Torsion
    double T = Se(5);
    theVector(9) =  T;
    theVector(3) = -T;
    
    // Moments about z and shears along y
    double M1 = Se(1);
    double M2 = Se(2);
    theVector(5)  = M1;
    theVector(11) = M2;
    double L = crdTransf->getInitialLength();
    double V = (M1+M2)/L;
    theVector(1) =  V+p0[1];
    theVector(7) = -V+p0[2];
    
    // Moments about y and shears along z
    M1 = Se(3);
    M2 = Se(4);
    theVector(4)  = M1;
    theVector(10) = M2;
    V = -(M1+M2)/L;
    theVector(2) = -V+p0[3];
    theVector(8) =  V+p0[4];
      
    return eleInfo.setVector(theVector);

  }
      
  // Chord rotation
  else if (responseID == 3) {
    vp = crdTransf->getBasicTrialDisp();
    return eleInfo.setVector(vp);
  }

  // Plastic rotation
  else if (responseID == 4) {
    this->getInitialFlexibility(fe);
    vp = crdTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, Se, -1.0);
    return eleInfo.setVector(vp);
  }

  // Point of inflection
  else if (responseID == 5) {
    static Vector LI(2);
    LI(0) = 0.0;
    LI(1) = 0.0;

    double L = crdTransf->getInitialLength();

    if (fabs(Se(1)+Se(2)) > DBL_EPSILON)
      LI(0) = Se(1)/(Se(1)+Se(2))*L;

    if (fabs(Se(3)+Se(4)) > DBL_EPSILON)
      LI(1) = Se(3)/(Se(3)+Se(4))*L;

    return eleInfo.setVector(LI);
  }

  // Tangent drift
  else if (responseID == 6) {
    double d2z = 0.0;
    double d2y = 0.0;
    double d3z = 0.0;
    double d3y = 0.0;

    double L = crdTransf->getInitialLength();

    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);

    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);

    // Location of inflection point from node I
    double LIz = 0.0;
    if (fabs(Se(1)+Se(2)) > DBL_EPSILON)
      LIz = Se(1)/(Se(1)+Se(2))*L;

    double LIy = 0.0;
    if (fabs(Se(3)+Se(4)) > DBL_EPSILON)
      LIy = Se(3)/(Se(3)+Se(4))*L;

    int i;
    for (i = 0; i < numSections; i++) {
      double x = pts[i]*L;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      if (x < LIz) {
	for (int j = 0; j < order; j++)
	  if (type(j) == SECTION_RESPONSE_MZ)
	    kappa += vs[i](j);
	double b = -LIz+x;
	d2z += (wts[i]*L)*kappa*b;
      }
      kappa = 0.0;
      if (x < LIy) {
	for (int j = 0; j < order; j++)
	  if (type(j) == SECTION_RESPONSE_MY)
	    kappa += vs[i](j);
	double b = -LIy+x;
	d2y += (wts[i]*L)*kappa*b;
      }
    }

    d2z += beamIntegr->getTangentDriftI(L, LIz, Se(1), Se(2));
    d2y += beamIntegr->getTangentDriftI(L, LIy, Se(3), Se(4), true);

    for (i = numSections-1; i >= 0; i--) {
      double x = pts[i]*L;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      if (x > LIz) {
	for (int j = 0; j < order; j++)
	  if (type(j) == SECTION_RESPONSE_MZ)
	    kappa += vs[i](j);
	double b = x-LIz;
	d3z += (wts[i]*L)*kappa*b;
      }
      kappa = 0.0;
      if (x > LIy) {
	for (int j = 0; j < order; j++)
	  if (type(j) == SECTION_RESPONSE_MY)
	    kappa += vs[i](j);
	double b = x-LIy;
	d3y += (wts[i]*L)*kappa*b;
      }
    }

    d3z += beamIntegr->getTangentDriftJ(L, LIz, Se(1), Se(2));
    d3y += beamIntegr->getTangentDriftJ(L, LIy, Se(3), Se(4), true);

    static Vector d(4);
    d(0) = d2z;
    d(1) = d3z;
    d(2) = d2y;
    d(3) = d3y;

    return eleInfo.setVector(d);
  }

  else
    return -1;
}

int
ForceBeamColumn3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // If the parameter belongs to the element itself
  if (strcmp(argv[0],"rho") == 0)
    return param.addObject(1, this);
  
  // If the parameter belongs to a section or lower
  if (strstr(argv[0],"section") != 0) {
    
    if (argc < 3)
      return -1;
    
    // Get section and material tag numbers from user input
    int paramSectionTag = atoi(argv[1]);
    
    // Find the right section and call its setParameter method
    int ok = 0;
    for (int i = 0; i < numSections; i++)
      if (paramSectionTag == sections[i]->getTag())
	ok += sections[i]->setParameter(&argv[2], argc-2, param);

    return ok;
  }
  
  else if (strstr(argv[0],"integration") != 0) {
    
    if (argc < 2)
      return -1;

    return beamIntegr->setParameter(&argv[1], argc-1, param);
  }
  else {
    return -1;
  }
}

int
ForceBeamColumn3d::updateParameter (int parameterID, Information &info)
{
  // If the parameterID value is not equal to 1 it belongs 
  // to section or material further down in the hierarchy. 
  
  return 0;

  /*
  if (parameterID == 1) {
    
    this->rho = info.theDouble;
    return 0;
    
  }
  else if (parameterID > 0 ) {
    
    // Extract the section number
    int sectionNumber = (int)( floor((double)parameterID) / (100000) );
    
    int ok = -1;
    for (int i=0; i<numSections; i++) {
      if (sectionNumber == sections[i]->getTag()) {
	ok = sections[i]->updateParameter(parameterID, info);
      }
    }
    
    if (ok < 0) {
      opserr << "ForceBeamColumn3d::updateParameter() - could not update parameter. " << endln;
      return ok;
    }
    else {
      return ok;
    }
  }
  else {
    opserr << "ForceBeamColumn3d::updateParameter() - could not update parameter. " << endln;
    return -1;
  }      
  */ 
}

void
ForceBeamColumn3d::setSectionPointers(int numSec, SectionForceDeformation **secPtrs)
{
  if (numSec > maxNumSections) {
    opserr << "Error: ForceBeamColumn3d::setSectionPointers -- max number of sections exceeded";
  }
  
  numSections = numSec;
  
  if (secPtrs == 0) {
    opserr << "Error: ForceBeamColumn3d::setSectionPointers -- invalid section pointer";
  }	  
  
  sections = new SectionForceDeformation *[numSections];
  if (sections == 0) {
    opserr << "Error: ForceBeamColumn3d::setSectionPointers -- could not allocate section pointers";
  }  
  
  for (int i = 0; i < numSections; i++) {
    
    if (secPtrs[i] == 0) {
      opserr << "Error: ForceBeamColumn3d::setSectionPointers -- null section pointer " << i << endln;
    }
    
    sections[i] = secPtrs[i]->getCopy();
    
    if (sections[i] == 0) {
      opserr << "Error: ForceBeamColumn3d::setSectionPointers -- could not create copy of section " << i << endln;
    }

    int order = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    for (int j = 0; j < order; j++) {
      if (code(j) == SECTION_RESPONSE_T)
	isTorsion = true;
    }
  }
  
  if (!isTorsion)
    opserr << "ForceBeamColumn3d::ForceBeamColumn3d -- no torsion detected in sections, " <<
      "continuing with element torsional stiffness GJ/L = " << 1.0/DefaultLoverGJ;

  // allocate section flexibility matrices and section deformation vectors
  fs  = new Matrix [numSections];
  if (fs == 0) {
    opserr << "ForceBeamColumn3d::setSectionPointers -- failed to allocate fs array";
  }
  
  vs = new Vector [numSections];
  if (vs == 0) {
    opserr << "ForceBeamColumn3d::setSectionPointers -- failed to allocate vs array";
  }
  
  Ssr  = new Vector [numSections];
  if (Ssr == 0) {
    opserr << "ForceBeamColumn3d::setSectionPointers -- failed to allocate Ssr array";
  }
  
  vscommit = new Vector [numSections];
  if (vscommit == 0) {
    opserr << "ForceBeamColumn3d::setSectionPointers -- failed to allocate vscommit array";   
  }
  
}

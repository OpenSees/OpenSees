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
                                                                        
// $Revision: 1.11 $
// $Date: 2003-05-12 23:44:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ForceBeamColumn2d.cpp,v $

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>


#include <Information.h>
#include <ForceBeamColumn2d.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>

#include <ElementResponse.h>
#include <ElementalLoad.h>

#define  NDM   2         // dimension of the problem (2d)
#define  NND   3         // number of nodal dof's
#define  NEGD  6         // number of element global dof's
#define  NEBD  3         // number of element dof's in the basic system

Matrix ForceBeamColumn2d::theMatrix(6,6);
Vector ForceBeamColumn2d::theVector(6);
double ForceBeamColumn2d::workArea[100];

Vector *ForceBeamColumn2d::vsSubdivide = 0;
Matrix *ForceBeamColumn2d::fsSubdivide = 0;
Vector *ForceBeamColumn2d::SsrSubdivide = 0;

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
ForceBeamColumn2d::ForceBeamColumn2d(): 
  Element(0,ELE_TAG_ForceBeamColumn2d), connectedExternalNodes(2), 
  beamIntegr(0), numSections(0), sections(0), crdTransf(0),
  rho(0.0), maxIters(0), tol(0.0),
  initialFlag(0),
  kv(NEBD,NEBD), Se(NEBD),
  kvcommit(NEBD,NEBD), Secommit(NEBD),
  fs(0), vs(0), Ssr(0), vscommit(0), sp(0), Ki(0)
{
  theNodes[0] = 0;  
  theNodes[1] = 0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  
  v0[0] = 0.0;
  v0[1] = 0.0;
  v0[2] = 0.0;

  if (vsSubdivide == 0)
    vsSubdivide  = new Vector [maxNumSections];
  if (fsSubdivide == 0)
    fsSubdivide  = new Matrix [maxNumSections];
  if (SsrSubdivide == 0)
    SsrSubdivide  = new Vector [maxNumSections];
  if (!vsSubdivide || !fsSubdivide || !SsrSubdivide) {
    opserr << "ForceBeamColumn2d::ForceBeamColumn2d() -- failed to allocate Subdivide arrays";   
    exit(-1);
  }
}

// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points. 
// allocates the necessary space needed by each object
ForceBeamColumn2d::ForceBeamColumn2d (int tag, int nodeI, int nodeJ,
					  int numSec, SectionForceDeformation **sec,
					  BeamIntegration &bi,
					  CrdTransf2d &coordTransf, double massDensPerUnitLength,
					  int maxNumIters, double tolerance):
  Element(tag,ELE_TAG_ForceBeamColumn2d), connectedExternalNodes(2),
  beamIntegr(0), numSections(0), sections(0), crdTransf(0),
  rho(massDensPerUnitLength),maxIters(maxNumIters), tol(tolerance), 
  initialFlag(0),
  kv(NEBD,NEBD), Se(NEBD), 
  kvcommit(NEBD,NEBD), Secommit(NEBD),
  fs(0), vs(0),Ssr(0), vscommit(0), sp(0), Ki(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;    
  
  beamIntegr = bi.getCopy();
  if (beamIntegr == 0) {
    opserr << "Error: ForceBeamColumn2d::ForceBeamColumn2d: could not create copy of beam integration object" << endln;
    exit(-1);
  }
  
  // get copy of the transformation object   
  crdTransf = coordTransf.getCopy(); 
  if (crdTransf == 0) {
    opserr << "Error: ForceBeamColumn2d::ForceBeamColumn2d: could not create copy of coordinate transformation object" << endln;
    exit(-1);
  }

  this->setSectionPointers(numSec, sec);
  
  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  
  v0[0] = 0.0;
  v0[1] = 0.0;
  v0[2] = 0.0;

  if (vsSubdivide == 0)
    vsSubdivide  = new Vector [maxNumSections];
  if (fsSubdivide == 0)
    fsSubdivide  = new Matrix [maxNumSections];
  if (SsrSubdivide == 0)
    SsrSubdivide  = new Vector [maxNumSections];
  if (!vsSubdivide || !fsSubdivide || !SsrSubdivide) {
    opserr << "ForceBeamColumn2d::ForceBeamColumn2d() -- failed to allocate Subdivide arrays";   
    exit(-1);
  }
}

// ~ForceBeamColumn2d():
// 	destructor
//      delete must be invoked on any objects created by the object
ForceBeamColumn2d::~ForceBeamColumn2d()
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
ForceBeamColumn2d::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
ForceBeamColumn2d::getExternalNodes(void) 
{
  return connectedExternalNodes;
}

Node **
ForceBeamColumn2d::getNodePtrs()
{
  return theNodes;
}

int
ForceBeamColumn2d::getNumDOF(void) 
{
  return NEGD;
}

void
ForceBeamColumn2d::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    opserr << "ForceBeamColumn2d::setDomain:  theDomain = 0 ";
    exit(0); 
  }

  // get pointers to the nodes
  
  int Nd1 = connectedExternalNodes(0);  
  int Nd2 = connectedExternalNodes(1);
  
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);  
  
  if (theNodes[0] == 0) {
    opserr << "ForceBeamColumn2d::setDomain: Nd1: ";
    opserr << Nd1 << "does not exist in model\n";
    exit(0);
  }
  
  if (theNodes[1] == 0) {
    opserr << "ForceBeamColumn2d::setDomain: Nd2: ";
    opserr << Nd2 << "does not exist in model\n";
    exit(0);
  }
  
  // call the DomainComponent class method 
  this->DomainComponent::setDomain(theDomain);
  
  // ensure connected nodes have correct number of dof's
  int dofNode1 = theNodes[0]->getNumberDOF();
  int dofNode2 = theNodes[1]->getNumberDOF();
  
  if ((dofNode1 != NND) || (dofNode2 != NND)) {
    opserr << "ForceBeamColumn2d::setDomain(): Nd2 or Nd1 incorrect dof ";
    exit(0);
  }
   
  // initialize the transformation
  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    opserr << "ForceBeamColumn2d::setDomain(): Error initializing coordinate transformation";  
    exit(0);
  }
    
  // get element length
  double L = crdTransf->getInitialLength();
  if (L == 0.0) {
    opserr << "ForceBeamColumn2d::setDomain(): Zero element length:" << this->getTag();  
    exit(0);
  }

  if (initialFlag == 0) 
    this->initializeSectionHistoryVariables();
}

int
ForceBeamColumn2d::commitState()
{
  int err = 0;
  int i = 0;

  // call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0) {
    opserr << "ForceBeamColumn2d::commitState () - failed in base class";
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

int ForceBeamColumn2d::revertToLastCommit()
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

int ForceBeamColumn2d::revertToStart()
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
ForceBeamColumn2d::getInitialStiff(void)
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
    opserr << "ForceBeamColumn2d::getInitialStiff() -- could not invert flexibility\n";
      

  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kvInit));
  
  return *Ki;
}

const Matrix &
ForceBeamColumn2d::getTangentStiff(void)
{
  crdTransf->update();	// Will remove once we clean up the corotational 2d transformation -- MHS
  return crdTransf->getGlobalStiffMatrix(kv, Se);
}
    
const Vector &
ForceBeamColumn2d::getResistingForce(void)
{
  // Will remove once we clean up the corotational 2d transformation -- MHS
  crdTransf->update();
  
  Vector p0Vec(p0, 3);
  
  return crdTransf->getGlobalResistingForce(Se, p0Vec);
}

void
ForceBeamColumn2d::initializeSectionHistoryVariables (void)
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
ForceBeamColumn2d::update()
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
    // initial tangent on first iteration then regular newton (if l==1), or 
    // initial tangent iterations (if l==2)

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
	  }

	  // Add effects of element loads
	  vr(0) += v0[0];
	  vr(1) += v0[1];
	  vr(2) += v0[2];
	  
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
		Ss(ii) =  xL1*SeTrial(1) + xL*SeTrial(2);
		break;
	      case SECTION_RESPONSE_VY:
		Ss(ii) = oneOverL*(SeTrial(1)+SeTrial(2));
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
	    fsSubdivide[i] = sections[i]->getSectionFlexibility();
	    
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
	      default:
		break;
	      }
	    }
	  }
	  
	  // calculate element stiffness matrix
	  // invert3by3Matrix(f, kv);	  
	  if (f.Solve(I, kvTrial) < 0)
	    opserr << "ForceBeamColumn2d::update() -- could not invert flexibility\n";
				    

	  // dv = vin + dvTrial  - vr
	  dv = vin;
	  dv += dvTrial;
	  dv -= vr;
	  
	  // dv.addVector(1.0, vr, -1.0);
	  
	  // dSe = kv * dv;
	  dSe.addMatrixVector(0.0, kvTrial, dv, 1.0);

	  dW = dv ^ dSe; 
	  
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
    opserr << "WARNING - ForceBeamColumn2d::update - failed to get compatible ";
    opserr << "element forces & deformations for element: ";
    opserr << this->getTag() << "(dW: << " << dW << ")\n";
    return -1;
  }

  initialFlag = 1;

  return 0;
}

void ForceBeamColumn2d::getForceInterpolatMatrix(double xi, Matrix &b, const ID &code)
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
    default:
      break;
    }
  }
}

void ForceBeamColumn2d::getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code)
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
    default:
      break;
    }
  }
}

const Matrix &
ForceBeamColumn2d::getDamp(void)
{
  theMatrix.Zero();
  
  return theMatrix; // zero matrix still
}

const Matrix &
ForceBeamColumn2d::getMass(void)
{ 
  theMatrix.Zero();
  
  double L = crdTransf->getInitialLength();
  if (rho != 0.0)
    theMatrix(0,0) = theMatrix(1,1) = theMatrix(3,3) = theMatrix(4,4) = 0.5*L*rho;
  
  return theMatrix;
}

void 
ForceBeamColumn2d::zeroLoad(void)
{
  if (sp != 0)
    sp->Zero();

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  v0[0] = 0.0;
  v0[1] = 0.0;
  v0[2] = 0.0;
}

int
ForceBeamColumn2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  
  if (sp == 0) {
    sp = new Matrix(3,numSections);
    if (sp == 0)
      opserr << "ForceBeamColumn2d::addLoad -- out of memory\n";
			    
  }

  double L = crdTransf->getInitialLength();

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wa = data(1)*loadFactor;  // Axial
    double wy = data(0)*loadFactor;  // Transverse

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < numSections; i++) {
      double x = xi[i]*L;
      // Axial
      s_p(0,i) += wa*(L-x);
      // Moment
      s_p(1,i) += wy*0.5*x*(x-L);
      // Shear
      s_p(2,i) += wy*(x-0.5*L);
    }

    // Accumulate reactions in basic system
    p0[0] -= wa*L;
    double V = 0.5*wy*L;
    p0[1] -= V;
    p0[2] -= V;

    // Accumulate initial deformations in basic system
    beamIntegr->addElasticDeformations(theLoad, loadFactor, L, v0);
  }
  else if (type == LOAD_TAG_Beam2dPointLoad) {
    double P = data(0)*loadFactor;
    double N = data(1)*loadFactor;
    double aOverL = data(2);
    double a = aOverL*L;

    double V1 = P*(1.0-aOverL);
    double V2 = P*aOverL;

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < numSections; i++) {
      double x = xi[i]*L;
      if (x <= a) {
	s_p(0,i) += N;
	s_p(1,i) -= x*V1;
	s_p(2,i) -= V1;
      }
      else {
	s_p(1,i) -= (L-x)*V2;
	s_p(2,i) += V2;
      }
    }

    // Accumulate reactions in basic system
    p0[0] -= N;
    p0[1] -= V1;
    p0[2] -= V2;

    // Accumulate initial deformations in basic system
    beamIntegr->addElasticDeformations(theLoad, loadFactor, L, v0);
  }

  else {
    opserr << "ForceBeamColumn2d::addLoad -- load type unknown for element with tag: " <<
      this->getTag() << endln;
    return -1;
  }

  return 0;
}

int 
ForceBeamColumn2d::addInertiaLoadToUnbalance(const Vector &accel)
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
  load(3) -= m*Raccel2(0);
  load(4) -= m*Raccel2(1);
  */

  return 0;
}

const Vector &
ForceBeamColumn2d::getResistingForceIncInertia()
{	
  // Check for a quick return
  if (rho == 0.0)
    return this->getResistingForce();
  
  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();
  
  // Compute the current resisting force
  theVector = this->getResistingForce();
  
  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;
  
  theVector(0) += m*accel1(0);
  theVector(1) += m*accel1(1);
  theVector(3) += m*accel2(0);
  theVector(4) += m*accel2(1);
  
  return theVector;
}

int
ForceBeamColumn2d::sendSelf(int commitTag, Channel &theChannel)
{  
  // place the integer data into an ID
  int dbTag = this->getDbTag();
  int i, j , k;
  int loc = 0;

  static ID idData(9);  // one bigger than needed so no clash later
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);
  idData(3) = numSections;
  idData(4) = maxIters;
  idData(5) = initialFlag;
  idData(6) = crdTransf->getClassTag();
  int crdTransfDbTag  = crdTransf->getDbTag();
  if (crdTransfDbTag  == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag  != 0) 
      crdTransf->setDbTag(crdTransfDbTag);
  }
  idData(7) = crdTransfDbTag;
  

  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "ForceBeamColumn2d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  // send the coordinate transformation
  
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "ForceBeamColumn2d::sendSelf() - failed to send crdTranf\n";
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
    opserr << "ForceBeamColumn2d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  //
  // send the sections
  //
  
  for (j = 0; j<numSections; j++) {
    if (sections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "ForceBeamColumn2d::sendSelf() - section " << 
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

  Vector dData(1+1+NEBD+NEBD*NEBD+secDefSize); 
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
  
  if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
    opserr << "ForceBeamColumn2d::sendSelf() - failed to send Vector data\n";
     return -1;
  }    

  return 0;
}    

int
ForceBeamColumn2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i,j,k;
  
  static ID idData(9); // one bigger than needed 

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    opserr << "ForceBeamColumn2d::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);
  maxIters = idData(4);
  initialFlag = idData(5);
  
  int crdTransfClassTag = idData(6);
  int crdTransfDbTag = idData(7);

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;

      crdTransf = theBroker.getNewCrdTransf2d(crdTransfClassTag);

      if (crdTransf == 0) {
	opserr << "ForceBeamColumn2d::recvSelf() - failed to obtain a CrdTrans object with classTag" <<
	  crdTransfClassTag << endln;
	exit(-1);
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf obkject
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     opserr << "ForceBeamColumn2d::sendSelf() - failed to recv crdTranf\n";
	     		     
     return -3;
  }      
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "ForceBeamColumn2d::recvSelf() - failed to recv ID data\n";
			    
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
      opserr << "ForceBeamColumn2d::recvSelf -- failed to allocate vscommit array\n";
			      
      return -1;
    }

    // Delete the old
    if (fs != 0)
      delete [] fs;

    // Allocate the right number
    fs = new Matrix[numSections];  
    if (fs == 0) {
      opserr << "ForceBeamColumn2d::recvSelf -- failed to allocate fs array\n";
      return -1;
    }

    // Delete the old
    if (vs != 0)
      delete [] vs;

    // Allocate the right number
    vs = new Vector[numSections];  
    if (vs == 0) {
      opserr << "ForceBeamColumn2d::recvSelf -- failed to allocate vs array\n";
      return -1;
    }

    // Delete the old
    if (Ssr != 0)
      delete [] Ssr;
    
    // Allocate the right number
    Ssr = new Vector[numSections];  
    if (Ssr == 0) {
      opserr << "ForceBeamColumn2d::recvSelf -- failed to allocate Ssr array\n";
      return -1;
    }

    // create a new array to hold pointers
    sections = new SectionForceDeformation *[idData(3)];
    if (sections == 0) {
      opserr << "ForceBeamColumn2d::recvSelf() - " << 
	"out of memory creating sections array of size" << idData(3) << endln;
      exit(-1);
    }    

    loc = 0;
    
    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      sections[i] = theBroker.getNewSection(sectClassTag);
      if (sections[i] == 0) {
	opserr << "ForceBeamColumn2d::recvSelf() - " << 
	  "Broker could not create Section of class type " << sectClassTag << endln;
	exit(-1);
      }
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "ForceBeamColumn2d::recvSelf() - section " << 
	  i << "failed to recv itself\n";
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
	  opserr << "ForceBeamColumn2d::recvSelf() - Broker could not create Section of class type" <<
	    sectClassTag << endln;
	  return -1;
	}
      }

      // recvvSelf on it
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "ForceBeamColumn2d::recvSelf() - section " << 
	  i << "failed to recv itself\n";
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
  
  Vector dData(1+1+NEBD+NEBD*NEBD+secDefSize);   
  
  if (theChannel.recvVector(dbTag, commitTag, dData) < 0)  {
    opserr << "ForceBeamColumn2d::sendSelf() - failed to send Vector data\n";
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

  initialFlag = 2;  

  return 0;
}

int
ForceBeamColumn2d::getInitialFlexibility(Matrix &fe)
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
      default:
	break;
      }
    }
  }
  
  return 0;
}

void ForceBeamColumn2d::compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const
{
  return;	       
}

void
ForceBeamColumn2d::Print(OPS_Stream &s, int flag)
{
  s << "\nElement: " << this->getTag() << " Type: ForceBeamColumn2d ";
  s << "\tConnected Nodes: " << connectedExternalNodes ;
  s << "\tNumber of Sections: " << numSections;
  s << "\tMass density: " << rho << endln;
  double P  = Secommit(0);
  double M1 = Secommit(1);
  double M2 = Secommit(2);
  double L = crdTransf->getInitialLength();
  double V = (M1+M2)/L;
  theVector(1) = V;
  theVector(4) = -V;
  s << "\tEnd 1 Forces (P V M): " << -P+p0[0] << " " << V+p0[1] << " " << M1 << endln;
  s << "\tEnd 2 Forces (P V M): " << P << " " << -V+p0[2] << " " << M2 << endln;
  
  if (flag == 1) { 
    for (int i = 0; i < numSections; i++)
      s << "\numSections "<<i<<" :" << *sections[i];
  }
}

OPS_Stream &operator<<(OPS_Stream &s, ForceBeamColumn2d &E)
{
  E.Print(s);
  return s;
}

int
ForceBeamColumn2d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
  // first determine the end points of the beam based on
  // the display factor (a measure of the distorted image)
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	
  
  const Vector &end1Disp = theNodes[0]->getDisp();
  const Vector &end2Disp = theNodes[1]->getDisp();
  
  static Vector v1(3);
  static Vector v2(3);
  
  for (int i = 0; i < 2; i++) {
    v1(i) = end1Crd(i) + end1Disp(i)*fact;
    v2(i) = end2Crd(i) + end2Disp(i)*fact;    
  }
  
  return theViewer.drawLine (v1, v2, 1.0, 1.0);
}

Response*
ForceBeamColumn2d::setResponse(const char **argv, int argc, Information &eleInformation)
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
  
  // chord rotation -
  else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0
	   || strcmp(argv[0],"basicDeformation") == 0)
    return new ElementResponse(this, 3, Vector(3));
  
  // plastic rotation -
  else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0)
    return new ElementResponse(this, 4, Vector(3));

  // point of inflection
  else if (strcmp(argv[0],"inflectionPoint") == 0)
    return new ElementResponse(this, 5, 0.0);
  
  // tangent drift
  else if (strcmp(argv[0],"tangentDrift") == 0)
    return new ElementResponse(this, 6, Vector(2));

  // section response -
  else if (strcmp(argv[0],"section") ==0) {
    if (argc <= 2)
      return 0;
    
    int sectionNum = atoi(argv[1]);
    if (sectionNum > 0 && sectionNum <= numSections)
      return sections[sectionNum-1]->setResponse(&argv[2], argc-2, eleInformation);
    else
      return 0;
  }
  
  else
    return 0;
}

int 
ForceBeamColumn2d::getResponse(int responseID, Information &eleInfo)
{
  static Vector vp(3);
  static Matrix fe(3,3);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());
  
  else if (responseID == 2) {
    theVector(3) =  Se(0);
    theVector(0) = -Se(0)+p0[0];
    theVector(2) = Se(1);
    theVector(5) = Se(2);
    double V = (Se(1)+Se(2))/crdTransf->getInitialLength();
    theVector(1) =  V+p0[1];
    theVector(4) = -V+p0[2];
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
    double LI = 0.0;

    if (fabs(Se(1)+Se(2)) > DBL_EPSILON) {
      double L = crdTransf->getInitialLength();
      
      LI = Se(1)/(Se(1)+Se(2))*L;
    }

    return eleInfo.setDouble(LI);
  }

  // Tangent drift
  else if (responseID == 6) {
    double d2 = 0.0;
    double d3 = 0.0;

    double L = crdTransf->getInitialLength();
    
    // Location of inflection point from node I
    double LI = 0.0;
    if (fabs(Se(1)+Se(2)) > DBL_EPSILON)
      LI = Se(1)/(Se(1)+Se(2))*L;
      
    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);
    
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    
    int i;
    for (i = 0; i < numSections; i++) {
      double x = pts[i]*L;
      if (x > LI)
	continue;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      for (int j = 0; j < order; j++)
	if (type(j) == SECTION_RESPONSE_MZ)
	  kappa += vs[i](j);
      double b = -LI+x;
      d2 += (wts[i]*L)*kappa*b;
    }
    
    d2 += beamIntegr->getTangentDriftI(L, LI, Se(1), Se(2));
    
    for (i = numSections-1; i >= 0; i--) {
      double x = pts[i]*L;
      if (x < LI)
	continue;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      for (int j = 0; j < order; j++)
	if (type(j) == SECTION_RESPONSE_MZ)
	  kappa += vs[i](j);
      double b = x-LI;
      d3 += (wts[i]*L)*kappa*b;
    }
    
    d3 += beamIntegr->getTangentDriftJ(L, LI, Se(1), Se(2));

    static Vector d(2);
    d(0) = d2;
    d(1) = d3;

    return eleInfo.setVector(d);
  }

  else
    return -1;
}

int
ForceBeamColumn2d::setParameter (const char **argv, int argc, Information &info)
{
  //
  // From the parameterID value it should be possible to extract
  // information about:
  //  1) Which parameter is in question. The parameter could
  //     be at element, section, or material level. 
  //  2) Which section and material number (tag) it belongs to. 
  //
  // To accomplish this the parameterID is given the following value:
  //     parameterID = type + 1000*matrTag + 100000*sectionTag
  // ...where 'type' is an integer in the range (1-99) and added 100
  // for each level (from material to section to element). 
  //
  // Example:
  //    If 'E0' (case 2) is random in material #3 of section #5
  //    the value of the parameterID at this (element) level would be:
  //    parameterID = 2 + 1000*3 + 100000*5 = 503002
  //    As seen, all given information can be extracted from this number. 
  //
  
  // Initial declarations
  int parameterID;

  // If the parameter belongs to the element itself
  if (strcmp(argv[0],"rho") == 0) {
    info.theType = DoubleType;
    return 1;
  }
  
  // If the parameter is belonging to a section or lower
  else if (strcmp(argv[0],"section") == 0) {
    
    // For now, no parameters of the section itself:
    if (argc<5) {
      opserr << "For now: cannot handle parameters of the section itself." << endln;
      return -1;
    }
    
    // Get section and material tag numbers from user input
    int paramSectionTag = atoi(argv[1]);
    
    // Find the right section and call its setParameter method
    for (int i=0; i<numSections; i++) {
      if (paramSectionTag == sections[i]->getTag()) {
	parameterID = sections[i]->setParameter(&argv[2], argc-2, info);
      }
    }
    
    // Check if the parameterID is valid
    if (parameterID < 0) {
      opserr << "ForceBeamColumn2d::setParameter() - could not set parameter. " << endln;
      return -1;
    }
    else {
      // Return the parameterID value (according to the above comments)
      return parameterID;
    }
  }
  
  // Otherwise parameter is unknown for this class
  else {
    return -1;
  }
}

int
ForceBeamColumn2d::updateParameter (int parameterID, Information &info)
{
  // If the parameterID value is not equal to 1 it belongs 
  // to section or material further down in the hierarchy. 
  
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
      opserr << "ForceBeamColumn2d::updateParameter() - could not update parameter. " << endln;
      return ok;
    }
    else {
      return ok;
    }
  }
  else {
    opserr << "ForceBeamColumn2d::updateParameter() - could not update parameter. " << endln;
    return -1;
  }       
}

void
ForceBeamColumn2d::setSectionPointers(int numSec, SectionForceDeformation **secPtrs)
{
  if (numSec > maxNumSections) {
    opserr << "Error: ForceBeamColumn2d::setSectionPointers -- max number of sections exceeded";
  }
  
  numSections = numSec;
  
  if (secPtrs == 0) {
    opserr << "Error: ForceBeamColumn2d::setSectionPointers -- invalid section pointer";
  }	  
  
  sections = new SectionForceDeformation *[numSections];
  if (sections == 0) {
    opserr << "Error: ForceBeamColumn2d::setSectionPointers -- could not allocate section pointers";
  }  
  
  for (int i = 0; i < numSections; i++) {
    
    if (secPtrs[i] == 0) {
      opserr << "Error: ForceBeamColumn2d::setSectionPointers -- null section pointer " << i << endln;
    }
    
    sections[i] = secPtrs[i]->getCopy();
    
    if (sections[i] == 0) {
      opserr << "Error: ForceBeamColumn2d::setSectionPointers -- could not create copy of section " << i << endln;
    }
  }
  
  // allocate section flexibility matrices and section deformation vectors
  fs  = new Matrix [numSections];
  if (fs == 0) {
    opserr << "ForceBeamColumn2d::setSectionPointers -- failed to allocate fs array";
  }
  
  vs = new Vector [numSections];
  if (vs == 0) {
    opserr << "ForceBeamColumn2d::setSectionPointers -- failed to allocate vs array";
  }
  
  Ssr  = new Vector [numSections];
  if (Ssr == 0) {
    opserr << "ForceBeamColumn2d::setSectionPointers -- failed to allocate Ssr array";
  }
  
  vscommit = new Vector [numSections];
  if (vscommit == 0) {
    opserr << "ForceBeamColumn2d::setSectionPointers -- failed to allocate vscommit array";   
  }
  
}

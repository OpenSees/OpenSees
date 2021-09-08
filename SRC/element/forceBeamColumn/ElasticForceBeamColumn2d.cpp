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
                                                                        
// $Revision: 1.1 $
// $Date: 2007-10-13 01:21:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ElasticForceBeamColumn2d.cpp,v $

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
#include <ElasticForceBeamColumn2d.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>
#include <elementAPI.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

Matrix ElasticForceBeamColumn2d::theMatrix(6,6);
Vector ElasticForceBeamColumn2d::theVector(6);
double ElasticForceBeamColumn2d::workArea[200];

void* OPS_ElasticForceBeamColumn2d()
{
    if(OPS_GetNumRemainingInputArgs() < 5) {
	opserr<<"insufficient arguments:eleTag,iNode,jNode,transfTag,integrationTag <-mass mass> <-cmass>\n";
	return 0;
    }

    // inputs: 
    int iData[5];
    int numData = 5;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	opserr<<"WARNING: invalid integer inputs\n";
	return 0;
    }

    // options
    double mass = 0.0;
    numData = 1;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type,"-mass") == 0) {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&mass) < 0) {
		    opserr<<"WARNING: invalid mass\n";
		    return 0;
		}
	    }
	}
    }

    // check transf
    CrdTransf* theTransf = OPS_getCrdTransf(iData[3]);
    if(theTransf == 0) {
	opserr<<"coord transfomration not found\n";
	return 0;
    }

    // check beam integrataion
    BeamIntegrationRule* theRule = OPS_getBeamIntegrationRule(iData[4]);
    if(theRule == 0) {
	opserr<<"beam integration not found\n";
	return 0;
    }
    BeamIntegration* bi = theRule->getBeamIntegration();
    if(bi == 0) {
	opserr<<"beam integration is null\n";
	return 0;
    }

    // check sections
    const ID& secTags = theRule->getSectionTags();
    SectionForceDeformation** sections = new SectionForceDeformation *[secTags.Size()];
    for(int i=0; i<secTags.Size(); i++) {
	sections[i] = OPS_getSectionForceDeformation(secTags(i));
	if(sections[i] == 0) {
	    opserr<<"section "<<secTags(i)<<"not found\n";
		delete [] sections;
	    return 0;
	}
    }
    
    Element *theEle =  new ElasticForceBeamColumn2d(iData[0],iData[1],iData[2],secTags.Size(),sections,*bi,*theTransf,mass);
    delete [] sections;
    return theEle;
}

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
ElasticForceBeamColumn2d::ElasticForceBeamColumn2d(): 
  Element(0,ELE_TAG_ElasticForceBeamColumn2d), connectedExternalNodes(2), 
  beamIntegr(0), numSections(0), crdTransf(0),
  rho(0.0), initialFlag(0),
  numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0),
  parameterID(0)
{
  theNodes[0] = 0;  
  theNodes[1] = 0;

  for (int i = 0; i < maxNumSections; i++)
    sections[i] = 0;
}

// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points. 
// allocates the necessary space needed by each object
ElasticForceBeamColumn2d::ElasticForceBeamColumn2d (int tag, 
						    int nodeI, int nodeJ,
						    int numSec, 
						    SectionForceDeformation **sec,
						    BeamIntegration &bi,
						    CrdTransf &coordTransf,
						    double massDensPerUnitLength):
  Element(tag,ELE_TAG_ElasticForceBeamColumn2d), connectedExternalNodes(2),
  beamIntegr(0), numSections(numSec), crdTransf(0),
  rho(massDensPerUnitLength),
  initialFlag(0),
  numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0),
  parameterID(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;    
  
  beamIntegr = bi.getCopy();
  if (beamIntegr == 0) {
    opserr << "Error: ElasticForceBeamColumn2d::ElasticForceBeamColumn2d: could not create copy of beam integration object" << endln;
  }
  
  // get copy of the transformation object   
  crdTransf = coordTransf.getCopy2d(); 
  if (crdTransf == 0) {
    opserr << "Error: ElasticForceBeamColumn2d::ElasticForceBeamColumn2d: could not create copy of coordinate transformation object" << endln;
  }

  if (numSections > maxNumSections) {
    opserr << "Error: ElasticForceBeamColumn2d::ElasticForceBeamColumn2d: numSections " << numSections << " exceeds max allowed, " << maxNumSections << endln;
    numSections = maxNumSections;
  }  

  int i;
  for (i = 0; i < numSections; i++) {
    sections[i] = sec[i]->getCopy();
    if (sections[i] == 0) {
      opserr << "Error: ElasticForceBeamColumn2d::ElasticForceBeamColumn2d: could not create copy of section object " << i << endln;
    }
  }
  for ( ; i < maxNumSections; i++)
    sections[i] = 0;
}

// ~ElasticForceBeamColumn2d():
// 	destructor
//      delete must be invoked on any objects created by the object
ElasticForceBeamColumn2d::~ElasticForceBeamColumn2d()
{
  for (int i=0; i < numSections; i++)
    if (sections[i] != 0)
      delete sections[i];
  
  if (sizeEleLoads != 0) {
      if (eleLoads != 0)
          delete[] eleLoads;

      if (eleLoadFactors != 0)
          delete[] eleLoadFactors;
  }

  if (crdTransf != 0)
    delete crdTransf;

  if (beamIntegr != 0)
    delete beamIntegr;
}

int
ElasticForceBeamColumn2d::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
ElasticForceBeamColumn2d::getExternalNodes(void) 
{
  return connectedExternalNodes;
}

Node **
ElasticForceBeamColumn2d::getNodePtrs()
{
  return theNodes;
}

int
ElasticForceBeamColumn2d::getNumDOF(void) 
{
  return NEGD;
}

void
ElasticForceBeamColumn2d::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    opserr << "ElasticForceBeamColumn2d::setDomain:  theDomain = 0 ";
  }

  // get pointers to the nodes
  
  int Nd1 = connectedExternalNodes(0);  
  int Nd2 = connectedExternalNodes(1);
  
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);  
  
  if (theNodes[0] == 0) {
    opserr << "ElasticForceBeamColumn2d::setDomain: Nd1: ";
    opserr << Nd1 << "does not exist in model\n";
  }
  
  if (theNodes[1] == 0) {
    opserr << "ElasticForceBeamColumn2d::setDomain: Nd2: ";
    opserr << Nd2 << "does not exist in model\n";
  }
  
  // call the DomainComponent class method 
  this->DomainComponent::setDomain(theDomain);
  
  // ensure connected nodes have correct number of dof's
  int dofNode1 = theNodes[0]->getNumberDOF();
  int dofNode2 = theNodes[1]->getNumberDOF();
  
  if ((dofNode1 != NND) || (dofNode2 != NND)) {
    opserr << "ElasticForceBeamColumn2d::setDomain(): Nd2 or Nd1 incorrect dof ";
  }
   
  // initialize the transformation
  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    opserr << "ElasticForceBeamColumn2d::setDomain(): Error initializing coordinate transformation";  
  }
    
  // get element length
  double L = crdTransf->getInitialLength();
  if (L == 0.0) {
    opserr << "ElasticForceBeamColumn2d::setDomain(): Zero element length:" << this->getTag();  
  }
}

int
ElasticForceBeamColumn2d::commitState()
{
  return crdTransf->commitState();
}

int ElasticForceBeamColumn2d::revertToLastCommit()
{
  return crdTransf->revertToLastCommit();
}

int ElasticForceBeamColumn2d::revertToStart()
{
  return crdTransf->revertToStart();
}


const Matrix &
ElasticForceBeamColumn2d::getInitialStiff(void)
{
  static Matrix f(NEBD, NEBD);   // element flexibility matrix  
  this->getInitialFlexibility(f);

  /*
  static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse  
  I.Zero();
  for (int i=0; i<NEBD; i++)
    I(i,i) = 1.0;
  
  // calculate element stiffness matrix
  // invert3by3Matrix(f, kv);

  static Matrix kvInit(NEBD, NEBD);
  if (f.Solve(I, kvInit) < 0)
    opserr << "ElasticForceBeamColumn2d::getInitialStiff() -- could not invert flexibility\n";
  */

  static Matrix kvInit(NEBD, NEBD);
  f.Invert(kvInit);

  static Vector SeInit(NEBD);
  SeInit.Zero();
  return crdTransf->getGlobalStiffMatrix(kvInit, SeInit);
}

const Matrix &
ElasticForceBeamColumn2d::getTangentStiff(void)
{
  crdTransf->update();	// Will remove once we clean up the corotational 2d transformation -- MHS
  static Matrix f(NEBD, NEBD);   // element flexibility matrix  
  this->getInitialFlexibility(f);

  /*
  static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse  
  I.Zero();
  for (int i=0; i<NEBD; i++)
    I(i,i) = 1.0;
  
  // calculate element stiffness matrix
  // invert3by3Matrix(f, kv);

  static Matrix kvInit(NEBD, NEBD);
  if (f.Solve(I, kvInit) < 0)
    opserr << "ElasticForceBeamColumn2d::getInitialStiff() -- could not invert flexibility\n";
  */

  static Matrix kvInit(NEBD, NEBD);
  f.Invert(kvInit);

  static Vector q(NEBD);
  q.Zero();
  this->computeBasicForces(q);

  return crdTransf->getGlobalStiffMatrix(kvInit, q);
}
    
void
ElasticForceBeamColumn2d::computeReactions(double *p0)
{
  int type;
  double L = crdTransf->getInitialLength();
  
  for (int i = 0; i < numEleLoads; i++) {
    
    double loadFactor = eleLoadFactors[i];
    const Vector &data = eleLoads[i]->getData(type, loadFactor);
    
    if (type == LOAD_TAG_Beam2dUniformLoad) {
      double wa = data(1)*loadFactor;  // Axial
      double wy = data(0)*loadFactor;  // Transverse

      p0[0] -= wa*L;
      double V = 0.5*wy*L;
      p0[1] -= V;
      p0[2] -= V;
    }
    else if (type == LOAD_TAG_Beam2dPointLoad) {
      double P = data(0)*loadFactor;
      double N = data(1)*loadFactor;
      double aOverL = data(2);
      
      if (aOverL < 0.0 || aOverL > 1.0)
	continue;
      
      
      double V1 = P*(1.0-aOverL);
      double V2 = P*aOverL;
      
      p0[0] -= N;
      p0[1] -= V1;
      p0[2] -= V2;
    }
  }
}

const Vector &
ElasticForceBeamColumn2d::getResistingForce(void)
{
  // Will remove once we clean up the corotational 2d transformation -- MHS
  crdTransf->update();
  
  double p0[3];
  Vector p0Vec(p0, 3);
  p0Vec.Zero();

  if (numEleLoads > 0)
    this->computeReactions(p0);

  static Matrix f(NEBD, NEBD);   // element flexibility matrix  
  this->getInitialFlexibility(f);

  static Vector Se(NEBD);
  this->computeBasicForces(Se);

  return crdTransf->getGlobalResistingForce(Se, p0Vec);
}

void
ElasticForceBeamColumn2d::computeBasicForces(Vector &q)
{
  if (q.Size() != NEBD) {
    opserr << "ElasticFBC2d::computeBasicForces -- q size not 3" << endln;
    return;
  }

  static Matrix f(NEBD, NEBD);   // element flexibility matrix  
  this->getInitialFlexibility(f);

  const Vector &v = crdTransf->getBasicTrialDisp();
  f.Solve(v, q);
}

/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS ********************
 */
int
ElasticForceBeamColumn2d::update()
{
  int ok = crdTransf->update();
	
  static Vector q(NEBD);
  q.Zero();
  this->computeBasicForces(q);
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  
  for (int i = 0; i < numSections; i++) {
    
    int order      = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    
    double xL  = xi[i];
    double xL1 = xL-1.0;
    
    static Vector s;
    s.setData(workArea, order);
    static Vector e;
    e.setData(&workArea[order], order);
    
    int ii;
    for (ii = 0; ii < order; ii++) {
      switch(code(ii)) {
      case SECTION_RESPONSE_P:
	s(ii) = q(0);
	break;
      case SECTION_RESPONSE_MZ:
	s(ii) =  xL1*q(1) + xL*q(2);
	break;
      case SECTION_RESPONSE_VY:
	s(ii) = oneOverL*(q(1)+q(2));
	break;
      default:
	s(ii) = 0.0;
	break;
      }
    }
    
    // Add the effects of element loads, if present
    // s = b*q + sp
    if (numEleLoads > 0)
      this->computeSectionForces(s, i);
    
    const Matrix &fs = sections[i]->getInitialFlexibility();
    e.addMatrixVector(0.0, fs, s, 1.0);
    
    ok += sections[i]->setTrialSectionDeformation(e);
  }
  
  return ok;
}

const Matrix &
ElasticForceBeamColumn2d::getMass(void)
{ 
  theMatrix.Zero();
  
  double L = crdTransf->getInitialLength();
  if (rho != 0.0)
    theMatrix(0,0) = theMatrix(1,1) = theMatrix(3,3) = theMatrix(4,4) = 0.5*L*rho;
  
  return theMatrix;
}

void 
ElasticForceBeamColumn2d::zeroLoad(void)
{
  // This is a semi-hack -- MHS
  numEleLoads = 0;

  return;
}

int
ElasticForceBeamColumn2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    if (numEleLoads == sizeEleLoads) {

        //
        // create larger arrays, copy old, delete old & set as new
        //

        ElementalLoad ** theNextEleLoads = new ElementalLoad *[sizeEleLoads + 1];
        double *theNextEleLoadFactors = new double[sizeEleLoads + 1];
        for (int i = 0; i<numEleLoads; i++) {
            theNextEleLoads[i] = eleLoads[i];
            theNextEleLoadFactors[i] = eleLoadFactors[i];
        }
        delete[] eleLoads;
        delete[] eleLoadFactors;
        eleLoads = theNextEleLoads;
        eleLoadFactors = theNextEleLoadFactors;

        // increment array size
        sizeEleLoads += 1;
    }

    eleLoadFactors[numEleLoads] = loadFactor;
    eleLoads[numEleLoads] = theLoad;
    numEleLoads++;

    return 0;
}

void
ElasticForceBeamColumn2d::computeSectionForces(Vector &sp, int isec)
{
  int type;

  double L = crdTransf->getInitialLength();

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  double x = xi[isec]*L;

  int order = sections[isec]->getOrder();
  const ID &code = sections[isec]->getType();

  for (int i = 0; i < numEleLoads; i++) {

    double loadFactor = eleLoadFactors[i];
    const Vector &data = eleLoads[i]->getData(type, loadFactor);
    
    if (type == LOAD_TAG_Beam2dUniformLoad) {
      double wa = data(1)*loadFactor;  // Axial
      double wy = data(0)*loadFactor;  // Transverse

      for (int ii = 0; ii < order; ii++) {
	
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  sp(ii) += wa*(L-x);
	  break;
	case SECTION_RESPONSE_MZ:
	  sp(ii) += wy*0.5*x*(x-L);
	  break;
	case SECTION_RESPONSE_VY:
	  sp(ii) += wy*(x-0.5*L);
	  break;
	default:
	  break;
	}
      }
    }
    else if (type == LOAD_TAG_Beam2dPointLoad) {
      double P = data(0)*loadFactor;
      double N = data(1)*loadFactor;
      double aOverL = data(2);
      
      if (aOverL < 0.0 || aOverL > 1.0)
	continue;
      
      double a = aOverL*L;
      
      double V1 = P*(1.0-aOverL);
      double V2 = P*aOverL;
      
      for (int ii = 0; ii < order; ii++) {
	
	if (x <= a) {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_P:
	    sp(ii) += N;
	    break;
	  case SECTION_RESPONSE_MZ:
	    sp(ii) -= x*V1;
	    break;
	  case SECTION_RESPONSE_VY:
	    sp(ii) -= V1;
	    break;
	  default:
	    break;
	  }
	}
	else {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_MZ:
	    sp(ii) -= (L-x)*V2;
	    break;
	  case SECTION_RESPONSE_VY:
	    sp(ii) += V2;
	    break;
	  default:
	    break;
	  }
	}
      }
    }
    else {
      opserr << "ElasticForceBeamColumn2d::addLoad -- load type unknown for element with tag: " <<
	this->getTag() << endln;
    }
  }
  
  // Don't think we need to do this anymore -- MHS
  //this->update(); // quick fix -- MHS
}

int 
ElasticForceBeamColumn2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  //const Vector &Raccel1 = theNodes[0]->getRV(accel);
  //  const Vector &Raccel2 = theNodes[1]->getRV(accel);    

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
ElasticForceBeamColumn2d::getResistingForceIncInertia()
{	
  // Compute the current resisting force
  theVector = this->getResistingForce();

  // Check for a quick return
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;
    
    theVector(0) += m*accel1(0);
    theVector(1) += m*accel1(1);
    theVector(3) += m*accel2(0);
    theVector(4) += m*accel2(1);
    
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
ElasticForceBeamColumn2d::sendSelf(int commitTag, Channel &theChannel)
{  
  return -1;
}    

int
ElasticForceBeamColumn2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}

int
ElasticForceBeamColumn2d::getInitialFlexibility(Matrix &fe)
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

void ElasticForceBeamColumn2d::compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const
{
   // get basic displacements and increments
   static Vector ub(NEBD);
   ub = crdTransf->getBasicTrialDisp();    

   double L = crdTransf->getInitialLength();
  
   // get integration point positions and weights
   //   const Matrix &xi_pt  = quadRule.getIntegrPointCoords(numSections);
   // get integration point positions and weights
   static double xi_pts[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi_pts);

   // setup Vandermode and CBDI influence matrices
   int i;
   double xi;
 
   // get CBDI influence matrix
   Matrix ls(numSections, numSections);
   getCBDIinfluenceMatrix(numSections, xi_pts, L, ls);

   // get section curvatures
   Vector kappa(numSections);  // curvature
   static Vector vs;              // section deformations 

   for (i=0; i<numSections; i++)
   {
       // THIS IS VERY INEFFICIENT ... CAN CHANGE LATER
       int sectionKey = 0;
       const ID &code = sections[i]->getType();
       int ii;
       for (ii = 0; ii < code.Size(); ii++)
	   if (code(ii) == SECTION_RESPONSE_MZ)
	   {
	       sectionKey = ii;
	       break;
	   }

       if (ii == code.Size()) {
	 opserr << "FATAL NLBeamColumn2d::compSectionDispls - section does not provide Mz response\n";
	 exit(-1);
       }
			
       // get section deformations
       vs = sections[i]->getSectionDeformation();
       kappa(i) = vs(sectionKey);
   }

   Vector w(numSections);
   static Vector xl(NDM), uxb(NDM);
   static Vector xg(NDM), uxg(NDM); 

   // w = ls * kappa;  
   w.addMatrixVector (0.0, ls, kappa, 1.0);
   
   for (i=0; i<numSections; i++)
   {
      xi = xi_pts[i];

      xl(0) = xi * L;
      xl(1) = 0;

      // get section global coordinates
      sectionCoords[i] = crdTransf->getPointGlobalCoordFromLocal(xl);

      // compute section displacements
      uxb(0) = xi * ub(0); // consider linear variation for axial displacement. CHANGE LATER!!!!!!!!!!
      uxb(1) = w(i);
             
      // get section displacements in global system 
      sectionDispls[i] = crdTransf->getPointGlobalDisplFromBasic(xi, uxb);
   }	       
   return;	       
}

void
ElasticForceBeamColumn2d::Print(OPS_Stream &s, int flag)
{
  static Vector Se(NEBD);
  static Vector vp(NEBD);
  static Matrix fe(NEBD,NEBD);

  if (flag == 2) {

    s << "#ElasticForceBeamColumn2D\n";

    const Vector &node1Crd = theNodes[0]->getCrds();
    const Vector &node2Crd = theNodes[1]->getCrds();	
    const Vector &node1Disp = theNodes[0]->getDisp();
    const Vector &node2Disp = theNodes[1]->getDisp();    
    
    s << "#NODE " << node1Crd(0) << " " << node1Crd(1) 
      << " " << node1Disp(0) << " " << node1Disp(1) << " " << node1Disp(2) << endln;
    
    s << "#NODE " << node2Crd(0) << " " << node2Crd(1) 
      << " " << node2Disp(0) << " " << node2Disp(1) << " " << node2Disp(2) << endln;

    this->computeBasicForces(Se);
    double P  = Se(0);
    double M1 = Se(1);
    double M2 = Se(2);
    double L = crdTransf->getInitialLength();
    double V = (M1+M2)/L;
    
    double p0[3]; p0[0] = 0.0; p0[1] = 0.0; p0[2] = 0.0;
    if (numEleLoads > 0)
      this->computeReactions(p0);

    s << "#END_FORCES " << -P+p0[0] << " " << V+p0[1] << " " << M1 << endln;
    s << "#END_FORCES " << P << " " << -V+p0[2] << " " << M2 << endln;

    // plastic hinge rotation
    this->getInitialFlexibility(fe);
    vp = crdTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, Se, -1.0);
    s << "#PLASTIC_HINGE_ROTATION " << vp[1] << " " << vp[2] << " " << 0.1*L << " " << 0.1*L << endln;
/*
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
	opserr << "NLBeamColumn3d::Print() -- failed to allocate coords array";   
	exit(-1);
      }
      
      int i;
      for (i = 0; i < numSections; i++)
	coords[i] = Vector(NDM);
      
      if (!displs) {
	opserr << "NLBeamColumn3d::Print() -- failed to allocate coords array";   
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
      s << "#SECTION " << (coords[i])(0) << " " << (coords[i])(1);       
      s << " " << (displs[i])(0) << " " << (displs[i])(1) << endln;
      sections[i]->Print(s, flag); 
    }
    */
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nElement: " << this->getTag() << " Type: ElasticForceBeamColumn2d ";
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho << endln;
    beamIntegr->Print(s, flag);
    crdTransf->Print(s, flag);

    this->computeBasicForces(Se);
    double P  = Se(0);
    double M1 = Se(1);
    double M2 = Se(2);
    double L = crdTransf->getInitialLength();
    double V = (M1+M2)/L;
    theVector(1) = V;
    theVector(4) = -V;
    double p0[3]; p0[0] = 0.0; p0[1] = 0.0; p0[2] = 0.0;
    if (numEleLoads > 0)
      this->computeReactions(p0);

    s << "\tEnd 1 Forces (P V M): " << -P+p0[0] << " " << V+p0[1] << " " << M1 << endln;
    s << "\tEnd 2 Forces (P V M): " << P << " " << -V+p0[2] << " " << M2 << endln;
    
    if (flag == 1) { 
      for (int i = 0; i < numSections; i++)
	s << "\numSections "<<i<<" :" << *sections[i];
    }
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
	  s << "\t\t\t{";
	  s << "\"name\": " << this->getTag() << ", ";
	  s << "\"type\": \"ElasticForceBeamColumn2d\", ";
	  s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
	  s << "\"sections\": [";
	  for (int i = 0; i < numSections - 1; i++)
		  s << "\"" << sections[i]->getTag() << "\", ";
	  s << "\"" << sections[numSections - 1]->getTag() << "\"], ";
	  s << "\"integration\": ";
	  beamIntegr->Print(s, flag);
	  s << ", \"massperlength\": " << rho << ", ";
	  s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"}";
  }
}

OPS_Stream &operator<<(OPS_Stream &s, ElasticForceBeamColumn2d &E)
{
  E.Print(s);
  return s;
}

int
ElasticForceBeamColumn2d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response*
ElasticForceBeamColumn2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","ElasticForceBeamColumn2d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

  // global force - 
  if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
      || strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 1, theVector);
  
  
  // local force -
  } else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_1");
    output.tag("ResponseType","V_1");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","V_2");
    output.tag("ResponseType","M_2");

    theResponse =  new ElementResponse(this, 2, theVector);
  

  // basic force -
  } else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","M_2");

    theResponse =  new ElementResponse(this, 7, Vector(3));

  // chord rotation -
  } else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0 
	     || strcmp(argv[0],"basicDeformation") == 0) {

    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta_1");
    output.tag("ResponseType","theta_2");

    theResponse =  new ElementResponse(this, 3, Vector(3));
  
  // plastic rotation -
  } else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {

    output.tag("ResponseType","epsP");
    output.tag("ResponseType","thetaP_1");
    output.tag("ResponseType","thetaP_2");

    theResponse =  new ElementResponse(this, 4, Vector(3));

  // point of inflection
  } else if (strcmp(argv[0],"inflectionPoint") == 0) {
    
    output.tag("ResponseType","inflectionPoint");

    theResponse =  new ElementResponse(this, 5, 0.0);
  
  // tangent drift
  } else if (strcmp(argv[0],"tangentDrift") == 0) {
    theResponse =  new ElementResponse(this, 6, Vector(2));

  // basic forces
  } else if (strcmp(argv[0],"basicForce") == 0)
    theResponse = new ElementResponse(this, 7, Vector(3));

  else if (strcmp(argv[0],"integrationPoints") == 0)
    theResponse = new ElementResponse(this, 10, Vector(numSections));

  else if (strcmp(argv[0],"integrationWeights") == 0)
    theResponse = new ElementResponse(this, 11, Vector(numSections));

  else if (strcmp(argv[0],"basicStiffness") == 0)
    theResponse = new ElementResponse(this, 12, Matrix(3,3));

  // section response -
  else if (strstr(argv[0],"sectionX") != 0) {
    if (argc > 2) {
      double sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamIntegr->getSectionLocations(numSections, L, xi);
      
      sectionLoc /= L;

      double minDistance = fabs(xi[0]-sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
	if (fabs(xi[i]-sectionLoc) < minDistance) {
	  minDistance = fabs(xi[i]-sectionLoc);
	  sectionNum = i;
	}
	  }

      output.tag("GaussPointOutput");
      output.attr("number",sectionNum+1);
      output.attr("eta",xi[sectionNum]*L);

      theResponse = sections[sectionNum]->setResponse(&argv[2], argc-2, output);
	}
  }

  // section response -
  else if (strstr(argv[0],"section") != 0) {
    if (argc > 2) {
      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections) {

	double xi[maxNumSections];
	double L = crdTransf->getInitialLength();
	beamIntegr->getSectionLocations(numSections, L, xi);

	output.tag("GaussPointOutput");
	output.attr("number",sectionNum);
	output.attr("eta",xi[sectionNum-1]*L);
	
	theResponse = sections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
	
      }
    }
  }
  
  output.endTag(); // ElementOutput

  return theResponse;
}

int 
ElasticForceBeamColumn2d::getResponse(int responseID, Information &eleInfo)
{
  static Vector Se(NEBD);
  static Vector vp(NEBD);
  static Matrix fe(NEBD,NEBD);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());
  
  else if (responseID == 2) {
    double p0[3]; p0[0] = 0.0; p0[1] = 0.0; p0[2] = 0.0;
    if (numEleLoads > 0)
      this->computeReactions(p0);
    this->computeBasicForces(Se);
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
  else if (responseID == 7) {
    this->computeBasicForces(Se);
    return eleInfo.setVector(Se);
  }
      
  // Chord rotation
  else if (responseID == 3) {
    vp = crdTransf->getBasicTrialDisp();
    return eleInfo.setVector(vp);
  }

  // Plastic rotation
  else if (responseID == 4) {
    this->computeBasicForces(Se);
    this->getInitialFlexibility(fe);
    vp = crdTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, Se, -1.0);
    return eleInfo.setVector(vp);
  }

  // Point of inflection
  else if (responseID == 5) {
    double LI = 0.0;
    this->computeBasicForces(Se);
    if (fabs(Se(1)+Se(2)) > DBL_EPSILON) {
      double L = crdTransf->getInitialLength();
      
      LI = Se(1)/(Se(1)+Se(2))*L;
    }

    return eleInfo.setDouble(LI);
  }

  else if (responseID == 7) {
    this->computeBasicForces(Se);
    return eleInfo.setVector(Se);
  }
  else if (responseID == 10) {
    double L = crdTransf->getInitialLength();
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i]*L;
    return eleInfo.setVector(locs);
  }

  else if (responseID == 11) {
    double L = crdTransf->getInitialLength();
    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wts[i]*L;
    return eleInfo.setVector(weights);
  }

  else if (responseID == 12) {
    static Matrix fb(NEBD,NEBD);
    this->getInitialFlexibility(fb);
    static Matrix kb(NEBD,NEBD);
    fb.Invert(kb);
    return eleInfo.setMatrix(kb);
  }

  else
    return -1;
}

int
ElasticForceBeamColumn2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return 0;

  if (strcmp(argv[0],"rho") == 0)
    return param.addObject(1, this);
  
  // section response -
  else if (strstr(argv[0],"sectionX") != 0) {
    if (argc > 2) {
      double sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamIntegr->getSectionLocations(numSections, L, xi);
      
      sectionLoc /= L;

      double minDistance = fabs(xi[0]-sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
	if (fabs(xi[i]-sectionLoc) < minDistance) {
	  minDistance = fabs(xi[i]-sectionLoc);
	  sectionNum = i;
	}
      }

      return sections[sectionNum]->setParameter(&argv[2], argc-2, param);
    }
  }

  // If the parameter belongs to a section or lower
  else if (strstr(argv[0],"section") != 0) {
    
    if (argc < 3)
      return 0;

    // Get section number: 1...Np
    int sectionNum = atoi(argv[1]);
    
    if (sectionNum > 0 && sectionNum <= numSections)
      return sections[sectionNum-1]->setParameter(&argv[2], argc-2, param);

    else
      return 0;

    /*
    // Find the right section and call its setParameter method
    int ok = 0;
    for (int i = 0; i < numSections; i++)
      if (paramSectionTag == sections[i]->getTag())
	ok += sections[i]->setParameter(&argv[2], argc-2, param);

    return ok;
    */
  }
  
  else if (strstr(argv[0],"integration") != 0) {
    
    if (argc < 2)
      return 0;

    return beamIntegr->setParameter(&argv[1], argc-1, param);
  }

  // Default, send to everything
  else {
    int ok = 0;
    for (int i = 0; i < numSections; i++)
      ok += sections[i]->setParameter(argv, argc, param);
    ok += beamIntegr->setParameter(argv, argc, param);
    return ok;
  }
}

int
ElasticForceBeamColumn2d::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  else
    return -1;
}

int
ElasticForceBeamColumn2d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;

  return 0;  
}

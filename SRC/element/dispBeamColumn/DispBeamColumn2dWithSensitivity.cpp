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

// $Revision: 1.2 $
// $Date: 2009-08-19 18:11:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/DispBeamColumn2dWithSensitivity.cpp,v $

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for DispBeamColumn2dWithSensitivity.

#include <DispBeamColumn2dWithSensitivity.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <elementAPI.h>

#define ELE_TAG_DispBeamColumn2dWithSensitivity 102030

Matrix DispBeamColumn2dWithSensitivity::K(6,6);
Vector DispBeamColumn2dWithSensitivity::P(6);
double DispBeamColumn2dWithSensitivity::workArea[100];

void* OPS_DispBeamColumn2dWithSensitivity()
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
    
    Element *theEle =  new DispBeamColumn2dWithSensitivity(iData[0],iData[1],iData[2],secTags.Size(),sections,*bi,*theTransf,mass);
    delete [] sections;
    return theEle;
}

DispBeamColumn2dWithSensitivity::DispBeamColumn2dWithSensitivity(int tag, int nd1, int nd2,
				   int numSec, SectionForceDeformation **s,
				   BeamIntegration& bi,
				   CrdTransf &coordTransf, double r)
:Element (tag, ELE_TAG_DispBeamColumn2dWithSensitivity),
 numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
  connectedExternalNodes(2),
  Q(6), q(3), rho(r)
{
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];

  if (theSections == 0) {
    opserr << "DispBeamColumn2dWithSensitivity::DispBeamColumn2dWithSensitivity - failed to allocate section model pointer\n";
    exit(-1);
  }

  for (int i = 0; i < numSections; i++) {

    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getCopy();

    // Check allocation
    if (theSections[i] == 0) {
      opserr << "DispBeamColumn2dWithSensitivity::DispBeamColumn2dWithSensitivity -- failed to get a copy of section model\n";
      exit(-1);
    }
  }

  beamInt = bi.getCopy();
  
  if (beamInt == 0) {
    opserr << "DispBeamColumn2d::DispBeamColumn2d - failed to copy beam integration\n";
    exit(-1);
  }

  crdTransf = coordTransf.getCopy2d();

  if (crdTransf == 0) {
    opserr << "DispBeamColumn2dWithSensitivity::DispBeamColumn2dWithSensitivity - failed to copy coordinate transformation\n";
    exit(-1);
  }

  // Set connected external node IDs
  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;

  theNodes[0] = 0;
  theNodes[1] = 0;

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END //////////////////////////////////////
}

DispBeamColumn2dWithSensitivity::DispBeamColumn2dWithSensitivity()
:Element (0, ELE_TAG_DispBeamColumn2dWithSensitivity),
 numSections(0), theSections(0), crdTransf(0), beamInt(0),
 connectedExternalNodes(2),
  Q(6), q(3), rho(0.0)
{
    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;

    theNodes[0] = 0;
    theNodes[1] = 0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END //////////////////////////////////////
}

DispBeamColumn2dWithSensitivity::~DispBeamColumn2dWithSensitivity()
{
  for (int i = 0; i < numSections; i++) {
    if (theSections[i])
      delete theSections[i];
  }
  
  // Delete the array of pointers to SectionForceDeformation pointer arrays
  /*if (theSections)
    delete [] theSections;
  */
  if (crdTransf)
    delete crdTransf;

  if (beamInt != 0)
    delete beamInt;

}

int
DispBeamColumn2dWithSensitivity::getNumExternalNodes() const
{
    return 2;
}

const ID&
DispBeamColumn2dWithSensitivity::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
DispBeamColumn2dWithSensitivity::getNodePtrs()
{
    return theNodes;
}

int
DispBeamColumn2dWithSensitivity::getNumDOF()
{
    return 6;
}

void
DispBeamColumn2dWithSensitivity::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);

    if (theNodes[0] == 0 || theNodes[1] == 0) {
	//opserr << "FATAL ERROR DispBeamColumn2dWithSensitivity (tag: %d), node not found in domain",
	//	this->getTag());

	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    if (dofNd1 != 3 || dofNd2 != 3) {
	//opserr << "FATAL ERROR DispBeamColumn2dWithSensitivity (tag: %d), has differing number of DOFs at its nodes",
	//	this->getTag());

	return;
    }

	if (crdTransf->initialize(theNodes[0], theNodes[1])) {
		// Add some error check
	}

	double L = crdTransf->getInitialLength();

	if (L == 0.0) {
		// Add some error check
	}

    this->DomainComponent::setDomain(theDomain);

	this->update();
}

int
DispBeamColumn2dWithSensitivity::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "DispBeamColumn2dWithSensitivity::commitState () - failed in base class";
    }

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    return retVal;
}

int
DispBeamColumn2dWithSensitivity::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    return retVal;
}

int
DispBeamColumn2dWithSensitivity::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    return retVal;
}

int
DispBeamColumn2dWithSensitivity::update(void)
{
  int err = 0;
  // Update the transformation
  crdTransf->update();

  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    Vector e(workArea, order);

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];

    int j;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	e(j) = oneOverL*v(0); break;
      case SECTION_RESPONSE_MZ:
	e(j) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2)); break;
      default:
	e(j) = 0.0; break;
      }
    }

//	opserr<<"element:"<<this->getTag()<<endln;
    // Set the section deformations
    err += theSections[i]->setTrialSectionDeformation(e);
  }
  if (err != 0) {
    opserr << "DispBeamColumn2dWithSensitivity::update() - failed setTrialSectionDeformations()\n";
    return err;
  }
  return 0;
}

const Matrix&
DispBeamColumn2dWithSensitivity::getTangentStiff()
{
  static Matrix kb(3,3);

  // Zero for integral
  kb.Zero();
  q.Zero();

 
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    Matrix ka(workArea, order, 3);
    ka.Zero();

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];

    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getSectionTangent();
    const Vector &s = theSections[i]->getStressResultant();

    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    //double wti = wts(i)*oneOverL;
    double wti = wt[i]*oneOverL;

    double tmp;
    int j, k;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < order; k++)
	  ka(k,0) += ks(k,j)*wti;
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < order; k++) {
	  tmp = ks(k,j)*wti;
	  ka(k,1) += (xi6-4.0)*tmp;
	  ka(k,2) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < 3; k++)
	  kb(0,k) += ka(j,k);
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < 3; k++) {
	  tmp = ka(j,k);
	  kb(1,k) += (xi6-4.0)*tmp;
	  kb(2,k) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }

    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));
    double si;
    for (j = 0; j < order; j++) {
      //si = s(j)*wts(i);
      si = s(j)*wt[i];
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	q(0) += si; break;
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si; break;
      default:
	break;
      }
    }

  }

  // Add effects of element loads, q = q(v) + q0
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];

  // Transform to global stiffness
  K = crdTransf->getGlobalStiffMatrix(kb, q);

  //opserr<<"in DispBeamColumn2dWithSensitivity: K is"<<K<<endln;
  return K;
}

const Matrix&
DispBeamColumn2dWithSensitivity::getInitialBasicStiff()
{
  static Matrix kb(3,3);

  // Zero for integral
  kb.Zero();

//  const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
//  const Vector &wts = quadRule.getIntegrPointWeights(numSections);

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    Matrix ka(workArea, order, 3);
    ka.Zero();

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];

    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getInitialTangent();

    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    //double wti = wts(i)*oneOverL;
    double wti = wt[i]*oneOverL;
    double tmp;
    int j, k;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < order; k++)
	  ka(k,0) += ks(k,j)*wti;
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < order; k++) {
	  tmp = ks(k,j)*wti;
	  ka(k,1) += (xi6-4.0)*tmp;
	  ka(k,2) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < 3; k++)
	  kb(0,k) += ka(j,k);
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < 3; k++) {
	  tmp = ka(j,k);
	  kb(1,k) += (xi6-4.0)*tmp;
	  kb(2,k) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }

  }

  return kb;
}

const Matrix&
DispBeamColumn2dWithSensitivity::getInitialStiff()
{
  const Matrix &kb = this->getInitialBasicStiff();

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);



  return K;
}

const Matrix&
DispBeamColumn2dWithSensitivity::getMass()
{
  K.Zero();

  if (rho == 0.0)
    return K;

  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;

  K(0,0) = K(1,1) = K(3,3) = K(4,4) = m;

  return K;
}

void
DispBeamColumn2dWithSensitivity::zeroLoad(void)
{
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  return;
}

int
DispBeamColumn2dWithSensitivity::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = crdTransf->getInitialLength();

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wt = data(0)*loadFactor;  // Transverse (+ve upward)
    double wa = data(1)*loadFactor;  // Axial (+ve from node I to J)

    double V = 0.5*wt*L;
    double M = V*L/6.0; // wt*L*L/12
    double P = wa*L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= V;
    p0[2] -= V;

    // Fixed end forces in basic system
    q0[0] -= 0.5*P;
    q0[1] -= M;
    q0[2] += M;
  }
  else if (type == LOAD_TAG_Beam2dPointLoad) {
    double P = data(0)*loadFactor;
    double N = data(1)*loadFactor;
    double aOverL = data(2);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;
    double b = L-a;

    // Reactions in basic system
    p0[0] -= N;
    double V1 = P*(1.0-aOverL);
    double V2 = P*aOverL;
    p0[1] -= V1;
    p0[2] -= V2;

    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    double M1 = -a * b2 * P * L2;
    double M2 = a2 * b * P * L2;
    q0[1] += M1;
    q0[2] += M2;
  }
  else {
    opserr << "DispBeamColumn2dWithSensitivity::DispBeamColumn2dWithSensitivity -- load type unknown for element with tag: "
	   << this->getTag() << "DispBeamColumn2dWithSensitivity::addLoad()\n";

    return -1;
  }

  return 0;
}

int
DispBeamColumn2dWithSensitivity::addInertiaLoadToUnbalance(const Vector &accel)
{
	// Check for a quick return
	if (rho == 0.0)
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = theNodes[0]->getRV(accel);
	const Vector &Raccel2 = theNodes[1]->getRV(accel);

    if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
      opserr << "DispBeamColumn2dWithSensitivity::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
      return -1;
    }

	double L = crdTransf->getInitialLength();
	double m = 0.5*rho*L;

    // Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	Q(0) -= m*Raccel1(0);
	Q(1) -= m*Raccel1(1);
	Q(3) -= m*Raccel2(0);
	Q(4) -= m*Raccel2(1);

    return 0;
}

const Vector&
DispBeamColumn2dWithSensitivity::getResistingForce()
{
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);  
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Zero for integration
  q.Zero();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];

    // Get section stress resultant
    const Vector &s = theSections[i]->getStressResultant();

    // Perform numerical integration on internal force
    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));

    double si;
    for (int j = 0; j < order; j++) {
      //si = s(j)*wts(i);
      si = s(j)*wt[i];
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	q(0) += si; break;
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si; break;
      default:
	break;
      }
    }

  }

  // Add effects of element loads, q = q(v) + q0
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];

  // Vector for reactions in basic system
  Vector p0Vec(p0, 3);

  P = crdTransf->getGlobalResistingForce(q, p0Vec);

  return P;
}

const Vector&
DispBeamColumn2dWithSensitivity::getResistingForceIncInertia()
{
  P = this->getResistingForce();
  
  // Subtract other external nodal loads ... P_res = P_int - P_ext
  P.addVector(1.0, Q, -1.0);
  
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();

    // Compute the current resisting force
    this->getResistingForce();

    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;

    P(0) += m*accel1(0);
    P(1) += m*accel1(1);
    P(3) += m*accel2(0);
    P(4) += m*accel2(1);

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P.addVector(1.0, this->getRayleighDampingForces(), 1.0);

  } else {

    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
  }

  return P;
}

int
DispBeamColumn2dWithSensitivity::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j;
  int loc = 0;

  static ID idData(7);  // one bigger than needed so no clash later
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);
  idData(3) = numSections;
  idData(4) = crdTransf->getClassTag();
  int crdTransfDbTag  = crdTransf->getDbTag();
  if (crdTransfDbTag  == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag  != 0)
      crdTransf->setDbTag(crdTransfDbTag);
  }
  idData(5) = crdTransfDbTag;
  if (alphaM != 0 || betaK != 0 || betaK0 != 0 || betaKc != 0) 
    idData(6) = 1;
  else
    idData(6) = 0;
  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "DispBeamColumn2dWithSensitivity::sendSelf() - failed to send ID data\n";
     return -1;
  }
  if (idData(6) == 1) {
    // send damping coefficients
    static Vector dData(4);
    dData(0) = alphaM;
    dData(1) = betaK;
    dData(2) = betaK0;
    dData(3) = betaKc;
    if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
      opserr << "DispBeamColumn2d::sendSelf() - failed to send double data\n";
      return -1;
    } 
  }

  // send the coordinate transformation

  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "DispBeamColumn2dWithSensitivity::sendSelf() - failed to send crdTranf\n";
     return -1;
  }


  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //

  ID idSections(2*numSections);
  loc = 0;
  for (i = 0; i<numSections; i++) {
    int sectClassTag = theSections[i]->getClassTag();
    int sectDbTag = theSections[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      theSections[i]->setDbTag(sectDbTag);
    }

    idSections(loc) = sectClassTag;
    idSections(loc+1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumn2dWithSensitivity::sendSelf() - failed to send ID data\n";
    return -1;
  }

  //
  // send the sections
  //

  for (j = 0; j<numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "DispBeamColumn2dWithSensitivity::sendSelf() - section " <<
	j << "failed to send itself\n";
      return -1;
    }
  }

  return 0;
}

int
DispBeamColumn2dWithSensitivity::recvSelf(int commitTag, Channel &theChannel,
			   FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i;

  static ID idData(7); // one bigger than needed so no clash with section ID

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    opserr << "DispBeamColumn2dWithSensitivity::recvSelf() - failed to recv ID data\n";
    return -1;
  }

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);

  int crdTransfClassTag = idData(4);
  int crdTransfDbTag = idData(5);

  if (idData(6) == 1) {
    // recv damping coefficients
    static Vector dData(4);
    if (theChannel.recvVector(dbTag, commitTag, dData) < 0) {
      opserr << "DispBeamColumn2d::sendSelf() - failed to recv double data\n";
      return -1;
    }    
    alphaM = dData(0);
    betaK = dData(1);
    betaK0 = dData(2);
    betaKc = dData(3);
  }

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;

      crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

      if (crdTransf == 0) {
	opserr << "DispBeamColumn2dWithSensitivity::recvSelf() - failed to obtain a CrdTrans object with classTag " <<
	  crdTransfClassTag << endln;
	  return -2;
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "DispBeamColumn2dWithSensitivity::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }

  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumn2dWithSensitivity::recvSelf() - failed to recv ID data\n";
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
	delete theSections[i];
      delete [] theSections;
    }

    // create a new array to hold pointers
    theSections = new SectionForceDeformation *[idData(3)];
    if (theSections == 0) {
opserr << "DispBeamColumn2dWithSensitivity::recvSelf() - out of memory creating sections array of size " <<
  idData(3) << endln;
      return -1;
    }

    // create a section and recvSelf on it
    numSections = idData(3);
    loc = 0;

    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      theSections[i] = theBroker.getNewSection(sectClassTag);
      if (theSections[i] == 0) {
	opserr << "DispBeamColumn2dWithSensitivity::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumn2dWithSensitivity::recvSelf() - section " << i << " failed to recv itself\n";
	return -1;
      }
    }

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
      if (theSections[i]->getClassTag() !=  sectClassTag) {
	// delete the old section[i] and create a new one
	delete theSections[i];
	theSections[i] = theBroker.getNewSection(sectClassTag);
	if (theSections[i] == 0) {
	opserr << "DispBeamColumn2dWithSensitivity::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
	}
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumn2dWithSensitivity::recvSelf() - section " << i << " failed to recv itself\n";
	return -1;
      }
    }
  }

  return 0;
}

void
DispBeamColumn2dWithSensitivity::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nDispBeamColumn2dWithSensitivity, element id:  " << this->getTag() << endln;
        s << "\tConnected external nodes:  " << connectedExternalNodes;
        s << "\tCoordTransf: " << crdTransf->getTag() << endln;
        s << "\tmass density:  " << rho << endln;

        double L = crdTransf->getInitialLength();
        double P = q(0);
        double M1 = q(1);
        double M2 = q(2);
        double V = (M1 + M2) / L;
        s << "\tEnd 1 Forces (P V M): " << -P + p0[0]
            << " " << V + p0[1] << " " << M1 << endln;
        s << "\tEnd 2 Forces (P V M): " << P
            << " " << -V + p0[2] << " " << M2 << endln;

        beamInt->Print(s, flag);

        for (int i = 0; i < numSections; i++)
            theSections[i]->Print(s, flag);
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"DispBeamColumn2dWithSensitivity\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"sections\": [";
        for (int i = 0; i < numSections - 1; i++)
            s << "\"" << theSections[i]->getTag() << "\", ";
        s << "\"" << theSections[numSections - 1]->getTag() << "\"], ";
        s << "\"integration\": ";
        beamInt->Print(s, flag);
        s << ", \"massperlength\": " << rho << ", ";
        s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"}";
    }
}


int
DispBeamColumn2dWithSensitivity::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

  static Vector v1(3);
  static Vector v2(3);

  if (displayMode >= 0) {
    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();
    
    for (int i = 0; i < 2; i++) {
      v1(i) = end1Crd(i) + end1Disp(i)*fact;
      v2(i) = end2Crd(i) + end2Disp(i)*fact;    
    }
  } else {
    int mode = displayMode  *  -1;
    const Matrix &eigen1 = theNodes[0]->getEigenvectors();
    const Matrix &eigen2 = theNodes[1]->getEigenvectors();
    if (eigen1.noCols() >= mode) {
      for (int i = 0; i < 2; i++) {
	v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	v2(i) = end2Crd(i) + eigen2(i,mode-1)*fact;    
      }    
    } else {
      for (int i = 0; i < 2; i++) {
	v1(i) = end1Crd(i);
	v2(i) = end2Crd(i);
      }    
    }
  }
	
  return theViewer.drawLine (v1, v2, 1.0, 1.0);
}

Response*
DispBeamColumn2dWithSensitivity::setResponse(const char **argv, int argc,
			      OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","DispBeamColumn2d");
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

    theResponse =  new ElementResponse(this, 1, P);
  
  
  // local force -
  } else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N1");
    output.tag("ResponseType","V1");
    output.tag("ResponseType","M1");
    output.tag("ResponseType","N2");
    output.tag("ResponseType","V2");
    output.tag("ResponseType","M2");

    theResponse =  new ElementResponse(this, 2, P);
  

  // basic force -
  } else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M1");
    output.tag("ResponseType","M2");

    theResponse =  new ElementResponse(this, 9, Vector(3));

  // chord rotation -
  } else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0 
	     || strcmp(argv[0],"basicDeformation") == 0) {

    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta1");
    output.tag("ResponseType","theta2");

    theResponse =  new ElementResponse(this, 3, Vector(3));
  
  // plastic rotation -
  } else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {

    output.tag("ResponseType","epsP");
    output.tag("ResponseType","theta1P");
    output.tag("ResponseType","theta2P");

    theResponse =  new ElementResponse(this, 4, Vector(3));

  // section response -
  } 
    // section response -
  
  else if (strstr(argv[0],"section") != 0) {
    if (argc > 2) {
      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections) {



	theResponse = theSections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
	
	output.endTag();
      }
    }
  }
  

  output.endTag();
  return theResponse;
}

int
DispBeamColumn2dWithSensitivity::getResponse(int responseID, Information &eleInfo)
{
  double V;
  double L = crdTransf->getInitialLength();

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2) {
      P(3) =  q(0);
      P(0) = -q(0)+p0[0];
      P(2) = q(1);
      P(5) = q(2);
      V = (q(1)+q(2))/L;
      P(1) =  V+p0[1];
      P(4) = -V+p0[2];
      return eleInfo.setVector(P);
  }

  else if (responseID == 9) {
    return eleInfo.setVector(q);
  }

  // Chord rotation
  else if (responseID == 3)
    return eleInfo.setVector(crdTransf->getBasicTrialDisp());

  // Plastic rotation
  else if (responseID == 4) {
    static Vector vp(3);
    static Vector ve(3);
    const Matrix &kb = this->getInitialBasicStiff();
    kb.Solve(q, ve);
    vp = crdTransf->getBasicTrialDisp();
    vp -= ve;
    return eleInfo.setVector(vp);
  }

 
  else
    return -1;
}



// AddingSensitivity:BEGIN ///////////////////////////////////
int
DispBeamColumn2dWithSensitivity::setParameter(const char **argv, int argc, Parameter &param)
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
/*  Michael    
    // Get section number: 1...Np
    int sectionNum = atoi(argv[1]);
    
    if (sectionNum > 0 && sectionNum <= numSections)
      return theSections[sectionNum-1]->setParameter(&argv[2], argc-2, param);
    
    else
      return -1;
*/

// Get section and material tag numbers from user input
		int paramSectionTag = atoi(argv[1]);
		int ok=0;
		// Find the right section and call its setParameter method
		for (int i=0; i<numSections; i++) {
			if (paramSectionTag == theSections[i]->getTag()) {
				ok += theSections[i]->setParameter(&argv[2], argc-2, param);
			}
		}
		
		// Check if the ok is valid
		if (ok < 0) {
			opserr << "DispBeamColumn2d::setParameter() - could not set parameter. " << endln;
			return -1;
		}
		else {
			// Return the ok value (according to the above comments)
			return ok;
		}
  }
  
 

  // Default, send to every object
  int ok = 0;
  for (int i = 0; i < numSections; i++)
    ok += theSections[i]->setParameter(argv, argc, param);
  // ok += beamInt->setParameter(argv, argc, param);
  return ok;
}

int
DispBeamColumn2dWithSensitivity::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  else
    return -1; 
}




int
DispBeamColumn2dWithSensitivity::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}





const Matrix &
DispBeamColumn2dWithSensitivity::getKiSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}

const Matrix &
DispBeamColumn2dWithSensitivity::getMassSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}



const Vector &
DispBeamColumn2dWithSensitivity::getResistingForceSensitivity(int gradNumber)
{
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

 
	// Zero for integration
	q.Zero();
	static Vector qsens(3);
	qsens.Zero();

	// Some extra declarations
	static Matrix kbmine(3,3);
	kbmine.Zero();

	int j, k;
	double d1oLdh = 0.0;

	// Check if a nodal coordinate is random
	bool randomNodeCoordinate = false;
	static ID nodeParameterID(2);
	nodeParameterID(0) = theNodes[0]->getCrdsSensitivity();
	nodeParameterID(1) = theNodes[1]->getCrdsSensitivity();
	if (nodeParameterID(0) != 0 || nodeParameterID(1) != 0) {

		randomNodeCoordinate = true;

		const Vector &ndICoords = theNodes[0]->getCrds();
		const Vector &ndJCoords = theNodes[1]->getCrds();

		double dx = ndJCoords(0) - ndICoords(0);
		double dy = ndJCoords(1) - ndICoords(1);

		if (nodeParameterID(0) == 1) // here x1 is random
		  d1oLdh = dx/(L*L*L);
		if (nodeParameterID(0) == 2) // here y1 is random
		  d1oLdh = dy/(L*L*L);

		if (nodeParameterID(1) == 1) // here x2 is random
		  d1oLdh = -dx/(L*L*L);
		if (nodeParameterID(1) == 2) // here y2 is random
		  d1oLdh = -dy/(L*L*L);
	}

	// Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];
    //double wti = wts(i);
    double wti = wt[i];

		// Get section stress resultant gradient
		const Vector &s = theSections[i]->getStressResultant();
		const Vector &sens = theSections[i]->getStressResultantSensitivity(gradNumber,true);

		// Perform numerical integration on internal force gradient
		//q.addMatrixTransposeVector(1.0, *B, s, wts(i));

    double si;
    double sensi;
    for (j = 0; j < order; j++) {
      si = s(j)*wti;
      sensi = sens(j)*wti;
      switch(code(j)) {
	case SECTION_RESPONSE_P:
	q(0) += si;
	qsens(0) += sensi;
	break;
	case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si;
	q(2) += (xi6-2.0)*si;
	qsens(1) += (xi6-4.0)*sensi;
	qsens(2) += (xi6-2.0)*sensi;
	break;
      default:
	break;
  }
  }

  if (randomNodeCoordinate) {


			// Perform numerical integration to obtain basic stiffness matrix
			//kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
			double tmp;

			const Matrix &ks = theSections[i]->getSectionTangent();
			Matrix ka(workArea, order, 3);
			ka.Zero();

			for (j = 0; j < order; j++) {
				switch(code(j)) {
				case SECTION_RESPONSE_P:
					for (k = 0; k < order; k++) {
						ka(k,0) += ks(k,j)*wti;
					}
					break;
				case SECTION_RESPONSE_MZ:
					for (k = 0; k < order; k++) {
						tmp = ks(k,j)*wti;
						ka(k,1) += (xi6-4.0)*tmp;
						ka(k,2) += (xi6-2.0)*tmp;
					}
					break;
				default:
					break;
				}
			}
			for (j = 0; j < order; j++) {
				switch (code(j)) {
				case SECTION_RESPONSE_P:
					for (k = 0; k < 3; k++) {
						kbmine(0,k) += ka(j,k);
					}
					break;
				case SECTION_RESPONSE_MZ:
					for (k = 0; k < 3; k++) {
						tmp = ka(j,k);
						kbmine(1,k) += (xi6-4.0)*tmp;
						kbmine(2,k) += (xi6-2.0)*tmp;
					}
					break;
				default:
					break;
				}
			}
		}

	}

	static Vector dqdh(3);
	const Vector &dAdh_u = crdTransf->getBasicTrialDispShapeSensitivity();
	//dqdh = (1.0/L) * (kbmine * dAdh_u);
	dqdh.addMatrixVector(0.0, kbmine, dAdh_u, oneOverL);

	static Vector dkbdh_v(3);
	const Vector &A_u = crdTransf->getBasicTrialDisp();
	//dkbdh_v = (d1oLdh) * (kbmine * A_u);
	dkbdh_v.addMatrixVector(0.0, kbmine, A_u, d1oLdh);

	// Transform forces
	static Vector dummy(3);		// No distributed loads

	// Term 5
	P = crdTransf->getGlobalResistingForce(qsens,dummy);

	if (randomNodeCoordinate) {

		opserr<< "fatal: DispBeamColumn2dWithSensitivity::getResistingForceSensitivity() random node may not work!"<<endln;

		// Term 1  . following code is not tested yet!

		P += crdTransf->getGlobalResistingForceShapeSensitivity(q,dummy);

		// Term 2
		P += crdTransf->getGlobalResistingForce(dqdh,dummy);

		// Term 4
		P += crdTransf->getGlobalResistingForce(dkbdh_v,dummy);  //*/
	}

//	opserr<<"DispBeamColumn2dWithSensitivity: sens RHS P is:"<<P<<endln;
	return P;
}




// NEW METHOD
int
DispBeamColumn2dWithSensitivity::commitSensitivity(int gradNumber, int numGrads)
{

    // Get basic deformation and sensitivities
	const Vector &v = crdTransf->getBasicTrialDisp();

	static Vector vsens(3);
	vsens = crdTransf->getBasicDisplSensitivity(gradNumber,1);

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);

	// Some extra declarations
	double d1oLdh = 0.0;

	// Check if a nodal coordinate is random
	bool randomNodeCoordinate = false;
	static ID nodeParameterID(2);
	nodeParameterID(0) = theNodes[0]->getCrdsSensitivity();
	nodeParameterID(1) = theNodes[1]->getCrdsSensitivity();
	if (nodeParameterID(0) != 0 || nodeParameterID(1) != 0) {

		vsens += crdTransf->getBasicTrialDispShapeSensitivity();

		randomNodeCoordinate = true;

		const Vector &ndICoords = theNodes[0]->getCrds();
		const Vector &ndJCoords = theNodes[1]->getCrds();

		double dx = ndJCoords(0) - ndICoords(0);
		double dy = ndJCoords(1) - ndICoords(1);

		if (nodeParameterID(0) == 1) // here x1 is random
		  d1oLdh = dx/(L*L*L);
		if (nodeParameterID(0) == 2) // here y1 is random
		  d1oLdh = dy/(L*L*L);

		if (nodeParameterID(1) == 1) // here x2 is random
		  d1oLdh = -dx/(L*L*L);
		if (nodeParameterID(1) == 2) // here y2 is random
		  d1oLdh = -dy/(L*L*L);
	}
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

	

    Vector e(workArea, order);
    
    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];

    for (int j = 0; j < order; j++) {
			switch(code(j)) {
			case SECTION_RESPONSE_P:
				e(j) = oneOverL*vsens(0)
				+ d1oLdh*v(0);
				break;
			case SECTION_RESPONSE_MZ:
				e(j) = oneOverL*((xi6-4.0)*vsens(1) + (xi6-2.0)*vsens(2))
				+ d1oLdh*((xi6-4.0)*v(1) + (xi6-2.0)*v(2));
				break;
			default:
				e(j) = 0.0;
				break;
			}
		}

		// Set the section deformations
		theSections[i]->commitSensitivity(e,gradNumber,numGrads);
	}

	return 0;
}


const Matrix &
DispBeamColumn2dWithSensitivity::getInitialStiffSensitivity(int gradNumber)
{
	//K.Zero();
	//return K;

	// Michele Barbato February 2007
	  const Matrix &kb = this->getInitialBasicStiffSensitivity(gradNumber);

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);

  return K;
  //end addition Barbato
}


const Matrix&
DispBeamColumn2dWithSensitivity::getInitialBasicStiffSensitivity(int gradNumber)
{
  static Matrix kb(3,3);

  // Zero for integral
  kb.Zero();

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  
//  const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
//  const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);



  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    Matrix ka(workArea, order, 3);
    ka.Zero();

//    double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];
    //double wti = wts(i);
    //double wti = wt[i];


    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getInitialTangentSensitivity(gradNumber);

    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    double wti = wt[i];
    double tmp;
    int j, k;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < order; k++)
	  ka(k,0) += ks(k,j)*wti;
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < order; k++) {
	  tmp = ks(k,j)*wti;
	  ka(k,1) += (xi6-4.0)*tmp;
	  ka(k,2) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < 3; k++)
	  kb(0,k) += ka(j,k);
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < 3; k++) {
	  tmp = ka(j,k);
	  kb(1,k) += (xi6-4.0)*tmp;
	  kb(2,k) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }

  }
 //opserr<<"kb:"<<kb<<endln;
  return kb;
}



// AddingSensitivity:END /////////////////////////////////////////////


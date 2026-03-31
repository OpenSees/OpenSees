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
// $Date: 2010-05-13 00:16:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/DispBeamColumn3dThermal.cpp,v $

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for DispBeamColumn3dThermal.

#include <DispBeamColumn3dThermal.h>
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
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <BeamIntegration.h>
#include <Parameter.h>
#include <math.h>
#include <NodalThermalAction.h>
#include <ThermalActionWrapper.h>
#include <FiberSectionGJThermal.h>
#include <FiberSection3dThermal.h>
#include <elementAPI.h>

Matrix DispBeamColumn3dThermal::K(12,12);
Vector DispBeamColumn3dThermal::P(12);
double DispBeamColumn3dThermal::workArea[200];

void* OPS_DispBeamColumn3dThermal()
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

    Element *theEle =  new DispBeamColumn3dThermal(iData[0],iData[1],iData[2],secTags.Size(),sections,
					    *bi,*theTransf,mass);
    delete [] sections;
    return theEle;
}

DispBeamColumn3dThermal::DispBeamColumn3dThermal(int tag, int nd1, int nd2,
				   int numSec, SectionForceDeformation **s,
				   BeamIntegration &bi,
				   CrdTransf &coordTransf, double r)
:Element (tag, ELE_TAG_DispBeamColumn3dThermal),
numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
connectedExternalNodes(2),
Q(12), q(6), rho(r), parameterID(0)
{
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];

  if (theSections == 0) {
    opserr << "DispBeamColumn3dThermal::DispBeamColumn3dThermal - failed to allocate section model pointer\n";
    exit(-1);
  }

  for (int i = 0; i < numSections; i++) {

    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getCopy();

    // Check allocation
    if (theSections[i] == 0) {
      opserr << "DispBeamColumn3dThermal::DispBeamColumn3dThermal -- failed to get a copy of section model\n";
      exit(-1);
    }
  }

  beamInt = bi.getCopy();

  if (beamInt == 0) {
    opserr << "DispBeamColumn3dThermal::DispBeamColumn3dThermal - failed to copy beam integration\n";
    exit(-1);
  }

  crdTransf = coordTransf.getCopy3d();

  if (crdTransf == 0) {
    opserr << "DispBeamColumn3dThermal::DispBeamColumn3dThermal - failed to copy coordinate transformation\n";
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
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  //add by J.Jiang




  residThermal[0] = 0.0;
  residThermal[1] = 0.0;
  residThermal[2] = 0.0;
  residThermal[3] = 0.0;
  residThermal[4] = 0.0;


  counterTemperature = 0;
  AverageThermalElong =0;
  for(int i=0; i<numSections; i++){
   SectionThermalElong[i]=0;
  }




}

DispBeamColumn3dThermal::DispBeamColumn3dThermal()
:Element (0, ELE_TAG_DispBeamColumn3dThermal),
numSections(0), theSections(0), crdTransf(0), beamInt(0),
connectedExternalNodes(2),
Q(12), q(6), rho(0.0), parameterID(0)
{
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  theNodes[0] = 0;
  theNodes[1] = 0;





  residThermal[0] = 0.0;
  residThermal[1] = 0.0;
  residThermal[2] = 0.0;
  residThermal[3] = 0.0;
  residThermal[4] = 0.0;


  counterTemperature = 0;
  AverageThermalElong =0;
  for(int i=0; i<numSections; i++){
   SectionThermalElong[i]=0;
  }


}

DispBeamColumn3dThermal::~DispBeamColumn3dThermal()
{
  for (int i = 0; i < numSections; i++) {
    if (theSections[i])
      delete theSections[i];
  }

  // Delete the array of pointers to SectionForceDeformation pointer arrays
  if (theSections)
    delete [] theSections;

  if (crdTransf)
    delete crdTransf;

  if (beamInt != 0)
    delete beamInt;
  //if (dataMix != 0)//J.Jiang
    //delete [] dataMix;

  //opserr<<"the element is deleted"<<endln;
}

int
DispBeamColumn3dThermal::getNumExternalNodes() const
{
    return 2;
}

const ID&
DispBeamColumn3dThermal::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
DispBeamColumn3dThermal::getNodePtrs()
{

    return theNodes;
}

int
DispBeamColumn3dThermal::getNumDOF()
{
    return 12;
}

void
DispBeamColumn3dThermal::setDomain(Domain *theDomain)
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
	//opserr << "FATAL ERROR DispBeamColumn3dThermal (tag: %d), node not found in domain",
	//	this->getTag());

	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    if (dofNd1 != 6 || dofNd2 != 6) {
	//opserr << "FATAL ERROR DispBeamColumn3dThermal (tag: %d), has differing number of DOFs at its nodes",
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
DispBeamColumn3dThermal::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "DispBeamColumn3dThermal::commitState () - failed in base class";
    }

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    return retVal;
}

int
DispBeamColumn3dThermal::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    return retVal;
}

int
DispBeamColumn3dThermal::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    return retVal;
}

int
DispBeamColumn3dThermal::update(void)
{
  int err = 0;

  // Update the transformation
  crdTransf->update();

  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();
 #ifdef _DEBUG
  if(this->getTag()==601|| this->getTag() == 602)
opserr<<"disp3d Tag: "<<this->getTag()<<" , Displacement: "<<v<<endln;
#endif

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

    double xi6 = 6.0*xi[i];

    int j;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	e(j) = oneOverL*v(0)-AverageThermalElong+SectionThermalElong[i];
	break;
      case SECTION_RESPONSE_MZ:
	e(j) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2));
	break;
      case SECTION_RESPONSE_MY:
	e(j) = oneOverL*((xi6-4.0)*v(3) + (xi6-2.0)*v(4));
	break;
      case SECTION_RESPONSE_T:
	e(j) = oneOverL*v(5);
	break;
      default:
	e(j) = 0.0;
	break;
      }
    }

    // Set the section deformations
#ifdef _DEBUG
	if (this->getTag() == 601 || this->getTag() == 602) {
		opserr << "disp3d Tag: " << this->getTag() << " SECTION "<<i<<"  deformation: " << e << endln;
	}
#endif
    err += theSections[i]->setTrialSectionDeformation(e);
  }

  if (err != 0) {
    opserr << "DispBeamColumn3dThermal::update() - failed setTrialSectionDeformations()\n";
    return err;
  }
  return 0;
}

const Matrix&
DispBeamColumn3dThermal::getTangentStiff()
{
  static Matrix kb(6,6);

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

	  Matrix ka(workArea, order, 6);
	  ka.Zero();

	  double xi6 = 6.0*xi[i];

	  // Get the section tangent stiffness and stress resultant
	  const Matrix &ks = theSections[i]->getSectionTangent();
	  const Vector &s = theSections[i]->getStressResultant();


#ifdef _DEBUG
	  if (this->getTag() == 601 || this->getTag() == 602) {
	  opserr << "DispBeam3d: " << this->getTag() << " Tangent stiff: " << ks << endln;
	  opserr << "DispBeam3d: " << this->getTag() << " Section Stress: " << s << endln;
  }
    #endif
    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
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
      case SECTION_RESPONSE_MY:
	for (k = 0; k < order; k++) {
	  tmp = ks(k,j)*wti;
	  ka(k,3) += (xi6-4.0)*tmp;
	  ka(k,4) += (xi6-2.0)*tmp;
	}
	break;
      case SECTION_RESPONSE_T:
	for (k = 0; k < order; k++)
	  ka(k,5) += ks(k,j)*wti;
	break;
      default:
	break;
      }
    }
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < 6; k++)
	  kb(0,k) += ka(j,k);
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < 6; k++) {
	  tmp = ka(j,k);
	  kb(1,k) += (xi6-4.0)*tmp;
	  kb(2,k) += (xi6-2.0)*tmp;
	}
	break;
      case SECTION_RESPONSE_MY:
	for (k = 0; k < 6; k++) {
	  tmp = ka(j,k);
	  kb(3,k) += (xi6-4.0)*tmp;
	  kb(4,k) += (xi6-2.0)*tmp;
	}
	break;
      case SECTION_RESPONSE_T:
	for (k = 0; k < 6; k++)
	  kb(5,k) += ka(j,k);
	break;
      default:
	break;
      }
    }

    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));
    double si;
    for (j = 0; j < order; j++) {
      si = s(j)*wt[i];
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	q(0) += si;
	break;
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si;
	break;
      case SECTION_RESPONSE_MY:
	q(3) += (xi6-4.0)*si; q(4) += (xi6-2.0)*si;
	break;
      case SECTION_RESPONSE_T:
	q(5) += si;
		  //q(5) += 0.;//J.Jiang add
	break;
      default:
	break;
      }
    }

  }

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

  // Transform to global stiffness
  K = crdTransf->getGlobalStiffMatrix(kb, q);

  #ifdef _bDEBUG
  if(this->getTag() == 601 || this->getTag() == 602){
  opserr<<"DispBeam3d "<<this->getTag()<<" q after tangent: "<<q<<endln;
	  }
  #endif

  return K;
}

const Matrix&
DispBeamColumn3dThermal::getInitialBasicStiff()
{
  static Matrix kb(6,6);

  // Zero for integral
  kb.Zero();

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

    Matrix ka(workArea, order, 6);
    ka.Zero();

    double xi6 = 6.0*xi[i];

    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getInitialTangent();

    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
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
      case SECTION_RESPONSE_MY:
	for (k = 0; k < order; k++) {
	  tmp = ks(k,j)*wti;
	  ka(k,3) += (xi6-4.0)*tmp;
	  ka(k,4) += (xi6-2.0)*tmp;
	}
	break;
      case SECTION_RESPONSE_T:
	for (k = 0; k < order; k++)
	  ka(k,5) += ks(k,j)*wti;
	break;
      default:
	break;
      }
    }
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < 6; k++)
	  kb(0,k) += ka(j,k);
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < 6; k++) {
	  tmp = ka(j,k);
	  kb(1,k) += (xi6-4.0)*tmp;
	  kb(2,k) += (xi6-2.0)*tmp;
	}
	break;
      case SECTION_RESPONSE_MY:
	for (k = 0; k < 6; k++) {
	  tmp = ka(j,k);
	  kb(3,k) += (xi6-4.0)*tmp;
	  kb(4,k) += (xi6-2.0)*tmp;
	}
	break;
      case SECTION_RESPONSE_T:
	for (k = 0; k < 6; k++)
	  kb(5,k) += ka(j,k);
	break;
      default:
	break;
      }
    }

  }

  return kb;
}

const Matrix&
DispBeamColumn3dThermal::getInitialStiff()
{
  const Matrix &kb = this->getInitialBasicStiff();

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);

  return K;
}

const Matrix&
DispBeamColumn3dThermal::getMass()
{
  K.Zero();

  if (rho == 0.0)
    return K;

  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;

  K(0,0) = K(1,1) = K(2,2) = K(6,6) = K(7,7) = K(8,8) = m;

  return K;
}

void
DispBeamColumn3dThermal::zeroLoad(void)
{
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  return;
}

int
DispBeamColumn3dThermal::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = crdTransf->getInitialLength();

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0)*loadFactor;  // Transverse
    double wz = data(1)*loadFactor;  // Transverse
    double wx = data(2)*loadFactor;  // Axial (+ve from node I to J)

    double Vy = 0.5*wy*L;
    double Mz = Vy*L/6.0; // wy*L*L/12
    double Vz = 0.5*wz*L;
    double My = Vz*L/6.0; // wz*L*L/12
    double P = wx*L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= Vy;
    p0[2] -= Vy;
    p0[3] -= Vz;
    p0[4] -= Vz;

    // Fixed end forces in basic system
    q0[0] -= 0.5*P;
    q0[1] -= Mz;
    q0[2] += Mz;
    q0[3] += My;
    q0[4] -= My;
  }
  else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py = data(0)*loadFactor;
    double Pz = data(1)*loadFactor;
    double N  = data(2)*loadFactor;
    double aOverL = data(3);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;
    double b = L-a;

    // Reactions in basic system
    p0[0] -= N;
    double V1, V2;
    V1 = Py*(1.0-aOverL);
    V2 = Py*aOverL;
    p0[1] -= V1;
    p0[2] -= V2;
    V1 = Pz*(1.0-aOverL);
    V2 = Pz*aOverL;
    p0[3] -= V1;
    p0[4] -= V2;

    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    double M1, M2;
    M1 = -a * b2 * Py * L2;
    M2 = a2 * b * Py * L2;
    q0[1] += M1;
    q0[2] += M2;
    M1 = -a * b2 * Pz * L2;
    M2 = a2 * b * Pz * L2;
    q0[3] -= M1;
    q0[4] -= M2;
  }

//J.Jiang add to consider thermal load
else if (type == LOAD_TAG_Beam3dThermalAction) {

	// load not inside fire load pattern
	 //static Vector factors(9);
	 //factors.Zero();
	 //factors += loadFactor;
	// return this->addLoad(theLoad, factors);
   /// This code block is added by LJ and copied from DispBeamColumn2d(Modified edition) for 'FireLoadPattern'--08-May-2012--//[END]
    counterTemperature = 1;
	for(int i=0; i<5; i++){
		residThermal[i]=0;
		}
	double ThermalN, ThermalMz, ThermalMy;
    ThermalN = 0.;
    ThermalMz = 0;
	ThermalMy = 0;

	double xi[maxNumSections];
	beamInt->getSectionLocations(numSections, L, xi);
	double wt[maxNumSections];
	beamInt->getSectionWeights(numSections, L, wt);

	// Zero for integration
	//q.Zero();
	  Vector* dataMixV;
      dataMixV = new Vector(data.Size());

      *dataMixV=data;

#ifdef _DEBUG
	  //if(this->getTag()==1){
	// opserr<<"dataMix in Element"<<this->getTag()<< "___dataMix "<<*dataMixV<<endln;
		//  }
#endif
	// Loop over the integration points
	for (int i = 0; i < numSections; i++) {
		// Get section stress resultant
		const Vector &s = theSections[i]->getTemperatureStress(*dataMixV);
		//opserr<< "Temperature  Stress "<<s<<endln;

		//apply temp along y
		residThermal[0] = -s(0);
		residThermal[1] = -s(1);
		residThermal[2] = s(1);
		residThermal[3] = -s(2);
		residThermal[4] = s(2);
   // opserr<<residThermal[1]<<" "<<residThermal[2]<<"  "<<residThermal[3]<<" "<<residThermal[4]<<endln;
  //apply temp along z
	//residThermal[0] = -s(0);
	//residThermal[1] = -0.;
	//residThermal[2] = -0.;
	//residThermal[3] = -s(1);
	//residThermal[4] = s(1);
	SectionThermalElong[i]=0;
  }
	AverageThermalElong=0;
  }
//Added by Liming for implementation of NodalThermalAction
else if (type == LOAD_TAG_NodalThermalAction) {

	// load not inside fire load pattern
  //static Vector factors(9);
  //factors.Zero();
  //factors += loadFactor;

	AverageThermalElong=0.0;
	for(int i=0; i<5; i++){
		residThermal[i]=0;
		}

	NodalThermalAction* theNodalThermal0 = theNodes[0]->getNodalThermalActionPtr();
	NodalThermalAction* theNodalThermal1 = theNodes[1]->getNodalThermalActionPtr();
	int type;
	const Vector &data0 = theNodalThermal0->getData(type);
	const Vector &data1 = theNodalThermal1->getData(type);
	Vector* Loc;Vector* NodalT0;Vector* NodalT1;

	if(data0.Size()==9){
	 Loc = new Vector(9);
	 NodalT0 = new Vector(9);
	 NodalT1 = new Vector(9);
	 #ifdef _BDEBUG
	 if(this->getTag()==1){
 	 //opserr<<"NodalT0: "<<data0<<endln<<"NodalT1: "<<data1;
		 }
	 #endif
	 for(int i =0; i<9;i++){
		 if(data0(2*i+1)-data1(2*i+1)>1e-8||data0(2*i+1)-data1(2*i+1)<-1e-8){
			 opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
			      << "incompatible loc input for datapoint "<< i << endln;
			 }
		 else{
			 (*Loc)(i)=data0(2*i+1);
			 (*NodalT0)(i)=data0(2*i);
			 (*NodalT1)(i)=data1(2*i);
			 }
	 }
	}
	//for 9 data POINTS
	else {
	 Loc = new Vector(10);
	 NodalT0 = new Vector(15);
	 NodalT1 = new Vector(15);

	for(int i =0; i<5;i++){
	 if(data0(2*i+1)-data1(2*i+1)>1e-8||data0(2*i+1)-data1(2*i+1)<-1e-8){
		opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
			  << "incompatible loc input for datapoint "<< i << endln;
	 }else{
		 //for loc
		(*Loc)(i)=data0(2*i+1);
		(*Loc)(i+5)=data0(3*i+12);
		//for NodalT0
		(*NodalT0)(i)=data0(2*i);
		(*NodalT0)(i+5)=data0(3*i+10);
		(*NodalT0)(i+10)=data0(3*i+11);
		//for NodalT1
		(*NodalT1)(i)=data1(2*i);
		(*NodalT1)(i+5)=data1(3*i+10);
		(*NodalT1)(i+10)=data1(3*i+11);
	 }
	}

	}
	//for 15 data points
 // opserr<<"data0: "<<data0<<endln<<"NodalT0: "<<NodalT0;
	double ThermalN, ThermalMz, ThermalMy;
	ThermalN =  0;
	ThermalMz = 0;
	ThermalMy = 0;

	double xi[maxNumSections];
	beamInt->getSectionLocations(numSections, L, xi);
	double wt[maxNumSections];
	beamInt->getSectionWeights(numSections, L, wt);

	// Loop over the integration points
	for (int i = 0; i < numSections; i++) {
    // Get section stress resultant
	Vector* dataMixV;
	if(NodalT0->Size()==9){
	  for(int m=0;m<9;m++){
		  (*dataMixV)(2*m)=(*NodalT0)(m)+xi[i]*((*NodalT1)(m)-(*NodalT0)(m)); //Linear temperature interpolation
		  (*dataMixV)(2*m+1)=(*Loc)(m);
		  (*dataMixV)(18+m)=1000;
	  }
	}
	else if(NodalT0->Size()==15){
	 for(int m=0;m<5;m++){
		(*dataMixV)(2*m)=(*NodalT0)(m)+xi[i]*((*NodalT1)(m)-(*NodalT0)(m)); ////5 temps through y, Linear temperature interpolation
		(*dataMixV)(3*m+10)=(*NodalT0)(m+5)+xi[i]*((*NodalT1)(m+5)-(*NodalT0)(m+5));///////5 temps through Z in bottom flange
		(*dataMixV)(3*m+11)=(*NodalT0)(m+10)+xi[i]*((*NodalT1)(m+10)-(*NodalT0)(m+10));/////5 temps through Z in top flange
		(*dataMixV)(2*m+1)=(*Loc)(m);////5 (*Loc)s through y
		(*dataMixV)(3*m+12)=(*Loc)(m+5);///////5 (*Loc)s through Z
	 }
	}
	const Vector &s = theSections[i]->getTemperatureStress(*dataMixV);
#ifdef _BDEBUG
	if(this->getTag()==1){
	opserr<<"Tag: "<<this->getTag()<< "Beam3D::dataMixv "<<endln<<dataMixV<<endln;
		}
#endif
	SectionThermalElong[i]= theSections[i]->getThermalElong()(0);
	//opserr<< "Beam3D::Temperature  Stress "<<s<<endln;
	//apply temp along y
	residThermal[0] = -s(0);
	residThermal[1] = -s(1);
	residThermal[2] = s(1);
	residThermal[3] = -s(2);
	residThermal[4] = s(2);
    // opserr<<residThermal[1]<<" "<<residThermal[2]<<"  "<<residThermal[3]<<" "<<residThermal[4]<<endln;
    //apply temp along z
    //residThermal[0] = -s(0);
    //residThermal[1] = -0.;
    //residThermal[2] = -0.;
    //residThermal[3] = -s(1);
    //residThermal[4] = s(1);
	double ThermalEloni = SectionThermalElong[i]*wt[i];
	AverageThermalElong+=ThermalEloni;

  }
	counterTemperature = 1;
}
//----------------------------------------------
  else if(type == LOAD_TAG_ThermalActionWrapper) {
    counterTemperature = 1;
	double ThermalN, ThermalMz, ThermalMy;
	ThermalN =  0;
	ThermalMz = 0;
	ThermalMy = 0;

	AverageThermalElong=0.0;

	double xi[maxNumSections];
	beamInt->getSectionLocations(numSections, L, xi);
	double wt[maxNumSections];
	beamInt->getSectionWeights(numSections, L, wt);

	// Loop over the integration points
    // Get section stress resultant
	 Vector theNode0Crds = theNodes[0]->getCrds();
	 Vector theNode1Crds = theNodes[1]->getCrds();
	 int ndm = theNode0Crds.Size();
	 Vector theIntCrds = Vector(ndm);

	for (int i = 0; i < numSections; i++) {

      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();

      double xi6 = 6.0*xi[i];
      theIntCrds.Zero();
	  for(int m = 0; m<ndm; m++){
			theIntCrds(m) =  theNode0Crds(m)+xi[i]*(theNode1Crds(m)- theNode0Crds(m));
	  }
#ifdef _BDEBUG
	  if(this->getTag()==1){
	 // opserr<<"DispBeamColumn2dThermal::addLoad at the section "<<i<<" --theIntCrds"<<theIntCrds<<endln;
		  }
#endif
	  //NodalThermalActions attached with the wrapper would have been updated by pattern;
	  Vector dataMixV = ((ThermalActionWrapper*) theLoad)->getIntData(theIntCrds);

      const Vector &s = theSections[i]->getTemperatureStress(dataMixV);    //contribuited by ThermalElongation
	  if(theSections[i]->getClassTag()==SEC_TAG_FiberSection3dThermal)
		  SectionThermalElong[i]= ((FiberSection3dThermal*)theSections[i])->getThermalElong()(0);
	  else if(theSections[i]->getClassTag()==SEC_TAG_FiberSectionGJThermal)
		  SectionThermalElong[i]= ((FiberSectionGJThermal*)theSections[i])->getThermalElong()(0);

	  //opserr<<"Thermal Stress "<<s<<endln;
     //apply temp along y
	  residThermal[0] = -s(0);
	  residThermal[1] = -s(1);
	  residThermal[2] = s(1);
	  residThermal[3] = -s(2);
	  residThermal[4] = s(2);
    // opserr<<residThermal[1]<<" "<<residThermal[2]<<"  "<<residThermal[3]<<" "<<residThermal[4]<<endln;
    //apply temp along z
    //residThermal[0] = -s(0);
    //residThermal[1] = -0.;
    //residThermal[2] = -0.;
    //residThermal[3] = -s(1);
    //residThermal[4] = s(1);
	double ThermalEloni = SectionThermalElong[i]*wt[i];
	AverageThermalElong+=ThermalEloni;
    }
	//end of for loop for section
	//counterTemperature = 1;
  }
  else {
    opserr << "DispBeamColumn3dThermal::addLoad() -- load type unknown for element with tag: " <<
      this->getTag() << endln;
    return -1;
  }

  return 0;
}



int
DispBeamColumn3dThermal::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0)
    return 0;

  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);

  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "DispBeamColumn3dThermal::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }

  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;

  // Want to add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  Q(0) -= m*Raccel1(0);
  Q(1) -= m*Raccel1(1);
  Q(2) -= m*Raccel1(2);
  Q(6) -= m*Raccel2(0);
  Q(7) -= m*Raccel2(1);
  Q(8) -= m*Raccel2(2);

  return 0;
}

const Vector&
DispBeamColumn3dThermal::getResistingForce()
{
  double L = crdTransf->getInitialLength();

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Zero for integration
  q.Zero();

 if (counterTemperature == 1) {
	 // to consider material softening
  this->update();
  }

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    double xi6 = 6.0*xi[i];

    // Get section stress resultant
    const Vector &s = theSections[i]->getStressResultant();
#ifdef _DEBUG
	if(this->getTag()==601){
		opserr<<"Resisting force s "<<s<<endln;
		}
#endif
    // Perform numerical integration on internal force
    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));

    double si;
    for (int j = 0; j < order; j++) {
      si = s(j)*wt[i];
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	q(0) += si;
	break;
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si;
	break;
      case SECTION_RESPONSE_MY:
	q(3) += (xi6-4.0)*si; q(4) += (xi6-2.0)*si;
	break;
      case SECTION_RESPONSE_T:
	q(5) += si;
	//q(5) += 0.;//J.Jiang add
	break;
      default:
	break;
      }
    }

  }
  if (counterTemperature == 1)//thermal load at first step
	{
	for ( int i=0; i<5; i++)
    q(i) +=residThermal[i];

	counterTemperature++;
	}


  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

#ifdef _DEBUG
  if (this->getTag() == 601 || this->getTag() == 602) {
	  opserr << "Tag: " << this->getTag() << " DispBeam3d::q  " << q << endln;
	  opserr << " ResidThermal  " << residThermal[0]<<" , "<< residThermal[1]<<" ," << residThermal[2] << endln;
  }


#endif
  // Transform forces
  Vector p0Vec(p0, 5);
  P = crdTransf->getGlobalResistingForce(q, p0Vec);
  //opserr<<"Beam "<<q<<endln;

  // Subtract other external nodal loads ... P_res = P_int - P_ext
   //opserr <<" Ele: "<< this->getTag()<< ", P "<<P<<endln;
  P.addVector(1.0, Q, -1.0);
  return P;
}

const Vector&
DispBeamColumn3dThermal::getResistingForceIncInertia()
{
  this->getResistingForce();

  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();

    // Compute the current resisting force
    this->getResistingForce();

    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;

    P(0) += m*accel1(0);
    P(1) += m*accel1(1);
    P(2) += m*accel1(2);
    P(6) += m*accel2(0);
    P(7) += m*accel2(1);
    P(8) += m*accel2(2);

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P += this->getRayleighDampingForces();

  } else {

    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P += this->getRayleighDampingForces();
  }

  return P;
}

int
DispBeamColumn3dThermal::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j;
  int loc = 0;

  static ID idData(9);  // one bigger than needed so no clash later
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

  idData(7) = beamInt->getClassTag();
  int beamIntDbTag  = beamInt->getDbTag();
  if (beamIntDbTag  == 0) {
    beamIntDbTag = theChannel.getDbTag();
    if (beamIntDbTag  != 0)
      beamInt->setDbTag(beamIntDbTag);
  }
  idData(8) = beamIntDbTag;

  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "DispBeamColumn3dThermal::sendSelf() - failed to send ID data\n";
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
      opserr << "DispBeamColumn3dThermal::sendSelf() - failed to send double data\n";
      return -1;
    }
  }

  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "DispBeamColumn3dThermal::sendSelf() - failed to send crdTranf\n";
     return -1;
  }

  // send the beam integration
  if (beamInt->sendSelf(commitTag, theChannel) < 0) {
    opserr << "DispBeamColumn3dThermal::sendSelf() - failed to send beamInt\n";
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
    opserr << "DispBeamColumn3dThermal::sendSelf() - failed to send ID data\n";
    return -1;
  }

  //
  // send the sections
  //

  for (j = 0; j<numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "DispBeamColumn3dThermal::sendSelf() - section " << j << "failed to send itself\n";
      return -1;
    }
  }

  return 0;
}

int
DispBeamColumn3dThermal::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i;

  static ID idData(9); // one bigger than needed so no clash with section ID

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    opserr << "DispBeamColumn3dThermal::recvSelf() - failed to recv ID data\n";
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
      opserr << "DispBeamColumn3dThermal::sendSelf() - failed to recv double data\n";
      return -1;
    }
    alphaM = dData(0);
    betaK = dData(1);
    betaK0 = dData(2);
    betaKc = dData(3);
  }

  int beamIntClassTag = idData(7);
  int beamIntDbTag = idData(8);

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;

      crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

      if (crdTransf == 0) {
	opserr << "DispBeamColumn3dThermal::recvSelf() - " <<
	  "failed to obtain a CrdTrans object with classTag" <<
	  crdTransfClassTag << endln;
	return -2;
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "DispBeamColumn3dThermal::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }

  // create a new beamInt object if one needed
  if (beamInt == 0 || beamInt->getClassTag() != beamIntClassTag) {
      if (beamInt != 0)
	  delete beamInt;

      beamInt = theBroker.getNewBeamIntegration(beamIntClassTag);

      if (beamInt == 0) {
	opserr << "DispBeamColumn3dThermal::recvSelf() - failed to obtain the beam integration object with classTag" <<
	  beamIntClassTag << endln;
	exit(-1);
      }
  }

  beamInt->setDbTag(beamIntDbTag);

  // invoke recvSelf on the beamInt object
  if (beamInt->recvSelf(commitTag, theChannel, theBroker) < 0)
  {
     opserr << "DispBeamColumn3dThermal::sendSelf() - failed to recv beam integration\n";
     return -3;
  }

  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumn3dThermal::recvSelf() - failed to recv ID data\n";
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
      opserr << "DispBeamColumn3dThermal::recvSelf() - out of memory creating sections array of size" <<
	idData(3) << endln;
      exit(-1);
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
	opserr << "DispBeamColumn3dThermal::recvSelf() - Broker could not create Section of class type" <<
	  sectClassTag << endln;
	exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumn3dThermal::recvSelf() - section " <<
	  i << "failed to recv itself\n";
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
	  opserr << "DispBeamColumn3dThermal::recvSelf() - Broker could not create Section of class type" <<
	    sectClassTag << endln;
	  exit(-1);
	}
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumn3dThermal::recvSelf() - section " <<
	  i << "failed to recv itself\n";
	return -1;
      }
    }
  }

  return 0;
}

void
DispBeamColumn3dThermal::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nDispBeamColumn3dThermal, element id:  " << this->getTag() << endln;
        s << "\tConnected external nodes:  " << connectedExternalNodes;
        s << "\tmass density:  " << rho << endln;

        double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
        double L = crdTransf->getInitialLength();
        double oneOverL = 1.0 / L;

        N = q(0);
        Mz1 = q(1);
        Mz2 = q(2);
        Vy = (Mz1 + Mz2)*oneOverL;
        My1 = q(3);
        My2 = q(4);
        Vz = -(My1 + My2)*oneOverL;
        T = q(5);



        s << "\tEnd 1 Forces (P Mz Vy My Vz T): "
            << -N + p0[0] << ' ' << Mz1 << ' ' << Vy + p0[1] << ' ' << My1 << ' ' << Vz + p0[3] << ' ' << -T << endln;
        s << "\tEnd 2 Forces (P Mz Vy My Vz T): "
            << N << ' ' << Mz2 << ' ' << -Vy + p0[2] << ' ' << My2 << ' ' << -Vz + p0[4] << ' ' << T << endln;

        //for (int i = 0; i < numSections; i++)
        //theSections[i]->Print(s,flag);
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"DispBeamColumn3dThermal\", ";
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
DispBeamColumn3dThermal::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numModes)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response*
DispBeamColumn3dThermal::setResponse(const char **argv, int argc, OPS_Stream &output)
{

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","DispBeamColumn3dThermal");
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


      theResponse = new ElementResponse(this, 1, P);

    // local force -
    }  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

      output.tag("ResponseType","N_1");
      output.tag("ResponseType","Vy_1");
      output.tag("ResponseType","Vz_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Vy_2");
      output.tag("ResponseType","Vz_2");
      output.tag("ResponseType","T_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");

      theResponse = new ElementResponse(this, 2, P);

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

  // section response -
  }
    else if (strstr(argv[0],"sectionX") != 0) {
      if (argc > 2) {
	float sectionLoc = atof(argv[1]);

	double xi[maxNumSections];
	double L = crdTransf->getInitialLength();
	beamInt->getSectionLocations(numSections, L, xi);

	sectionLoc /= L;

	float minDistance = fabs(xi[0]-sectionLoc);
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

	theResponse = theSections[sectionNum]->setResponse(&argv[2], argc-2, output);
      }
    }

    else if (strcmp(argv[0],"section") ==0) {
      if (argc > 1) {

	int sectionNum = atoi(argv[1]);

	if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {

	  double xi[maxNumSections];
	  double L = crdTransf->getInitialLength();
	  beamInt->getSectionLocations(numSections, L, xi);

	  output.tag("GaussPointOutput");
	  output.attr("number",sectionNum);
	  output.attr("eta",xi[sectionNum-1]*L);

	  theResponse =  theSections[sectionNum-1]->setResponse(&argv[2], argc-2, output);

	  output.endTag();
	} else if (sectionNum == 0) { // argv[1] was not an int, we want all sections,

	  CompositeResponse *theCResponse = new CompositeResponse();
	  int numResponse = 0;
	  double xi[maxNumSections];
	  double L = crdTransf->getInitialLength();
	  beamInt->getSectionLocations(numSections, L, xi);

	  for (int i=0; i<numSections; i++) {

	    output.tag("GaussPointOutput");
	    output.attr("number",i+1);
	    output.attr("eta",xi[i]*L);

	    Response *theSectionResponse = theSections[i]->setResponse(&argv[1], argc-1, output);

	    output.endTag();

	    if (theSectionResponse != 0) {
	      numResponse = theCResponse->addResponse(theSectionResponse);
	    }
	  }

	  if (numResponse == 0) // no valid responses found
	    delete theCResponse;
	  else
	    theResponse = theCResponse;
	}
      }
    }

  output.endTag();
  return theResponse;
}

int
DispBeamColumn3dThermal::getResponse(int responseID, Information &eleInfo)
{
  double N, V, M1, M2, T;
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2) {
    // Axial
    N = q(0);
    P(6) =  N;
    P(0) = -N+p0[0];

    // Torsion
    T = q(5);
    P(9) =  T;
    P(3) = -T;

    // Moments about z and shears along y
    M1 = q(1);
    M2 = q(2);
    P(5)  = M1;
    P(11) = M2;
    V = (M1+M2)*oneOverL;
    P(1) =  V+p0[1];
    P(7) = -V+p0[2];

    // Moments about y and shears along z
    M1 = q(3);
    M2 = q(4);
    P(4)  = M1;
    P(10) = M2;
    V = -(M1+M2)*oneOverL;
    P(2) = -V+p0[3];
    P(8) =  V+p0[4];

    return eleInfo.setVector(P);
  }

  // Chord rotation
  else if (responseID == 3)
    return eleInfo.setVector(crdTransf->getBasicTrialDisp());

  // Plastic rotation
  else if (responseID == 4) {
    static Vector vp(6);
    static Vector ve(6);
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
DispBeamColumn3dThermal::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // If the parameter belongs to the element itself
  if (strcmp(argv[0],"rho") == 0) {
    param.setValue(rho);
    return param.addObject(1, this);
  }
  
  if (strstr(argv[0],"sectionX") != 0) {
    if (argc < 3)
		return -1;

	float sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamInt->getSectionLocations(numSections, L, xi);

      sectionLoc /= L;

      float minDistance = fabs(xi[0]-sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
	if (fabs(xi[i]-sectionLoc) < minDistance) {
	  minDistance = fabs(xi[i]-sectionLoc);
	  sectionNum = i;
	}
	  }
	return theSections[sectionNum]->setParameter(&argv[2], argc-2, param);
  }
  // If the parameter belongs to a section or lower
  if (strstr(argv[0],"section") != 0) {

    if (argc < 3)
      return -1;

    // Get section and material tag numbers from user input
    int paramSectionTag = atoi(argv[1]);

    // Find the right section and call its setParameter method
    int ok = 0;
    for (int i = 0; i < numSections; i++)
      if (paramSectionTag == theSections[i]->getTag())
	ok += theSections[i]->setParameter(&argv[2], argc-2, param);

    return ok;
  }

  if (strstr(argv[0],"integration") != 0) {

    if (argc < 2)
      return -1;

    return beamInt->setParameter(&argv[1], argc-1, param);
  }

  // Default, send to every object
  int ok = 0;
  int result = -1;

  for (int i = 0; i < numSections; i++) {
    ok = theSections[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }
  
  ok = beamInt->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;

  return result;
}

int
DispBeamColumn3dThermal::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  else
    return -1;
}




int
DispBeamColumn3dThermal::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;

  return 0;
}

const Matrix &
DispBeamColumn3dThermal::getKiSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}

const Matrix &
DispBeamColumn3dThermal::getMassSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}



const Vector &
DispBeamColumn3dThermal::getResistingForceSensitivity(int gradNumber)
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
  static Vector dqdh(6);
  dqdh.Zero();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];
    //double wti = wts(i);
    double wti = wt[i];

    // Get section stress resultant gradient
    const Vector &dsdh = theSections[i]->getStressResultantSensitivity(gradNumber,true);

    // Perform numerical integration on internal force gradient
    double sensi;
    for (int j = 0; j < order; j++) {
      sensi = dsdh(j)*wti;
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	dqdh(0) += sensi;
	break;
      case SECTION_RESPONSE_MZ:
	dqdh(1) += (xi6-4.0)*sensi;
	dqdh(2) += (xi6-2.0)*sensi;
	break;
      case SECTION_RESPONSE_MY:
	dqdh(3) += (xi6-4.0)*sensi;
	dqdh(4) += (xi6-2.0)*sensi;
	break;
      case SECTION_RESPONSE_T:
	dqdh(5) += sensi;
	break;
      default:
	break;
      }
    }
  }

  // Transform forces
  static Vector dp0dh(6);		// No distributed loads

  P.Zero();

  //////////////////////////////////////////////////////////////

  if (crdTransf->isShapeSensitivity()) {

    // Perform numerical integration to obtain basic stiffness matrix
    // Some extra declarations
    static Matrix kbmine(6,6);
    kbmine.Zero();
    q.Zero();

    double tmp;

    int j, k;

    for (int i = 0; i < numSections; i++) {

      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();

      //double xi6 = 6.0*pts(i,0);
      double xi6 = 6.0*xi[i];
      //double wti = wts(i);
      double wti = wt[i];

      const Vector &s = theSections[i]->getStressResultant();
      const Matrix &ks = theSections[i]->getSectionTangent();

      Matrix ka(workArea, order, 6);
      ka.Zero();

      double si;
      for (j = 0; j < order; j++) {
	si = s(j)*wti;
	switch(code(j)) {
	case SECTION_RESPONSE_P:
	  q(0) += si;
	  for (k = 0; k < order; k++) {
	    ka(k,0) += ks(k,j)*wti;
	  }
	  break;
	case SECTION_RESPONSE_MZ:
	  q(1) += (xi6-4.0)*si;
	  q(2) += (xi6-2.0)*si;
	  for (k = 0; k < order; k++) {
	    tmp = ks(k,j)*wti;
	    ka(k,1) += (xi6-4.0)*tmp;
	    ka(k,2) += (xi6-2.0)*tmp;
	  }
	  break;
	case SECTION_RESPONSE_MY:
	  q(3) += (xi6-4.0)*si;
	  q(4) += (xi6-2.0)*si;
	  for (k = 0; k < order; k++) {
	    tmp = ks(k,j)*wti;
	    ka(k,3) += (xi6-4.0)*tmp;
	    ka(k,4) += (xi6-2.0)*tmp;
	  }
	  break;
	case SECTION_RESPONSE_T:
	  q(5) += si;
	  for (k = 0; k < order; k++) {
	    ka(k,5) += ks(k,j)*wti;
	  }
	  break;
	default:
	  break;
	}
      }
      for (j = 0; j < order; j++) {
	switch (code(j)) {
	case SECTION_RESPONSE_P:
	  for (k = 0; k < 6; k++) {
	    kbmine(0,k) += ka(j,k);
	  }
	  break;
	case SECTION_RESPONSE_MZ:
	  for (k = 0; k < 6; k++) {
	    tmp = ka(j,k);
	    kbmine(1,k) += (xi6-4.0)*tmp;
	    kbmine(2,k) += (xi6-2.0)*tmp;
	  }
	  break;
	case SECTION_RESPONSE_MY:
	  for (k = 0; k < 6; k++) {
	    tmp = ka(j,k);
	    kbmine(3,k) += (xi6-4.0)*tmp;
	    kbmine(4,k) += (xi6-2.0)*tmp;
	  }
	  break;
	case SECTION_RESPONSE_T:
	  for (k = 0; k < 6; k++) {
	    kbmine(5,k) += ka(j,k);
	  }
	  break;
	default:
	  break;
	}
      }
    }

    const Vector &A_u = crdTransf->getBasicTrialDisp();
    double dLdh = crdTransf->getdLdh();
    double d1overLdh = -dLdh/(L*L);
    // a^T k_s dadh v
    dqdh.addMatrixVector(1.0, kbmine, A_u, d1overLdh);

    // k dAdh u
    const Vector &dAdh_u = crdTransf->getBasicTrialDispShapeSensitivity();
    dqdh.addMatrixVector(1.0, kbmine, dAdh_u, oneOverL);

    // dAdh^T q
    P += crdTransf->getGlobalResistingForceShapeSensitivity(q, dp0dh, gradNumber);
  }

  // A^T (dqdh + k dAdh u)
  P += crdTransf->getGlobalResistingForce(dqdh, dp0dh);

  return P;
}



// NEW METHOD
int
DispBeamColumn3dThermal::commitSensitivity(int gradNumber, int numGrads)
{
  // Get basic deformation and sensitivities
  const Vector &v = crdTransf->getBasicTrialDisp();

  static Vector dvdh(6);
  dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);

  // Some extra declarations
  double d1oLdh = crdTransf->getd1overLdh();

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
	e(j) = oneOverL*dvdh(0)
	  + d1oLdh*v(0);
	break;
      case SECTION_RESPONSE_MZ:
	e(j) = oneOverL*((xi6-4.0)*dvdh(1) + (xi6-2.0)*dvdh(2))
	  + d1oLdh*((xi6-4.0)*v(1) + (xi6-2.0)*v(2));
	break;
      case SECTION_RESPONSE_MY:
	e(j) = oneOverL*((xi6-4.0)*dvdh(3) + (xi6-2.0)*dvdh(4))
	  + d1oLdh*((xi6-4.0)*v(3) + (xi6-2.0)*v(4));
	break;
      case SECTION_RESPONSE_T:
	e(j) = oneOverL*dvdh(5)
	  + d1oLdh*v(5);
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


// AddingSensitivity:END /////////////////////////////////////////////

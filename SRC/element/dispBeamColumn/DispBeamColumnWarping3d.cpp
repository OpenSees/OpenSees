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
// $Date: 2008/11/04 21:32:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/DispBeamColumnWarping3d.cpp,v $

// Modified by Xi Zhang from University of Sydney, Australia (include warping degrees of freedom). Refer to 
// Formulation and Implementation of Three-dimensional Doubly Symmetric Beam-Column Analyses with Warping Effects in OpenSees
// Research Report R917, School of Civil Engineering, University of Sydney.
// Description: This file contains the class definition for DispBeamColumnWarping3d which include warping degrees of freedom.
// Modified by Xi Zhang from University of Sydney, Australia (include warping degrees of freedom). Refer to 
// Formulation and Implementation of Three-dimensional Doubly Symmetric Beam-Column Analyses with Warping Effects in OpenSees
// Research Report R917, School of Civil Engineering, University of Sydney.

#include <DispBeamColumnWarping3d.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <elementAPI.h>
#include <string.h>
#include <math.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <BeamIntegration.h>
#include <Parameter.h>
using std::string;
using namespace std;


Matrix DispBeamColumnWarping3d::K(14,14);
Vector DispBeamColumnWarping3d::P(14);
double DispBeamColumnWarping3d::workArea[200];

void* OPS_DispBeamColumnWarping3d()
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
    int cmass = 0;
    numData = 1;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type,"-cMass") == 0) {
	    cmass = 1;
	} else if(strcmp(type,"-mass") == 0) {
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
    
    Element *theEle =  new DispBeamColumnWarping3d(iData[0],iData[1],iData[2],secTags.Size(),sections,
					    *bi,*theTransf,mass);
    delete [] sections;
    return theEle;
}

DispBeamColumnWarping3d::DispBeamColumnWarping3d(int tag, int nd1, int nd2,
				   int numSec, SectionForceDeformation **s,
				   BeamIntegration &bi,
				   CrdTransf &coordTransf, double r)
:Element (tag, ELE_TAG_DispBeamColumnWarping3d),
numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
connectedExternalNodes(2), 
Q(14), q(9), rho(r), parameterID(0)
{
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];
  
  if (theSections == 0) {
    opserr << "DispBeamColumnWarping3d::DispBeamColumnWarping3d - failed to allocate section model pointer\n";
    exit(-1);
  }
  
  for (int i = 0; i < numSections; i++) {
    
    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getCopy();
    
    // Check allocation
    if (theSections[i] == 0) {
      opserr << "DispBeamColumnWarping3d::DispBeamColumnWarping3d -- failed to get a copy of section model\n";
      exit(-1);
    }
  }
  
  beamInt = bi.getCopy();
  
  if (beamInt == 0) {
    opserr << "DispBeamColumnWarping3d::DispBeamColumnWarping3d - failed to copy beam integration\n";
    exit(-1);
  }

  crdTransf = coordTransf.getCopy3d();
  
  if (crdTransf == 0) {
    opserr << "DispBeamColumnWarping3d::DispBeamColumnWarping3d - failed to copy coordinate transformation\n";
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
  q0[5] = 0.0;
  q0[6] = 0.0;
  q0[7] = 0.0;
  q0[8] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;
}

DispBeamColumnWarping3d::DispBeamColumnWarping3d()
:Element (0, ELE_TAG_DispBeamColumnWarping3d),
numSections(0), theSections(0), crdTransf(0), beamInt(0),
connectedExternalNodes(2), 
Q(14), q(9), rho(0.0), parameterID(0)
{
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;
  q0[5] = 0.0;
  q0[6] = 0.0;
  q0[7] = 0.0;
  q0[8] = 0.0;


  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  theNodes[0] = 0;
  theNodes[1] = 0;
}

DispBeamColumnWarping3d::~DispBeamColumnWarping3d()
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
}

int
DispBeamColumnWarping3d::getNumExternalNodes() const
{
    return 2;
}

const ID&
DispBeamColumnWarping3d::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
DispBeamColumnWarping3d::getNodePtrs()
{

    return theNodes;
}

int
DispBeamColumnWarping3d::getNumDOF()
{
    return 14;
}

void
DispBeamColumnWarping3d::setDomain(Domain *theDomain)
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
      opserr << "FATAL ERROR DispBeamColumnWarping3d (tag: %d), node not found in domain",
	this->getTag();
	
	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    // release this restriction
	/*if (dofNd1 != 6 || dofNd2 != 6) {
	//opserr << "FATAL ERROR DispBeamColumnWarping3d (tag: %d), has differing number of DOFs at its nodes",
	//	this->getTag());
	
	return;
    }*/

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
DispBeamColumnWarping3d::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "DispBeamColumnWarping3d::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
      retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    return retVal;
}

int
DispBeamColumnWarping3d::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
      retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    return retVal;
}

int
DispBeamColumnWarping3d::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
      retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    return retVal;
}

int
DispBeamColumnWarping3d::update(void)
{
  int err = 0;

  // Update the transformation
  crdTransf->update();
  
  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();
 
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  double oneOverLsquare = oneOverL*oneOverL;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    Vector e(workArea, 8);
      
    double xi6 = 6.0*xi[i];
    double xi12 = 12.0*xi[i];
    double xi1 = xi[i];
    
    int j;

    // here total strain (not incremental) is used
    e(0) = oneOverL*v(8);
    e(1) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(5));
    e(2) = oneOverL*((xi6-4.0)*v(2) + (xi6-2.0)*v(6));
    e(4) =oneOverL*(6.0*xi1*xi1-6.0*xi1)*v(0)+(1.0-4.0*xi1+3*xi1*xi1)*v(3)+oneOverL*(6.0*xi1-6.0*xi1*xi1)*v(4)+(3.0*xi1*xi1-2.0*xi1)*v(7);
    e(3) = -oneOverLsquare*(xi12-6.0)*v(0)-oneOverL*(xi6-4.0)*v(3)-oneOverLsquare*(6.0-xi12)*v(4)-oneOverL*(xi6-2.0)*v(7);
    e(5) = (1.0+3.0*xi1*xi1-4.0*xi1)*v(1)+(3.0*xi1*xi1-2.0*xi1)*v(5);
    e(6) = (1.0+3.0*xi1*xi1-4.0*xi1)*v(2)+(3.0*xi1*xi1-2.0*xi1)*v(6);
    e(7) = (1.0-3.0*xi1*xi1+2.0*xi1*xi1*xi1)*v(0)+L*xi1*(1.0-xi1)*(1.0-xi1)*v(3)+(3.0*xi1*xi1-2.0*xi1*xi1*xi1)*v(4)+L*xi1*xi1*(xi1-1.0)*v(7);
	
    // Set the section deformations
    err += theSections[i]->setTrialSectionDeformation(e);
  }

  if (err != 0) {
    opserr << "DispBeamColumnWarping3d::update() - failed setTrialSectionDeformations()\n";
    return err;
  }
  
  return 0;
}

const Matrix&
DispBeamColumnWarping3d::getTangentStiff()
{
  static Matrix kb(9,9);
  static Matrix N1(6,8);
  static Matrix N2(8,9);
  static Matrix N3(8,8);
  static Matrix kbPart1(9,9);
  static Matrix Gmax(8,8);
  static Matrix kbPart2(9,9);
  const Vector &v = crdTransf->getBasicTrialDisp();
  
  // Zero for integral
  kb.Zero();

  q.Zero();
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  double oneOverLsquare = oneOverL*oneOverL;
  double oneOverLcube = oneOverLsquare*oneOverL;

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

    N1.Zero();
    N2.Zero();
    N3.Zero();
    kbPart1.Zero();
    Gmax.Zero();
    kbPart2.Zero();
    
    double xi6 = 6.0*xi[i];
    double xi12 = 12.0*xi[i];
    double xi1 = xi[i];
    double dNv1 = 1.0+3.0*xi1*xi1-4.0*xi1;
    double dNv2 = 3.0*xi1*xi1-2.0*xi1;
    double dNw1 = dNv1;
    double dNw2 = dNv2;
    double ddNv1 = 6.0*xi1*oneOverL-4.0*oneOverL;
    double ddNv2 = 6.0*xi1*oneOverL-2.0*oneOverL;
    double ddNw1 = ddNv1;
    double ddNw2 = ddNv2;
    double Nf1 = 1.0-3.0*xi1*xi1+2.0*xi1*xi1*xi1;
    double Nf2 = xi1*L*(1.0-xi1)*(1.0-xi1);
    double Nf3 = 3.0*xi1*xi1-2.0*xi1*xi1*xi1;
    double Nf4 = xi1*xi1*L*(xi1-1.0);
    double dNf1 = 6.0*xi1*xi1*oneOverL-6.0*xi1*oneOverL;
    double dNf2 = dNv1;
    double dNf3 = 6.0*xi1*oneOverL-6.0*xi1*xi1*oneOverL;
    double dNf4 = dNv2;
    double ddNf1 = 12.0*xi1*oneOverLsquare-6.0*oneOverLsquare;
    double ddNf2 = ddNv1;
    double ddNf3 = 6.0*oneOverLsquare-12.0*xi1*oneOverLsquare;
    double ddNf4 = ddNv2;
    N1(0,0) = 1.0;
    N1(0,1) = dNv1*v(1)+dNv2*v(5);
    N1(0,2) = dNw1*v(2)+dNw2*v(6);
    N1(1,3) = 1.0;
    N1(1,4) = Nf1*v(0)+Nf2*v(3)+Nf3*v(4)+Nf4*v(7);
    N1(1,5) = ddNw1*v(2)+ddNw2*v(6);
    N1(2,3) = -N1(1,4);
    N1(2,4) = 1.0;
    N1(2,5) = -ddNv1*v(1)-ddNv2*v(5);
    N1(3,6) = dNf1*v(0)+dNf2*v(3)+dNf3*v(4)+dNf4*v(7);
    N1(4,7) = -1.0;
    N1(5,6) = 1.0;
    N2(0,8) = oneOverL;
    N2(1,1) = dNv1;
    N2(1,5) = dNv2;
    N2(2,2) = dNw1;
    N2(2,6) = dNw2;
    N2(3,1) = ddNv1;
    N2(3,5) = ddNv2;
    N2(4,2) = ddNw1;
    N2(4,6) = ddNw2;
    N2(5,0) = Nf1;
    N2(5,3) = Nf2;
    N2(5,4) = Nf3;
    N2(5,7) = Nf4;
    N2(6,0) = dNf1;
    N2(6,3) = dNf2;
    N2(6,4) = dNf3;
    N2(6,7) = dNf4;
    N2(7,0) = ddNf1;
    N2(7,3) = ddNf2;
    N2(7,4) = ddNf3;
    N2(7,7) = ddNf4;
    
    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getSectionTangent();

    //calculate kb, refer to Alemdar

    N3.addMatrixTripleProduct(0.0, N1, ks, 1.0);
    kbPart1.addMatrixTripleProduct(0.0, N2, N3, 1.0);
    const Vector &s = theSections[i]->getStressResultant();
    
    Gmax(1,1) = s(0);
    Gmax(2,2) = s(0);
    Gmax(3,5) = s(1);
    Gmax(4,5) = s(2);
    Gmax(5,3) = s(2);
    Gmax(5,4) = s(1);
    Gmax(6,6) = s(3);
    kbPart2.addMatrixTripleProduct(0.0, N2, Gmax, 1.0);
    
    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    double wti = wt[i];
	
    for (int j=0; j<9; j++)
      {
	for (int k=0; k<9; k++)
	  {
	    kb(j,k) +=kbPart1(j,k)*L*wti + kbPart2(j,k)*L*wti;
	  }
      }
   
    static Vector qProduct1(8);
    static Vector qProduct2(9);
    qProduct1.Zero();
    qProduct2.Zero();
    qProduct1.addMatrixTransposeVector(0.0, N1, s, 1.0);
    qProduct2.addMatrixTransposeVector(0.0, N2, qProduct1, 1.0);
    
    for (int j=0; j<9; j++)
      {
	q(j) += qProduct2(j)*L*wti;
      }
  }
  
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

  // Transform to global stiffness
  K = crdTransf->getGlobalStiffMatrix(kb, q);
  
  return K;
}

const Matrix&
DispBeamColumnWarping3d::getInitialBasicStiff()
{
  static Matrix kb(9,9);
  static Matrix N1(6,8);
  static Matrix N2(8,9);
  static Matrix N3(8,8);
  static Matrix kbPart1(9,9);
  static Matrix Gmax(8,8);
  static Matrix kbPart2(9,9);
  
  // Zero for integral
  kb.Zero();
  const Vector &v = crdTransf->getBasicTrialDisp();
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  double oneOverLsquare = oneOverL*oneOverL;
  double oneOverLcube = oneOverLsquare*oneOverL;

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
    
    N1.Zero();
    N2.Zero();
    N3.Zero();
    kbPart1.Zero();
    Gmax.Zero();
    kbPart2.Zero();
    
    
    double xi6 = 6.0*xi[i];
    double xi12=12*xi[i];
    double xi1=xi[i];
    double dNv1 = 1.0+3.0*xi1*xi1-4.0*xi1;
    double dNv2 = 3.0*xi1*xi1-2.0*xi1;
    double dNw1 = dNv1;
    double dNw2 = dNv2;
    double ddNv1 = 6.0*xi1*oneOverL-4.0*oneOverL;
    double ddNv2 = 6.0*xi1*oneOverL-2.0*oneOverL;
    double ddNw1 = ddNv1;
    double ddNw2 = ddNv2;
    double Nf1 = 1.0-3.0*xi1*xi1+2.0*xi1*xi1*xi1;
    double Nf2 = xi1*L*(1.0-xi1)*(1.0-xi1);
    double Nf3 = 3.0*xi1*xi1-2.0*xi1*xi1*xi1;
    double Nf4 = xi1*xi1*L*(xi1-1.0);
    double dNf1 = 6.0*xi1*xi1*oneOverL-6.0*xi1*oneOverL;
    double dNf2 = dNv1;
    double dNf3 = 6.0*xi1*oneOverL-6.0*xi1*xi1*oneOverL;
    double dNf4 = dNv2;
    double ddNf1 = 12.0*xi1*oneOverLsquare-6.0*oneOverLsquare;
    double ddNf2 = ddNv1;
    double ddNf3 = 6.0*oneOverLsquare-12.0*xi1*oneOverLsquare;
    double ddNf4 = ddNv2;
    N1(0,0) = 1.0;
    N1(0,1) = dNv1*v(1)+dNv2*v(5);
    N1(0,2) = dNw1*v(2)+dNw2*v(6);
    N1(1,3) = 1.0;
    N1(1,4) = Nf1*v(0)+Nf2*v(3)+Nf3*v(4)+Nf4*v(7);
    N1(1,5) = ddNw1*v(2)+ddNw2*v(6);
    N1(2,3) = -N1(1,4);
    N1(2,4) = 1.0;
    N1(2,5) = -ddNv1*v(1)-ddNv2*v(5);
    N1(3,6) = dNf1*v(0)+dNf2*v(3)+dNf3*v(4)+dNf4*v(7);
    N1(4,7) = -1.0;
    N1(5,6) = 1.0;
    N2(0,8) = oneOverL;
    N2(1,1) = dNv1;
    N2(1,5) = dNv2;
    N2(2,2) = dNw1;
    N2(2,6) = dNw2;
    N2(3,1) = ddNv1;
    N2(3,5) = ddNv2;
    N2(4,2) = ddNw1;
    N2(4,6) = ddNw2;
    N2(5,0) = Nf1;
    N2(5,3) = Nf2;
    N2(5,4) = Nf3;
    N2(5,7) = Nf4;
    N2(6,0) = dNf1;
    N2(6,3) = dNf2;
    N2(6,4) = dNf3;
    N2(6,7) = dNf4;
    N2(7,0) = ddNf1;
    N2(7,3) = ddNf2;
    N2(7,4) = ddNf3;
    N2(7,7) = ddNf4;
    
    
    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getInitialTangent();
    N3.addMatrixTripleProduct(0.0, N1, ks, 1.0);
    kbPart1.addMatrixTripleProduct(0.0, N2, N3, 1.0);
    const Vector &s = theSections[i]->getStressResultant();
    Gmax(1,1) = s(0);
    Gmax(2,2) = s(0);
    Gmax(3,5) = s(1);
    Gmax(4,5) = s(2);
    Gmax(5,3) = s(2);
    Gmax(5,4) = s(1);
    Gmax(6,6) = s(3);
    kbPart2.addMatrixTripleProduct(0.0, N2, Gmax, 1.0);
    //opserr<<"s"<<s;
    
    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    double wti = wt[i];
    
    for (int j=0; j<9; j++)
      {
	for (int k=0; k<9; k++)
	  {
	    kb(j,k) +=kbPart1(j,k)*L*wti + kbPart2(j,k)*L*wti;
	  }
      }
    
	/*kb(0,0) += (ks(3,3)*(xi12-6.0)*(xi12-6.0)*oneOverLcube+ks(4,4)*oneOverL*(6*xi1*xi1-xi6)*(6*xi1*xi1-xi6))*wti;
	kb(0,3) += (ks(3,3)*(4.0-xi6)*(6.0-xi12)*oneOverLsquare+ks(4,4)*(6.0*xi1*xi1-xi6)*(1.0-4.0*xi1+3.0*xi1*xi1))*wti;
	kb(3,0) += (ks(3,3)*(4.0-xi6)*(6.0-xi12)*oneOverLsquare+ks(4,4)*(6.0*xi1*xi1-xi6)*(1.0-4.0*xi1+3.0*xi1*xi1))*wti;
	kb(0,4) += -(ks(3,3)*(xi12-6.0)*(xi12-6.0)*oneOverLcube+ks(4,4)*oneOverL*(6*xi1*xi1-xi6)*(6*xi1*xi1-xi6))*wti;
	kb(4,0) += -(ks(3,3)*(xi12-6.0)*(xi12-6.0)*oneOverLcube+ks(4,4)*oneOverL*(6*xi1*xi1-xi6)*(6*xi1*xi1-xi6))*wti;
	kb(0,7) += (ks(3,3)*(xi12-6.0)*(xi6-2.0)*oneOverLsquare+ks(4,4)*(6.0*xi1*xi1-xi6)*(3.0*xi1*xi1-2.0*xi1))*wti;
	kb(7,0) += (ks(3,3)*(xi12-6.0)*(xi6-2.0)*oneOverLsquare+ks(4,4)*(6.0*xi1*xi1-xi6)*(3.0*xi1*xi1-2.0*xi1))*wti;
	kb(1,1) += ks(1,1)*(xi6-4.0)*(xi6-4.0)*oneOverL*wti;
	kb(1,5) += ks(1,1)*(xi6-4.0)*(xi6-2.0)*oneOverL*wti;
	kb(5,1) += ks(1,1)*(xi6-4.0)*(xi6-2.0)*oneOverL*wti;
	kb(2,2) += ks(2,2)*(xi6-4.0)*(xi6-4.0)*oneOverL*wti;
	kb(2,6) += ks(2,2)*(xi6-4.0)*(xi6-2.0)*oneOverL*wti;
	kb(6,2) += ks(2,2)*(xi6-4.0)*(xi6-2.0)*oneOverL*wti;
	kb(3,3) += (ks(3,3)*(xi6-4.0)*(xi6-4.0)*oneOverL+ks(4,4)*L*(1.0-4.0*xi1+3.0*xi1*xi1)*(1.0-4.0*xi1+3.0*xi1*xi1))*wti;
	kb(3,4) += -(ks(3,3)*(4.0-xi6)*(6.0-xi12)*oneOverLsquare+ks(4,4)*(6.0*xi1*xi1-xi6)*(1.0-4.0*xi1+3.0*xi1*xi1))*wti;
	kb(4,3) += -(ks(3,3)*(4.0-xi6)*(6.0-xi12)*oneOverLsquare+ks(4,4)*(6.0*xi1*xi1-xi6)*(1.0-4.0*xi1+3.0*xi1*xi1))*wti;
	kb(3,7) += (ks(3,3)*(xi6-4.0)*(xi6-2.0)*oneOverL+ks(4,4)*L*(1.0-4.0*xi1+3.0*xi1*xi1)*(3.0*xi1*xi1-2.0*xi1))*wti;
	kb(4,4) += (ks(3,3)*(xi12-6.0)*(xi12-6.0)*oneOverLcube+ks(4,4)*(xi6-6.0*xi1*xi1)*(xi6-6.0*xi1*xi1)*oneOverL)*wti;
	kb(4,7) += (ks(3,3)*(6.0-xi12)*(xi6-2.0)*oneOverLsquare+ks(4,4)*(xi6-6.0*xi1*xi1)*(3.0*xi1*xi1-2.0*xi1))*wti;
	kb(7,4) += (ks(3,3)*(6.0-xi12)*(xi6-2.0)*oneOverLsquare+ks(4,4)*(xi6-6.0*xi1*xi1)*(3.0*xi1*xi1-2.0*xi1))*wti;
	kb(5,5) += ks(1,1)*(xi6-2.0)*(xi6-2.0)*oneOverL*wti;
	kb(6,6) += ks(2,2)*(xi6-2.0)*(xi6-2.0)*oneOverL*wti;
	kb(7,3) += (ks(3,3)*(xi6-4.0)*(xi6-2.0)*oneOverL+ks(4,4)*L*(1.0-4.0*xi1+3.0*xi1*xi1)*(3.0*xi1*xi1-2.0*xi1))*wti;
	kb(7,7) += (ks(3,3)*(xi6-2.0)*(xi6-2.0)*oneOverL+ks(4,4)*L*(2.0*xi1-3.0*xi1*xi1)*(2.0*xi1-3.0*xi1*xi1))*wti;
	kb(8,8) += ks(0,0)*oneOverL*wti;*/
  }
  
  return kb;
}

const Matrix&
DispBeamColumnWarping3d::getInitialStiff()
{
  const Matrix &kb = this->getInitialBasicStiff();

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);
  
  return K;
}

const Matrix&
DispBeamColumnWarping3d::getMass()
{
  K.Zero();
  
  if (rho == 0.0)
    return K;
  
  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;
  
  K(0,0) = K(1,1) = K(2,2) = K(7,7) = K(8,8) = K(9,9) = m;
  
  return K;
}

void
DispBeamColumnWarping3d::zeroLoad(void)
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
DispBeamColumnWarping3d::addLoad(ElementalLoad *theLoad, double loadFactor)
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
  else {
    opserr << "DispBeamColumnWarping3d::addLoad() -- load type unknown for element with tag: " << 
      this->getTag() << endln;
    return -1;
  }

  return 0;
}

int 
DispBeamColumnWarping3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0) 
    return 0;
  
  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  if (7 != Raccel1.Size() || 7 != Raccel2.Size()) {
    opserr << "DispBeamColumnWarping3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
    return -1;
  }
  
  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;
  
  // Want to add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  Q(0) -= m*Raccel1(0);
  Q(1) -= m*Raccel1(1);
  Q(2) -= m*Raccel1(2);
  Q(7) -= m*Raccel2(0);
  Q(8) -= m*Raccel2(1);
  Q(9) -= m*Raccel2(2);
  
  return 0;
}

const Vector&
DispBeamColumnWarping3d::getResistingForce()
{
  double L = crdTransf->getInitialLength();
  double oneOverL=1.0/L;
  double oneOverLsquare = oneOverL*oneOverL;
  double oneOverLcube = oneOverLsquare*oneOverL;
  static Matrix N1(6,8);
  static Matrix N2(8,9);
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);
  const Vector &v = crdTransf->getBasicTrialDisp();
  // Zero for integration
  q.Zero();
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    N1.Zero();
    N2.Zero();
    
    double xi6 = 6.0*xi[i];
    double xi12=12*xi[i];
    double xi1=xi[i];
    double dNv1 = 1.0+3.0*xi1*xi1-4.0*xi1;
    double dNv2 = 3.0*xi1*xi1-2.0*xi1;
    double dNw1 = dNv1;
    double dNw2 = dNv2;
    double ddNv1 = 6.0*xi1*oneOverL-4.0*oneOverL;
    double ddNv2 = 6.0*xi1*oneOverL-2.0*oneOverL;
    double ddNw1 = ddNv1;
    double ddNw2 = ddNv2;
    double Nf1 = 1.0-3.0*xi1*xi1+2.0*xi1*xi1*xi1;
    double Nf2 = xi1*L*(1.0-xi1)*(1.0-xi1);
    double Nf3 = 3.0*xi1*xi1-2.0*xi1*xi1*xi1;
    double Nf4 = xi1*xi1*L*(xi1-1.0);
    double dNf1 = 6.0*xi1*xi1*oneOverL-6.0*xi1*oneOverL;
    double dNf2 = dNv1;
    double dNf3 = 6.0*xi1*oneOverL-6.0*xi1*xi1*oneOverL;
    double dNf4 = dNv2;
    double ddNf1 = 12.0*xi1*oneOverLsquare-6.0*oneOverLsquare;
    double ddNf2 = ddNv1;
    double ddNf3 = 6.0*oneOverLsquare-12.0*xi1*oneOverLsquare;
    double ddNf4 = ddNv2;
    N1(0,0) = 1.0;
    N1(0,1) = dNv1*v(1)+dNv2*v(5);
    N1(0,2) = dNw1*v(2)+dNw2*v(6);
    N1(1,3) = 1.0;
    N1(1,4) = Nf1*v(0)+Nf2*v(3)+Nf3*v(4)+Nf4*v(7);
    N1(1,5) = ddNw1*v(2)+ddNw2*v(6);
    N1(2,3) = -N1(1,4);
    N1(2,4) = 1.0;
    N1(2,5) = -ddNv1*v(1)-ddNv2*v(5);
    N1(3,6) = dNf1*v(0)+dNf2*v(3)+dNf3*v(4)+dNf4*v(7);
    N1(4,7) = -1.0;
    N1(5,6) = 1.0;
    N2(0,8) = oneOverL;
    N2(1,1) = dNv1;
    N2(1,5) = dNv2;
    N2(2,2) = dNw1;
    N2(2,6) = dNw2;
    N2(3,1) = ddNv1;
    N2(3,5) = ddNv2;
    N2(4,2) = ddNw1;
    N2(4,6) = ddNw2;
    N2(5,0) = Nf1;
    N2(5,3) = Nf2;
    N2(5,4) = Nf3;
    N2(5,7) = Nf4;
    N2(6,0) = dNf1;
    N2(6,3) = dNf2;
    N2(6,4) = dNf3;
    N2(6,7) = dNf4;
    N2(7,0) = ddNf1;
    N2(7,3) = ddNf2;
    N2(7,4) = ddNf3;
    N2(7,7) = ddNf4;
    
    // Get section stress resultant
    const Vector &s = theSections[i]->getStressResultant();
  
    double wti=wt[i];
  /*q(0) += ((6.0-xi12)*oneOverL*s(3)+(6.0*xi1*xi1-xi6)*s(4))*wti;
  q(1) += (xi6-4.0)*s(1)*wti;
  q(2) += (xi6-4.0)*s(2)*wti;
  q(3) += ((4.0-xi6)*s(3)+L*(1.0-4.0*xi1+3.0*xi1*xi1)*s(4))*wti;
  q(4) += ((xi12-6.0)*oneOverL*s(3)+(xi6-6.0*xi1*xi1)*s(4))*wti;
  q(5) += (xi6-2.0)*s(1)*wti;
  q(6) += (xi6-2.0)*s(2)*wti;
  q(7) += ((2.0-xi6)*s(3)+(3.0*xi1*xi1-2.0*xi1)*L*s(4))*wti;
  q(8) += s(0)*wti;*/

  	
    static Vector qProduct1(8);
    static Vector qProduct2(9);
    qProduct1.Zero();
    qProduct2.Zero();
    qProduct1.addMatrixTransposeVector(0.0, N1, s, 1.0);
    qProduct2.addMatrixTransposeVector(0.0, N2, qProduct1, 1.0);
    
    for (int j=0; j<9; j++)
      {
	q(j) += qProduct2(j)*L*wti;
      }
  }
  
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

  // Transform forces
  Vector p0Vec(p0, 5);
  P = crdTransf->getGlobalResistingForce(q, p0Vec);
  
  // Subtract other external nodal loads ... P_res = P_int - P_ext
  P.addVector(1.0, Q, -1.0);
  
  return P;
}

const Vector&
DispBeamColumnWarping3d::getResistingForceIncInertia()
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
    P(7) += m*accel2(0);
    P(8) += m*accel2(1);
    P(9) += m*accel2(2);

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
DispBeamColumnWarping3d::sendSelf(int commitTag, Channel &theChannel)
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
    opserr << "DispBeamColumnWarping3d::sendSelf() - failed to send ID data\n";
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
      opserr << "DispBeamColumnWarping3d::sendSelf() - failed to send double data\n";
      return -1;
    }    
  }

  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "DispBeamColumnWarping3d::sendSelf() - failed to send crdTranf\n";
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
    opserr << "DispBeamColumnWarping3d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  //
  // send the sections
  //
  
  for (j = 0; j<numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "DispBeamColumnWarping3d::sendSelf() - section " << j << "failed to send itself\n";
      return -1;
    }
  }

  return 0;
}

int
DispBeamColumnWarping3d::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i;
  
  static ID idData(7); // one bigger than needed so no clash with section ID

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    opserr << "DispBeamColumnWarping3d::recvSelf() - failed to recv ID data\n";
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
      opserr << "DispBeamColumnWarping3d::sendSelf() - failed to recv double data\n";
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
	opserr << "DispBeamColumnWarping3d::recvSelf() - " <<
	  "failed to obtain a CrdTrans object with classTag" <<
	  crdTransfClassTag << endln;
	return -2;	  
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "DispBeamColumnWarping3d::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }      
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumnWarping3d::recvSelf() - failed to recv ID data\n";
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
      opserr << "DispBeamColumnWarping3d::recvSelf() - out of memory creating sections array of size" <<
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
	opserr << "DispBeamColumnWarping3d::recvSelf() - Broker could not create Section of class type" <<
	  sectClassTag << endln;
	exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumnWarping3d::recvSelf() - section " <<
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
	  opserr << "DispBeamColumnWarping3d::recvSelf() - Broker could not create Section of class type" <<
	    sectClassTag << endln;
	  exit(-1);
	}
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumnWarping3d::recvSelf() - section " << 
	  i << "failed to recv itself\n";
	return -1;
      }     
    }
  }

  return 0;
}

void
DispBeamColumnWarping3d::Print(OPS_Stream &s, int flag)
{
  s << "\nDispBeamColumnWarping3d, element id:  " << this->getTag() << endln;
  s << "\tConnected external nodes:  " << connectedExternalNodes;
  s << "\tmass density:  " << rho << endln;

  double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  N   = q(0);
  Mz1 = q(1);
  Mz2 = q(2);
  Vy  = (Mz1+Mz2)*oneOverL;
  My1 = q(3);
  My2 = q(4);
  Vz  = -(My1+My2)*oneOverL;
  T   = q(5);

  s << "\tEnd 1 Forces (P Mz Vy My Vz T): "
    << -N+p0[0] << ' ' << Mz1 << ' ' <<  Vy+p0[1] << ' ' << My1 << ' ' <<  Vz+p0[3] << ' ' << -T << endln;
  s << "\tEnd 2 Forces (P Mz Vy My Vz T): "
    <<  N << ' ' << Mz2 << ' ' << -Vy+p0[2] << ' ' << My2 << ' ' << -Vz+p0[4] << ' ' <<  T << endln;
  
  //for (int i = 0; i < numSections; i++)
  //theSections[i]->Print(s,flag);
}


int
DispBeamColumnWarping3d::displaySelf(Renderer &theViewer, int displayMode, float fact)
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
DispBeamColumnWarping3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","DispBeamColumnWarping3d");
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
    if (argc > 2) {
    
      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections) {
	
      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamInt->getSectionLocations(numSections, L, xi);
  
	output.tag("GaussPointOutput");
	output.attr("number",sectionNum);
	output.attr("eta",xi[sectionNum-1]*L);

	theResponse =  theSections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
	
	output.endTag();
      }
    }
  }
  
  output.endTag();
  return theResponse;
}

int 
DispBeamColumnWarping3d::getResponse(int responseID, Information &eleInfo)
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
DispBeamColumnWarping3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;
  
  // If the parameter belongs to the element itself
  if (strcmp(argv[0],"rho") == 0)
    return param.addObject(1, this);
  
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
  
  else if (strstr(argv[0],"integration") != 0) {
    
    if (argc < 2)
      return -1;

    return beamInt->setParameter(&argv[1], argc-1, param);
  }

  // Default, send to every object
  int ok = 0;
  for (int i = 0; i < numSections; i++)
    ok += theSections[i]->setParameter(argv, argc, param);
  ok += beamInt->setParameter(argv, argc, param);
  return ok;
}

int
DispBeamColumnWarping3d::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  else
    return -1;  
}




int
DispBeamColumnWarping3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}

const Matrix &
DispBeamColumnWarping3d::getKiSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}

const Matrix &
DispBeamColumnWarping3d::getMassSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}



const Vector &
DispBeamColumnWarping3d::getResistingForceSensitivity(int gradNumber)
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
DispBeamColumnWarping3d::commitSensitivity(int gradNumber, int numGrads)
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


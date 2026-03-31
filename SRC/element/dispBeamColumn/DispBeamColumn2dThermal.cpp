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
// $Date: 2011-07-18 10:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/DispBeamColumn2dThermal.cpp,v $
                                                                        
// Written: MHS
// Created: Feb 2001
// Modified: Jian Zhang[University of Edinburgh]
// Modified: Panagiotis Kotsovinos[University of Edinburgh]
// Modified: Jian Jiang[University of Edinburgh]
// Modified: Liming Jiang[University of Edinburgh,2014



#include <DispBeamColumn2dThermal.h>
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
#include <NodalThermalAction.h>
#include <ThermalActionWrapper.h>
#include <elementAPI.h>

#include <math.h>
#include <StandardStream.h>
#include <fstream>
#include <stdlib.h>
#include <FiberSection2dThermal.h>

Matrix DispBeamColumn2dThermal::K(6,6);
Vector DispBeamColumn2dThermal::P(6);
double DispBeamColumn2dThermal::workArea[100];

void* OPS_DispBeamColumn2dThermal()
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
    
    Element *theEle =  new DispBeamColumn2dThermal(iData[0],iData[1],iData[2],secTags.Size(),sections,
					    *bi,*theTransf,mass);
    delete [] sections;
    return theEle;
}

DispBeamColumn2dThermal::DispBeamColumn2dThermal(int tag, int nd1, int nd2,
				   int numSec, SectionForceDeformation **s,
				   BeamIntegration& bi,
				   CrdTransf &coordTransf, double r)
:Element (tag, ELE_TAG_DispBeamColumn2dThermal), 
 numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
  connectedExternalNodes(2),
  Q(6), q(3), rho(r), parameterID(0)
{
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];
    
  if (theSections == 0) {
    opserr << "DispBeamColumn2dThermal::DispBeamColumn2dThermal - failed to allocate section model pointer\n";
    exit(-1);
  }

  for (int i = 0; i < numSections; i++) {
    
    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getCopy();
    
    // Check allocation
    if (theSections[i] == 0) {
      opserr << "DispBeamColumn2dThermal::DispBeamColumn2dThermal -- failed to get a copy of section model\n";
      exit(-1);
    }
  }
  
  beamInt = bi.getCopy();
  
  if (beamInt == 0) {
    opserr << "DispBeamColumn2dThermal::DispBeamColumn2dThermal - failed to copy beam integration\n";
    exit(-1);
  }

  crdTransf = coordTransf.getCopy2d();
  
  if (crdTransf == 0) {
    opserr << "DispBeamColumn2dThermal::DispBeamColumn2dThermal - failed to copy coordinate transformation\n";
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

  dataMix = new double [27];
  for(int i=0; i<27 ;i++){
		dataMix[i] = 0;
  }
  for(int i=0; i<10 ;i++){
		SectionThermalElong[i] = 0;
  }
   
  
  //PK

  q0Temperature[0] = 0.0; 
  q0Temperature[1] = 0.0; 
  q0Temperature[2] = 0.0; 
  q0TemperatureP[0] = 0.0; 
  q0TemperatureP[1] = 0.0; 
  q0TemperatureP[2] = 0.0; 
  counterTemperature = 0;
  AverageThermalElong =0.0;

  
//PK  loadfactors
 loadFactor2=0;
loadFactor3=0;
loadFactor4=0;
loadFactor5=0;
loadFactor6=0;
loadFactor7=0;
loadFactor8=0;
loadFactor9=0;
//pk end
}

DispBeamColumn2dThermal::DispBeamColumn2dThermal()
:Element (0, ELE_TAG_DispBeamColumn2dThermal),
 numSections(0), theSections(0), crdTransf(0), beamInt(0),
 connectedExternalNodes(2),
  Q(6), q(3), rho(0.0), parameterID(0)
{
    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;

    theNodes[0] = 0;
    theNodes[1] = 0;

  dataMix = new double [27];
  for(int i=0; i<27 ;i++){
		dataMix[i] = 0;
  }
   for(int i=0; i<10 ;i++){
		SectionThermalElong[i] = 0;
  }

  q0Temperature[0] = 0.0; 
  q0Temperature[1] = 0.0; 
  q0Temperature[2] = 0.0; 
  q0TemperatureP[0] = 0.0; 
  q0TemperatureP[1] = 0.0; 
  q0TemperatureP[2] = 0.0; 
  counterTemperature = 0;
  AverageThermalElong =0.0;

  //PK  loadfactors
 loadFactor2=0;
loadFactor3=0;
loadFactor4=0;
loadFactor5=0;
loadFactor6=0;
loadFactor7=0;
loadFactor8=0;
loadFactor9=0;
//pk end
}

DispBeamColumn2dThermal::~DispBeamColumn2dThermal()
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

  if (dataMix != 0)
    delete [] dataMix;
}

int
DispBeamColumn2dThermal::getNumExternalNodes() const
{
    return 2;
}

const ID&
DispBeamColumn2dThermal::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
DispBeamColumn2dThermal::getNodePtrs()
{
    return theNodes;
}

int
DispBeamColumn2dThermal::getNumDOF()
{
    return 6;
}

void
DispBeamColumn2dThermal::setDomain(Domain *theDomain)
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
	//opserr << "FATAL ERROR DispBeamColumn2dThermal (tag: %d), node not found in domain",
	//	this->getTag());
	
	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    if (dofNd1 != 3 || dofNd2 != 3) {
	//opserr << "FATAL ERROR DispBeamColumn2dThermal (tag: %d), has differing number of DOFs at its nodes",
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
DispBeamColumn2dThermal::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "DispBeamColumn2dThermal::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    return retVal;
}

int
DispBeamColumn2dThermal::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    return retVal;
}

int
DispBeamColumn2dThermal::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    return retVal;
}

int
DispBeamColumn2dThermal::update(void)
{
  //opserr<<"**Now Update: Ele "<< this->getTag()<<" **"<<endln;
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
  
  //opserr<< "Basic Deformation: "<<v<<endln;
  // Loop over the integration points
  #ifdef _DEBUG
  //opserr << "BeamNUT_" << this->getTag() << " ,  v  " << v << endln;
	  //<<" , Average epsT "<<AverageThermalElong<<endln;
  #endif
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();   // return 2
    const ID &code = theSections[i]->getType();   // code: Section_Response_P, Section_Response_MZ
    
    Vector e(workArea, order);
    
    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];
    int j;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	e(j) = oneOverL*v(0)-AverageThermalElong+SectionThermalElong[i]; break;
      case SECTION_RESPONSE_MZ:
	e(j) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2)); 
	break;
      default:
	e(j) = 0.0; break;
      }
    }
    
    // Set the section deformations	
   // err += theSections[i]->setTrialSectionDeformation(e);
    //J.Z

     //FMK err += theSections[i]->setTrialSectionDeformationTemperature(e,dataMix); 
#ifdef _BDEBUG
	//opserr << "Beam_" << this->getTag() << " Section " << i << endln;
#endif
	Vector dataMixV(dataMix,27);
    err += theSections[i]->setTrialSectionDeformation(e);  

  }
  
  if (err != 0) {
    opserr << "DispBeamColumn2dThermal::update() - failed setTrialSectionDeformations()\n";
    return err;
  }

  return 0;
}

const Matrix&
DispBeamColumn2dThermal::getTangentStiff()
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
#ifdef _DEBUG
 // opserr<<"Beam_"<<this->getTag()<<" Section "<<i<<endln<<" ,  K  "<<ks<<endln;
#endif
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

  return K;
}

const Matrix&
DispBeamColumn2dThermal::getInitialBasicStiff()
{
  static Matrix kb(3,3);

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
DispBeamColumn2dThermal::getInitialStiff()
{
  const Matrix &kb = this->getInitialBasicStiff();

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);

  return K;
}

const Matrix&
DispBeamColumn2dThermal::getMass()
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
DispBeamColumn2dThermal::zeroLoad(void)
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
DispBeamColumn2dThermal::addLoad(ElementalLoad *theLoad, double loadFactor)
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

  else if (type == LOAD_TAG_Beam2dThermalAction) {
     // load not inside fire load pattern
	 //const Vector &data = theLoad->getData(type, loadFactor);

	 /// This code block is added by LJ and copied from DispBeamColumn2d(Modified edition) for 'FireLoadPattern'--08-May-2012--//[END]
    counterTemperature = 1;

    q0Temperature[0] = 0.0; 
    q0Temperature[1] = 0.0; 
    q0Temperature[2] = 0.0; 
    
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
      
      double xi6 = 6.0*xi[i];
      
      // Get section stress resultant 
      // FMk const Vector &s = theSections[i]->getTemperatureStress(dataMix);

      Vector dataMixV(27);
	  for(int m=0;m<9;m++){
		  dataMixV(2*m)=data(2*m); //Linear temperature interpolation
		  dataMixV(2*m+1)=data(2*m+1);
		  dataMixV(18+m)=1000;
		  }
      const Vector &s = theSections[i]->getTemperatureStress(dataMixV);    //contribuited by ThermalElongation
#ifdef _BDEBUG
	  if(this->getTag()==1)
	  opserr<<"Thermal Stress "<<s<<endln;
#endif
	  if(s!=0){
      double si;
      for (int j = 0; j < order; j++) {
		si = s(j)*wt[i];
		switch(code(j)) {
		case SECTION_RESPONSE_P:
		q0Temperature[0] += si; break;
		case SECTION_RESPONSE_MZ:
		q0Temperature[1] += (xi6-4.0)*si; 
		q0Temperature[2] += (xi6-2.0)*si; break;
		default:
		break;
		}
      }
	 }

    }
    
    // q0[0] -= q0Temperature[0];
    // q0[1] -= q0Temperature[1];
    // q0[2] -= q0Temperature[2];

  }
  else if(type == LOAD_TAG_NodalThermalAction) {
	  //NodalLoad* theNodalThermal0,theNodalThermal1;	 
	 NodalThermalAction* theNodalThermal0 = theNodes[0]->getNodalThermalActionPtr();
	 NodalThermalAction* theNodalThermal1 = theNodes[1]->getNodalThermalActionPtr();	
	 int type;
	 const Vector &data0 = theNodalThermal0->getData(type);
	 const Vector &data1 = theNodalThermal1->getData(type);
	 Vector Loc(9);
	 Vector NodalT0(9);
	 Vector NodalT1(9);
	 #ifdef _DEBUG
 	 opserr<<"NodalT0: "<<data0<<endln<<"NodalT1: "<<data1;
	 #endif
	 for(int i =0; i<9;i++){
		 if(data0(2*i+1)-data1(2*i+1)>1e-8||data0(2*i+1)-data1(2*i+1)<-1e-8){
			 opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
			      << "incompatible loc input for datapoint "<< i << endln;
			 }
		 else{
			 Loc(i)=data0(2*i+1);
			 NodalT0(i)=data0(2*i);
			 NodalT1(i)=data1(2*i);
			 }
		 }

	 double L = crdTransf->getInitialLength();
	  counterTemperature = 1;
	AverageThermalElong =0.0;
    q0Temperature[0] = 0.0;
    q0Temperature[1] = 0.0; 
    q0Temperature[2] = 0.0; 

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
      
      double xi6 = 6.0*xi[i];
      
      // Get section stress resultant 
      // FMk const Vector &s = theSections[i]->getTemperatureStress(dataMix);
	  Vector dataMixV(27);
	  for(int m=0;m<9;m++){
		  dataMixV(2*m)=NodalT0(m)+xi[i]*(NodalT1(m)-NodalT0(m)); //Linear temperature interpolation
		  dataMixV(2*m+1)=Loc(m);
		  dataMixV(18+m)=1000;
		  }
      const Vector &s = theSections[i]->getTemperatureStress(dataMixV);    //contribuited by ThermalElongation
	  if(theSections[i]->getClassTag()==SEC_TAG_FiberSection2dThermal){
		SectionThermalElong[i]= ((FiberSection2dThermal*)theSections[i])->getThermalElong()(0);
	  }
	  //opserr<<"Thermal Stress "<<s<<endln;
      double si;
	  double ThermalEloni;
	  if(s!=0){
      for (int j = 0; j < order; j++) {
		si = s(j)*wt[i];
		ThermalEloni = SectionThermalElong[i]*wt[i];
		switch(code(j)) {
		case SECTION_RESPONSE_P:
		q0Temperature[0] += si; 
		AverageThermalElong += ThermalEloni; break;
		case SECTION_RESPONSE_MZ:
		q0Temperature[1] += (xi6-4.0)*si; 
		q0Temperature[2] += (xi6-2.0)*si; break;
		default:
		break;
		}
      }  
	  }
    }
    
    // q0[0] -= q0Temperature[0];
    // q0[1] -= q0Temperature[1];
    // q0[2] -= q0Temperature[2];
}
//----------------------------------------------
  else if(type == LOAD_TAG_ThermalActionWrapper) {
	  //NodalLoad* theNodalThermal0,theNodalThermal1;	 	
	 counterTemperature = 1;
	AverageThermalElong =0.0;
    q0Temperature[0] = 0.0;
    q0Temperature[1] = 0.0; 
    q0Temperature[2] = 0.0; 

    double oneOverL = 1.0/L;
    
    //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
    //const Vector &wts = quadRule.getIntegrPointWeights(numSections);  
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);
	//xi: 0.0~1.0;
    double wt[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wt);
    
    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {
      
      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();
      
      double xi6 = 6.0*xi[i];
      
      // Get section stress resultant 
      // FMk const Vector &s = theSections[i]->getTemperatureStress(dataMix);
	  Vector theNode0Crds = theNodes[0]->getCrds();
	  Vector theNode1Crds = theNodes[1]->getCrds();
	  
	  int ndm = theNode0Crds.Size();
	  Vector theIntCrds = Vector(ndm);
	  for(int m = 0; m<ndm; m++){
			theIntCrds(m) =  theNode0Crds(m)+xi[i]*(theNode1Crds(m)- theNode0Crds(m));
	  }
#ifdef _BDEBUG
	  opserr<<"DispBeamColumn2dThermal::addLoad at the ele "<<this->getTag()<< " section "<<i<<" --theIntCrds"<<theIntCrds<<endln;
#endif
	  //NodalThermalActions attached with the wrapper would have been updated by pattern;
	  Vector dataMixV = ((ThermalActionWrapper*) theLoad)->getIntData(theIntCrds);
	  
      const Vector &s = theSections[i]->getTemperatureStress(dataMixV);    //contribuited by ThermalElongation
	  if(theSections[i]->getClassTag()==SEC_TAG_FiberSection2dThermal){
		SectionThermalElong[i]= ((FiberSection2dThermal*)theSections[i])->getThermalElong()(0);
	  }
#ifdef _BDEBUG
	  //opserr<<"Thermal Stress "<<s<<endln;
#endif
      double si;
	  double ThermalEloni;
	  if(s!=0){
      for (int j = 0; j < order; j++) {
		si = s(j)*wt[i];
		ThermalEloni = SectionThermalElong[i]*wt[i];
		switch(code(j)) {
		case SECTION_RESPONSE_P:
		q0Temperature[0] += si; 
		AverageThermalElong += ThermalEloni; break;
		case SECTION_RESPONSE_MZ:
		q0Temperature[1] += (xi6-4.0)*si; 
		q0Temperature[2] += (xi6-2.0)*si; break;
		default:
		break;
		}
      }  
	  //end of for loop for j order
	  }
	  //end of if s!=0;
    }
	//end of for loop
    
    // q0[0] -= q0Temperature[0];
    // q0[1] -= q0Temperature[1];
    // q0[2] -= q0Temperature[2];
}
 //----------------------------------------
else {
	   opserr << "DispBeamColumn2dThermal::addLoad(double) -- load type " << theLoad->getClassType() 
		  << "unknown for element with tag: " << this->getTag() << "\n"; 
    
    return -1;
  }
  
  return 0;
}


int 
DispBeamColumn2dThermal::addLoad(ElementalLoad *theLoad, const Vector &factors)
{
  int type;
  const Vector &data = theLoad->getData(type, factors(0));
  double L = crdTransf->getInitialLength();

  /// This code block is added by LJ and copied from DispBeamColumn2d(Modified edition) for 'FireLoadPattern'--08-May-2012--//[BEGIN]
 //JZ 07/10 /////////////////////////////////////////////////////////////start
  if (type == LOAD_TAG_Beam2dThermalAction) {
  //PK
    //  Vector &factors = theLoad->getfactors();
    
    double loadFactor = factors(0);
    loadFactor2=factors(1);
    loadFactor3=factors(2);
    loadFactor4=factors(3);
    loadFactor5=factors(4);
    loadFactor6=factors(5);
    loadFactor7=factors(6);
    loadFactor8=factors(7);
    loadFactor9=factors(8);
    
    dataMix[0] = data(0)*loadFactor;
    dataMix[2] = data(2)*loadFactor2;
    dataMix[4] = data(4)*loadFactor3;
    dataMix[6] = data(6)*loadFactor4;
    dataMix[8] = data(8)*loadFactor5;
    dataMix[10] = data(10)*loadFactor6;
    dataMix[12] = data(12)*loadFactor7;
    dataMix[14] = data(14)*loadFactor8;
    dataMix[16] = data(16)*loadFactor9;
    
    dataMix[1] = data(1);
    dataMix[3] = data(3);
    dataMix[5] = data(5);
    dataMix[7] = data(7);
    dataMix[9] = data(9);
    dataMix[11] = data(11);
    dataMix[13] = data(13);
    dataMix[15] = data(15);
    dataMix[17] = data(17);
    
    //PK add the maximum temperatures to be passed along with the factor ones in the dataMix
    //18-26
    dataMix[18] = data(0);
    dataMix[19] = data(2);
    dataMix[20] = data(4);
    dataMix[21] = data(6);
    dataMix[22] = data(8);
    dataMix[23] = data(10);
    dataMix[24] = data(12);
    dataMix[25] = data(14);
    dataMix[26] = data(16);
    //PK end of change
    
    /// This code block is added by LJ and copied from DispBeamColumn2d(Modified edition) for 'FireLoadPattern'--08-May-2012--//[END]
    counterTemperature = 1;

    q0Temperature[0] = 0.0; 
    q0Temperature[1] = 0.0; 
    q0Temperature[2] = 0.0; 
    
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
      
      double xi6 = 6.0*xi[i];
      
      // Get section stress resultant 
      // FMk const Vector &s = theSections[i]->getTemperatureStress(dataMix);
      Vector dataMixV(dataMix, 27);
      const Vector &s = theSections[i]->getTemperatureStress(dataMixV);    //contribuited by ThermalElongation
	  //opserr<<"Thermal Stress "<<s<<endln;
      double si;
      for (int j = 0; j < order; j++) {
		si = s(j)*wt[i];
		switch(code(j)) {
		case SECTION_RESPONSE_P:
		q0Temperature[0] += si; break;
		case SECTION_RESPONSE_MZ:
		q0Temperature[1] += (xi6-4.0)*si; 
		q0Temperature[2] += (xi6-2.0)*si; break;
		default:
		break;
		}
      }   
    }
    
    // q0[0] -= 0;
    // q0[1] -= 0;
    // q0[2] -= 0;
	
  }  
  else {
	  opserr << "DispBeamColumn2dThermal::addLoad(Vector) -- load type " << theLoad->getClassType() 
		  << "unknown for element with tag: " << this->getTag() << "\n"; 
    return -1;
  }
  
  return 0;
}

int 
DispBeamColumn2dThermal::addInertiaLoadToUnbalance(const Vector &accel)
{
	// Check for a quick return
	if (rho == 0.0) 
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = theNodes[0]->getRV(accel);
	const Vector &Raccel2 = theNodes[1]->getRV(accel);

    if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
      opserr << "DispBeamColumn2dThermal::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
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
DispBeamColumn2dThermal::getResistingForce()
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
  
  if (counterTemperature ==1 ) {
     // this->update();
	 // counterTemperature++;
  }
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
#ifdef _BDEBUG
	if(this->getTag()==1)
	opserr<<"ele 1"<<s<<endln;
#endif
    double si;
    for (int j = 0; j < order; j++) {
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


	if (counterTemperature == 1) { 
		#ifdef _BDEBUG
		   opserr <<" Ele: "<< this->getTag()<<" ,q before Thermal force"<<q<<endln;
		#endif
		
		q(0) -= (q0Temperature[0]);
		q(1) -= (q0Temperature[1]);
		q(2) -= (q0Temperature[2]);
		//q0TemperatureP[0]= q0Temperature[0];
		//q0TemperatureP[1]= q0Temperature[1];
		//q0TemperatureP[2]= q0Temperature[2];
		counterTemperature++;
		  }

#ifdef _BDEBUG
 opserr <<" Ele: "<< this->getTag()<< ", q after Thermal force"<<q<<endln;
#endif




  // Vector for reactions in basic system
  Vector p0Vec(p0, 3);

  P = crdTransf->getGlobalResistingForce(q, p0Vec);

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  //P.addVector(1.0, Q, -1.0);
// opserr <<" Ele: "<< this->getTag()<< ", P "<<P<<endln;
  P(0) -= Q(0);
  P(1) -= Q(1);
  P(2) -= Q(2);
  P(3) -= Q(3);
  P(4) -= Q(4);
  P(5) -= Q(5);
  return P;

   
   
}

const Vector&
DispBeamColumn2dThermal::getResistingForceIncInertia()
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
    P(3) += m*accel2(0);
    P(4) += m*accel2(1);
    
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
DispBeamColumn2dThermal::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j;
  int loc = 0;
  
  static Vector data(14);
  data(0) = this->getTag();
  data(1) = connectedExternalNodes(0);
  data(2) = connectedExternalNodes(1);
  data(3) = numSections;
  data(4) = crdTransf->getClassTag();
  int crdTransfDbTag  = crdTransf->getDbTag();
  if (crdTransfDbTag  == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag  != 0) 
      crdTransf->setDbTag(crdTransfDbTag);
  }
  data(5) = crdTransfDbTag;
  data(6) = beamInt->getClassTag();
  int beamIntDbTag  = beamInt->getDbTag();
  if (beamIntDbTag  == 0) {
    beamIntDbTag = theChannel.getDbTag();
    if (beamIntDbTag  != 0) 
      beamInt->setDbTag(beamIntDbTag);
  }
  data(7) = beamIntDbTag;
  data(8) = rho;
  //data(9) = cMass;
  data(10) = alphaM;
  data(11) = betaK;
  data(12) = betaK0;
  data(13) = betaKc;
  
  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "DispBeamColumn2d::sendSelf() - failed to send data Vector\n";
     return -1;
  }

  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "DispBeamColumn2dThermal::sendSelf() - failed to send crdTranf\n";
     return -1;
  }      

  // send the beam integration
  if (beamInt->sendSelf(commitTag, theChannel) < 0) {
    opserr << "DispBeamColumn2d::sendSelf() - failed to send beamInt\n";
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
    opserr << "DispBeamColumn2dThermal::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  //
  // send the sections
  //
  
  for (j = 0; j<numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "DispBeamColumn2dThermal::sendSelf() - section " << 
	j << "failed to send itself\n";
      return -1;
    }
  }

  return 0;
}

int
DispBeamColumn2dThermal::recvSelf(int commitTag, Channel &theChannel,
			   FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i;
  
  static Vector data(14);

  if (theChannel.recvVector(dbTag, commitTag, data) < 0)  {
    opserr << "DispBeamColumn2dThermal::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  this->setTag((int)data(0));
  connectedExternalNodes(0) = (int)data(1);
  connectedExternalNodes(1) = (int)data(2);
  int nSect = (int)data(3);
  int crdTransfClassTag = (int)data(4);
  int crdTransfDbTag = (int)data(5);

  int beamIntClassTag = (int)data(6);
  int beamIntDbTag = (int)data(7);
  
  rho = data(8);
  //cMass = (int)data(9);
  
  alphaM = data(10);
  betaK = data(11);
  betaK0 = data(12);
  betaKc = data(13);  

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;

      crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

      if (crdTransf == 0) {
	opserr << "DispBeamColumn2dThermal::recvSelf() - failed to obtain a CrdTrans object with classTag " <<
	  crdTransfClassTag << endln;
	  return -2;	  
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "DispBeamColumn2dThermal::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }      

  // create a new beamInt object if one needed
  if (beamInt == 0 || beamInt->getClassTag() != beamIntClassTag) {
      if (beamInt != 0)
	  delete beamInt;

      beamInt = theBroker.getNewBeamIntegration(beamIntClassTag);

      if (beamInt == 0) {
	opserr << "DispBeamColumn2d::recvSelf() - failed to obtain the beam integration object with classTag" <<
	  beamIntClassTag << endln;
	exit(-1);
      }
  }

  beamInt->setDbTag(beamIntDbTag);

  // invoke recvSelf on the beamInt object
  if (beamInt->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     opserr << "DispBeamColumn2d::sendSelf() - failed to recv beam integration\n";
     return -3;
  }
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*nSect);
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumn2dThermal::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  //
  // now receive the sections
  //
  
  if (numSections != nSect) {

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
    theSections = new SectionForceDeformation *[nSect];
    if (theSections == 0) {
opserr << "DispBeamColumn2dThermal::recvSelf() - out of memory creating sections array of size " <<
  nSect << endln;
      return -1;
    }    

    // create a section and recvSelf on it
    numSections = nSect;
    loc = 0;
    
    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      theSections[i] = theBroker.getNewSection(sectClassTag);
      if (theSections[i] == 0) {
	opserr << "DispBeamColumn2dThermal::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumn2dThermal::recvSelf() - section " << i << " failed to recv itself\n";
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
	opserr << "DispBeamColumn2dThermal::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
	}
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumn2dThermal::recvSelf() - section " << i << " failed to recv itself\n";
	return -1;
      }     
    }
  }

  return 0;
}

void
DispBeamColumn2dThermal::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nDispBeamColumn2dThermal, element id:  " << this->getTag() << endln;
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
        s << "\"type\": \"DispBeamColumn2dThermal\", ";
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
DispBeamColumn2dThermal::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes, int numModes)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response*
DispBeamColumn2dThermal::setResponse(const char **argv, int argc,
			      OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","DispBeamColumn2dThermal");
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
  else if (strstr(argv[0],"section") != 0) {
    if (argc > 2) {
      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections) {

	output.tag("GaussPointOutput");
	output.attr("number",sectionNum);
	double xi[maxNumSections];
	double L = crdTransf->getInitialLength();
	beamInt->getSectionLocations(numSections, L, xi);
	output.attr("eta",xi[sectionNum-1]*L);

	theResponse = theSections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
	
	output.endTag();
      }
    }
  }
  
  // curvature sensitivity along element length
  else if (strcmp(argv[0],"dcurvdh") == 0)
    return new ElementResponse(this, 5, Vector(numSections));
  
  // basic deformation sensitivity
  else if (strcmp(argv[0],"dvdh") == 0)
    return new ElementResponse(this, 6, Vector(3));
  
  else if (strcmp(argv[0],"integrationPoints") == 0)
    return new ElementResponse(this, 7, Vector(numSections));
  
  else if (strcmp(argv[0],"integrationWeights") == 0)
    return new ElementResponse(this, 8, Vector(numSections));
  
  output.endTag();
  return theResponse;
}

int 
DispBeamColumn2dThermal::getResponse(int responseID, Information &eleInfo)
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

  // Curvature sensitivity
  else if (responseID == 5) {
    /*
      Vector curv(numSections);
      const Vector &v = crdTransf->getBasicDispGradient(1);
      
      double L = crdTransf->getInitialLength();
      double oneOverL = 1.0/L;
      //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
      double pts[2];
      pts[0] = 0.0;
      pts[1] = 1.0;
      
      // Loop over the integration points
      for (int i = 0; i < numSections; i++) {
	int order = theSections[i]->getOrder();
	const ID &code = theSections[i]->getType();
	//double xi6 = 6.0*pts(i,0);
	double xi6 = 6.0*pts[i];
	curv(i) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2));
      }
      
      return eleInfo.setVector(curv);
    */

    Vector curv(numSections);

    /*
    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {
      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();
      const Vector &dedh = theSections[i]->getdedh();
      for (int j = 0; j < order; j++) {
	if (code(j) == SECTION_RESPONSE_MZ)
	  curv(i) = dedh(j);
      }
    }
    */

    return eleInfo.setVector(curv);
  }

  // Basic deformation sensitivity
  else if (responseID == 6) {  
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(1);
    return eleInfo.setVector(dvdh);
  }

  else if (responseID == 7) {
    //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = xi[i]*L;
    return eleInfo.setVector(locs);
  }

  else if (responseID == 8) {
    //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
    double wt[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wt);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wt[i]*L;
    return eleInfo.setVector(weights);
  }

  else
    return -1;
}


// AddingSensitivity:BEGIN ///////////////////////////////////
int
DispBeamColumn2dThermal::setParameter(const char **argv, int argc, Parameter &param)
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
    
    // Get section number: 1...Np
    int sectionNum = atoi(argv[1]);
    
    if (sectionNum > 0 && sectionNum <= numSections)
      return theSections[sectionNum-1]->setParameter(&argv[2], argc-2, param);
    
    else
      return -1;
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
DispBeamColumn2dThermal::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  else
    return -1;  
}




int
DispBeamColumn2dThermal::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}



const Matrix &
DispBeamColumn2dThermal::getKiSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}

const Matrix &
DispBeamColumn2dThermal::getMassSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}



const Vector &
DispBeamColumn2dThermal::getResistingForceSensitivity(int gradNumber)
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
  static Vector dqdh(3);
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
      default:
	break;
      }
    }
  }
  
  // Transform forces
  static Vector dp0dh(3);		// No distributed loads

  P.Zero();

  //////////////////////////////////////////////////////////////

  if (crdTransf->isShapeSensitivity()) {
    
    // Perform numerical integration to obtain basic stiffness matrix
    // Some extra declarations
    static Matrix kbmine(3,3);
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
      
      Matrix ka(workArea, order, 3);
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
DispBeamColumn2dThermal::commitSensitivity(int gradNumber, int numGrads)
{
  // Get basic deformation and sensitivities
  const Vector &v = crdTransf->getBasicTrialDisp();
  
  static Vector dvdh(3);
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


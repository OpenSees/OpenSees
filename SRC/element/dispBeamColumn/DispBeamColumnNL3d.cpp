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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for DispBeamColumnNL3d.

#include <DispBeamColumnNL3d.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <MatrixUtil.h>
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
#include <CompositeResponse.h>
#include <math.h>
#include <ElementalLoad.h>
#include <elementAPI.h>
#include <map>

Matrix DispBeamColumnNL3d::K(12,12);
Vector DispBeamColumnNL3d::P(12);
double DispBeamColumnNL3d::workArea[200];

void* OPS_DispBeamColumnNL3d()
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
	if(strcmp(type, "-cMass") == 0) {
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
    
    Element *theEle =  new DispBeamColumnNL3d(iData[0],iData[1],iData[2],secTags.Size(),sections,
					      *bi,*theTransf,mass);
    delete [] sections;
    return theEle;
}

void* OPS_DispBeamColumnNL3d(const ID &info)
{
    // data
    int iData[5];
    int numData;
    double mass = 0.0;
    int cmass = 0;

    // regular element, not in a mesh, get tags
    if (info.Size() == 0) {
	if(OPS_GetNumRemainingInputArgs() < 5) {
	    opserr<<"insufficient arguments:eleTag,iNode,jNode,transfTag,integrationTag <-mass mass> <-cmass>\n";
	    return 0;
	}

	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();
	if(ndm != 2 || ndf != 3) {
	    opserr<<"ndm must be 2 and ndf must be 3\n";
	    return 0;
	}

	// inputs:
	numData = 3;
	if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	    opserr<<"WARNING: invalid integer inputs\n";
	    return 0;
	}
    }

    // regular element, or in a mesh
    if (info.Size()==0 || info(0)==1) {
	if(OPS_GetNumRemainingInputArgs() < 2) {
	    opserr<<"insufficient arguments: transfTag,integrationTag\n";
	    return 0;
	}

	numData = 2;
	if(OPS_GetIntInput(&numData,&iData[3]) < 0) {
	    opserr << "WARNING invalid int inputs\n";
	    return 0;
	}

	// options
	numData = 1;
	while(OPS_GetNumRemainingInputArgs() > 0) {
	    const char* type = OPS_GetString();
	    if(strcmp(type, "-cMass") == 0) {
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
    }

    // store data for different mesh
    static std::map<int, Vector> meshdata;
    if (info.Size()>0 && info(0)==1) {
	if (info.Size() < 2) {
	    opserr << "WARNING: need info -- inmesh, meshtag\n";
	    return 0;
	}

	// save the data for a mesh
	Vector& mdata = meshdata[info(1)];
	mdata.resize(4);
	mdata(0) = iData[3];
	mdata(1) = iData[4];
	mdata(2) = mass;
	mdata(3) = cmass;
	return &meshdata;

    } else if (info.Size()>0 && info(0)==2) {
	if (info.Size() < 5) {
	    opserr << "WARNING: need info -- inmesh, meshtag, eleTag, nd1, nd2\n";
	    return 0;
	}

	// get the data for a mesh
	Vector& mdata = meshdata[info(1)];
	if (mdata.Size() < 4) return 0;

	iData[0] = info(2);
	iData[1] = info(3);
	iData[2] = info(4);
	iData[3] = mdata(0);
	iData[4] = mdata(1);
	mass = mdata(2);
	cmass = mdata(3);
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
    
    Element *theEle =  new DispBeamColumnNL3d(iData[0],iData[1],iData[2],secTags.Size(),sections,
					      *bi,*theTransf,mass);
    delete [] sections;
    return theEle;
}

DispBeamColumnNL3d::DispBeamColumnNL3d(int tag, int nd1, int nd2,
				   int numSec, SectionForceDeformation **s,
				   BeamIntegration& bi,
				   CrdTransf &coordTransf, double r)
:Element (tag, ELE_TAG_DispBeamColumnNL3d), 
 numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
  connectedExternalNodes(2),
  Q(12), q(6), rho(r), parameterID(0)
{
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];
    
  if (theSections == 0) {
    opserr << "DispBeamColumnNL3d::DispBeamColumnNL3d - failed to allocate section model pointer\n";
    exit(-1);
  }

  for (int i = 0; i < numSections; i++) {
    
    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getCopy();
    
    // Check allocation
    if (theSections[i] == 0) {
      opserr << "DispBeamColumnNL3d::DispBeamColumnNL3d -- failed to get a copy of section model\n";
      exit(-1);
    }
  }
  
  beamInt = bi.getCopy();
  
  if (beamInt == 0) {
    opserr << "DispBeamColumnNL3d::DispBeamColumnNL3d - failed to copy beam integration\n";
    exit(-1);
  }

  crdTransf = coordTransf.getCopy3d();
  
  if (crdTransf == 0) {
    opserr << "DispBeamColumnNL3d::DispBeamColumnNL3d - failed to copy coordinate transformation\n";
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
}

DispBeamColumnNL3d::DispBeamColumnNL3d()
:Element (0, ELE_TAG_DispBeamColumnNL3d),
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
}

DispBeamColumnNL3d::~DispBeamColumnNL3d()
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
DispBeamColumnNL3d::getNumExternalNodes() const
{
    return 2;
}

const ID&
DispBeamColumnNL3d::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
DispBeamColumnNL3d::getNodePtrs()
{
    return theNodes;
}

int
DispBeamColumnNL3d::getNumDOF()
{
    return 12;
}

void
DispBeamColumnNL3d::setDomain(Domain *theDomain)
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
      opserr << "WARNING DispBeamColumnNL3d (tag: %d), node not found in domain" << this->getTag() << endln;;
      return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    if (dofNd1 != 6 || dofNd2 != 6) {
	//opserr << "FATAL ERROR DispBeamColumnNL3d (tag: %d), has differing number of DOFs at its nodes",
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
DispBeamColumnNL3d::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "DispBeamColumnNL3d::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    return retVal;
}

int
DispBeamColumnNL3d::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    return retVal;
}

int
DispBeamColumnNL3d::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    return retVal;
}

int
DispBeamColumnNL3d::update(void)
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
    double zeta = xi[i];

    double thetaz = (3*zeta*zeta-4*zeta+1)*v(1) + (3*zeta*zeta-2*zeta)*v(2);
    double thetay = (3*zeta*zeta-4*zeta+1)*v(3) + (3*zeta*zeta-2*zeta)*v(4);    

    int j;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	e(j) = oneOverL*v(0) + 0.5*thetaz*thetaz + 0.5*thetay*thetay; break;
      case SECTION_RESPONSE_MZ:
	e(j) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2)); break;
      case SECTION_RESPONSE_MY:
	e(j) = oneOverL*((xi6-4.0)*v(3) + (xi6-2.0)*v(4)); break;
      case SECTION_RESPONSE_T:
	e(j) = oneOverL*v(5); break;	
      default:
	e(j) = 0.0; break;
      }
    }
    
    // Set the section deformations
    err += theSections[i]->setTrialSectionDeformation(e);
  }
  
  if (err != 0) {
    opserr << "DispBeamColumnNL3d::update() - failed setTrialSectionDeformations()\n";
    return err;
  }

  return 0;
}

void
DispBeamColumnNL3d::getBasicStiff(Matrix &kb, int initial)
{
  // Zero for integral
  kb.Zero();
  
  // Update the transformation
  crdTransf->update();
  
  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
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
    double zeta = xi[i];

    double c1 = 3*zeta*zeta-4*zeta+1;
    double c2 = 3*zeta*zeta-2*zeta;
    double thetaz = c1*v(1) + c2*v(2);
    double thetay = c1*v(3) + c2*v(4);    

    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getSectionTangent();
    const Vector &s = theSections[i]->getStressResultant();

    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    //double wti = wts(i)*oneOverL;
    double wti = wt[i]*oneOverL;
    double tmp;
    int j, k;

    // Axial force term: int_0^L N*C'*C dx
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	tmp = s(j)*wt[i]*L;
	kb(1,1) += tmp*c1*c1;
	kb(2,1) += tmp*c2*c1;
	kb(1,2) += tmp*c1*c2;
	kb(2,2) += tmp*c2*c2;

	kb(3,3) += tmp*c1*c1;
	kb(4,3) += tmp*c2*c1;
	kb(3,4) += tmp*c1*c2;
	kb(4,4) += tmp*c2*c2;	
	break;
      }
    }

    Matrix B(order,6);
    Matrix Cz(order,6);
    Matrix Cy(order,6);    
    static Matrix Cz1(1,6);
    static Matrix Cy1(1,6);    
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	B(j,0) = 1.0;

	Cz(j,1) = c1;
	Cz(j,2) = c2;
	Cz1(0,1) = c1;
	Cz1(0,2) = c2;

	Cy(j,3) = c1;
	Cy(j,4) = c2;
	Cy1(0,3) = c1;
	Cy1(0,4) = c2;	
	break;
      case SECTION_RESPONSE_T:
	B(j,5) = 1.0;
	break;
      case SECTION_RESPONSE_MZ:
	B(j,1) = xi6-4.0;
	B(j,2) = xi6-2.0;
	break;
      case SECTION_RESPONSE_MY:
	B(j,3) = xi6-4.0;
	B(j,4) = xi6-2.0;	
	break;
      default:
	break;
      }
    }

    // B'*ks*B
    kb.addMatrixTripleProduct(1.0, B, ks, wti);
    
    Matrix kC(order,6);

    // B'*ks*C*theta
    kC.addMatrixProduct(0.0, ks, Cz, 1.0);
    kb.addMatrixTransposeProduct(1.0, B, kC, thetaz*wt[i]);
    kC.addMatrixProduct(0.0, ks, Cy, 1.0);
    kb.addMatrixTransposeProduct(1.0, B, kC, thetay*wt[i]);

    Matrix ks1(1,order);
    static Matrix ksB(1,6);

    for (j = 0; j < order; j++) {
      if (code(j) == SECTION_RESPONSE_P) {
	for (int jj = 0; jj < order; jj++)
	  ks1(0,jj) = ks(j,jj);

	ksB.addMatrixProduct(0.0, ks1, B, 1.0);
	kb.addMatrixTransposeProduct(1.0, Cz1, ksB, thetaz*wt[i]);
	ksB.addMatrixProduct(0.0, ks1, Cz, 1.0);
	kb.addMatrixTransposeProduct(1.0, Cz1, ksB, thetaz*thetaz*wt[i]*L);
	ksB.addMatrixProduct(0.0, ks1, Cz, 1.0);
	kb.addMatrixTransposeProduct(1.0, Cy1, ksB, thetaz*thetay*wt[i]*L);	

	ksB.addMatrixProduct(0.0, ks1, B, 1.0);
	kb.addMatrixTransposeProduct(1.0, Cy1, ksB, thetay*wt[i]);
	ksB.addMatrixProduct(0.0, ks1, Cy, 1.0);
	kb.addMatrixTransposeProduct(1.0, Cy1, ksB, thetay*thetay*wt[i]*L);
	ksB.addMatrixProduct(0.0, ks1, Cy, 1.0);
	kb.addMatrixTransposeProduct(1.0, Cz1, ksB, thetay*thetaz*wt[i]*L);		
      }
    }

    continue;
    
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
	for (int jj = 0; jj < order; jj++) {
	  switch(code(jj)) {
	  case SECTION_RESPONSE_P:
	    tmp = ks(k,jj)*wt[i]*thetaz;
	    ka(k,1) += tmp*c1;
	    ka(k,2) += tmp*c2;
	    break;
	  }
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
}

const Matrix&
DispBeamColumnNL3d::getTangentStiff()
{
  static Matrix kb(6,6);

  this->getBasicStiff(kb);

  // Zero for integral
  q.Zero();
  
  // Update the transformation
  crdTransf->update();
  
  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();

  double L = crdTransf->getInitialLength();
  
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

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];
    double zeta = xi[i];

    double thetaz = (3*zeta*zeta-4*zeta+1)*v(1) + (3*zeta*zeta-2*zeta)*v(2);
    double thetay = (3*zeta*zeta-4*zeta+1)*v(3) + (3*zeta*zeta-2*zeta)*v(4);    

    // Get the section tangent stiffness and stress resultant
    const Vector &s = theSections[i]->getStressResultant();
        
    // Perform numerical integration
    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));
    double si;
    for (int j = 0; j < order; j++) {
      //si = s(j)*wts(i);
      si = s(j)*wt[i];
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	q(0) += si; break;
      case SECTION_RESPONSE_T:
	q(5) += si; break;	
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si;
	q(2) += (xi6-2.0)*si; 
	for (int jj = 0; jj < order; jj++) {
	  switch(code(jj)) {
	  case SECTION_RESPONSE_P:
	    q(1) += (3*zeta*zeta-4*zeta+1)*thetaz*s(jj)*wt[i]*L;
	    q(2) += (3*zeta*zeta-2*zeta)*thetaz*s(jj)*wt[i]*L;
	    break;
	  }
	}
	break;
      case SECTION_RESPONSE_MY:
	q(3) += (xi6-4.0)*si;
	q(4) += (xi6-2.0)*si; 
	for (int jj = 0; jj < order; jj++) {
	  switch(code(jj)) {
	  case SECTION_RESPONSE_P:
	    q(3) += (3*zeta*zeta-4*zeta+1)*thetay*s(jj)*wt[i]*L;
	    q(4) += (3*zeta*zeta-2*zeta)*thetay*s(jj)*wt[i]*L;
	    break;
	  }
	}
	break;	
      default:
	break;
      }
    }
    
  }
  
  // Add effects of element loads, q = q(v) + q0
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
DispBeamColumnNL3d::getInitialBasicStiff()
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
      default:
	break;
      }
    }
    
  }

  return kb;
}

const Matrix&
DispBeamColumnNL3d::getInitialStiff()
{
  const Matrix &kb = this->getInitialBasicStiff();

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);

  return K;
}

const Matrix&
DispBeamColumnNL3d::getMass()
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
DispBeamColumnNL3d::zeroLoad(void)
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
DispBeamColumnNL3d::addLoad(ElementalLoad *theLoad, double loadFactor)
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
    opserr << "DispBeamColumnNL3d::addLoad() -- load type unknown for element with tag: " << 
      this->getTag() << endln;
    return -1;
  }

  return 0;
}

int 
DispBeamColumnNL3d::addInertiaLoadToUnbalance(const Vector &accel)
{
	// Check for a quick return
	if (rho == 0.0) 
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = theNodes[0]->getRV(accel);
	const Vector &Raccel2 = theNodes[1]->getRV(accel);

    if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
      opserr << "DispBeamColumnNL3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
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
DispBeamColumnNL3d::getResistingForce()
{
  // Update the transformation
  crdTransf->update();
  
  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();

  double L = crdTransf->getInitialLength();
  
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
    double zeta = xi[i];

    double c1 = 3*zeta*zeta-4*zeta+1;
    double c2 = 3*zeta*zeta-2*zeta;
    double thetaz = c1*v(1) + c2*v(2);
    double thetay = c1*v(3) + c2*v(4);    
    
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
      case SECTION_RESPONSE_T:
	q(5) += si; break;
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si; 
	q(2) += (xi6-2.0)*si; 
	for (int jj = 0; jj < order; jj++) {
	  switch(code(jj)) {
	  case SECTION_RESPONSE_P:
	    q(1) += c1*thetaz*s(jj)*wt[i]*L;
	    q(2) += c2*thetaz*s(jj)*wt[i]*L;
	    break;
	  }
	}
	break;
      case SECTION_RESPONSE_MY:
	q(3) += (xi6-4.0)*si; 
	q(4) += (xi6-2.0)*si; 
	for (int jj = 0; jj < order; jj++) {
	  switch(code(jj)) {
	  case SECTION_RESPONSE_P:
	    q(3) += c1*thetay*s(jj)*wt[i]*L;
	    q(4) += c2*thetay*s(jj)*wt[i]*L;
	    break;
	  }
	}
	break;
      default:
	break;
      }
    }
    
  }
  
  // Add effects of element loads, q = q(v) + q0
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];    

  // Vector for reactions in basic system
  Vector p0Vec(p0, 5);

  P = crdTransf->getGlobalResistingForce(q, p0Vec);
  
  // Subtract other external nodal loads ... P_res = P_int - P_ext
  P.addVector(1.0, Q, -1.0);

  return P;
}

const Vector&
DispBeamColumnNL3d::getResistingForceIncInertia()
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
DispBeamColumnNL3d::sendSelf(int commitTag, Channel &theChannel)
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
    opserr << "DispBeamColumnNL3d::sendSelf() - failed to send ID data\n";
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
      opserr << "DispBeamColumnNL3d::sendSelf() - failed to send double data\n";
      return -1;
    }    
  }

  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "DispBeamColumnNL3d::sendSelf() - failed to send crdTranf\n";
     return -1;
  }      

  // send the beam integration
  if (beamInt->sendSelf(commitTag, theChannel) < 0) {
    opserr << "DispBeamColumnNL3d::sendSelf() - failed to send beamInt\n";
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
    opserr << "DispBeamColumnNL3d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  //
  // send the sections
  //
  
  for (j = 0; j<numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "DispBeamColumnNL3d::sendSelf() - section " << 
	j << "failed to send itself\n";
      return -1;
    }
  }

  return 0;
}

int
DispBeamColumnNL3d::recvSelf(int commitTag, Channel &theChannel,
			   FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i;
  
  static ID idData(9); // one bigger than needed so no clash with section ID

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    opserr << "DispBeamColumnNL3d::recvSelf() - failed to recv ID data\n";
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
      opserr << "DispBeamColumnNL3d::sendSelf() - failed to recv double data\n";
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
	opserr << "DispBeamColumnNL3d::recvSelf() - failed to obtain a CrdTrans object with classTag " <<
	  crdTransfClassTag << endln;
	  return -2;	  
      }
  }
  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "DispBeamColumnNL3d::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }      

  // create a new beamInt object if one needed
  if (beamInt == 0 || beamInt->getClassTag() != beamIntClassTag) {
      if (beamInt != 0)
	  delete beamInt;

      beamInt = theBroker.getNewBeamIntegration(beamIntClassTag);

      if (beamInt == 0) {
	opserr << "DispBeamColumnNL3d::recvSelf() - failed to obtain the beam integration object with classTag" <<
	  beamIntClassTag << endln;
	exit(-1);
      }
  }

  beamInt->setDbTag(beamIntDbTag);

  // invoke recvSelf on the beamInt object
  if (beamInt->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     opserr << "DispBeamColumnNL3d::sendSelf() - failed to recv beam integration\n";
     return -3;
  }      

  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumnNL3d::recvSelf() - failed to recv ID data\n";
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
opserr << "DispBeamColumnNL3d::recvSelf() - out of memory creating sections array of size " <<
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
	opserr << "DispBeamColumnNL3d::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumnNL3d::recvSelf() - section " << i << " failed to recv itself\n";
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
	opserr << "DispBeamColumnNL3d::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
	}
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumnNL3d::recvSelf() - section " << i << " failed to recv itself\n";
	return -1;
      }     
    }
  }

  return 0;
}

void
DispBeamColumnNL3d::Print(OPS_Stream &s, int flag)
{
  s << "\nDispBeamColumnNL3d, element id:  " << this->getTag() << endln;
  s << "\tConnected external nodes:  " << connectedExternalNodes;
  s << "\tCoordTransf: " << crdTransf->getTag() << endln;
  s << "\tmass density:  " << rho << endln;
  s << "\tNum sections:  " << numSections << endln;
  
  double L = crdTransf->getInitialLength();
  double P  = q(0);
  double M1 = q(1);
  double M2 = q(2);
  double V = (M1+M2)/L;
  s << "\tEnd 1 Forces (P V M): " << -P+p0[0]
    << " " << V+p0[1] << " " << M1 << endln;
  s << "\tEnd 2 Forces (P V M): " << P
    << " " << -V+p0[2] << " " << M2 << endln;

  beamInt->Print(s, flag);

  for (int i = 0; i < numSections; i++)
    theSections[i]->Print(s,flag);
}


int
DispBeamColumnNL3d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response*
DispBeamColumnNL3d::setResponse(const char **argv, int argc,
			      OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","DispBeamColumnNL3d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

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


    theResponse =  new ElementResponse(this, 1, P);
  
  
  // local force -
  } else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

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

    theResponse =  new ElementResponse(this, 2, P);
  

  // basic force -
  } else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M1");
    output.tag("ResponseType","M2");

    theResponse =  new ElementResponse(this, 9, Vector(6));

  // basic stiffness -
  } else if (strcmp(argv[0],"basicStiffness") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M1");
    output.tag("ResponseType","M2");

    theResponse =  new ElementResponse(this, 19, Matrix(6,6));

  // chord rotation -
  } else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0 
	     || strcmp(argv[0],"basicDeformation") == 0) {

    output.tag("ResponseType","eps");
    output.tag("ResponseType","thetaZ_1");
    output.tag("ResponseType","thetaZ_2");
    output.tag("ResponseType","thetaY_1");
    output.tag("ResponseType","thetaY_2");
    output.tag("ResponseType","thetaX");

    theResponse =  new ElementResponse(this, 3, Vector(6));
  
  // plastic rotation -
  } else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {

    output.tag("ResponseType","eps");
    output.tag("ResponseType","thetaZ_1");
    output.tag("ResponseType","thetaZ_2");
    output.tag("ResponseType","thetaY_1");
    output.tag("ResponseType","thetaY_2");
    output.tag("ResponseType","thetaX");

    theResponse =  new ElementResponse(this, 4, Vector(6));

  } else if (strcmp(argv[0],"RayleighForces") == 0 || strcmp(argv[0],"rayleighForces") == 0) {

    theResponse =  new ElementResponse(this, 12, P);
  }

  else if (strcmp(argv[0],"sectionDisplacements") == 0) {
    if (argc > 1 && strcmp(argv[1],"local") == 0)
      theResponse = new ElementResponse(this, 1111, Matrix(numSections,3));
    else
      theResponse = new ElementResponse(this, 111, Matrix(numSections,3));
  }
  
  // section response -
  else if (strcmp(argv[0],"sectionX") == 0) {
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
  else if (strcmp(argv[0],"section") == 0) {

    if (argc > 1) {
      
      int sectionNum = atoi(argv[1]);
      
      if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {

	output.tag("GaussPointOutput");
	output.attr("number",sectionNum);
	double xi[maxNumSections];
	double L = crdTransf->getInitialLength();
	beamInt->getSectionLocations(numSections, L, xi);
	output.attr("eta",xi[sectionNum-1]*L);

	if (strcmp(argv[2],"dsdh") != 0) {
	  theResponse = theSections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
	} else {
	  int order = theSections[sectionNum-1]->getOrder();
	  theResponse = new ElementResponse(this, 76, Vector(order));
	  Information &info = theResponse->getInformation();
	  info.theInt = sectionNum;
	}
	
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

  else if (strcmp(argv[0],"sectionTags") == 0)
    theResponse = new ElementResponse(this, 110, ID(numSections));

  if (theResponse == 0)
    theResponse = crdTransf->setResponse(argv, argc, output);
  
  output.endTag();
  return theResponse;
}

int 
DispBeamColumnNL3d::getResponse(int responseID, Information &eleInfo)
{
  double V;
  double L = crdTransf->getInitialLength();

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 12)
    return eleInfo.setVector(this->getRayleighDampingForces());

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

  else if (responseID == 19) {
    static Matrix kb(6,6);
    this->getBasicStiff(kb);
    return eleInfo.setMatrix(kb);
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

  else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = theSections[i]->getTag();
    return eleInfo.setID(tags);
  }
  
  else if (responseID == 111 || responseID == 1111) {
    double L = crdTransf->getInitialLength();
    double pts[maxNumSections];
    beamInt->getSectionLocations(numSections, L, pts);
    // CBDI influence matrix
    Matrix ls(numSections, numSections);
    getCBDIinfluenceMatrix(numSections, pts, L, ls);
    // Curvature vector
    Vector kappaz(numSections); // about section z
    Vector kappay(numSections); // about section y
    for (int i = 0; i < numSections; i++) {
      const ID &code = theSections[i]->getType();
      const Vector &e = theSections[i]->getSectionDeformation();
      int order = theSections[i]->getOrder();
      for (int j = 0; j < order; j++) {
	if (code(j) == SECTION_RESPONSE_MZ)
	  kappaz(i) += e(j);
	if (code(j) == SECTION_RESPONSE_MY)
	  kappay(i) += e(j);
      }
    }
    // Displacement vector
    Vector dispsy(numSections); // along local y
    Vector dispsz(numSections); // along local z    
    dispsy.addMatrixVector(0.0, ls, kappaz,  1.0);
    dispsz.addMatrixVector(0.0, ls, kappay, -1.0);    
    beamInt->getSectionLocations(numSections, L, pts);
    static Vector uxb(3);
    static Vector uxg(3);
    Matrix disps(numSections,3);
    static Vector vp(6);
    vp = crdTransf->getBasicTrialDisp();
    for (int i = 0; i < numSections; i++) {
      uxb(0) = pts[i]*vp(0); // linear shape function
      uxb(1) = dispsy(i);
      uxb(2) = dispsz(i);
      if (responseID == 111)
	uxg = crdTransf->getPointGlobalDisplFromBasic(pts[i],uxb);
      else
	uxg = crdTransf->getPointLocalDisplFromBasic(pts[i],uxb);
      disps(i,0) = uxg(0);
      disps(i,1) = uxg(1);
      disps(i,2) = uxg(2);            
    }
    return eleInfo.setMatrix(disps);
  }
  
  else
    return Element::getResponse(responseID, eleInfo);
}

int 
DispBeamColumnNL3d::getResponseSensitivity(int responseID, int gradNumber,
					   Information &eleInfo)
{
  // Basic deformation sensitivity
  if (responseID == 3) {  
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
    return eleInfo.setVector(dvdh);
  }

  // Basic force sensitivity
  else if (responseID == 9) {
    static Vector dqdh(6);

    dqdh.Zero();

    return eleInfo.setVector(dqdh);
  }

  // dsdh
  else if (responseID == 76) {

    int sectionNum = eleInfo.theInt;
    int order = theSections[sectionNum-1]->getOrder();
    const ID &code = theSections[sectionNum-1]->getType();

    Vector dsdh(order);
    dsdh = theSections[sectionNum-1]->getStressResultantSensitivity(gradNumber, true);

    const Vector &v = crdTransf->getBasicTrialDisp();
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0/L;

    const Matrix &ks = theSections[sectionNum-1]->getSectionTangent();

    Vector dedh(order);

    //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);

    double x = xi[sectionNum-1];

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*x;
    double zeta = x;

    double theta = (3*zeta*zeta-4*zeta+1)*v(1) + (3*zeta*zeta-2*zeta)*v(2);
    double dthetadh = (3*zeta*zeta-4*zeta+1)*dvdh(1) + (3*zeta*zeta-2*zeta)*dvdh(2);

    int j;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	dedh(j) = oneOverL*dvdh(0) + theta*dthetadh; break;
      case SECTION_RESPONSE_MZ:
	dedh(j) = oneOverL*((xi6-4.0)*dvdh(1) + (xi6-2.0)*dvdh(2)); break;
      default:
	dedh(j) = 0.0; break;
      }
    }

    dsdh.addMatrixVector(1.0, ks, dedh, 1.0);

    return eleInfo.setVector(dsdh);
  }

  else
    return -1;
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
DispBeamColumnNL3d::setParameter(const char **argv, int argc, Parameter &param)
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
    
    // Get section number: 1...Np
    int sectionNum = atoi(argv[1]);
    
    if (sectionNum > 0 && sectionNum <= numSections)
      return theSections[sectionNum-1]->setParameter(&argv[2], argc-2, param);
    
    else
      return -1;
  }
  
  else if (strstr(argv[0],"integration") != 0) {
    
    if (argc < 2)
      return -1;

    return beamInt->setParameter(&argv[1], argc-1, param);
  }
  int result =-1;
  // Default, send to every object
  int ok = 0;
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
DispBeamColumnNL3d::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  else
    return -1;  
}




int
DispBeamColumnNL3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}



const Matrix &
DispBeamColumnNL3d::getInitialStiffSensitivity(int gradNumber)
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

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];
    
    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getInitialTangentSensitivity(gradNumber);
    
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

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);
  
  return K;
}

const Matrix &
DispBeamColumnNL3d::getMassSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}



const Vector &
DispBeamColumnNL3d::getResistingForceSensitivity(int gradNumber)
{
  // Update the transformation
  crdTransf->update();
  
  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();

  double L = crdTransf->getInitialLength();
  
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
    double zeta = xi[i];

    double c1 = 3*zeta*zeta-4*zeta+1;
    double c2 = 3*zeta*zeta-2*zeta;
    double theta = c1*v(1) + c2*v(2);

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
	dqdh(0) += sensi; // 1
	dqdh(1) += c1*theta*sensi*L; // 2
	dqdh(2) += c2*theta*sensi*L; // 2
	break;
      case SECTION_RESPONSE_MZ:
	dqdh(1) += (xi6-4.0)*sensi; // 1
	dqdh(2) += (xi6-2.0)*sensi; // 1
	break;
      default:
	break;
      }
    }
  }
  


  if (crdTransf->isShapeSensitivity()) {
    double dLdh = crdTransf->getdLdh();

    double dptsdh[maxNumSections];
    beamInt->getLocationsDeriv(numSections, L, dLdh, dptsdh);

    double dwtsdh[maxNumSections];
    beamInt->getWeightsDeriv(numSections, L, dLdh, dwtsdh);

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

      double dwtLdh = wt[i]*dLdh + dwtsdh[i]*L;
      double dptLdh = xi[i]*dLdh + dptsdh[i]*L;

      int order = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();

      double xi6 = 6.0*xi[i];
      double zeta = xi[i];

      double c1 = 3*zeta*zeta-4*zeta+1;
      double c2 = 3*zeta*zeta-2*zeta;
      double theta = c1*v(1) + c2*v(2);

      Matrix B(order,3);
      Matrix dBdh(order,3);
      Matrix C(order,3);
      Matrix dCdh(order,3);
      Matrix C1(1,3);
      Matrix dC1dh(1,3);

      const Vector &s = theSections[i]->getStressResultant();
      const Matrix &ks = theSections[i]->getSectionTangent();

      double N = 0.0;

      for(int j = 0; j < order; j++) {
	switch(code(j)) {
	case SECTION_RESPONSE_P:
	  N += s(j);

	  B(j,0) = 1.0/L;
	  dBdh(j,0) = -dLdh/(L*L);

	  C(j,1) = c1;
	  C(j,2) = c2;
	  dCdh(j,1) = -(xi6-4.0)*xi[i]/L*dLdh + (xi6-4.0)/L*dptLdh;
	  dCdh(j,2) = -(xi6-2.0)*xi[i]/L*dLdh + (xi6-2.0)/L*dptLdh;

	  C1(0,1) = c1;
	  C1(0,2) = c2;
	  dC1dh(0,1) = -(xi6-4.0)*xi[i]/L*dLdh + (xi6-4.0)/L*dptLdh;
	  dC1dh(0,2) = -(xi6-2.0)*xi[i]/L*dLdh + (xi6-2.0)/L*dptLdh;
	  break;
	case SECTION_RESPONSE_MZ:
	  B(j,1) = (xi6-4.0)/L;
	  B(j,2) = (xi6-2.0)/L;
	  dBdh(j,1) = -(12*xi[i]-4.0)/(L*L)*dLdh + 6/(L*L)*dptLdh;
	  dBdh(j,2) = -(12*xi[i]-2.0)/(L*L)*dLdh + 6/(L*L)*dptLdh;
	  break;
	}
      }

      dqdh.addMatrixTransposeVector(1.0, dBdh, s, wt[i]*L); // 3

      dqdh(1) += dC1dh(0,1)*theta*N*wt[i]*L; // 4
      dqdh(2) += dC1dh(0,2)*theta*N*wt[i]*L; // 4

      double dCdhv = dC1dh(0,1)*v(1) + dC1dh(0,2)*v(2);
      dqdh(1) += C1(0,1)*dCdhv*N*wt[i]*L; // 5
      dqdh(2) += C1(0,2)*dCdhv*N*wt[i]*L; // 5

      Matrix tmp(order,3);
      tmp = dBdh;
      tmp.addMatrix(1.0, dCdh, theta);

      Matrix tmp2(3,3);
      tmp2.addMatrixTripleProduct(0.0, C, ks, tmp, 1.0);
      dqdh.addMatrixVector(1.0, tmp2, v, theta*wt[i]*L); // 7

      tmp2.addMatrixTripleProduct(0.0, B, ks, tmp, 1.0);
      dqdh.addMatrixVector(1.0, tmp2, v, wt[i]*L); // 6

      dqdh.addMatrixTransposeVector(1.0, B, s, dwtLdh); // 8

      dqdh(1) += C1(0,1)*theta*N*dwtLdh; // 9
      dqdh(2) += C1(0,2)*theta*N*dwtLdh; // 9

      continue;

      double si;
      for(int j = 0; j < order; j++) {
	switch(code(j)) {
	  si = s(j)*wt[i];
	case SECTION_RESPONSE_P:
	  dqdh(0) += -dLdh/(L*L);
	  dqdh(1) += (-6*xi[i]+4.0)*xi[i]/L*theta*si*L*dLdh;
	  dqdh(2) += (-6*xi[i]+2.0)*xi[i]/L*theta*si*L*dLdh;
	  break;
	case SECTION_RESPONSE_MZ:
	  dqdh(1) += (-12*xi[i]+4.0)/(L*L)*si*L*dLdh; 
	  dqdh(2) += (-12*xi[i]+2.0)/(L*L)*si*L*dLdh; 
	  break;
	default:
	  break;
	}
      }

    }

  }


  // Transform forces
  static Vector dp0dh(5);		// No distributed loads
  dp0dh.Zero();

  P.Zero();

  //////////////////////////////////////////////////////////////

  if (crdTransf->isShapeSensitivity()) {
    
    // Perform numerical integration to obtain basic stiffness matrix
    // Some extra declarations
    static Matrix kbmine(6,6);
    this->getBasicStiff(kbmine);
    
    // k dAdh u
    const Vector &dAdh_u = crdTransf->getBasicTrialDispShapeSensitivity();
    //dqdh.addMatrixVector(1.0, kbmine, dAdh_u, oneOverL);
    dqdh.addMatrixVector(1.0, kbmine, dAdh_u, 1.0);

    // dAdh^T q
    P += crdTransf->getGlobalResistingForceShapeSensitivity(q, dp0dh, gradNumber);
  }
  
  // A^T (dqdh + k dAdh u)
  P += crdTransf->getGlobalResistingForce(dqdh, dp0dh);
  
  return P;
}



// NEW METHOD
int
DispBeamColumnNL3d::commitSensitivity(int gradNumber, int numGrads)
{
  // Get basic deformation and sensitivities
  const Vector &v = crdTransf->getBasicTrialDisp();
  
  static Vector dvdh(6);
  dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);

  // Some extra declarations
  double d1oLdh = crdTransf->getd1overLdh();
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    
    Vector e(workArea, order);
    
    double xi6 = 6.0*xi[i];
    
    double zeta = xi[i];

    double c1 = 3*zeta*zeta-4*zeta+1;
    double c2 = 3*zeta*zeta-2*zeta;
    double theta = c1*v(1) + c2*v(2);

    for (int j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	e(j) = oneOverL*dvdh(0)
	  + d1oLdh*v(0)
	  + theta*(c1*dvdh(1)+c2*dvdh(2)); 
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


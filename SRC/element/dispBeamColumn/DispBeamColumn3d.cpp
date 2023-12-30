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
// $URL$

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for DispBeamColumn3d.

#include <DispBeamColumn3d.h>
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
#include <elementAPI.h>
#include <string>

Matrix DispBeamColumn3d::K(12,12);
Vector DispBeamColumn3d::P(12);
double DispBeamColumn3d::workArea[200];

void* OPS_DispBeamColumn3d()
{
    int dampingTag = 0;
    Damping* theDamping = 0;
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
    else if (strcmp(type, "-damp") == 0) {

        if (OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
            theDamping = OPS_getDamping(dampingTag);
            if (theDamping == 0) {
                opserr << "damping not found\n";
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
    
    Element *theEle =  new DispBeamColumn3d(iData[0],iData[1],iData[2],secTags.Size(),sections,
					    *bi,*theTransf,mass,cmass, theDamping);
    delete [] sections;
    return theEle;
}

void *OPS_DispBeamColumn3d(const ID &info) {
    // data needed
    int iData[5];
    double mass = 0.0;
    int cmass = 0;
    int numData;
    int dampingTag = 0;
    Damping* theDamping = 0;
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if (ndm != 3 || ndf != 6) {
        opserr << "ndm must be 3 and ndf must be 6\n";
        return 0;
    }

    // 1. regular elements
    if (info.Size() == 0) {
        numData = 3;
        if (OPS_GetNumRemainingInputArgs() < numData) {
            opserr << "insufficient "
                      "arguments:eleTag,iNode,jNode\n";
            return 0;
        }
        if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
            opserr << "WARNING invalid int inputs\n";
            return 0;
        }
    }

    // 2. regular elements or save data
    if (info.Size() == 0 || info(0) == 1) {
        numData = 2;
        if (OPS_GetNumRemainingInputArgs() < numData) {
            opserr << "insufficient "
                      "arguments:transfTag,integrationTag\n";
            return 0;
        }
        if (OPS_GetIntInput(&numData, &iData[3]) < 0) {
            opserr << "WARNING invalid int inputs\n";
            return 0;
        }

        numData = 1;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char *type = OPS_GetString();
            if (strcmp(type, "-cMass") == 0) {
                cmass = 1;
            } else if (strcmp(type, "-mass") == 0) {
                if (OPS_GetNumRemainingInputArgs() > 0) {
                    if (OPS_GetDoubleInput(&numData, &mass) < 0) {
                        opserr << "WARNING: invalid mass\n";
                        return 0;
                    }
                }
            }
            else if (strcmp(type, "-damp") == 0) {

                if (OPS_GetNumRemainingInputArgs() > 0) {
                    if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
                    theDamping = OPS_getDamping(dampingTag);
                    if (theDamping == 0) {
                        opserr << "damping not found\n";
                        return 0;
                    }
                }
            }
        }
    }

    // 3: save data
    static std::map<int, Vector> meshdata;
    if (info.Size() > 0 && info(0) == 1) {
        if (info.Size() < 2) {
            opserr << "WARNING: need info -- inmesh, meshtag\n";
            return 0;
        }

        // save the data for a mesh
        Vector &mdata = meshdata[info(1)];
        mdata.resize(4);
        mdata(0) = iData[3];
        mdata(1) = iData[4];
        mdata(2) = mass;
        mdata(3) = cmass;
        return &meshdata;
    }

    // 4: load data
    if (info.Size() > 0 && info(0) == 2) {
        if (info.Size() < 5) {
            opserr << "WARNING: need info -- inmesh, meshtag, "
                      "eleTag, nd1, nd2\n";
            return 0;
        }

        // get the data for a mesh
        Vector &mdata = meshdata[info(1)];
        if (mdata.Size() < 5) return 0;

        iData[0] = info(2);
        iData[1] = info(3);
        iData[2] = info(4);
        iData[3] = mdata(0);
        iData[4] = mdata(1);
        mass = mdata(2);
        cmass = mdata(3);
    }

    // 5: create element
    CrdTransf *theTransf = OPS_getCrdTransf(iData[3]);
    if (theTransf == 0) {
        opserr << "coord transfomration not found\n";
        return 0;
    }

    // check beam integrataion
    BeamIntegrationRule *theRule =
        OPS_getBeamIntegrationRule(iData[4]);
    if (theRule == 0) {
        opserr << "beam integration not found\n";
        return 0;
    }
    BeamIntegration *bi = theRule->getBeamIntegration();
    if (bi == 0) {
        opserr << "beam integration is null\n";
        return 0;
    }

    // check sections
    const ID &secTags = theRule->getSectionTags();
    SectionForceDeformation **sections =
        new SectionForceDeformation *[secTags.Size()];
    for (int i = 0; i < secTags.Size(); i++) {
        sections[i] = OPS_getSectionForceDeformation(secTags(i));
        if (sections[i] == 0) {
            opserr << "section " << secTags(i) << "not found\n";
            delete[] sections;
            return 0;
        }
    }

    Element *theEle = new DispBeamColumn3d(
        iData[0], iData[1], iData[2], secTags.Size(), sections, *bi,
        *theTransf, mass, cmass, theDamping);
    delete[] sections;
    return theEle;
}

DispBeamColumn3d::DispBeamColumn3d(int tag, int nd1, int nd2,
				   int numSec, SectionForceDeformation **s,
				   BeamIntegration &bi,
				   CrdTransf &coordTransf, double r, int cm,
				   Damping *damping)
:Element (tag, ELE_TAG_DispBeamColumn3d),
numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
connectedExternalNodes(2), 
Q(12), q(6), rho(r), cMass(cm), parameterID(0), theDamping(0)
{
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];
  
  if (theSections == 0) {
    opserr << "DispBeamColumn3d::DispBeamColumn3d - failed to allocate section model pointer\n";
    exit(-1);
  }
  
  for (int i = 0; i < numSections; i++) {
    
    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getCopy();
    
    // Check allocation
    if (theSections[i] == 0) {
      opserr << "DispBeamColumn3d::DispBeamColumn3d -- failed to get a copy of section model\n";
      exit(-1);
    }
  }
  
  beamInt = bi.getCopy();
  
  if (beamInt == 0) {
    opserr << "DispBeamColumn3d::DispBeamColumn3d - failed to copy beam integration\n";
    exit(-1);
  }

  crdTransf = coordTransf.getCopy3d();
  
  if (crdTransf == 0) {
    opserr << "DispBeamColumn3d::DispBeamColumn3d - failed to copy coordinate transformation\n";
    exit(-1);
  }
  
  if (damping)
  {
    theDamping =(*damping).getCopy();
    
    if (!theDamping) {
      opserr << "DispBeamColumn3d::DispBeamColumn3d - failed to copy damping\n";
      exit(-1);
    }
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

DispBeamColumn3d::DispBeamColumn3d()
:Element (0, ELE_TAG_DispBeamColumn3d),
numSections(0), theSections(0), crdTransf(0), beamInt(0),
connectedExternalNodes(2), 
Q(12), q(6), rho(0.0), cMass(0), parameterID(0),
theDamping(0)
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

DispBeamColumn3d::~DispBeamColumn3d()
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

	if (theDamping) delete theDamping;
}

int
DispBeamColumn3d::getNumExternalNodes() const
{
    return 2;
}

const ID&
DispBeamColumn3d::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
DispBeamColumn3d::getNodePtrs()
{

    return theNodes;
}

int
DispBeamColumn3d::getNumDOF()
{
    return 12;
}

void
DispBeamColumn3d::setDomain(Domain *theDomain)
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
	//opserr << "FATAL ERROR DispBeamColumn3d (tag: %d), node not found in domain",
	//	this->getTag());
	
	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    if (dofNd1 != 6 || dofNd2 != 6) {
	//opserr << "FATAL ERROR DispBeamColumn3d (tag: %d), has differing number of DOFs at its nodes",
	//	this->getTag());
	
	return;
    }

	if (crdTransf->initialize(theNodes[0], theNodes[1])) {
		// Add some error check
	}

  // initialize the damping
  if (theDamping && theDamping->setDomain(theDomain, 6)) {
    opserr << "DispBeamColumn3d::setDomain(): Error initializing damping";  
    exit(0);
  }

	double L = crdTransf->getInitialLength();

	if (L == 0.0) {
		// Add some error check
	}

    this->DomainComponent::setDomain(theDomain);

	this->update();
}

int
DispBeamColumn3d::setDamping(Domain *theDomain, Damping *damping)
{
  if (theDomain && damping)
  {
    if (theDamping) delete theDamping;

    theDamping =(*damping).getCopy();
    
    if (!theDamping) {
      opserr << "DispBeamColumn3d::setDamping -- failed to get copy of damping\n";
      return -1;
    }
    if (theDamping->setDomain(theDomain, 6)) {
      opserr << "DispBeamColumn3d::setDamping -- Error initializing damping\n";
      return -2;
    }
  }
  
  return 0;
}

int
DispBeamColumn3d::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "DispBeamColumn3d::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    if (theDamping) retVal += theDamping->commitState();

    return retVal;
}

int
DispBeamColumn3d::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    if (theDamping) retVal += theDamping->revertToLastCommit();

    return retVal;
}

int
DispBeamColumn3d::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    if (theDamping) retVal += theDamping->revertToStart();

    return retVal;
}

int
DispBeamColumn3d::update(void)
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
      
    double xi6 = 6.0*xi[i];
    
    int j;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	e(j) = oneOverL*v(0);
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
    err += theSections[i]->setTrialSectionDeformation(e);
  }

  if (err != 0) {
    opserr << "DispBeamColumn3d::update() - failed setTrialSectionDeformations()\n";
    return err;
  }
  return 0;
}

const Matrix&
DispBeamColumn3d::getTangentStiff()
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
  //   opserr << this->getTag() << " " << K;
  return K;
}

void
DispBeamColumn3d::getBasicStiff(Matrix &kb, int initial)
{
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
    const Matrix &ks = (initial) ? theSections[i]->getInitialTangent() : theSections[i]->getSectionTangent();
    
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
  
  if(theDamping) kb *= theDamping->getStiffnessMultiplier();
}

const Matrix&
DispBeamColumn3d::getInitialStiff()
{
  static Matrix kb(6,6);

  this->getBasicStiff(kb, 1);

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);
  
  return K;
}

const Matrix&
DispBeamColumn3d::getMass()
{
  K.Zero();
  
  if (rho == 0.0)
    return K;
  
  double L = crdTransf->getInitialLength();
  if (cMass == 0)  {
    // lumped mass matrix
    double m = 0.5*rho*L;
    K(0,0) = K(1,1) = K(2,2) = K(6,6) = K(7,7) = K(8,8) = m;
  } else  {
    // consistent mass matrix
    static Matrix ml(12,12);
    double m = rho*L/420.0;
    ml(0,0) = ml(6,6) = m*140.0;
    ml(0,6) = ml(6,0) = m*70.0;
    //ml(3,3) = ml(9,9) = m*(Jx/A)*140.0;  // CURRENTLY NO TORSIONAL MASS 
    //ml(3,9) = ml(9,3) = m*(Jx/A)*70.0;   // CURRENTLY NO TORSIONAL MASS
    
    ml(2,2) = ml(8,8) = m*156.0;
    ml(2,8) = ml(8,2) = m*54.0;
    ml(4,4) = ml(10,10) = m*4.0*L*L;
    ml(4,10) = ml(10,4) = -m*3.0*L*L;
    ml(2,4) = ml(4,2) = -m*22.0*L;
    ml(8,10) = ml(10,8) = -ml(2,4);
    ml(2,10) = ml(10,2) = m*13.0*L;
    ml(4,8) = ml(8,4) = -ml(2,10);
    
    ml(1,1) = ml(7,7) = m*156.0;
    ml(1,7) = ml(7,1) = m*54.0;
    ml(5,5) = ml(11,11) = m*4.0*L*L;
    ml(5,11) = ml(11,5) = -m*3.0*L*L;
    ml(1,5) = ml(5,1) = m*22.0*L;
    ml(7,11) = ml(11,7) = -ml(1,5);
    ml(1,11) = ml(11,1) = -m*13.0*L;
    ml(5,7) = ml(7,5) = -ml(1,11);
    
    // transform local mass matrix to global system
    K = crdTransf->getGlobalMatrixFromLocal(ml);
  }
  
  return K;
}

void
DispBeamColumn3d::zeroLoad(void)
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
DispBeamColumn3d::addLoad(ElementalLoad *theLoad, double loadFactor)
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
    opserr << "DispBeamColumn3d::addLoad() -- load type unknown for element with tag: " << 
      this->getTag() << endln;
    return -1;
  }

  return 0;
}

int 
DispBeamColumn3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0) 
    return 0;
  
  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "DispBeamColumn3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }
  
  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;

    Q(0) -= m*Raccel1(0);
    Q(1) -= m*Raccel1(1);
    Q(2) -= m*Raccel1(2);
    Q(6) -= m*Raccel2(0);
    Q(7) -= m*Raccel2(1);
    Q(8) -= m*Raccel2(2);

  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector Raccel(12);
    for (int i=0; i<6; i++)  {
      Raccel(i)   = Raccel1(i);
      Raccel(i+6) = Raccel2(i);
    }
    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }
  
  return 0;
}

const Vector&
DispBeamColumn3d::getResistingForce()
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
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    double xi6 = 6.0*xi[i];
    
    // Get section stress resultant
    const Vector &s = theSections[i]->getStressResultant();
    
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

  if (theDamping) theDamping->update(q);

  // Transform forces
  Vector p0Vec(p0, 5);
  P = crdTransf->getGlobalResistingForce(q, p0Vec);

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  if (rho != 0)
    P.addVector(1.0, Q, -1.0);
  
  return P;
}

const Vector &
DispBeamColumn3d::getDampingForce(void)
{
  crdTransf->update();

  return crdTransf->getGlobalResistingForce(theDamping->getDampingForce(), Vector(5));
}

const Vector&
DispBeamColumn3d::getResistingForceIncInertia()
{
  P = this->getResistingForce();
  
  if (theDamping) P += this->getDampingForce();
  
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;
  
    P(0) += m*accel1(0);
    P(1) += m*accel1(1);
    P(2) += m*accel1(2);
    P(6) += m*accel2(0);
    P(7) += m*accel2(1);
    P(8) += m*accel2(2);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector accel(12);
    for (int i=0; i<6; i++)  {
      accel(i)   = accel1(i);
      accel(i+6) = accel2(i);
    }
    P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }
    
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
DispBeamColumn3d::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j;
  int loc = 0;
  
  static Vector data(16);
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
  data(9) = cMass;
  data(10) = alphaM;
  data(11) = betaK;
  data(12) = betaK0;
  data(13) = betaKc;
  
  data(14) = 0;
  data(15) = 0;
  if (theDamping) {
    data(14) = theDamping->getClassTag();
    int dbTag = theDamping->getDbTag();
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	      theDamping->setDbTag(dbTag);
	  }
    data(15) = dbTag;
  }

  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "DispBeamColumn3d::sendSelf() - failed to send data Vector\n";
     return -1;
  }    
  
  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
     opserr << "DispBeamColumn3d::sendSelf() - failed to send crdTranf\n";
     return -1;
  }      

  // send the beam integration
  if (beamInt->sendSelf(commitTag, theChannel) < 0) {
    opserr << "DispBeamColumn3d::sendSelf() - failed to send beamInt\n";
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
    opserr << "DispBeamColumn3d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  //
  // send the sections
  //
  
  for (j = 0; j<numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "DispBeamColumn3d::sendSelf() - section " << j << "failed to send itself\n";
      return -1;
    }
  }

  // Ask the Damping to send itself
  if (theDamping && theDamping->sendSelf(commitTag, theChannel) < 0) {
      opserr << "DispBeamColumn3d::sendSelf -- could not send Damping\n";
      return -1;
  }

  return 0;
}

int
DispBeamColumn3d::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i;
  
  static Vector data(16);

  if (theChannel.recvVector(dbTag, commitTag, data) < 0)  {
    opserr << "DispBeamColumn3d::recvSelf() - failed to recv data Vector\n";
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
  cMass = (int)data(9);
  
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
	opserr << "DispBeamColumn3d::recvSelf() - " <<
	  "failed to obtain a CrdTrans object with classTag" <<
	  crdTransfClassTag << endln;
	return -2;	  
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "DispBeamColumn3d::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }      

  // create a new beamInt object if one needed
  if (beamInt == 0 || beamInt->getClassTag() != beamIntClassTag) {
      if (beamInt != 0)
	  delete beamInt;

      beamInt = theBroker.getNewBeamIntegration(beamIntClassTag);

      if (beamInt == 0) {
	opserr << "DispBeamColumn3d::recvSelf() - failed to obtain the beam integration object with classTag" <<
	  beamIntClassTag << endln;
	exit(-1);
      }
  }

  beamInt->setDbTag(beamIntDbTag);

  // invoke recvSelf on the beamInt object
  if (beamInt->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     opserr << "DispBeamColumn3d::sendSelf() - failed to recv beam integration\n";
     return -3;
  }      
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*nSect);
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumn3d::recvSelf() - failed to recv ID data\n";
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
      opserr << "DispBeamColumn3d::recvSelf() - out of memory creating sections array of size" <<
	nSect << endln;
      exit(-1);
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
	opserr << "DispBeamColumn3d::recvSelf() - Broker could not create Section of class type" <<
	  sectClassTag << endln;
	exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumn3d::recvSelf() - section " <<
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
	  opserr << "DispBeamColumn3d::recvSelf() - Broker could not create Section of class type" <<
	    sectClassTag << endln;
	  exit(-1);
	}
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "DispBeamColumn3d::recvSelf() - section " << 
	  i << "failed to recv itself\n";
	return -1;
      }     
    }
  }


  // Check if the Damping is null; if so, get a new one
  int dmpTag = (int)data(14);
  if (dmpTag) {
    if (theDamping == 0) {
      theDamping = theBroker.getNewDamping(dmpTag);
      if (theDamping == 0) {
        opserr << "DispBeamColumn3d::recvSelf -- could not get a Damping\n";
        exit(-1);
      }
    }
  
    // Check that the Damping is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (theDamping->getClassTag() != dmpTag) {
      delete theDamping;
      theDamping = theBroker.getNewDamping(dmpTag);
      if (theDamping == 0) {
        opserr << "DispBeamColumn3d::recvSelf -- could not get a Damping\n";
        exit(-1);
      }
    }
  
    // Now, receive the Damping
    theDamping->setDbTag((int)data(15));
    if (theDamping->recvSelf(commitTag, theChannel, theBroker) < 0) {
      opserr << "DispBeamColumn3d::recvSelf -- could not receive Damping\n";
      exit(-1);
    }
  }
  else {
    if (theDamping) {
      delete theDamping;
      theDamping = 0;
    }
  }
    
  return 0;
}

void
DispBeamColumn3d::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "\nDispBeamColumn3d, element id:  " << this->getTag() << endln;
		s << "\tConnected external nodes:  " << connectedExternalNodes;
		s << "\tCoordTransf: " << crdTransf->getTag() << endln;
		s << "\tmass density:  " << rho << ", cMass: " << cMass << endln;

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
		s << "Number of sections: " << numSections << endln;
		beamInt->Print(s, flag);

		for (int i = 0; i < numSections; i++) {
		  //opserr << "Section Type: " << theSections[i]->getClassTag() << endln;
		  theSections[i]->Print(s,flag);
		}
		//  if (rho != 0)
		//    opserr << "Mass: \n" << this->getMass();
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"DispBeamColumn3d\", ";
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
DispBeamColumn3d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numModes)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response*
DispBeamColumn3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","DispBeamColumn3d");
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
    }
    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {
      output.tag("ResponseType","N");
      output.tag("ResponseType","M1");
      output.tag("ResponseType","M2");

      theResponse = new ElementResponse(this, 9, Vector(6));
    }
    else if (strcmp(argv[0],"basicStiffness") == 0) {
      output.tag("ResponseType","N");
      output.tag("ResponseType","M1");
      output.tag("ResponseType","M2");

      theResponse = new ElementResponse(this, 19, Matrix(6,6));
    // global damping force - 
    } else if (theDamping && (strcmp(argv[0],"globalDampingForce") == 0 || strcmp(argv[0],"globalDampingForces") == 0)) {

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


      theResponse = new ElementResponse(this, 21, P);

    // local damping force -
    } else if (theDamping && (strcmp(argv[0],"localDampingForce") == 0 || strcmp(argv[0],"localDampingForces") == 0)) {

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

      theResponse = new ElementResponse(this, 22, P);

    } else if (theDamping && (strcmp(argv[0],"basicDampingForce") == 0 || strcmp(argv[0],"basicDampingForces") == 0)) {

      output.tag("ResponseType","N");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Mz_2");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","T");
    
      theResponse = new ElementResponse(this, 23, Vector(6));

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
  

  } else if (strcmp(argv[0],"RayleighForces") == 0 || strcmp(argv[0],"rayleighForces") == 0) {

    theResponse =  new ElementResponse(this, 12, P);

  }   
    else if (strcmp(argv[0],"integrationPoints") == 0)
      theResponse = new ElementResponse(this, 10, Vector(numSections));

    else if (strcmp(argv[0],"integrationWeights") == 0)
      theResponse = new ElementResponse(this, 11, Vector(numSections));

    else if (strcmp(argv[0],"sectionTags") == 0)
      theResponse = new ElementResponse(this, 110, ID(numSections));

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
    // by SAJalali
    else if (strcmp(argv[0], "energy") == 0) {
      theResponse = new ElementResponse(this, 13, 0.0);
    }
    
    if (theResponse == 0)
      theResponse = crdTransf->setResponse(argv, argc, output);
    
  output.endTag();
  return theResponse;
}

int 
DispBeamColumn3d::getResponse(int responseID, Information &eleInfo)
{
  double N, V, M1, M2, T;
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 12)
    return eleInfo.setVector(this->getRayleighDampingForces());
    
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
    V = (M1+M2)*oneOverL;
    P(2) = -V+p0[3];
    P(8) =  V+p0[4];

    return eleInfo.setVector(P);
  }

  else if (responseID == 21)
    return eleInfo.setVector(this->getDampingForce());

  else if (responseID == 22) {
    Vector Sd(6);
    Sd = theDamping->getDampingForce();
    // Axial
    N = Sd(0);
    P(6) =  N;
    P(0) = -N;
    
    // Torsion
    T = Sd(5);
    P(9) =  T;
    P(3) = -T;
    
    // Moments about z and shears along y
    M1 = Sd(1);
    M2 = Sd(2);
    P(5)  = M1;
    P(11) = M2;
    V = (M1+M2)*oneOverL;
    P(1) =  V;
    P(7) = -V;
    
    // Moments about y and shears along z
    M1 = Sd(3);
    M2 = Sd(4);
    P(4)  = M1;
    P(10) = M2;
    V = (M1+M2)*oneOverL;
    P(2) = -V;
    P(8) =  V;
    return eleInfo.setVector(P);
  }

  else if (responseID == 23)
    return eleInfo.setVector(theDamping->getDampingForce());

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
    static Matrix kb(6,6);
    this->getBasicStiff(kb,1);
    kb.Solve(q, ve);
    vp = crdTransf->getBasicTrialDisp();
    vp -= ve;
    return eleInfo.setVector(vp);
  }

  else if (responseID == 10) {
    double L = crdTransf->getInitialLength();
    double pts[maxNumSections];
    beamInt->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i]*L;
    return eleInfo.setVector(locs);
  }

  else if (responseID == 11) {
    double L = crdTransf->getInitialLength();
    double wts[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wts);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wts[i]*L;
    return eleInfo.setVector(weights);
  }

  else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = theSections[i]->getTag();
    return eleInfo.setID(tags);
  }
  
  //by SAJalali
  else if (responseID == 13) {
	  double xi[maxNumSections];
	  double L = crdTransf->getInitialLength();
	  beamInt->getSectionWeights(numSections, L, xi);
	  double energy = 0;
	  for (int i = 0; i < numSections; i++) {
		  energy += theSections[i]->getEnergy()*xi[i] * L;
	  }
	  return eleInfo.setDouble(energy);
  }

  else
    return -1;
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
DispBeamColumn3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // don't do anything if MaterialStageParameter calls this element
  if (strcmp(argv[0],"updateMaterialStage") == 0) {
      return -1;
  }
  
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
    
    // Get section number
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
DispBeamColumn3d::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  else
    return -1;  
}




int
DispBeamColumn3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}

const Matrix &
DispBeamColumn3d::getInitialStiffSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}

const Matrix &
DispBeamColumn3d::getMassSensitivity(int gradNumber)
{
  K.Zero();
  
  if (rho == 0.0 || parameterID != 1)
    return K;
  
  double L = crdTransf->getInitialLength();
  if (cMass == 0)  {
    // lumped mass matrix
    //double m = 0.5*rho*L;
    double m = 0.5*L;
    K(0,0) = K(1,1) = K(2,2) = K(6,6) = K(7,7) = K(8,8) = m;
  } else  {
    // consistent mass matrix
    static Matrix ml(12,12);
    //double m = rho*L/420.0;
    double m = L/420.0;
    ml(0,0) = ml(6,6) = m*140.0;
    ml(0,6) = ml(6,0) = m*70.0;
    //ml(3,3) = ml(9,9) = m*(Jx/A)*140.0;  // CURRENTLY NO TORSIONAL MASS 
    //ml(3,9) = ml(9,3) = m*(Jx/A)*70.0;   // CURRENTLY NO TORSIONAL MASS
    
    ml(2,2) = ml(8,8) = m*156.0;
    ml(2,8) = ml(8,2) = m*54.0;
    ml(4,4) = ml(10,10) = m*4.0*L*L;
    ml(4,10) = ml(10,4) = -m*3.0*L*L;
    ml(2,4) = ml(4,2) = -m*22.0*L;
    ml(8,10) = ml(10,8) = -ml(2,4);
    ml(2,10) = ml(10,2) = m*13.0*L;
    ml(4,8) = ml(8,4) = -ml(2,10);
    
    ml(1,1) = ml(7,7) = m*156.0;
    ml(1,7) = ml(7,1) = m*54.0;
    ml(5,5) = ml(11,11) = m*4.0*L*L;
    ml(5,11) = ml(11,5) = -m*3.0*L*L;
    ml(1,5) = ml(5,1) = m*22.0*L;
    ml(7,11) = ml(11,7) = -ml(1,5);
    ml(1,11) = ml(11,1) = -m*13.0*L;
    ml(5,7) = ml(7,5) = -ml(1,11);
    
    // transform local mass matrix to global system
    K = crdTransf->getGlobalMatrixFromLocal(ml);
  }
  
  return K;
}



const Vector &
DispBeamColumn3d::getResistingForceSensitivity(int gradNumber)
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
DispBeamColumn3d::commitSensitivity(int gradNumber, int numGrads)
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


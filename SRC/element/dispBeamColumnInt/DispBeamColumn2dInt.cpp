// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumnInt/DispBeamColumn2dInt.cpp,v $

// $Version$

// $Date: 2008-04-14 21:22:20 $



// Created: 07/04

// Modified by: LMS 

// Description: This file contains the class implementation of DispBeamColumn2dInt.Based on DispBeamColumn2d.cpp.





#include "DispBeamColumn2dInt.h"

#include <Node.h>

#include "FiberSection2dInt.h"

#include "LinearCrdTransf2dInt.h"

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

#include <ElementalLoad.h>

#include <LegendreBeamIntegration.h>

#include <elementAPI.h>

#include <string>

#include <vector>

void* OPS_DispBeamColumn2dInt()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    int ok = 0;
    if (ndm == 2 && ndf == 3)
	ok = 1;

    if (ok == 0) {
	opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
	       << " not compatible with dispBeamColumn element" << endln;
	return 0;
    }

    if (OPS_GetNumRemainingInputArgs() < 7) {			//8
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag? C1? t1? NStrip1? t2? NStrip2? t3? NStrip3?\n";
	return 0;
    }

    // get the id and end nodes 
    int eleTag, iNode, jNode, nIP, transfTag;
    double C1;
    int secTag[10]; // Max size of integration rule ... can change if needed

    int idata[4];
    int numdata = 4;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid dispBeamColumn int inputs" << endln;
	return 0;
    }
    eleTag = idata[0];
    iNode = idata[1];
    jNode = idata[2];
    nIP = idata[3];

    const char* type = OPS_GetString();
    if (strcmp(type, "-sections") == 0) {
	if (nIP > OPS_GetNumRemainingInputArgs()) {
	    opserr << "WARNING insufficient number of section tags - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	    return 0;
	}
	int section;
	numdata = 1;
	for (int i = 0; i < nIP; i++) {
	    if (OPS_GetIntInput(&numdata, &section) < 0) {
		opserr << "WARNING invalid secTag - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
		return 0;
	    }
	    secTag[i] = section;
	}
    }
	
    else {
	OPS_ResetCurrentInputArg(-1);
	int section;
	numdata = 1;
	if (OPS_GetIntInput(&numdata, &section) < 0) {
	    opserr << "WARNING invalid secTag - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	    return 0;
	}
	for (int i = 0; i < nIP; i++)
	    secTag[i] = section;
    }
	
    if (OPS_GetNumRemainingInputArgs() > 0) {
	numdata = 1;
	if (OPS_GetIntInput(&numdata, &transfTag) < 0) {
	    opserr << "WARNING invalid transfTag? - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	    return 0;
	}
    }

    numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &C1) < 0) {
	opserr << "WARNING invalid dispBeamColumn C1" << endln;
	return 0;
    }


    double massDens = 0.0;

    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* massarg = OPS_GetString();
	if (strcmp(massarg,"-mass") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
	    numdata = 1;
	    if (OPS_GetDoubleInput(&numdata, &massDens) < 0) {
		opserr << "WARNING invalid massDens - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag? C1? t? NStrip?\n";
		return 0;
	    }
	}
    }

    SectionForceDeformation **sections = new SectionForceDeformation* [nIP];
	
    if (!sections) {
	opserr << "WARNING TclElmtBuilder - addFrameElement - Insufficient memory to create sections\n";
	return 0;
    }
	
    for (int j=0; j<nIP; j++) {
	SectionForceDeformation *theSection = OPS_getSectionForceDeformation(secTag[j]);
	  
	if (theSection == 0) {
	    opserr << "WARNING TclElmtBuilder - frameElement - no Section found with tag ";
	    opserr << secTag[j] << endln;
	    delete [] sections;
	    return 0;
	}

	sections[j] = theSection;
    }
	
    Element *theElement = 0;

    if (ndm == 2) {

	CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
      
	if (theTransf == 0) {
	    opserr << "WARNING transformation not found\n";
	    opserr << "transformation: " << transfTag;
	    opserr << "\ndispBeamColumn element: " << eleTag << endln;
	    return 0;
	}

	// now create the DispBeamColumn and add it to the Domain
	theElement = new DispBeamColumn2dInt(eleTag,iNode,jNode,nIP,sections,*theTransf,C1,massDens);

	delete [] sections;
    }

    // if get here we have successfully created the element and added it to the domain
    return theElement;
}



Matrix DispBeamColumn2dInt::K(6,6);

Vector DispBeamColumn2dInt::P(6);

double DispBeamColumn2dInt::workArea[100];

LegendreBeamIntegration DispBeamColumn2dInt::quadRule;



DispBeamColumn2dInt::DispBeamColumn2dInt(int tag, 

					 int nd1, 

					 int nd2,	

					 int numSec, 

					 SectionForceDeformation **s,

					 CrdTransf &coordTransf, 

					 double C, double r)

  :Element (tag, ELE_TAG_DispBeamColumn2dInt), 

  numSections(numSec), theSections(0), crdTransf(0), C1(C),

  connectedExternalNodes(2),

  Q(6), q(6), rho(r)

{

  // Allocate arrays of pointers to SectionForceDeformations

  theSections = new FiberSection2dInt *[numSections];

    

  if (theSections == 0) {

    opserr << "DispBeamColumn2dInt::DispBeamColumn2dInt - failed to allocate section model pointer\n";

    exit(-1);

  }



  for (int i = 0; i < numSections; i++) {

    

    // Get copies of the material model for each integration point

    SectionForceDeformation *theSection = s[i]->getCopy();

    if (theSections == 0 || theSection->getClassTag() != SEC_TAG_FiberSection2dInt) {

	opserr << "DispBeamColumn2dInt::DispBeamColumn2dInt -- failed to get a copy of section model\n";

	exit(-1);

    } else 

      theSections[i] = (FiberSection2dInt *)theSection;;

  }



  CrdTransf *theCoord = coordTransf.getCopy2d();

  if (theCoord == 0 || theCoord->getClassTag() != CRDTR_TAG_LinearCrdTransf2dInt) {

    opserr << "DispBeamColumn2dInt::DispBeamColumn2dInt -- failed to get a copy of coordinate transformation\n";

    if (theCoord == 0)

      opserr << "COPY ZERO\n";

    else

      opserr << "COPY NON _ZERO CLASTAG " << theCoord->getClassTag() << endln;

    exit(-1);

  } 

  

  crdTransf = (LinearCrdTransf2dInt *)theCoord;

  

  if (crdTransf == 0) {

    opserr << "DispBeamColumn2dInt::DispBeamColumn2dInt - failed to copy coordinate transformation\n";

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



 

// AddingSensitivity:BEGIN /////////////////////////////////////

	parameterID = 0;

// AddingSensitivity:END //////////////////////////////////////

}



DispBeamColumn2dInt::DispBeamColumn2dInt()

:Element (0, ELE_TAG_DispBeamColumn2dInt),

 numSections(0), theSections(0), crdTransf(0), C1(0.0),

 connectedExternalNodes(2),

  Q(6), q(6), rho(0.0)

{

    q0[0] = 0.0;

    q0[1] = 0.0;

    q0[2] = 0.0;

    q0[3] = 0.0;

    q0[4] = 0.0;

    q0[5] = 0.0;



    theNodes[0] = 0;

    theNodes[1] = 0;



// AddingSensitivity:BEGIN /////////////////////////////////////

	parameterID = 0;

// AddingSensitivity:END //////////////////////////////////////

}



DispBeamColumn2dInt::~DispBeamColumn2dInt()

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

}



int

DispBeamColumn2dInt::getNumExternalNodes() const

{

    return 2;

}



const ID&

DispBeamColumn2dInt::getExternalNodes()

{

    return connectedExternalNodes;

}



Node **

DispBeamColumn2dInt::getNodePtrs()

{

    return theNodes;

}



int

DispBeamColumn2dInt::getNumDOF()

{

    return 6;

}



void

DispBeamColumn2dInt::setDomain(Domain *theDomain)

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



	return;

    }



    int dofNd1 = theNodes[0]->getNumberDOF();

    int dofNd2 = theNodes[1]->getNumberDOF();

    

    if (dofNd1 != 3 || dofNd2 != 3) {



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

DispBeamColumn2dInt::commitState()

{

    int retVal = 0;





    // call element commitState to do any base class stuff

    if ((retVal = this->Element::commitState()) != 0) {

      opserr << "DispBeamColumn2dInt::commitState () - failed in base class";

    }    



    // Loop over the integration points and commit the material states

    for (int i = 0; i < numSections; i++)

		retVal += theSections[i]->commitStateB();



    retVal += crdTransf->commitState();					



    return retVal;

}



int

DispBeamColumn2dInt::revertToLastCommit()

{

    int retVal = 0;



	double L = crdTransf->getInitialLength();



    // Loop over the integration points and revert to last committed state

    for (int i = 0; i < numSections; i++)

		retVal += theSections[i]->revertToLastCommitB(L);



    retVal += crdTransf->revertToLastCommit();



    return retVal;

}



int

DispBeamColumn2dInt::revertToStart()

{

    int retVal = 0;

	

    crdTransf->getInitialLength();



    // Loop over the integration points and revert states to start

    for (int i = 0; i < numSections; i++)

		retVal += theSections[i]->revertToStartB();



    retVal += crdTransf->revertToStart();



    return retVal;

}



int

DispBeamColumn2dInt::update(void)

{

  // Update the transformation

  crdTransf->update();

  

  // Get basic deformations

  const Vector &v = crdTransf->getBasicTrialDispInt();		

  double L = crdTransf->getInitialLength();

  double oneOverL = 1.0/L;

  double pts[20];

  quadRule.getSectionLocations(numSections, L, pts);



  // Loop over the integration points

  for (int i = 0; i < numSections; i++) {

    

    int order = theSections[i]->getOrder();				

    const ID &code = theSections[i]->getType();			

														

    Vector e(workArea, order);

    

    double xi = 2.0*pts[i]-1.0;						



    int j;

    for (j = 0; j < order; j++) {

		

      switch(code(j)) {

      case SECTION_RESPONSE_P:

	e(j) = oneOverL*(-v(0)+v(3)); break;   // axial strain

      case SECTION_RESPONSE_MZ:

	e(j) = oneOverL*(-1.0+3.0*(1.0-2.0*C1)*xi)*(v(2)- v(5)); break;	// curvature

      case SECTION_RESPONSE_VY:		       // shear strain

	e(j) = oneOverL*(-v(1)+v(4))-C1*v(2)+(C1-1.0)*v(5); break;

      default:

	break;

      }

    }

    

    // Set the section deformations

    theSections[i]->setTrialSectionDeformationB(e,L);



  }

  

  return 0;

}



const Matrix&

DispBeamColumn2dInt::getTangentStiff()

{

  static Matrix kb(6,6);



  // Zero for integral

  kb.Zero();

  q.Zero();

  

  double L = crdTransf->getInitialLength();

  double oneOverL = 1.0/L;

  

  double pts[20];

  quadRule.getSectionLocations(numSections, L, pts);

  double wts[20];

  quadRule.getSectionWeights(numSections, L, wts);

  

  // Loop over the integration points

  for (int i = 0; i < numSections; i++) {

    

    int order = theSections[i]->getOrder();

    const ID &code = theSections[i]->getType();



    double xi = 2.0*pts[i]-1.0;



    // Get the section tangent stiffness and stress resultant

    const Matrix &ks = theSections[i]->getSectionTangent();			

    const Vector &s = theSections[i]->getStressResultant();			

        

    // Perform numerical integration

    double wti = wts[i]*oneOverL;

 

	double d11, d12, d13, d21, d22, d23, d31, d32, d33;

	

	d11 = ks(0,0);	d12 = ks(0,1);	d13 = ks(0,2);	

	d21 = ks(1,0);	d22 = ks(1,1);	d23 = ks(1,2);	

	d31 = ks(2,0);	d32 = ks(2,1);	d33 = ks(2,2);	



	kb(0,0) += wti*(d11);

	kb(0,1) += wti*(d13);

	kb(0,2) += wti*(d21 + C1*d13*L - 3.0*d21*xi + 6.0*C1*d21*xi);

	kb(0,3) += wti*(-d11);

	kb(0,4) += wti*(-d13);

	kb(0,5) += wti*(-((-1 + C1)*d13*L) + d21*(-1.0 + (3.0 - 6.0*C1)*xi));



	kb(1,0) += wti*(d31);

	kb(1,1) += wti*(d33);

	kb(1,2) += wti*(d32 + C1*d33*L - 3.0*d32*xi + 6.0*C1*d32*xi);

	kb(1,3) += wti*(-d31);

	kb(1,4) += wti*(-d33);

	kb(1,5) += wti*(-((-1.0 + C1)*d33*L) + d32*(-1.0 + (3.0 - 6.0*C1)*xi));



	kb(2,0) +=wti*(d21 + C1*d31*L - 3.0*d21*xi + 6.0*C1*d21*xi);

	kb(2,1) +=wti*(d23 + C1*d33*L - 3.0*d23*xi + 6.0*C1*d23*xi);

	kb(2,2) +=wti*(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi) + C1*L*(d23 + d32 + C1*d33*L - 3.0*d23*xi + 6.0*C1*d23*xi - 3.0*d32*xi + 6.0*C1*d32*xi));

	kb(2,3) +=wti*(-d21 - C1*d31*L + 3.0*d21*xi - 6.0*C1*d21*xi);

	kb(2,4) +=wti*(-d23 - C1*d33*L + 3.0*d23*xi - 6.0*C1*d23*xi);

	kb(2,5) +=wti*(-(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi)) - L*((-1.0 + C1)*d23*(1.0 + (-3.0 + 6.0*C1)*xi) + C1*((-1.0 + C1)*d33*L + d32*(1.0 - 3.0*xi + 6.0*C1*xi))));



	kb(3,0) +=wti*(-d11);

	kb(3,1) +=wti*(-d13);

	kb(3,2) +=wti*(-d21 - C1*d13*L + 3.0*d21*xi - 6.0*C1*d21*xi);

	kb(3,3) +=wti*(d11);

	kb(3,4) +=wti*(d13);

	kb(3,5) +=wti*((-1.0 + C1)*d13*L + d21*(1.0 + (-3.0 + 6.0*C1)*xi));



	kb(4,0) +=wti*(-d31);

	kb(4,1) +=wti*(-d33);

	kb(4,2) +=wti*(-d32 - C1*d33*L + 3.0*d32*xi - 6.0*C1*d32*xi);

	kb(4,3) +=wti*(d31);

	kb(4,4) +=wti*(d33);

	kb(4,5) +=wti*((-1.0 + C1)*d33*L + d32*(1.0 + (-3.0 + 6.0*C1)*xi));



	kb(5,0) +=wti*(-((-1.0 + C1)*d31*L) + d21*(-1.0 + (3.0 - 6.0*C1)*xi));

	kb(5,1) +=wti*(-((-1.0 + C1)*d33*L) + d23*(-1.0 + (3.0 - 6.0*C1)*xi));

	kb(5,2) +=wti*(-(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi)) - L*(d32*(-1.0 + 3.0*xi) + C1*(d23 + d32 - d33*L - 3.0*d23*xi - 9.0*d32*xi) + C1*C1*(d33*L + 6.0*(d23 + d32)*xi)));

	kb(5,3) +=wti*((-1.0 + C1)*d31*L + d21*(1.0 + (-3.0 + 6.0*C1)*xi));

	kb(5,4) +=wti*((-1.0 + C1)*d33*L + d23*(1.0 + (-3.0 + 6.0*C1)*xi));

	kb(5,5) +=wti*(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi) + (-1.0 + C1)*L*((-1.0 + C1)*d33*L + d32*(1.0 - 3.0*xi + 6.0*C1*xi) + d23*(1.0 + (-3.0 + 6.0*C1)*xi)));

   

    double s1, s2, s3;

    double wto = wts[i];



	s1=s(0);	s2=s(1);	s3=s(2);



   q(0)+= wto*(-s1);

   q(1)+= wto*(-s3);

   q(2)+= wto*(-s2 - C1*L*s3 + 3.0*s2*xi - 6.0*C1*s2*xi);

   q(3)+= wto*(s1);

   q(4)+= wto*(s3);

   q(5)+= wto*((-1.0 + C1)*L*s3 + s2*(1.0 + (-3.0 + 6.0*C1)*xi));

    

  }

  

  // Add effects of element loads, q = q(v) + q0		

  q(0) += q0[0];

  q(1) += q0[1];

  q(2) += q0[2];

  q(3) += q0[3];

  q(4) += q0[4];

  q(5) += q0[5];



  // Transform to global stiffness

  K = crdTransf->getGlobalStiffMatrixInt(kb, q);



  return K;

}



const Matrix&

DispBeamColumn2dInt::getInitialBasicStiff()

{

  static Matrix kb(6,6);



  // Zero for integral

  kb.Zero();

  

  double L = crdTransf->getInitialLength();

  double oneOverL = 1.0/L;

  

  double pts[20];

  quadRule.getSectionLocations(numSections, L, pts);

  double wts[20];

  quadRule.getSectionWeights(numSections, L, wts);

  

  // Loop over the integration points

  for (int i = 0; i < numSections; i++) {

    

    int order = theSections[i]->getOrder();

    const ID &code = theSections[i]->getType();

  

    double xi = 2.0*pts[i]-1.0;

    

    // Get the section tangent stiffness and stress resultant

    const Matrix &ks = theSections[i]->getInitialTangent();

    

    // Perform numerical integration

    

    double wti = wts[i]*oneOverL;

  

	double d11, d12, d13, d21, d22, d23, d31, d32, d33;

	

	d11 = ks(0,0);	d12 = ks(0,1);	d13 = ks(0,2);	

	d21 = ks(1,0);	d22 = ks(1,1);	d23 = ks(1,2);	

	d31 = ks(2,0);	d32 = ks(2,1);	d33 = ks(2,2);	



	kb(0,0) += wti*(d11);

	kb(0,1) += wti*(d13);

	kb(0,2) += wti*(d21 + C1*d13*L - 3.0*d21*xi + 6.0*C1*d21*xi);

	kb(0,3) += wti*(-d11);

	kb(0,4) += wti*(-d13);

	kb(0,5) += wti*(-((-1.0 + C1)*d13*L) + d21*(-1.0 + (3.0 - 6.0*C1)*xi));



	kb(1,0) += wti*(d31);

	kb(1,1) += wti*(d33);

	kb(1,2) += wti*(d32 + C1*d33*L - 3.0*d32*xi + 6.0*C1*d32*xi);

	kb(1,3) += wti*(-d31);

	kb(1,4) += wti*(-d33);

	kb(1,5) += wti*(-((-1.0 + C1)*d33*L) + d32*(-1.0 + (3.0 - 6.0*C1)*xi));



	kb(2,0) +=wti*(d21 + C1*d31*L - 3.0*d21*xi + 6.0*C1*d21*xi);

	kb(2,1) +=wti*(d23 + C1*d33*L - 3.0*d23*xi + 6.0*C1*d23*xi);

	kb(2,2) +=wti*(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi) + C1*L*(d23 + d32 + C1*d33*L - 3.0*d23*xi + 6.0*C1*d23*xi - 3.0*d32*xi + 6.0*C1*d32*xi));

	kb(2,3) +=wti*(-d21 - C1*d31*L + 3.0*d21*xi - 6.0*C1*d21*xi);

	kb(2,4) +=wti*(-d23 - C1*d33*L + 3.0*d23*xi - 6.0*C1*d23*xi);

	kb(2,5) +=wti*(-(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi)) - L*((-1.0 + C1)*d23*(1.0 + (-3.0 + 6.0*C1)*xi) + C1*((-1.0 + C1)*d33*L + d32*(1.0 - 3.0*xi + 6.0*C1*xi))));



	kb(3,0) +=wti*(-d11);

	kb(3,1) +=wti*(-d13);

	kb(3,2) +=wti*(-d21 - C1*d13*L + 3.0*d21*xi - 6.0*C1*d21*xi);

	kb(3,3) +=wti*(d11);

	kb(3,4) +=wti*(d13);

	kb(3,5) +=wti*((-1.0 + C1)*d13*L + d21*(1.0 + (-3.0 + 6.0*C1)*xi));



	kb(4,0) +=wti*(-d31);

	kb(4,1) +=wti*(-d33);

	kb(4,2) +=wti*(-d32 - C1*d33*L + 3.0*d32*xi - 6.0*C1*d32*xi);

	kb(4,3) +=wti*(d31);

	kb(4,4) +=wti*(d33);

	kb(4,5) +=wti*((-1.0 + C1)*d33*L + d32*(1.0 + (-3.0 + 6.0*C1)*xi));



	kb(5,0) +=wti*(-((-1.0 + C1)*d31*L) + d21*(-1.0 + (3.0 - 6.0*C1)*xi));

	kb(5,1) +=wti*(-((-1.0 + C1)*d33*L) + d23*(-1.0 + (3.0 - 6.0*C1)*xi));

	kb(5,2) +=wti*(-(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi)) - L*(d32*(-1.0 + 3.0*xi) + C1*(d23 + d32 - d33*L - 3.0*d23*xi - 9.0*d32*xi) + C1*C1*(d33*L + 6.0*(d23 + d32)*xi)));

	kb(5,3) +=wti*((-1.0 + C1)*d31*L + d21*(1.0 + (-3.0 + 6.0*C1)*xi));

	kb(5,4) +=wti*((-1.0 + C1)*d33*L + d23*(1.0 + (-3.0 + 6.0*C1)*xi));

	kb(5,5) +=wti*(d22*(1.0 + (-3.0 + 6.0*C1)*xi)*(1.0 + (-3.0 + 6.0*C1)*xi) + (-1.0 + C1)*L*((-1.0 + C1)*d33*L + d32*(1.0 - 3.0*xi + 6.0*C1*xi) + d23*(1.0 + (-3.0 + 6.0*C1)*xi)));

  }



  return kb;

}



const Matrix&

DispBeamColumn2dInt::getInitialStiff()

{

  const Matrix &kb = this->getInitialBasicStiff();



  // Transform to global stiffness

  K = crdTransf->getInitialGlobalStiffMatrixInt(kb);



  return K;

}



const Matrix&

DispBeamColumn2dInt::getMass()

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

DispBeamColumn2dInt::zeroLoad(void)

{

  Q.Zero();



  q0[0] = 0.0;

  q0[1] = 0.0;

  q0[2] = 0.0;

  q0[3] = 0.0;

  q0[4] = 0.0;

  q0[5] = 0.0;

  

  return;

}



int 

DispBeamColumn2dInt::addLoad(ElementalLoad *theLoad, double loadFactor)

{

  int type;

  const Vector &data = theLoad->getData(type, loadFactor);

  double L = crdTransf->getInitialLength();

  

  if (type == LOAD_TAG_Beam2dUniformLoad) {

    double wt = data(0)*loadFactor;  // Transverse (+ve upward)

    double wa = data(1)*loadFactor;  // Axial (+ve from node I to J)



    // Fixed end forces in basic system

    q0[0] += wa*L*0.5;

    q0[1] += wt*L*0.5;

    q0[2] += wt*L*L/12.0;

    q0[3] += wa*L*0.5;

    q0[4] += wt*L*0.5;

    q0[5] += -wt*L*L/12.0;



  }

  else if (type == LOAD_TAG_Beam2dPointLoad) {

    double V = data(0)*loadFactor;

    double N = data(1)*loadFactor;

    double aOverL = data(2);					



    // Fixed end forces in basic system

   	

	double M1 = L*V*aOverL*(1.0-aOverL)*(1.0-C1-aOverL+aOverL*2.0*C1);

    q0[0] += N*(1.0-aOverL);

    q0[1] += V*(1.0-aOverL);

    q0[2] += M1;

    q0[3] += N*aOverL;

    q0[4] += V*aOverL;

    q0[5] += -M1;



  }

  else {

    opserr << "DispBeamColumn2dInt::DispBeamColumn2dInt -- load type unknown for element with tag: "

	   << this->getTag() << "DispBeamColumn2dInt::addLoad()\n"; 

			    

    return -1;

  }



  return 0;

}



int 

DispBeamColumn2dInt::addInertiaLoadToUnbalance(const Vector &accel)

{

	// Check for a quick return

	if (rho == 0.0) 

		return 0;



	// Get R * accel from the nodes

	const Vector &Raccel1 = theNodes[0]->getRV(accel);

	const Vector &Raccel2 = theNodes[1]->getRV(accel);



    if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {

      opserr << "DispBeamColumn2dInt::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";

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

DispBeamColumn2dInt::getResistingForce()

{

  double L = crdTransf->getInitialLength();



  double pts[20];

  quadRule.getSectionLocations(numSections, L, pts);

  double wts[20];

  quadRule.getSectionWeights(numSections, L, wts);

  

  // Zero for integration

  q.Zero();

  

  // Loop over the integration points

  for (int i = 0; i < numSections; i++) {

    

    int order = theSections[i]->getOrder();

    const ID &code = theSections[i]->getType();

  

    double xi = 2.0*pts[i]-1.0;



    // Get section stress resultant

    const Vector &s = theSections[i]->getStressResultant();

    

    // Perform numerical integration on internal force

       

    double s1, s2, s3;

    double wto = wts[i];

	double L = crdTransf->getInitialLength();

	

	s1=s(0);	s2=s(1);	s3=s(2);



   q(0)+= wto*(-s1);

   q(1)+= wto*(-s3);

   q(2)+= wto*(-s2 - C1*L*s3 + 3.0*s2*xi - 6.0*C1*s2*xi);

   q(3)+= wto*(s1);

   q(4)+= wto*(s3);

   q(5)+= wto*((-1.0 + C1)*L*s3 + s2*(1.0 + (-3.0 + 6.0*C1)*xi));



 

  }

  

  // Add effects of element loads, q = q(v) + q0		

  q(0) += q0[0];

  q(1) += q0[1];

  q(2) += q0[2];

  q(3) += q0[3];

  q(4) += q0[4];

  q(5) += q0[5];





  // Vector for reactions in basic system



  P = crdTransf->getGlobalResistingForceInt(q, q);	

  

  // Subtract other external nodal loads ... P_res = P_int - P_ext

  P(0) -= Q(0);

  P(1) -= Q(1);

  P(2) -= Q(2);

  P(3) -= Q(3);

  P(4) -= Q(4);

  P(5) -= Q(5);

  

  return P;

}



const Vector&

DispBeamColumn2dInt::getResistingForceIncInertia()

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

DispBeamColumn2dInt::sendSelf(int commitTag, Channel &theChannel)

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

  

  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {

    opserr << "DispBeamColumn2dInt::sendSelf() - failed to send ID data\n";

     return -1;

  }    



  // send the coordinate transformation

  

  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {

     opserr << "DispBeamColumn2dInt::sendSelf() - failed to send crdTranf\n";

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

    opserr << "DispBeamColumn2d::sendSelf() - failed to send ID data\n";

    return -1;

  }    



  // send the sections

  

  for (j = 0; j<numSections; j++) {

    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {

      opserr << "DispBeamColumn2dInt::sendSelf() - section " << 

	j << "failed to send itself\n";

      return -1;

    }

  }



  return 0;

}



int

DispBeamColumn2dInt::recvSelf(int commitTag, Channel &theChannel,

			   FEM_ObjectBroker &theBroker)

{

  // receive the integer data containing tag, numSections and coord transformation info

 

	int dbTag = this->getDbTag();

  int i;

  

  static ID idData(7); // one bigger than needed so no clash with section ID



  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {

    opserr << "DispBeamColumn2dInt::recvSelf() - failed to recv ID data\n";

    return -1;

  }    



  this->setTag(idData(0));

  connectedExternalNodes(0) = idData(1);

  connectedExternalNodes(1) = idData(2);

  

  int crdTransfClassTag = idData(4);

  int crdTransfDbTag = idData(5);



  // create a new crdTransf object if one needed

  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {

      if (crdTransf != 0)

	  delete crdTransf;



      crdTransf = new LinearCrdTransf2dInt();



      if (crdTransf == 0) {

	opserr << "DispBeamColumn2d::recvSelf() - failed to obtain a CrdTrans object with classTag " <<

	  crdTransfClassTag << endln;

	  return -2;	  

      }

  }



  crdTransf->setDbTag(crdTransfDbTag);



  // invoke recvSelf on the crdTransf object

  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {

    opserr << "DispBeamColumn2dInt::sendSelf() - failed to recv crdTranf\n";

    return -3;

  }      

  

  //

  // recv an ID for the sections containing each sections dbTag and classTag

  //



  ID idSections(2*idData(3));

  int loc = 0;



  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {

    opserr << "DispBeamColumn2dInt::recvSelf() - failed to recv ID data\n";

    return -1;

  }    



    // now receive the sections

   

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

    theSections = new FiberSection2dInt *[idData(3)];

    if (theSections == 0) {

opserr << "DispBeamColumn2dInt::recvSelf() - out of memory creating sections array of size " <<

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

      theSections[i] = new FiberSection2dInt();

      if (theSections[i] == 0) {

	opserr << "DispBeamColumn2dInt::recvSelf() - Broker could not create Section of class type " <<

	  sectClassTag << endln;

	exit(-1);

      }

      theSections[i]->setDbTag(sectDbTag);

      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {

	opserr << "DispBeamColumn2dInt::recvSelf() - section " << i << " failed to recv itself\n";

	return -1;

      }     

    }



  } else {



    // for each existing section, check it is of correct type

    // (if not delete old & create a new one) then recvSelf on it

      

    loc = 0;

    for (i=0; i<numSections; i++) {

      int sectClassTag = idSections(loc);

      int sectDbTag = idSections(loc+1);

      loc += 2;



      // check of correct type

      if (theSections[i]->getClassTag() !=  sectClassTag) {

	// delete the old section[i] and create a new one

	delete theSections[i];

	theSections[i] = new FiberSection2dInt();;

	if (theSections[i] == 0) {

	opserr << "DispBeamColumn2dInt::recvSelf() - Broker could not create Section of class type " <<

	  sectClassTag << endln;

	exit(-1);

	}

      }



      // recvSelf on it

      theSections[i]->setDbTag(sectDbTag);

      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {

	opserr << "DispBeamColumn2dInt::recvSelf() - section " << i << " failed to recv itself\n";

	return -1;

      }     

    }

  }



  return 0;

}



void

DispBeamColumn2dInt::Print(OPS_Stream &s, int flag)		

{

  s << "\nDispBeamColumn2dInt, element id:  " << this->getTag() << endln;

  s << "\tConnected external nodes:  " << connectedExternalNodes;

  s << "\tCoordTransf: " << crdTransf->getTag() << endln;

  s << "\tmass density:  " << rho << endln;





  s << "\tEnd 1 Forces (P V M): " << -q(0)

    << " " << q(1) << " " << q(2) << endln;

  s << "\tEnd 2 Forces (P V M): " << q(3)

    << " " << -q(4) << " " << q(5) << endln;



}





int

DispBeamColumn2dInt::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)

{

    // first determine the end points of the quad based on

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

DispBeamColumn2dInt::setResponse(const char **argv, int argc, OPS_Stream &s)

{

    // global force - 

    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0

		|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)

		return new ElementResponse(this, 1, P);



    // local force -

    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)

		return new ElementResponse(this, 2, P);



    // chord rotation -

    else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0

	     || strcmp(argv[0],"basicDeformation") == 0)

      return new ElementResponse(this, 3, Vector(3));

    

    // plastic rotation -

    else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0)

      return new ElementResponse(this, 4, Vector(3));

    

    // section response -

    else if (strcmp(argv[0],"section") == 0 || strcmp(argv[0],"-section") == 0) {

      if (argc <= 2)

	return 0;

      

      int sectionNum = atoi(argv[1]);

      if (sectionNum > 0 && sectionNum <= numSections)

	return theSections[sectionNum-1]->setResponse(&argv[2], argc-2, s);

      else

	return 0;

    }

    

    else

      return 0;

}



int 

DispBeamColumn2dInt::getResponse(int responseID, Information &eleInfo)		//L modify

{

  double V;

  double L = crdTransf->getInitialLength();



  if (responseID == 1)

    return eleInfo.setVector(this->getResistingForce());



  else if (responseID == 2) {

      P(3) =  q(3);

      P(0) = -q(0);

      P(2) = q(2);

      P(5) = q(5);

      P(1) =  q(1);

      P(4) = -q(4);

      return eleInfo.setVector(P);

  }



  // Chord rotation

  else if (responseID == 3)

    return eleInfo.setVector(crdTransf->getBasicTrialDispInt());



  // Plastic rotation

  else if (responseID == 4) {

    static Vector vp(6);

    static Vector ve(6);

    const Matrix &kb = this->getInitialBasicStiff();

    kb.Solve(q, ve);

    vp = crdTransf->getBasicTrialDispInt();

    vp -= ve;

    return eleInfo.setVector(vp);

  }



  else

    return -1;

}





// AddingSensitivity:BEGIN ///////////////////////////////////


const Matrix &

DispBeamColumn2dInt::getKiSensitivity(int gradNumber)

{

	K.Zero();

	return K;

}



const Matrix &

DispBeamColumn2dInt::getMassSensitivity(int gradNumber)

{

	K.Zero();

	return K;

}





const Vector &

DispBeamColumn2dInt::getResistingForceSensitivity(int gradNumber)

{

	static Vector dummy(3);		// No distributed loads

	return dummy;

}





// NEW METHOD

int

DispBeamColumn2dInt::commitSensitivity(int gradNumber, int numGrads)

{

 	return 0;

}





// AddingSensitivity:END /////////////////////////////////////////////




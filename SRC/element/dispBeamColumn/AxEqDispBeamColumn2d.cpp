/* *******************************************************************************
Copyright (c) 2016-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS 
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

********************************************************************************* */

// we specify what header files we need
#include "AxEqDispBeamColumn2d.h"
#include <elementAPI.h>
#include <G3Globals.h>

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>


#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <BeamIntegration.h>
#include <CrdTransf.h>
#include <ElementalLoad.h>
#include <CompositeResponse.h>
#include <Parameter.h>


// initialise the class wide variables
Matrix AxEqDispBeamColumn2d::K(6,6);
Vector AxEqDispBeamColumn2d::P(6);
double AxEqDispBeamColumn2d::workArea[100];


static int numMyDBEle = 0;

void *
OPS_AxEqDispBeamColumn2d(void)
{
  // print out a message about who wrote this element & any copyright info wanted
  if (numMyDBEle == 0) {
    opserr << "AxEqDispBeamColumn2d element - Written by Danilo Tarquini 2017 \n";
    numMyDBEle++;
  }
  
  // initializing a pointer to the element class
  Element *theEle = 0;
  
  // added by DANILO: initialization of the beam integration
  BeamIntegration *beamIntegr = 0;
  
  // creation of an empty element, required for parallel processing only
  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new AxEqDispBeamColumn2d();
    return theEle;
  }
  
  // check the required input parameters are given
  // for this element it is not possible to define different integration sections along the element, to specify element masses, and to form consistent mass matrix
  // Differently from a classical DB element, 
  // 1) a tolerance on the axial force unbalance $tol in the sections MUST be specified; It represents the unbalance in axial force at each IP that is deemed acceptable. It depends on the analysis performed and on the employed units
  // 2) a max num of internal element iterations can be specified. Default value is : $maxIters=20
  if (numRemainingArgs < 7) {
    opserr << "insufficient arguments: 1)eleTag? 2)iNode? 3)jNode? 4)numIntgrPts? 5)-$secTag? 6)$transfTag? 7)$tol optionals: <-integration $intType> <-iter $maxIters>";
    numMyDBEle++;
  }

  // get the input parameters 
  int iData[6];
  int numData = 6;
  
  // check the that 6 integer are given as input data
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }

  // variables storing the main quantities required to create the element
  int eleTag = iData[0];           // element tag
  int iNode = iData[1];            // initial node
  int jNode = iData[2];            // final node
  int nIP = iData[3];              // number of integration points
  int secTag = iData[4];           // section tag
  int transfTag = iData[5];        // transformation tag
  
  // added by DANILO: reading of the specified axial force tolerance
  double tolerance;
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &tolerance) != 0) {
    opserr << "WARNING error reading tolerance" << eleTag << endln;
    return 0;
  }
  
  // optional quantities that can be defined for the element
  double mass = 0.0;
  int cMass = 0;
  int numdata = 1;
  int maxNumIters = 20;
  
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* type = OPS_GetString();
    if (strcmp(type, "-cMass") == 0) {
      opserr << "WARNING: Consistent mass matrix not available for this element, Lumped mass matrix is used \n";
      cMass = 0; // changed by Danilo with respect to DB element
    }
    else if (strcmp(type, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetDoubleInput(&numData, &mass) < 0) {
	  opserr << "WARNING: invalid mass\n";
	  return 0;
	}
	opserr << "WARNING: Element mass cannot be defined for this element\n";
	return 0;
      }
    }
    // added by DANILO: optional parameter that can be specified is the integration type (LOBATTO or LEGENDRE) 
    else if (strcmp(type, "-integration") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
	const char* typeIntegration = OPS_GetString();
	if (strcmp(typeIntegration, "Lobatto") == 0)
	  beamIntegr = new LobattoBeamIntegration();
	else if (strcmp(typeIntegration, "Legendre") == 0)
	  beamIntegr = new LegendreBeamIntegration();
	else {
	  opserr << "WARNING: invalid integration type\n";
	  return 0;
	}
      }
    }
    // added by DANILO: optional parameter that can be specified is the number of internal element iterations 
    else if (strcmp(type, "-iter") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numData, &maxNumIters) < 0 && maxNumIters<=0) {
	  opserr << "WARNING: invalid max number of iterations\n";
	  return 0;
	}
      }
    }
    
  }
  
  // added by DANILO: if Beam Integration is not defined, use gauss legendere
  if (beamIntegr == 0) {
    beamIntegr = new LegendreBeamIntegration();
  }
  
  // get the beam transformation
  CrdTransf* theTransf = OPS_GetCrdTransf(transfTag);  // get the transformation
  if (theTransf == 0) {                                 // return error if the transformation is not found
    opserr << "coord transfomration not found\n";
    return 0;
  }

  // pointer to pointers of sections of type secTag
  SectionForceDeformation *theSection = OPS_getSectionForceDeformation(secTag);   // get the section of type sectionTag
  SectionForceDeformation** sections = new SectionForceDeformation *[nIP];        // pointer to vector of pointers to allocate nIP sections
  // check that the section tag exists
  if (theSection == 0) {
    opserr << "WARNING section not found\n";
    opserr << "Section: " << secTag;
    opserr << " element: " << eleTag << endln;
    return 0;
  }
  //definition of the vector of section pointers
  for (int i = 0; i < nIP; i++)
    sections[i] = theSection;
  
  // now creation of the element
  theEle = new AxEqDispBeamColumn2d(eleTag, iNode, jNode, nIP, sections, *beamIntegr, *theTransf, tolerance, mass, cMass, maxNumIters);

  delete[] sections;
  delete beamIntegr;

  return theEle;
}


// DEFAULT CONSTRUCTOR WITH GIVEN IPNUT PARAMETERS 
AxEqDispBeamColumn2d::AxEqDispBeamColumn2d(int tag, int nd1, int nd2,
	int numSec, SectionForceDeformation **s,
	BeamIntegration& bi,
	CrdTransf &coordTransf, double tolerance, double r, int cm, int maxNumIters)
	:Element(tag, ELE_TAG_AxEqDispBeamColumn2d), //0 is the class tag since it is an external element
	numSections(numSec), theSections(0), crdTransf(0), beamInt(0), // initialization of some element members
	connectedExternalNodes(2),
	Q(6), q(3), rho(r), cMass(cm), parameterID(0), vCommitted(3), eCommitted(2), e0Committed(numSec), curvCommitted(numSec), tol(tolerance), maxIters(maxNumIters)
{
	// Allocate arrays of pointers to SectionForceDeformations
	theSections = new SectionForceDeformation *[numSections];

	if (theSections == 0) {
		opserr << "AxEqDispBeamColumn2d::AxEqDispBeamColumn2d - failed to allocate section model pointer\n";
		exit(-1);
	}


	for (int i = 0; i < numSections; i++) {

		// Get copies of the material model for each integration point
		theSections[i] = s[i]->getCopy();

		// Check allocation
		if (theSections[i] == 0) {
			opserr << "AxEqDispBeamColumn2d::AxEqDispBeamColumn2d -- failed to get a copy of section model\n";
			exit(-1);
		}
	}

	beamInt = bi.getCopy();

	if (beamInt == 0) {
		opserr << "AxEqDispBeamColumn2d::AxEqDispBeamColumn2d - failed to copy beam integration\n";
		exit(-1);
	}

	crdTransf = coordTransf.getCopy2d();

	if (crdTransf == 0) {
		opserr << "AxEqDispBeamColumn2d::AxEqDispBeamColumn2d - failed to copy coordinate transformation\n";
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
}

// CONSTRUCTOR WITH NO INPUT PARAMETERS
AxEqDispBeamColumn2d::AxEqDispBeamColumn2d()
	:Element(0, ELE_TAG_AxEqDispBeamColumn2d),
	numSections(0), theSections(0), crdTransf(0), beamInt(0),
	connectedExternalNodes(2),
	Q(6), q(3), rho(0.0), cMass(0), parameterID(0), vCommitted(3), eCommitted(2), e0Committed(0), curvCommitted(0)
{
	q0[0] = 0.0;
	q0[1] = 0.0;
	q0[2] = 0.0;

	p0[0] = 0.0;
	p0[1] = 0.0;
	p0[2] = 0.0;

	theNodes[0] = 0;
	theNodes[1] = 0;

}

// DESTRUCTOR
AxEqDispBeamColumn2d::~AxEqDispBeamColumn2d()
{
	for (int i = 0; i < numSections; i++) {
		if (theSections[i])
			delete theSections[i];
	}

	// Delete the array of pointers to SectionForceDeformation pointer arrays
	if (theSections)
		delete[] theSections;

	if (crdTransf)
		delete crdTransf;

	if (beamInt != 0)
		delete beamInt;
}

// method returning the number of nodes to which the element is connected
int
AxEqDispBeamColumn2d::getNumExternalNodes() const
{
	return 2;
}

// method returning the tags of the nodes to which the element is connected
const ID&
AxEqDispBeamColumn2d::getExternalNodes()
{
	return connectedExternalNodes;
}

// method returning the pointers to the nodes to which the element is connected
Node **
AxEqDispBeamColumn2d::getNodePtrs()
{
	return theNodes;
}

// method returning the number of global DOFs belonging to the element
int
AxEqDispBeamColumn2d::getNumDOF()
{
	return 6;
}

// method adding the element to the DOMAIN, setting the right nodes to the element, initialising the transformation object
void
AxEqDispBeamColumn2d::setDomain(Domain *theDomain)
{
  opserr << "setDomain - START\n";
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
		opserr << "WARNING AxEqDispBeamColumn2d (tag: %d), node not found in domain" << this->getTag() << endln;;
		return;
	}

	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();

	if (dofNd1 != 3 || dofNd2 != 3) {
		//opserr << "FATAL ERROR AxEqDispBeamColumn2d (tag: %d), has differing number of DOFs at its nodes",
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
  opserr << "setDomain - 9\n";
	this->update();
  opserr << "setDomain - DONE\n";
}

// method committing the state of the object (called when we are at a convergence point)
int
AxEqDispBeamColumn2d::commitState()
{
	int retVal = 0;

	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "AxEqDispBeamColumn2d::commitState () - failed in base class";
	}

	// Loop over the integration points and commit the material states
	for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

	retVal += crdTransf->commitState();
	// added by Danilo: variables to store committed displacements and sectional deformations
	vCommitted = crdTransf->getBasicTrialDisp();

	// Loop over the integration points
	for (int i = 0; i < numSections; i++) {
		eCommitted = theSections[i]->getSectionDeformation();
		e0Committed(i) = eCommitted(0);
		curvCommitted(i) = eCommitted(1);
		//opserr << "eCommitted=" << eCommitted << endln;
	}
	
	//opserr << "e0Committed=" << e0Committed << endln;
	//opserr << "curvCommitted=" << curvCommitted << endln;
	return retVal;
}

// METHOD added by DANILO: calculation of the unbalance in the sectional axial force
double
AxEqDispBeamColumn2d::getSectionalAxialForceUnbalance() {
	
	double resN = 0;
	// loo over all integration sections
	for (int i = 0; i < numSections-1; i++) {
		const Vector &sect1 = theSections[i]->getStressResultant();
		double Ns1 = sect1[0];
		const Vector &sect2 = theSections[i+1]->getStressResultant();
		double Ns2 = sect2[0];
		resN += abs(Ns2 - Ns1);
	}
	return resN;
}

// method reverting to the last committed state of the object (called when we are not at a convergence point)
int
AxEqDispBeamColumn2d::revertToLastCommit()
{
	int retVal = 0;

	// Loop over the integration points and revert to last committed state
	for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

	retVal += crdTransf->revertToLastCommit();
	
	return retVal;
}

// method reverting to the initial state of the object (all set to 0)
int
AxEqDispBeamColumn2d::revertToStart()
{
	int retVal = 0;

	// Loop over the integration points and revert states to start
	for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();

	retVal += crdTransf->revertToStart();
	
	return retVal;
}


// modified by DANILO: method updating the state of the object 
int
AxEqDispBeamColumn2d::update(void)
{
	// CASE IN WHICH THE SECTIONAL AXIAL FORCES ARE IN EQUILIBRIUM: SAME PROCEDURE AS FOR THE CLASSICAL DB ELEMENT 
	// However some modifications were made by DT in order to have obtain same procedure described in Tarquini et al. 2017
	// namely the use of displacement increments as starting point computed as total displacements (trial) - committed (at convergence) displacements
	//opserr << "I AM HERE" << endln;
		flagDBae = 0;
		//opserr << "I am here NORMAL procedure" << endln;
		int err = 0;

		// Update the transformation
		crdTransf->update();

		// Get basic deformations
		const Vector &v = crdTransf->getBasicTrialDisp();
		Vector dv = v - vCommitted;
		//opserr << "ubsc_trial=" << v << "ubsc_comm="<< vCommitted << endln;

		double L = crdTransf->getInitialLength();
		double oneOverL = 1.0 / L;

		//const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
		double xi[maxNumSections];
		beamInt->getSectionLocations(numSections, L, xi);


		// Loop over the integration points
		for (int i = 0; i < numSections; i++) {

			int order = theSections[i]->getOrder();
			const ID &code = theSections[i]->getType();

			//Vector eCommitted(order); // eCommitted(workArea,order); I do not understand what workArea is needed for
			Vector de(order);
			Vector e(order);

			//double xi6 = 6.0*pts(i,0);
			double xi6 = 6.0*xi[i];
			
			int j;
			for (j = 0; j < order; j++) {
				switch (code(j)) {
				case SECTION_RESPONSE_P:
					//eCommitted(j) = oneOverL*vCommitted(0);
					//de(j) = oneOverL*dv(0); break;
					de(j) = oneOverL*dv(0); 
					e(j) = e0Committed(i) + de(j); break;
				case SECTION_RESPONSE_MZ:
					//eCommitted(j) = oneOverL*((xi6 - 4.0)*vCommitted(1) + (xi6 - 2.0)*vCommitted(2));
					//de(j) = oneOverL*((xi6 - 4.0)*dv(1) + (xi6 - 2.0)*dv(2)); break;
					de(j) = oneOverL*((xi6 - 4.0)*dv(1) + (xi6 - 2.0)*dv(2)); 
					e(j) = curvCommitted(i) + de(j); break;
				default:
					//eCommitted(j) = 0.0;
					//de(j) = 0.0; break;
					de(j) = 0.0; break;
					e(j) = 0.0; break;
				}
			}


			//e = eCommitted + de;
			// Set the section deformations
			err += theSections[i]->setTrialSectionDeformation(e);
			//opserr << "sectE" << e << endln;
		}

		// Check the unbalance in the axial forces
		if (getSectionalAxialForceUnbalance() < tol) {

			if (err != 0) {
				opserr << "AxEqDispBeamColumn2d::update() - failed setTrialSectionDeformations()\n";
				return err;
			}
		return 0;
		}
		// CASE IN WHICH THE SECTIONAL AXIAL FORCES ARE NOT IN EQUILIBRIUM: PROCEDURE DESCRIBED IN Tarquini et al. 2017 
		else {
			flagDBae = 1;
		//opserr << "I am here DB procedure" << endln;
			int err = 0;
			bool converged = false;
			int nIntIter = 0;
			while (converged == false && nIntIter < maxIters) {// In case of non equilibrium:
				// 1) get the increment of axial strains
				Vector deAE(numSections);
				deAE = getAxialStrainIncrement();
				// 2) get current sectional strains and add the axial equilibrium term to the axial strain			
				for (int i = 0; i < numSections; i++) {   // Loop over the integration points
					int order = theSections[i]->getOrder();
					Vector eCurrent(order);
					eCurrent = theSections[i]->getSectionDeformation();
					eCurrent(0) = eCurrent(0) + deAE(i);
					err += theSections[i]->setTrialSectionDeformation(eCurrent); //set new deformation as new trial and compute sectional responses
					//opserr << "sectdef" << eCurrent << endln;
				}

				Vector sectionalN(numSections);
				for (int i = 0; i < numSections; i++) {
					Vector sectP = theSections[i]->getStressResultant();
					sectionalN(i) = sectP(0);
				}
				//opserr << "sectN" << sectionalN << endln;

				nIntIter = nIntIter + 1;				
				if (getSectionalAxialForceUnbalance() < tol) {
					Vector sectionalN(numSections);
					for (int i = 0; i < numSections; i++) {
						Vector sectP = theSections[i]->getStressResultant();
						sectionalN(i) = sectP(0);
					}				
					//opserr << "sectNconv" << sectionalN << endln;
					//opserr << "CONVERGENCE DB_AE \n";
					//opserr << "unbalance" << getSectionalAxialForceUnbalance() << endln;
					//opserr << "nIntIter=" << nIntIter << endln;
					converged = true;
				}
				} // end while

			if (converged == false) {
				opserr << "WARNING - DB_AE NOT CONVERGING \n";
				opserr << "Disp=" << crdTransf->getBasicTrialDisp() << endln;
				//opserr << "unbalance" << getSectionalAxialForceUnbalance() << endln;



				return -1;
			}
			return 0;
			} 			
} // end method

// METHOD added by DANILO: calculation of the unbalance in the sectional axial force
Vector
AxEqDispBeamColumn2d::getAxialStrainIncrement(void) {
	
	// 1) definition of the element length
	double L = crdTransf->getInitialLength();

	// 2) definition of the vector for the integration weights
	double wi[maxNumSections];
	beamInt->getSectionWeights(numSections, L, wi);  // be careful that wi sums up to 1!!!			
	for (int i = 0; i < numSections; i++)  // set sum of integration weigths to 2
		wi[i] *= 2;
	

	// 3) get axial basic displacement
	const Vector &v = crdTransf->getBasicTrialDisp();
	double u1Bsc = v[0];
	//opserr << "ubsc" << v  << endln;
	

	// 4) get sectional axial forces
	Vector sectionalN(numSections);
	for (int i = 0; i < numSections; i++) {
		Vector sectP = theSections[i]->getStressResultant();
		sectionalN(i) = sectP(0);
	}
	//opserr << "sectN" << sectionalN << endln;

	// 5) get sectional axial stiffness component
	Vector sectionalK11(numSections);
	for (int i = 0; i < numSections; i++) {
		const Matrix &sectK = theSections[i]->getSectionTangent();
		sectionalK11(i) = sectK(0, 0);
	}
	//opserr << "sectionalK11" << sectionalK11 << endln;

	// 6) get current sectional axial strain
	Vector e0Current(numSections);
	for (int i = 0; i < numSections; i++) {
		int order = theSections[i]->getOrder();
		Vector e(order);
		e = theSections[i]->getSectionDeformation();
		e0Current(i) = e(0);
	}
	//opserr << "secte" << e0Current << endln;
	
	// 7) calculation of the increment of axial strains
	Vector deAxialEq(numSections);
	double sumNumerator = 0;
	double sumDenominator = 0;
	// 7.1) calculation of the summation oterms in the numerator and denominator
	for (int i = 0; i < numSections; i++) {
		sumNumerator = sumNumerator + wi[i] * (e0Current(i)+(sectionalN(0) - sectionalN(i)) / sectionalK11(i));
		sumDenominator = sumDenominator + wi[i] * (sectionalK11(0) / sectionalK11(i));
	}
	// 7.2) axial equilibrium strain increments for the first integration point
	deAxialEq(0) = (2 * u1Bsc / L - sumNumerator) / sumDenominator;
	
	// 7.3) axial equilibrium strain increments for all other integration points
	for (int i = 0; i < numSections; i++) {
		deAxialEq(i) = ((sectionalN(0) - sectionalN(i)) + sectionalK11(0)*deAxialEq(0)) / sectionalK11(i);
	}
	//opserr << "sectdeAX" << deAxialEq << endln;


	return deAxialEq;
}




// modified by DANILO: method computing the tangent stiffness matrix in the basic reference system
void
AxEqDispBeamColumn2d::getBasicStiff(Matrix &kb, int initial)
{

	// CASE IN WHICH THE SECTIONAL AXIAL FORCES ARE IN EQUILIBRIUM: SAME PROCEDURE AS FOR THE CLASSICAL DB ELEMENT 
	if (flagDBae==0) {
		//opserr << "I am here NORMAL Stiff" << endln;
	// Zero for integral
	kb.Zero();

	double L = crdTransf->getInitialLength();
	double oneOverL = 1.0 / L;

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

		double xi6 = 6.0*xi[i];

		// Get the section tangent stiffness and stress resultant
		const Matrix &ks = theSections[i]->getSectionTangent();

		// Perform numerical integration
		//kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
		//double wti = wts(i)*oneOverL;
		double wti = wt[i] * oneOverL;
		double tmp;
		int j, k;
		for (j = 0; j < order; j++) {
			switch (code(j)) {
			case SECTION_RESPONSE_P:
				for (k = 0; k < order; k++)
					ka(k, 0) += ks(k, j)*wti;
				break;
			case SECTION_RESPONSE_MZ:
				for (k = 0; k < order; k++) {
					tmp = ks(k, j)*wti;
					ka(k, 1) += (xi6 - 4.0)*tmp;
					ka(k, 2) += (xi6 - 2.0)*tmp;
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
					kb(0, k) += ka(j, k);
				break;
			case SECTION_RESPONSE_MZ:
				for (k = 0; k < 3; k++) {
					tmp = ka(j, k);
					kb(1, k) += (xi6 - 4.0)*tmp;
					kb(2, k) += (xi6 - 2.0)*tmp;
				}
				break;
			default:
				break;
			}
		}
	}
	}
	// CASE IN WHICH THE SECTIONAL AXIAL FORCES ARE NOT IN EQUILIBRIUM: PROCEDURE DESCRIBED IN Tarquini et al. 2017 

	else {
		//opserr << "I am here DBaeStiff" << endln;
		
		// initialization required variables
		double de01dtheta1_num = 0;  // numerator derivative axial strain first IP with respect to element basic displacement (rotation1)
		double de01dtheta2_num = 0; // numerator derivative axial strain first IP with respect to element basic displacement (rotation2)
		double de01du_den = 0;
		Vector de0du1(numSections);  // vector containing derivatives of axial strain with respect to element basic displacement (u1)
		Vector de0dtheta1(numSections); // vector containing derivatives of axial strain with respect to element basic displacement (rotation1)
		Vector de0dtheta2(numSections); // vector containing derivatives of axial strain with respect to element basic displacement (rotation2)


		// 1) definition of the element length and section weights and locations
		double L = crdTransf->getInitialLength();
		double oneOverL = 1.0 / L;
		double wi[maxNumSections];
		double xi[maxNumSections];	
		double xi6[maxNumSections];
		beamInt->getSectionWeights(numSections, L, wi);  // be careful that wi sums up to 1!!!
		beamInt->getSectionLocations(numSections, L, xi);
		for (int i = 0; i < numSections; i++) {  // set sum of integration weigths to 2
			wi[i] *= 2;
			xi6[i] = 6.0*xi[i];
		}

		

		// 2) get sectional axial stiffness component
		Vector sectionalK11(numSections);
		Vector sectionalK12(numSections);
		for (int i = 0; i < numSections; i++) {
			const Matrix &sectK = theSections[i]->getSectionTangent();
			sectionalK11(i) = sectK(0, 0);
			sectionalK12(i) = sectK(0, 1);
		}


		// calculation derivative of axial strain FIRST IP with respect to the basic displacement vector(U1, theta1, theta2)
		for (int i = 0; i < numSections; i++) {
			de01dtheta1_num = de01dtheta1_num - wi[i] * ((sectionalK12(0)*(oneOverL*(xi6[0] - 4.0))- sectionalK12(i)*(oneOverL*(xi6[i] - 4.0)))/ sectionalK11(i));
			de01dtheta2_num = de01dtheta2_num - wi[i] * ((sectionalK12(0)*(oneOverL*(xi6[0] - 2.0)) - sectionalK12(i)*(oneOverL*(xi6[i] - 2.0))) / sectionalK11(i));
			de01du_den = de01du_den + wi[i] * (sectionalK11(0) / sectionalK11(i));
		}

		de0du1(0) = 2 / (L * de01du_den);
		de0dtheta1(0) = de01dtheta1_num / de01du_den;
		de0dtheta2(0) = de01dtheta2_num / de01du_den;

		// calculation derivative of axial strain AT ALL IPs with respect to the basic displacement vector(U1, theta1, theta2)
		for (int i = 0; i < numSections; i++) {
			de0du1(i) = sectionalK11(0)*de0du1(0) / sectionalK11(i);
			de0dtheta1(i) = (sectionalK12(0)*(oneOverL*(xi6[0] - 4.0)) - sectionalK12(i)*(oneOverL*(xi6[i] - 4.0)) + sectionalK11(0)*de0dtheta1(0)) / sectionalK11(i);
			de0dtheta2(i) = (sectionalK12(0)*(oneOverL*(xi6[0] - 2.0)) - sectionalK12(i)*(oneOverL*(xi6[i] - 2.0)) + sectionalK11(0)*de0dtheta2(0)) / sectionalK11(i);
		}
		
		// compute the element end forces p_bsc and the element tangent stiffness matrix
		kb.Zero();


		for (int i = 0; i < numSections; i++) {
			const Matrix &sectK = theSections[i]->getSectionTangent();

			Matrix Bx(2, 3);
			Matrix BBx(2, 3);
			
			Bx.Zero();
			Bx(0, 0) = 1/L;
			Bx(0, 1) = 0;
			Bx(0, 2) = 0;
			Bx(1, 0) = 0;
			Bx(1, 1) = oneOverL*(xi6[i] - 4.0);		
			Bx(1, 2) = oneOverL*(xi6[i] - 2.0);
			
			
			BBx.Zero();
			BBx(0,0) = de0du1(i);
			BBx(0, 1) = de0dtheta1(i);
			BBx(0,2)= de0dtheta2(i);
			BBx(1, 0) = Bx(1, 0);
			BBx(1, 1) = Bx(1, 1);
			BBx(1, 2) = Bx(1, 2);

			kb.addMatrixTripleProduct(1.0, Bx, sectK, BBx, wi[i]*L/2);
		}
	} // end else
}// end method

// method transforming the basic stiffness matrix to global stiffness matrix
const Matrix&
AxEqDispBeamColumn2d::getTangentStiff()
{
	static Matrix kb(3, 3);

	//opserr << "I am here 2" << endln;

	this->getBasicStiff(kb);

	//opserr << "I am here 3" << endln;


	// Zero for integral
	q.Zero();

	double L = crdTransf->getInitialLength();
	double oneOverL = 1.0 / L;

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

		// Get the section tangent stiffness and stress resultant
		const Vector &s = theSections[i]->getStressResultant();

		// Perform numerical integration
		//q.addMatrixTransposeVector(1.0, *B, s, wts(i));
		double si;
		for (int j = 0; j < order; j++) {
			//si = s(j)*wts(i);
			si = s(j)*wt[i];
			switch (code(j)) {
			case SECTION_RESPONSE_P:
				q(0) += si; break;
			case SECTION_RESPONSE_MZ:
				q(1) += (xi6 - 4.0)*si; q(2) += (xi6 - 2.0)*si; break;
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

// method computing the initial stiffness matrix in the basic reference system
const Matrix&
AxEqDispBeamColumn2d::getInitialBasicStiff()
{
	static Matrix kb(3, 3);

	// Zero for integral
	kb.Zero();

	double L = crdTransf->getInitialLength();
	double oneOverL = 1.0 / L;

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
		double wti = wt[i] * oneOverL;
		double tmp;
		int j, k;
		for (j = 0; j < order; j++) {
			switch (code(j)) {
			case SECTION_RESPONSE_P:
				for (k = 0; k < order; k++)
					ka(k, 0) += ks(k, j)*wti;
				break;
			case SECTION_RESPONSE_MZ:
				for (k = 0; k < order; k++) {
					tmp = ks(k, j)*wti;
					ka(k, 1) += (xi6 - 4.0)*tmp;
					ka(k, 2) += (xi6 - 2.0)*tmp;
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
					kb(0, k) += ka(j, k);
				break;
			case SECTION_RESPONSE_MZ:
				for (k = 0; k < 3; k++) {
					tmp = ka(j, k);
					kb(1, k) += (xi6 - 4.0)*tmp;
					kb(2, k) += (xi6 - 2.0)*tmp;
				}
				break;
			default:
				break;
			}
		}

	}

	return kb;
}

// method returning the initial stiffness matrix in the basic reference system
const Matrix&
AxEqDispBeamColumn2d::getInitialStiff()
{
	const Matrix &kb = this->getInitialBasicStiff();

	// Transform to global stiffness
	K = crdTransf->getInitialGlobalStiffMatrix(kb);

	return K;
}

// method returning the element mass matrix (lumped or consistent depending on parameter cMass)
const Matrix&
AxEqDispBeamColumn2d::getMass()
{
	K.Zero();

	if (rho == 0.0)
		return K;

	double L = crdTransf->getInitialLength();
	if (cMass == 0) {
		// lumped mass matrix
		double m = 0.5*rho*L;
		K(0, 0) = K(1, 1) = K(3, 3) = K(4, 4) = m;
	}
	else {
		// consistent mass matrix
		static Matrix ml(6, 6);
		double m = rho*L / 420.0;
		ml(0, 0) = ml(3, 3) = m*140.0;
		ml(0, 3) = ml(3, 0) = m*70.0;

		ml(1, 1) = ml(4, 4) = m*156.0;
		ml(1, 4) = ml(4, 1) = m*54.0;
		ml(2, 2) = ml(5, 5) = m*4.0*L*L;
		ml(2, 5) = ml(5, 2) = -m*3.0*L*L;
		ml(1, 2) = ml(2, 1) = m*22.0*L;
		ml(4, 5) = ml(5, 4) = -ml(1, 2);
		ml(1, 5) = ml(5, 1) = -m*13.0*L;
		ml(2, 4) = ml(4, 2) = -ml(1, 5);

		// transform local mass matrix to global system
		K = crdTransf->getGlobalMatrixFromLocal(ml);
	}

	return K;
}

// method setting to zero the element load
void
AxEqDispBeamColumn2d::zeroLoad(void)
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


// method adding element loads (uniform or point) if defined
int
AxEqDispBeamColumn2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	int type;
	const Vector &data = theLoad->getData(type, loadFactor);
	double L = crdTransf->getInitialLength();

	if (type == LOAD_TAG_Beam2dUniformLoad) {
		double wt = data(0)*loadFactor;  // Transverse (+ve upward)
		double wa = data(1)*loadFactor;  // Axial (+ve from node I to J)

		double V = 0.5*wt*L;
		double M = V*L / 6.0; // wt*L*L/12
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
		double b = L - a;

		// Reactions in basic system
		p0[0] -= N;
		double V1 = P*(1.0 - aOverL);
		double V2 = P*aOverL;
		p0[1] -= V1;
		p0[2] -= V2;

		double L2 = 1.0 / (L*L);
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
		opserr << "AxEqDispBeamColumn2d::AxEqDispBeamColumn2d -- load type unknown for element with tag: "
			<< this->getTag() << "AxEqDispBeamColumn2d::addLoad()\n";

		return -1;
	}

	return 0;
}

// method adding inertia loads (for dynamic analysis) if element mass is defined
int
AxEqDispBeamColumn2d::addInertiaLoadToUnbalance(const Vector &accel)
{
	// Check for a quick return
	if (rho == 0.0)
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = theNodes[0]->getRV(accel);
	const Vector &Raccel2 = theNodes[1]->getRV(accel);

	if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
		opserr << "AxEqDispBeamColumn2d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
		return -1;
	}

	// want to add ( - fact * M R * accel ) to unbalance
	if (cMass == 0) {
		// take advantage of lumped mass matrix
		double L = crdTransf->getInitialLength();
		double m = 0.5*rho*L;

		Q(0) -= m*Raccel1(0);
		Q(1) -= m*Raccel1(1);
		Q(3) -= m*Raccel2(0);
		Q(4) -= m*Raccel2(1);
	}
	else {
		// use matrix vector multip. for consistent mass matrix
		static Vector Raccel(6);
		for (int i = 0; i<3; i++) {
			Raccel(i) = Raccel1(i);
			Raccel(i + 3) = Raccel2(i);
		}
		Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
	}

	return 0;
}

// method calculating and returning the resisting forces in the global reference system
const Vector&
AxEqDispBeamColumn2d::getResistingForce()
{
	double L = crdTransf->getInitialLength();

	double oneOverL = 1.0 / L;

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
			switch (code(j)) {
			case SECTION_RESPONSE_P:
				q(0) += si; break;
			case SECTION_RESPONSE_MZ:
				q(1) += (xi6 - 4.0)*si; q(2) += (xi6 - 2.0)*si; break;
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

// method calculating and returning the resisting forces including the inertia term (dynamic analysis if  element mass is defined)
const Vector&
AxEqDispBeamColumn2d::getResistingForceIncInertia()
{
	P = this->getResistingForce();

	// Subtract other external nodal loads ... P_res = P_int - P_ext
	P.addVector(1.0, Q, -1.0);

	if (rho != 0.0) {
		const Vector &accel1 = theNodes[0]->getTrialAccel();
		const Vector &accel2 = theNodes[1]->getTrialAccel();

		if (cMass == 0) {
			// take advantage of lumped mass matrix
			double L = crdTransf->getInitialLength();
			double m = 0.5*rho*L;

			P(0) += m*accel1(0);
			P(1) += m*accel1(1);
			P(3) += m*accel2(0);
			P(4) += m*accel2(1);
		}
		else {
			// use matrix vector multip. for consistent mass matrix
			static Vector accel(6);
			for (int i = 0; i<3; i++) {
				accel(i) = accel1(i);
				accel(i + 3) = accel2(i);
			}
			P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
		}

		// add the damping forces if rayleigh damping
		if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			P.addVector(1.0, this->getRayleighDampingForces(), 1.0);

	}
	else {

		// add the damping forces if rayleigh damping
		if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
	}

	return P;
}

// method sending element information (for output)
int
AxEqDispBeamColumn2d::sendSelf(int commitTag, Channel &theChannel)
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
	int crdTransfDbTag = crdTransf->getDbTag();
	if (crdTransfDbTag == 0) {
		crdTransfDbTag = theChannel.getDbTag();
		if (crdTransfDbTag != 0)
			crdTransf->setDbTag(crdTransfDbTag);
	}
	data(5) = crdTransfDbTag;
	data(6) = beamInt->getClassTag();
	int beamIntDbTag = beamInt->getDbTag();
	if (beamIntDbTag == 0) {
		beamIntDbTag = theChannel.getDbTag();
		if (beamIntDbTag != 0)
			beamInt->setDbTag(beamIntDbTag);
	}
	data(7) = beamIntDbTag;
	data(8) = rho;
	data(9) = cMass;
	data(10) = alphaM;
	data(11) = betaK;
	data(12) = betaK0;
	data(13) = betaKc;

	if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
		opserr << "AxEqDispBeamColumn2d::sendSelf() - failed to send data Vector\n";
		return -1;
	}

	// send the coordinate transformation
	if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
		opserr << "AxEqDispBeamColumn2d::sendSelf() - failed to send crdTranf\n";
		return -1;
	}

	// send the beam integration
	if (beamInt->sendSelf(commitTag, theChannel) < 0) {
		opserr << "AxEqDispBeamColumn2d::sendSelf() - failed to send beamInt\n";
		return -1;
	}

	//
	// send an ID for the sections containing each sections dbTag and classTag
	// if section ha no dbTag get one and assign it
	//
	ID idSections(2 * numSections);
	loc = 0;
	for (i = 0; i<numSections; i++) {
		int sectClassTag = theSections[i]->getClassTag();
		int sectDbTag = theSections[i]->getDbTag();
		if (sectDbTag == 0) {
			sectDbTag = theChannel.getDbTag();
			theSections[i]->setDbTag(sectDbTag);
		}

		idSections(loc) = sectClassTag;
		idSections(loc + 1) = sectDbTag;
		loc += 2;
	}

	if (theChannel.sendID(dbTag, commitTag, idSections) < 0) {
		opserr << "AxEqDispBeamColumn2d::sendSelf() - failed to send ID data\n";
		return -1;
	}

	//
	// send the sections
	//

	for (j = 0; j<numSections; j++) {
		if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
			opserr << "AxEqDispBeamColumn2d::sendSelf() - section " <<
				j << "failed to send itself\n";
			return -1;
		}
	}

	return 0;
}

// method receiving element information (for output)
int
AxEqDispBeamColumn2d::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	//
	// receive the integer data containing tag, numSections and coord transformation info
	//
	int dbTag = this->getDbTag();
	int i;

	static Vector data(14);

	if (theChannel.recvVector(dbTag, commitTag, data) < 0) {
		opserr << "AxEqDispBeamColumn2d::recvSelf() - failed to recv data Vector\n";
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
			opserr << "AxEqDispBeamColumn2d::recvSelf() - failed to obtain a CrdTrans object with classTag " <<
				crdTransfClassTag << endln;
			return -2;
		}
	}
	crdTransf->setDbTag(crdTransfDbTag);

	// invoke recvSelf on the crdTransf object
	if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
		opserr << "AxEqDispBeamColumn2d::sendSelf() - failed to recv crdTranf\n";
		return -3;
	}

	// create a new beamInt object if one needed
	if (beamInt == 0 || beamInt->getClassTag() != beamIntClassTag) {
		if (beamInt != 0)
			delete beamInt;

		beamInt = theBroker.getNewBeamIntegration(beamIntClassTag);

		if (beamInt == 0) {
			opserr << "AxEqDispBeamColumn2d::recvSelf() - failed to obtain the beam integration object with classTag" <<
				beamIntClassTag << endln;
			exit(-1);
		}
	}

	beamInt->setDbTag(beamIntDbTag);

	// invoke recvSelf on the beamInt object
	if (beamInt->recvSelf(commitTag, theChannel, theBroker) < 0)
	{
		opserr << "AxEqDispBeamColumn2d::sendSelf() - failed to recv beam integration\n";
		return -3;
	}


	//
	// recv an ID for the sections containing each sections dbTag and classTag
	//

	ID idSections(2 * nSect);
	int loc = 0;

	if (theChannel.recvID(dbTag, commitTag, idSections) < 0) {
		opserr << "AxEqDispBeamColumn2d::recvSelf() - failed to recv ID data\n";
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
			for (int i = 0; i<numSections; i++)
				delete theSections[i];
			delete[] theSections;
		}

		// create a new array to hold pointers
		theSections = new SectionForceDeformation *[nSect];
		if (theSections == 0) {
			opserr << "AxEqDispBeamColumn2d::recvSelf() - out of memory creating sections array of size " <<
				nSect << endln;
			return -1;
		}

		// create a section and recvSelf on it
		numSections = nSect;
		loc = 0;

		for (i = 0; i<numSections; i++) {
			int sectClassTag = idSections(loc);
			int sectDbTag = idSections(loc + 1);
			loc += 2;
			theSections[i] = theBroker.getNewSection(sectClassTag);
			if (theSections[i] == 0) {
				opserr << "AxEqDispBeamColumn2d::recvSelf() - Broker could not create Section of class type " <<
					sectClassTag << endln;
				exit(-1);
			}
			theSections[i]->setDbTag(sectDbTag);
			if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
				opserr << "AxEqDispBeamColumn2d::recvSelf() - section " << i << " failed to recv itself\n";
				return -1;
			}
		}

	}
	else {

		// 
		// for each existing section, check it is of correct type
		// (if not delete old & create a new one) then recvSelf on it
		//

		loc = 0;
		for (i = 0; i<numSections; i++) {
			int sectClassTag = idSections(loc);
			int sectDbTag = idSections(loc + 1);
			loc += 2;

			// check of correct type
			if (theSections[i]->getClassTag() != sectClassTag) {
				// delete the old section[i] and create a new one
				delete theSections[i];
				theSections[i] = theBroker.getNewSection(sectClassTag);
				if (theSections[i] == 0) {
					opserr << "AxEqDispBeamColumn2d::recvSelf() - Broker could not create Section of class type " <<
						sectClassTag << endln;
					exit(-1);
				}
			}

			// recvSelf on it
			theSections[i]->setDbTag(sectDbTag);
			if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
				opserr << "AxEqDispBeamColumn2d::recvSelf() - section " << i << " failed to recv itself\n";
				return -1;
			}
		}
	}
	return 0;
}

// method printing element information (for output)
void
AxEqDispBeamColumn2d::Print(OPS_Stream &s, int flag)
{
	s << "\nAxEqDispBeamColumn2d, element id:  " << this->getTag() << endln;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tCoordTransf: " << crdTransf->getTag() << endln;
	s << "\tmass density:  " << rho << ", cMass: " << cMass << endln;

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

// method displaynig element (for output) 
int
AxEqDispBeamColumn2d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes, int numModes)
{
	static Vector v1(3);
	static Vector v2(3);

	if (displayMode >= 0) {

		theNodes[0]->getDisplayCrds(v1, fact);
		theNodes[1]->getDisplayCrds(v2, fact);

	}
	else {

		theNodes[0]->getDisplayCrds(v1, 0.);
		theNodes[1]->getDisplayCrds(v2, 0.);

		// add eigenvector values
		int mode = displayMode  *  -1;
		const Matrix &eigen1 = theNodes[0]->getEigenvectors();
		const Matrix &eigen2 = theNodes[1]->getEigenvectors();
		if (eigen1.noCols() >= mode) {
			for (int i = 0; i < 2; i++) {
				v1(i) += eigen1(i, mode - 1)*fact;
				v2(i) += eigen2(i, mode - 1)*fact;
			}
		}
	}

	return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

// method for outputting the element response (set the element response)
Response*
AxEqDispBeamColumn2d::setResponse(const char **argv, int argc,
	OPS_Stream &output)
{
	Response *theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", "AxEqDispBeamColumn2d");
	output.attr("eleTag", this->getTag());
	output.attr("node1", connectedExternalNodes[0]);
	output.attr("node2", connectedExternalNodes[1]);

	// global force - 
	if (strcmp(argv[0], "forces") == 0 || strcmp(argv[0], "force") == 0
		|| strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

		output.tag("ResponseType", "Px_1");
		output.tag("ResponseType", "Py_1");
		output.tag("ResponseType", "Mz_1");
		output.tag("ResponseType", "Px_2");
		output.tag("ResponseType", "Py_2");
		output.tag("ResponseType", "Mz_2");

		theResponse = new ElementResponse(this, 1, P);


		// local force -
	}
	else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {

		output.tag("ResponseType", "N1");
		output.tag("ResponseType", "V1");
		output.tag("ResponseType", "M1");
		output.tag("ResponseType", "N2");
		output.tag("ResponseType", "V2");
		output.tag("ResponseType", "M2");

		theResponse = new ElementResponse(this, 2, P);


		// basic force -
	}
	else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) {

		output.tag("ResponseType", "N");
		output.tag("ResponseType", "M1");
		output.tag("ResponseType", "M2");

		theResponse = new ElementResponse(this, 9, Vector(3));

		// basic stiffness -
	}
	else if (strcmp(argv[0], "basicStiffness") == 0) {

		output.tag("ResponseType", "N");
		output.tag("ResponseType", "M1");
		output.tag("ResponseType", "M2");

		theResponse = new ElementResponse(this, 19, Matrix(3, 3));

		// chord rotation -
	}
	else if (strcmp(argv[0], "chordRotation") == 0 || strcmp(argv[0], "chordDeformation") == 0
		|| strcmp(argv[0], "basicDeformation") == 0) {

		output.tag("ResponseType", "eps");
		output.tag("ResponseType", "theta1");
		output.tag("ResponseType", "theta2");

		theResponse = new ElementResponse(this, 3, Vector(3));

		// plastic rotation -
	}
	else if (strcmp(argv[0], "plasticRotation") == 0 || strcmp(argv[0], "plasticDeformation") == 0) {

		output.tag("ResponseType", "epsP");
		output.tag("ResponseType", "theta1P");
		output.tag("ResponseType", "theta2P");

		theResponse = new ElementResponse(this, 4, Vector(3));

	}
	else if (strcmp(argv[0], "RayleighForces") == 0 ||
		strcmp(argv[0], "rayleighForces") == 0 ||
		strcmp(argv[0], "dampingForces") == 0) {

		theResponse = new ElementResponse(this, 12, P);
	}

	// section response -
	else if (strstr(argv[0], "sectionX") != 0) {
		if (argc > 2) {
			float sectionLoc = atof(argv[1]);

			double xi[maxNumSections];
			double L = crdTransf->getInitialLength();
			beamInt->getSectionLocations(numSections, L, xi);

			sectionLoc /= L;

			float minDistance = fabs(xi[0] - sectionLoc);
			int sectionNum = 0;
			for (int i = 1; i < numSections; i++) {
				if (fabs(xi[i] - sectionLoc) < minDistance) {
					minDistance = fabs(xi[i] - sectionLoc);
					sectionNum = i;
				}
			}

			output.tag("GaussPointOutput");
			output.attr("number", sectionNum + 1);
			output.attr("eta", xi[sectionNum] * L);

			theResponse = theSections[sectionNum]->setResponse(&argv[2], argc - 2, output);
		}
	}
	else if (strstr(argv[0], "section") != 0) {

		if (argc > 1) {

			int sectionNum = atoi(argv[1]);

			if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {

				output.tag("GaussPointOutput");
				output.attr("number", sectionNum);
				double xi[maxNumSections];
				double L = crdTransf->getInitialLength();
				beamInt->getSectionLocations(numSections, L, xi);
				output.attr("eta", xi[sectionNum - 1] * L);

				if (strcmp(argv[2], "dsdh") != 0) {
					theResponse = theSections[sectionNum - 1]->setResponse(&argv[2], argc - 2, output);
				}
				else {
					int order = theSections[sectionNum - 1]->getOrder();
					theResponse = new ElementResponse(this, 76, Vector(order));
					Information &info = theResponse->getInformation();
					info.theInt = sectionNum;
				}

				output.endTag();

			}
			else if (sectionNum == 0) { // argv[1] was not an int, we want all sections, 

				CompositeResponse *theCResponse = new CompositeResponse();
				int numResponse = 0;
				double xi[maxNumSections];
				double L = crdTransf->getInitialLength();
				beamInt->getSectionLocations(numSections, L, xi);

				for (int i = 0; i<numSections; i++) {

					output.tag("GaussPointOutput");
					output.attr("number", i + 1);
					output.attr("eta", xi[i] * L);

					Response *theSectionResponse = theSections[i]->setResponse(&argv[1], argc - 1, output);

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
	else if (strcmp(argv[0], "dcurvdh") == 0)
		return new ElementResponse(this, 5, Vector(numSections));

	// basic deformation sensitivity
	else if (strcmp(argv[0], "dvdh") == 0)
		return new ElementResponse(this, 6, Vector(3));

	else if (strcmp(argv[0], "integrationPoints") == 0)
		return new ElementResponse(this, 7, Vector(numSections));

	else if (strcmp(argv[0], "integrationWeights") == 0)
		return new ElementResponse(this, 8, Vector(numSections));

	output.endTag();

	if (theResponse == 0)
		return Element::setResponse(argv, argc, output);
	else
		return theResponse;
}

// method for outputting the element response (retrieve the element response)
int
AxEqDispBeamColumn2d::getResponse(int responseID, Information &eleInfo)
{
	double V;
	double L = crdTransf->getInitialLength();

	if (responseID == 1)
		return eleInfo.setVector(this->getResistingForce());

	else if (responseID == 12) {
		P.Zero();
		P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
		return eleInfo.setVector(P);

	}
	else if (responseID == 2) {
		P(3) = q(0);
		P(0) = -q(0) + p0[0];
		P(2) = q(1);
		P(5) = q(2);
		V = (q(1) + q(2)) / L;
		P(1) = V + p0[1];
		P(4) = -V + p0[2];
		return eleInfo.setVector(P);
	}

	else if (responseID == 9) {
		return eleInfo.setVector(q);
	}

	else if (responseID == 19) {
		static Matrix kb(3, 3);
		this->getBasicStiff(kb);
		return eleInfo.setMatrix(kb);
	}

	// Chord rotation
	else if (responseID == 3) {
		return eleInfo.setVector(crdTransf->getBasicTrialDisp());
	}

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
			locs(i) = xi[i] * L;
		return eleInfo.setVector(locs);
	}

	else if (responseID == 8) {
		//const Vector &wts = quadRule.getIntegrPointWeights(numSections);
		double wt[maxNumSections];
		beamInt->getSectionWeights(numSections, L, wt);
		Vector weights(numSections);
		for (int i = 0; i < numSections; i++)
			weights(i) = wt[i] * L;
		return eleInfo.setVector(weights);
	}

	else
		return Element::getResponse(responseID, eleInfo);
}
// method for outputting the element response (retrieve the element response if sensitivity is involved)
int
AxEqDispBeamColumn2d::getResponseSensitivity(int responseID, int gradNumber,
	Information &eleInfo)
{
	// Basic deformation sensitivity
	if (responseID == 3) {
		const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
		return eleInfo.setVector(dvdh);
	}

	// Basic force sensitivity
	else if (responseID == 9) {
		static Vector dqdh(3);

		dqdh.Zero();

		return eleInfo.setVector(dqdh);
	}

	// dsdh
	else if (responseID == 76) {

		int sectionNum = eleInfo.theInt;
		int order = theSections[sectionNum - 1]->getOrder();
		const ID &code = theSections[sectionNum - 1]->getType();

		Vector dsdh(order);
		dsdh = theSections[sectionNum - 1]->getStressResultantSensitivity(gradNumber, true);

		const Vector &v = crdTransf->getBasicTrialDisp();
		const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

		double L = crdTransf->getInitialLength();
		double oneOverL = 1.0 / L;

		const Matrix &ks = theSections[sectionNum - 1]->getSectionTangent();

		Vector dedh(order);

		//const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
		double xi[maxNumSections];
		beamInt->getSectionLocations(numSections, L, xi);

		double x = xi[sectionNum - 1];

		//double xi6 = 6.0*pts(i,0);
		double xi6 = 6.0*x;

		int j;
		for (j = 0; j < order; j++) {
			switch (code(j)) {
			case SECTION_RESPONSE_P:
				dedh(j) = oneOverL*dvdh(0); break;
			case SECTION_RESPONSE_MZ:
				dedh(j) = oneOverL*((xi6 - 4.0)*dvdh(1) + (xi6 - 2.0)*dvdh(2)); break;
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
AxEqDispBeamColumn2d::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 1)
		return -1;

	// If the parameter belongs to the element itself
	if (strcmp(argv[0], "rho") == 0)
		return param.addObject(1, this);

	if (strstr(argv[0], "sectionX") != 0) {
		if (argc < 3)
			return -1;

		float sectionLoc = atof(argv[1]);

		double xi[maxNumSections];
		double L = crdTransf->getInitialLength();
		beamInt->getSectionLocations(numSections, L, xi);

		sectionLoc /= L;

		float minDistance = fabs(xi[0] - sectionLoc);
		int sectionNum = 0;
		for (int i = 1; i < numSections; i++) {
			if (fabs(xi[i] - sectionLoc) < minDistance) {
				minDistance = fabs(xi[i] - sectionLoc);
				sectionNum = i;
			}
		}
		return theSections[sectionNum]->setParameter(&argv[2], argc - 2, param);
	}
	// If the parameter belongs to a section or lower
	if (strstr(argv[0], "section") != 0) {

		if (argc < 3)
			return -1;

		// Get section number: 1...Np
		int sectionNum = atoi(argv[1]);

		if (sectionNum > 0 && sectionNum <= numSections)
			return theSections[sectionNum - 1]->setParameter(&argv[2], argc - 2, param);

		else
			return -1;
	}

	else if (strstr(argv[0], "integration") != 0) {

		if (argc < 2)
			return -1;

		return beamInt->setParameter(&argv[1], argc - 1, param);
	}
	int result = -1;
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
AxEqDispBeamColumn2d::updateParameter(int parameterID, Information &info)
{
	if (parameterID == 1) {
		rho = info.theDouble;
		return 0;
	}
	else
		return -1;
}




int
AxEqDispBeamColumn2d::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}



const Matrix &
AxEqDispBeamColumn2d::getKiSensitivity(int gradNumber)
{
	K.Zero();
	return K;
}

const Matrix &
AxEqDispBeamColumn2d::getMassSensitivity(int gradNumber)
{
	K.Zero();

	if (rho == 0.0 || parameterID != 1)
		return K;

	double L = crdTransf->getInitialLength();
	if (cMass == 0) {
		// lumped mass matrix
		//double m = 0.5*rho*L;
		double m = 0.5*L;
		K(0, 0) = K(1, 1) = K(3, 3) = K(4, 4) = m;
	}
	else {
		// consistent mass matrix
		static Matrix ml(6, 6);
		//double m = rho*L/420.0;    
		double m = L / 420.0;
		ml(0, 0) = ml(3, 3) = m*140.0;
		ml(0, 3) = ml(3, 0) = m*70.0;

		ml(1, 1) = ml(4, 4) = m*156.0;
		ml(1, 4) = ml(4, 1) = m*54.0;
		ml(2, 2) = ml(5, 5) = m*4.0*L*L;
		ml(2, 5) = ml(5, 2) = -m*3.0*L*L;
		ml(1, 2) = ml(2, 1) = m*22.0*L;
		ml(4, 5) = ml(5, 4) = -ml(1, 2);
		ml(1, 5) = ml(5, 1) = -m*13.0*L;
		ml(2, 4) = ml(4, 2) = -ml(1, 5);

		// transform local mass matrix to global system
		K = crdTransf->getGlobalMatrixFromLocal(ml);
	}

	return K;
}



const Vector &
AxEqDispBeamColumn2d::getResistingForceSensitivity(int gradNumber)
{
	double L = crdTransf->getInitialLength();
	double oneOverL = 1.0 / L;

	//const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
	//const Vector &wts = quadRule.getIntegrPointWeights(numSections);
	double xi[maxNumSections];
	beamInt->getSectionLocations(numSections, L, xi);
	double wt[maxNumSections];
	beamInt->getSectionWeights(numSections, L, wt);

	double dLdh = crdTransf->getdLdh();
	double d1oLdh = crdTransf->getd1overLdh();

	double dptsdh[maxNumSections];
	beamInt->getLocationsDeriv(numSections, L, dLdh, dptsdh);

	double dwtsdh[maxNumSections];
	beamInt->getWeightsDeriv(numSections, L, dLdh, dwtsdh);

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
		const Vector &dsdh = theSections[i]->getStressResultantSensitivity(gradNumber, true);

		// Perform numerical integration on internal force gradient
		double sensi;
		for (int j = 0; j < order; j++) {
			sensi = dsdh(j)*wti;
			switch (code(j)) {
			case SECTION_RESPONSE_P:
				dqdh(0) += sensi;
				break;
			case SECTION_RESPONSE_MZ:
				dqdh(1) += (xi6 - 4.0)*sensi;
				dqdh(2) += (xi6 - 2.0)*sensi;
				break;
			default:
				break;
			}
		}

		const Vector &s = theSections[i]->getStressResultant();

		double dxi6dh = 6.0*dptsdh[i];
		double dwtLdh = wt[i] * dLdh + dwtsdh[i] * L;
		//dwtLdh = dwtsdh[i];

		// Perform numerical integration on internal force gradient
		for (int j = 0; j < order; j++) {
			switch (code(j)) {
			case SECTION_RESPONSE_P:
				//dqdh(0) += d1oLdh*s(j)*wti*L;
				break;
			case SECTION_RESPONSE_MZ:
				//dqdh(1) += (dxi6dh*oneOverL + d1oLdh*(xi6-4.0))*s(j)*wti*L;
				//dqdh(2) += (dxi6dh*oneOverL + d1oLdh*(xi6-2.0))*s(j)*wti*L;

				//dqdh(1) += oneOverL*(xi6-4.0)*s(j)*dwtLdh;
				//dqdh(2) += oneOverL*(xi6-2.0)*s(j)*dwtLdh;
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
		static Matrix kbmine(3, 3);
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
				switch (code(j)) {
				case SECTION_RESPONSE_P:
					q(0) += si;
					for (k = 0; k < order; k++) {
						ka(k, 0) += ks(k, j)*wti;
					}
					break;
				case SECTION_RESPONSE_MZ:
					q(1) += (xi6 - 4.0)*si;
					q(2) += (xi6 - 2.0)*si;
					for (k = 0; k < order; k++) {
						tmp = ks(k, j)*wti;
						ka(k, 1) += (xi6 - 4.0)*tmp;
						ka(k, 2) += (xi6 - 2.0)*tmp;
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
						kbmine(0, k) += ka(j, k);
					}
					break;
				case SECTION_RESPONSE_MZ:
					for (k = 0; k < 3; k++) {
						tmp = ka(j, k);
						kbmine(1, k) += (xi6 - 4.0)*tmp;
						kbmine(2, k) += (xi6 - 2.0)*tmp;
					}
					break;
				default:
					break;
				}
			}
		}

		const Vector &A_u = crdTransf->getBasicTrialDisp();
		double dLdh = crdTransf->getdLdh();
		double d1overLdh = -dLdh / (L*L);
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
AxEqDispBeamColumn2d::commitSensitivity(int gradNumber, int numGrads)
{
	// Get basic deformation and sensitivities
	const Vector &v = crdTransf->getBasicTrialDisp();

	static Vector dvdh(3);
	dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

	double L = crdTransf->getInitialLength();
	double oneOverL = 1.0 / L;
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
			switch (code(j)) {
			case SECTION_RESPONSE_P:
				e(j) = oneOverL*dvdh(0)
					+ d1oLdh*v(0);
				break;
			case SECTION_RESPONSE_MZ:
				e(j) = oneOverL*((xi6 - 4.0)*dvdh(1) + (xi6 - 2.0)*dvdh(2))
					+ d1oLdh*((xi6 - 4.0)*v(1) + (xi6 - 2.0)*v(2));
				break;
			default:
				e(j) = 0.0;
				break;
			}
		}

		// Set the section deformations
		theSections[i]->commitSensitivity(e, gradNumber, numGrads);
	}

	return 0;
}


// AddingSensitivity:END /////////////////////////////////////////////


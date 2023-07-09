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
                                                                        
// $Revision: 1.1
// $Date: 
// $Source$
                                                                        
// Written: Kathryn A. Petek
// Created: 02/04
//
// Revisions
//    02/04 created
//    11/10 F.Mckenna and C.McGann: changes for incorporation into main source code
//
// Description: This file contains the implementation for the SimpleContact2D class.
//

#include "SimpleContact2D.h"
#include <Information.h>
#include <ElementResponse.h>

#include <ID.h> 

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <G3Globals.h>
#include <ErrorHandler.h>
#include <Parameter.h>
#include <ContactMaterial2D.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h> 

#include <elementAPI.h>
#define OPS_Export 

static int num_SimpleContact2D = 0;

OPS_Export void *
OPS_SimpleContact2D(void)
{
  if (num_SimpleContact2D == 0) {
    num_SimpleContact2D++;
    //OPS_Error("SimpleContact2D element - Written: K.Petek, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n", 1);
    opserr << "SimpleContact2D element - Written: K.Petek, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 8) {
    opserr << "Invalid #args,  want: element SimpleContact2D eleTag? iNode? jNode? secondaryNode? lambdaNode? matTag? tolGap? tolForce?\n";
    return 0;
  }
    
  int    iData[6];
  double dData[2];

  int numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SimpleContact2DElement" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element SimpleContact2D " << iData[0] << endln;
    return 0;	
  }

  int matID = iData[5];
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element SimpleContact2D " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  // Parsing was successful, allocate the material
  theElement = new SimpleContact2D(iData[0], iData[1], iData[2], iData[3], iData[4], 
				   *theMaterial,
				   dData[0], dData[1]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SimpleContact2DElement\n";
    return 0;
  }

  return theElement;
}


// constructors:
SimpleContact2D::SimpleContact2D(int tag,  int Nd1, int Nd2, 
				 int NdS, int NdL, NDMaterial &theMat, double tolG, double tolF )
 : Element(tag,ELE_TAG_SimpleContact2D),     
   externalNodes(SC_NUM_NODE),
   tangentStiffness(SC_NUM_DOF, SC_NUM_DOF),
   internalForces(SC_NUM_DOF),
   n(SC_NUM_NDF), 
   T(SC_NUM_NDF),
   Bn(SC_NUM_DDOF), 
   Bs(SC_NUM_DDOF),
   dcrd1(SC_NUM_NDF),
   dcrd2(SC_NUM_NDF),
   dcrdS(SC_NUM_NDF),
   dispL(SC_NUM_NDF)
{
#ifdef DEBUG
	opserr << "SimpleContact2D::SimpleContact2D(...)" << endln;
#endif
	externalNodes(0) = Nd1;
	externalNodes(1) = Nd2;
	externalNodes(2) = NdS;
	externalNodes(3) = NdL;

        MyTag = tag;

	tolGap = tolG;
	tolForce = tolF;

	inContact          = true;
	was_inContact      = true;
	should_be_released = false;
	to_be_released     = false;
	in_bounds          = false;

	gap    = 0.0;
	slip   = 0.0;
	lambda = 0.0;
	
	xsi_n  = 0.0;	// this is the state variable, not xsi_nplus1

	theMaterial = 0;
	NDMaterial *theMatCopy = theMat.getCopy("ContactMaterial2D");
	if (theMatCopy != 0) {
	  theMaterial = (ContactMaterial2D *)theMatCopy;
	} else {
	  opserr << "SimpleContact2D::SimpleContact2D - material needs to be of type Contact2D for ele: " << this->getTag() << endln;
	}

	if (theMaterial == 0) {
	  opserr << "SimpleContact2D::SimpleContact2D - failed allocate material model pointer\n";
	  exit(-1);
	}
#ifdef DEBUG
        if (DEBUG_LEVEL > 1) {
			opserr << "   external nodes: " << externalNodes;
			opserr << "   T:  " << T ;
			opserr << "   n:  " << n ;
			opserr << "   Bn: " << Bn;
			opserr << "   Bs: " << Bs;
		}
#endif
}

SimpleContact2D::SimpleContact2D()
 :Element(0,ELE_TAG_SimpleContact2D),     
   externalNodes(SC_NUM_NODE),
   tangentStiffness(SC_NUM_DOF, SC_NUM_DOF),
   internalForces(SC_NUM_DOF),
   n(SC_NUM_NDF), 
   T(SC_NUM_NDF),
   Bn(SC_NUM_DDOF), 
   Bs(SC_NUM_DDOF),
   dcrd1(SC_NUM_NDF),
   dcrd2(SC_NUM_NDF),
   dcrdS(SC_NUM_NDF),
   dispL(SC_NUM_NDF)
{
    MyTag = 0;

#ifdef DEBUG
	opserr << "SimpleContact2D::SimpleContact2D()" << endln;
#endif
}

//  destructor:
SimpleContact2D::~SimpleContact2D()
{
#ifdef DEBUG
	opserr << "SimpleContact2D::~SimpleContact2D(): " << MyTag << endln;
#endif
	delete theMaterial;
}


int
SimpleContact2D::getNumExternalNodes(void) const
{
#ifdef DEBUG
	opserr << "SimpleContact2D::getNumExternalNodes(...): " << MyTag << endln;
#endif
    return SC_NUM_NODE;
}

const ID &
SimpleContact2D::getExternalNodes(void) 
{
#ifdef DEBUG
	opserr << "SimpleContact2D::getExternalNodes(...): " << MyTag << endln;
#endif
    return externalNodes;
}


Node **
SimpleContact2D::getNodePtrs(void)
{
#ifdef DEBUG
	opserr << "SimpleContact2D::getNodePtrs(...): " << MyTag << endln;
#endif
  return theNodes;			
}

int
SimpleContact2D::getNumDOF(void) 
{
#ifdef DEBUG
	opserr << "SimpleContact2D::getNumDOF(...): " << MyTag << endln;
#endif
    return SC_NUM_DOF;
}


void
SimpleContact2D::setDomain(Domain *theDomain)
{
#ifdef DEBUG
	opserr << "SimpleContact2D::setDomain(...): " << MyTag << endln;
#endif
	int Nd1 = externalNodes(0);
	int Nd2 = externalNodes(1);
	int NdS = externalNodes(2);
	int NdL = externalNodes(3);

	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);
	theNodes[2] = theDomain->getNode(NdS);
	theNodes[3] = theDomain->getNode(NdL);

	int i;
	for (i = 0; i < 4; i++) {
	   if (theNodes[i] == 0)
		return;  // don't go any further - otherwise segmentation fault
	}

	dcrd1 = theNodes[0]->getCrds();
	dcrd2 = theNodes[1]->getCrds();
	dcrdS = theNodes[2]->getCrds();
	dispL.Zero();

	// length of primary segment
	Vector L = (dcrd2 - dcrd1);
	Lprimary  = L.Norm();
	Lsquare  = L^L;

	// adjust cohesion force
	theMaterial->ScaleCohesion(Lprimary);
	theMaterial->ScaleTensileStrength(Lprimary);
	
	// error check that node 1 and node 2 are not in same location
	if (Lprimary < tolGap ) { 
	  opserr << "SimpleContact2D::SimpleContact2D - node 1 and node 2 share same coordinates\n";
	  opserr << "Program Terminated\n";
	  exit(-1);
	}
	
	
	// tangent vector
	T = L/Lprimary;
	
	// normal vector to primary surface
	n(0) = -T(1);
	n(1) =  T(0);
	// n(2) =  0.0;
	
	// initialize xsi_n
	xsi_n = (2*dcrdS - dcrd1 - dcrd2 )^T / Lprimary;
	
	
	// call the base class method
        this->DomainComponent::setDomain(theDomain);

	//this->update();

#ifdef DEBUG
		if (DEBUG_LEVEL > 1) {
			opserr << "   T: " << T ;
			opserr << "   n: " << n ;
		}
#endif

}   	 


int
SimpleContact2D::commitState()
{
#ifdef DEBUG
	opserr << "SimpleContact2D::commitState(): " << MyTag << endln;
#endif

	xsi_n          = xsi_nplus1;  
	was_inContact  = ( gap < tolGap );
	in_bounds      = ( fabs(xsi_n) <= 1 );
	to_be_released = (should_be_released || !in_bounds);
	inContact = ( was_inContact && !to_be_released && in_bounds );


	int retVal = 0;
	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "SimpleContact2D::commitState () - failed in base class";
		}    
	retVal = theMaterial->commitState();
	return retVal; 
}

int
SimpleContact2D::revertToLastCommit()
{
#ifdef DEBUG
    opserr << "SimpleContact2D:: revertToLastCommit\n";
#endif
	
	return theMaterial->revertToLastCommit();
}

int
SimpleContact2D::revertToStart()
{
#ifdef DEBUG
	  opserr << "SimpleContact2D:: revertToStart()\n";
#endif
	  
	  inContact = false;
	  was_inContact = false;
	  should_be_released = false;
	  to_be_released = false;
	  in_bounds = false;
	  
	  return theMaterial->revertToStart();
}

int
SimpleContact2D::update(void)
{
#ifdef DEBUG
	opserr << "SimpleContact2D::update(void): " << MyTag <<endln;
#endif

	double tensileStrength;

	dcrd1 = theNodes[0]->getCrds() + theNodes[0]->getTrialDisp();
	dcrd2 = theNodes[1]->getCrds() + theNodes[1]->getTrialDisp();
	dcrdS = theNodes[2]->getCrds() + theNodes[2]->getTrialDisp();
	dispL = theNodes[3]->getTrialDisp();


	// Vector t = T;
		
	N1 = 0.5*(1 - xsi_n);
	N2 = 0.5*(1 + xsi_n);

	gap = n ^ ( dcrdS - N1*dcrd1 - N2*dcrd2 );  

	xsi_nplus1 = (2*dcrdS - dcrd1 - dcrd2)^T / Lprimary;
		
	Bn(0) = -N1*n(0);
	Bn(1) = -N1*n(1);
	Bn(2) = -N2*n(0);
	Bn(3) = -N2*n(1);
	Bn(4) = n(0);
	Bn(5) = n(1);

	Bs(0) = -N1*T(0);
	Bs(1) = -N1*T(1);
	Bs(2) = -N2*T(0);
	Bs(3) = -N2*T(1);
	Bs(4) = T(0);
	Bs(5) = T(1);

	lambda = dispL(0);

	tensileStrength = theMaterial->getTensileStrength();

	// default: not in contact ( inContact = false)
	slip = 0;				
	
	// define state of the contact element
	should_be_released = ( lambda <= -(tensileStrength + tolForce ) );

#ifdef DEBUG
	if (DEBUG_LEVEL > 0) {
		opserr << "   CONTACT:            " << inContact << endln;
		opserr << "   should be released: " << should_be_released << endln;
	}
#endif

	if (inContact) {	

	   slip = 0.5 * (xsi_nplus1 - xsi_n) * Lprimary;
	
	   Vector strain(3);

	   strain(0) = gap; 
	   strain(1) = slip;
	   strain(2) = lambda;
	
       theMaterial->setTrialStrain(strain);
	}

	else if (to_be_released) {
	// prevents sliding & stabilizes behavior in lag step

		Vector strain(3);

        strain(0) = gap; 
        strain(1) = 0.0;
        strain(2) = lambda;
       
        theMaterial->setTrialStrain(strain);
        }

  return 0;
}


const Matrix &
SimpleContact2D::getTangentStiff(void)
{
#ifdef DEBUG
	opserr << "SimpleContact2D::getTangentStiff: " << MyTag << endln;
#endif

	// initialize Kt
	tangentStiffness.Zero();

	if (inContact) {		// in contact

            Matrix Cmat = theMaterial->getTangent();

            double Cnl = Cmat(0,2);
            double Css = Cmat(1,1);
            double Csl = Cmat(1,2);

#ifdef DEBUG
	if (DEBUG_LEVEL > 2) {
		opserr << "   C_nl = " << Cnl
               << "   C_ss = " << Css
               << "   C_sl = " << Csl
               << endln;
	}
#endif

            int i, j;

            // frictionless contact part

            if (Cnl != 0.0) {
#ifdef DEBUG
	if (DEBUG_LEVEL > 1) {
		opserr << "   ** tangent: normal" << endln;
	}
#endif
                // assume Cnl == 1 (!)

                for (i = 0; i < SC_NUM_DDOF; i++) {

                    tangentStiffness(i,SC_NUM_DDOF) -= Bn(i);
                    tangentStiffness(SC_NUM_DDOF,i) -= Bn(i);
                    }
		tangentStiffness(SC_NUM_DOF-1,SC_NUM_DOF-1) = 1.0;
                }

            // contribution from friction

            if (Css != 0.0) {

            // sticking
#ifdef DEBUG
	if (DEBUG_LEVEL > 1) {
        opserr << "   ** tangent: sticking" << endln;
	}
#endif

                for (i = 0; i < SC_NUM_DDOF; i++) {
                    for (j = 0; j < SC_NUM_DDOF; j++) {

                            tangentStiffness(i,j) += Bs(i)*Bs(j)*Css;
                        }
                    }
                }       // end sticking

            if (Csl != 0.0) {

            // sliding
#ifdef DEBUG
	if (DEBUG_LEVEL > 1) {
        opserr << "   ** tangent: sliding" << endln;
	}
#endif
                for (i = 0; i < SC_NUM_DDOF; i++) {				
                    // assume first is the row
                    tangentStiffness(i,SC_NUM_DDOF) += Bs(i)*Csl;
                    }

                }       // end sliding


	} else {				// not in contact

		tangentStiffness(SC_NUM_DOF-2,SC_NUM_DOF-2) = 1.0;
		tangentStiffness(SC_NUM_DOF-1,SC_NUM_DOF-1) = 1.0;
	}

#ifdef DEBUG
	if (DEBUG_LEVEL > 1) {
		opserr << "   K_T: " << tangentStiffness;
	}
#endif
	return tangentStiffness;
}

const Matrix &
SimpleContact2D::getInitialStiff(void)
{
#ifdef DEBUG
	opserr << "SimpleContact2D::getInitialStiff: " << MyTag << endln;
#endif

  return getTangentStiff();
}
    
void 
SimpleContact2D::zeroLoad(void)
{
#ifdef DEBUG
	opserr << "SimpleContact2D::zeroLoad(): " << MyTag << endln;
#endif

  return;
}

int 
SimpleContact2D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
#ifdef DEBUG
	opserr << "SimpleContact2D::addLoad(...): " << MyTag << endln;
#endif

  return 0;
}

int 
SimpleContact2D::addInertiaLoadToUnbalance(const Vector &accel)
{
#ifdef DEBUG
	opserr << "SimpleContact2D::addInertiaLoadToUnbalance(...): " << MyTag << endln;
#endif

  return 0;
}

const Vector &
SimpleContact2D::getResistingForce()
{
#ifdef DEBUG
	opserr << "SimpleContact2D::getResistingForce(): " << MyTag << endln;;
#endif

	double	t_n = 0.0;	// normal force
	double	t_s = 0.0;	// shear force

	int	i;	// loop counter
		
	// initialize Fi
	internalForces.Zero();

	if (inContact) {	// in contact

            // get contact stresse/forces
            Vector stress = theMaterial->getStress();

            t_n = stress(0);
            t_s = stress(1);

            for (i=0; i<6; i++) {
                internalForces(i) = -t_n*Bn(i) + t_s*Bs(i);
            }
            internalForces(6) = -gap;
	}
        else {
            internalForces(6) = lambda;
        }

#ifdef DEBUG
            opserr << "   Gap   = " << gap    << "   Slip = " << slip << endln;
            opserr << "   T_n   = " << t_n    << "   T_s  = " << t_s  << endln;
            opserr << "   F_int = " << internalForces;
#endif

	return internalForces;
}


const Vector &
SimpleContact2D::getResistingForceIncInertia()
{
#ifdef DEBUG
	opserr << "SimpleContact2D::getResistingForceIncInertia: " << MyTag << endln;
#endif

  return getResistingForce();
}


int
SimpleContact2D::sendSelf(int commitTag, Channel &theChannel)
{
  	int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // SimpleContact2D packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(6);
  data(0) = this->getTag();
  data(1) = SC_NUM_DOF;
  data(2) = tolGap;
  data(3) = tolForce;
  data(4) = theMaterial->getClassTag();

  int matDbTag = theMaterial->getDbTag();

  // NOTE: we do have to ensure that the material has a database
  // tag if we are sending to a database channel.
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  data(5) = matDbTag;

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING SimpleContact2D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }	      

  // SimpleContact2D then sends the tags of it's four nodes

  res = theChannel.sendID(dataTag, commitTag, externalNodes);
  if (res < 0) {
    opserr <<"WARNING SimpleContact2D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }

  // finally SimpleContact2D asks it's material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr <<"WARNING SimpleContact2D::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    return -3;
  }

  return 0;
}

int
SimpleContact2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  // NOTE: THIS HAS BEEN COPIED & MODIFIED FROM SimpleContact2D.cpp
  int res;
  int dataTag = this->getDbTag();

  // SimpleContact2D creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector

  static Vector data(6);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING SimpleContact2D::recvSelf() - failed to receive Vector\n";
    return -1;
  }	      

  this->setTag((int)data(0));
  // SC_NUM_DOF = (int)data(1);	// this must not be since SC_NUM_DOF is used to initialize variables and thus must not change
  tolGap = data(2);
  tolForce = data(3);
  
  // SimpleContact2D now receives the tags of it's four external nodes
  res = theChannel.recvID(dataTag, commitTag, externalNodes);
  if (res < 0) {
    opserr <<"WARNING SimpleContact2D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  // finally SimpleContact2D creates a material object of the correct type,
  // sets its database tag and asks this new object to recveive itself.

  int matClass = (int)data(4);
  int matDb = (int)data(5);

  // check if we have a material object already & if we do if of right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

    // if old one .. delete it
    if (theMaterial != 0)
      delete theMaterial;

    // create a new material object
    NDMaterial *theMatCopy = theBroker.getNewNDMaterial(matClass);
    theMaterial = (ContactMaterial2D *)theMatCopy;

    if (theMaterial == 0) {
      opserr <<"WARNING SimpleContact2D::recvSelf() - " << this->getTag() 
	<< " failed to get a blank Material of type " << matClass << endln;
      return -3;
    }
  }

  theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr <<"WARNING SimpleContact2D::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
    return -3;    
  }

  return 0;
}


int
SimpleContact2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  return 0;
}


void
SimpleContact2D::Print(OPS_Stream &s, int flag)
{
 	opserr << "SimpleContact2D, element id:  " << this->getTag() << endln;
	opserr << "   Connected external nodes:  " ; 
	for (int i = 0; i<SC_NUM_NODE; i++)
	{
		opserr << externalNodes(i)<< " ";
	}
	return;
}


Response*
SimpleContact2D::setResponse(const char **argv, int argc, OPS_Stream &eleInfo)
{
#ifdef DEBUG
	opserr << "SimpleContact2D::setResponse(...): " << MyTag << endln;
#endif

    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
    	return new ElementResponse(this, 1, Vector(2));
	
    } else if (strcmp(argv[0],"frictionforce") == 0 || strcmp(argv[0],"frictionforces") == 0) {
    	return new ElementResponse(this, 2, Vector(2));

    } else if (strcmp(argv[0],"forcescalar") == 0 || strcmp(argv[0],"forcescalars") == 0) {
    	return new ElementResponse(this, 3, Vector(2));
    
    // otherwise response quantity is unknown for the SimpleContact2D class
     } else {
		 opserr << "SimpleContact2D::setResponse(const char **argv, int argc, OPS_Stream &eleInfo): "
		  << argv[0] << " unknown request" << endln;
         return 0;
	 }
}


int 
SimpleContact2D::getResponse(int responseID, Information &eleInfo)
{
#ifdef DEBUG
	opserr << "SimpleContact2D::getResponse(...): " << MyTag << endln;
#endif

	Vector force(2);

	// get contact stresse/forces
    Vector stress = theMaterial->getStress();	

  	if (responseID == 1) {

	  	force = stress(0)*n + stress(1)*T;
    	return eleInfo.setVector(force);

  	} else if (responseID == 2) {

	  	force = stress(1)*T;
    	return eleInfo.setVector(force);

  	} else if (responseID == 3) {

	  	force(0) = stress(0);
      	force(1) = stress(1);
    	return eleInfo.setVector(force);
  
    } else {
    	return -1;
	}
}

int
SimpleContact2D::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"friction") == 0) {
		return param.addObject(1, this);
	}
	
	return -1;
}

int
SimpleContact2D::updateParameter(int parameterID, Information &info)
{
	int res = -1;
	int matRes =  theMaterial->updateParameter(parameterID, info);
	if (matRes != -1) {
		res = matRes;
	}
	return res;
}

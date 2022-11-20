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
                                                                        
// 'Grip n Grab' ratcheting, tension-only device material model
// - file edited from EPPGapMaterial file
// Jarrod Cook, University of Canterbury, Christchurch, New Zealand

//	^
//      |
//	|                ________(3)________
//      |               /                  /
//	F              /                  /
//	O             /                  /
//	R            /                  /
//	C          (2)                (4)
//	E          /                  /
//	|         /                  /
//	|        /                  /
//	|__(1)__/     _____(5)_____/
//  --------------DISPLACEMENT---------->
//
//	LOADING
// 		BELOW ENGAGEMENT THRESHOLD (1)
// 		ELASTIC REGION (2)
//		BEYOND YIELD (3)
//	UNLOADING
//		ELASTIC RECOVERY (4) 
// 		BELOW ENGAGEMENT THRESHOLD (5)

/* ************************************************************************** */

//Include Directives
//The first part of the file contains the list of includes. It is necessary to have an 
//#include directive for each class and api file that is used within the .cpp file and 
//is not included in the header.

#include <stdlib.h>

#include <GNGMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

//Added when trying to add args to eleResponse
#include <MaterialResponse.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Stream.h>

//External Procedure
//This is the all important external procedure that the interpreter will parse when it 
//comes across your element on the command line. You need to parse the command line,
//create a material using the command line arguments you parsed and then return this
//material. The name of the procedure must be OPS_YourClassName (no exceptions). If this
//procedure is missing or the name is incorrect, your material will fail to load.

//NOTE: parsing the command line is easy with some other procedures that are defined in 
//the elementAPI.h file. In the example we show how to get integer and double values from 
//the command line. Other options such as character strings and obtaining the number of 
//input arguments are also available.

static int numGNGMaterials = 0;

void* OPS_GNGMaterial()
{
	
	if (numGNGMaterials == 0) {
		numGNGMaterials++;
		//opserr << "Grip 'n' Grab device installed in this structure!" << endln;
    }
	
    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 4) {
      opserr << "Invalid #args,  want: uniaxialMaterial GNG tag E sigY P <eta>" << endln;
	return 0;
    }
  
    int tag;
    double dData[4];
    dData[3] = 0.0; // setting default eta to 0.

    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
	opserr << "WARNING invalid tag for uniaxialMaterial GNG" << endln;
	return 0;
    }

    numData = OPS_GetNumRemainingInputArgs();
    if(numData > 4) numData = 4;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxial GNG" << endln;
	return 0;	
    }

    // Parsing was successful, allocate the material
    theMaterial = new GNGMaterial(tag, dData[0], dData[1], dData[2], dData[3]);
    if (theMaterial == 0) {
      opserr << "WARNING could not create uniaxialMaterial of type GNG" << endln;
	return 0;
    }

    return theMaterial;
}

//full constructor
GNGMaterial::GNGMaterial(int tag, double e, double sigY0, double p, double eta0)//, int accum)
:UniaxialMaterial(tag,MAT_TAG_GNG),
 commitStrain(0.0), trialStrain(0.0), E(e), sigY(sigY0), P(p), eta(eta0), epsE(0.0),
 epsP(0.0), sigP(0.0), pdemand(0.0), nratchet(0)
{
	if (E == 0.0) {
	  opserr << "GNGMaterial::GNGMaterial -- E is zero, continuing with E = sigY/0.002" << endln;
	  if (sigY != 0.0)
	    E = fabs(sigY)/0.002;
	  else {
	    opserr << "GNGMaterial::GNGMaterial -- E and sigY are zero" << endln;
	    exit(-1);
	  }
	}
	else
		
	  epsY = epsE + sigY/E;

	if (sigY*P<0) { // To Remove...
	  opserr << "GNGMaterial::GNGMaterial -- Alternate signs on sigY and E encountered, continuing anyway" << endln;
	}
        
        if ( (eta >= 1) || (eta <= -1) ) {
          opserr << "GNGMaterial::GNGMaterial -- value of eta must be -1 <= eta <= 1, setting eta to 0" << endln;
          eta = 0;
        }
        
}

//null constructor
GNGMaterial::GNGMaterial()
:UniaxialMaterial(0,MAT_TAG_GNG),
 commitStrain(0.0), trialStrain(0.0),
 E(0.0), sigY(0.0), P(0.0), eta(0.0), epsY(0.0), epsE(0.0),
 epsP(0.0), sigP(0.0), pdemand(0.0), nratchet(0.0)
{

}

//Destructor
//Then we provide the destructor. In the destructor all memory that the the object created 
//or was passed to it in the constructor must be destroyed. For this example we have no 
//such memory. We could have left the destructor out entirely. Hoowever, it is good 
//practice to leave it in your source code.

//destructor
GNGMaterial::~GNGMaterial()
{
  // does nothing
}

//setTrialStrain() Method
//This, as mentioned, is the method called when the element has computed a new strain 
//for the element. The element will make subsequent calls to getStress() and getTangent() 
//to obtain new values of these for the new strain. This is typically the most complicated 
//method to write and to determine the theory for before you even write the code. All 
//subsequent methods are trivial.

int 
GNGMaterial::setTrialStrain(double strain, double strainRate)
{
	
  // set the trial strain
  trialStrain = strain;

  // determine trial stress and tangent
  
  // LOADING
  if (trialStrain > epsP) {
    
    // BEYOND YIELD (3)
    if (trialStrain >= epsY) {
        trialStress = sigY + eta*E*(trialStrain-epsY); // PLASTIC
		trialTangent = eta*E;
        
    // BELOW ENGAGEMENT THRESHOLD (1)
    } else if (trialStrain <= epsE) {
        trialStress = 0; // NO STRESS
		trialTangent = 0;
        
    }
	// ELASTIC REGION (2)
	else {
        trialStress = E*(trialStrain-epsE); // ELASTIC
		trialTangent = E;
    }
	
  }
  // UNLOADING 
  else {
    
    // BELOW ENGAGEMENT THRESHOLD (5)
    if (trialStrain <= epsE) {
        trialStress = 0; // NO STRESS
		trialTangent = 0;
        
    // ELASTIC RECOVERY (4)
    } else {
        trialStress = E*(trialStrain-epsE); // ELASTIC
		trialTangent = E;
    }
  }
  
  if (trialStrain < 0) {
	  
	  trialTangent = 0;
	  
  }
  
  return 0;
  
}

//Trivial Methods
//Next comes 3 rather simple methods that return basic information computed in the 
//setTrialStrain(). You do of course have the option to ignore the setTrialStrain() 
//method and compute the stress and tangent quantities again in the interests of 
//saving memory.

//send back the strain
double 
GNGMaterial::getStrain(void)
{
    return trialStrain;
}

//send back the stress
double 
GNGMaterial::getStress(void)
{
  return trialStress;

}

//send back the tangent
double 
GNGMaterial::getTangent(void)
{
  return trialTangent;
}

//send back the tangent
double 
GNGMaterial::getInitialTangent(void)
{
  if (epsE > 0.0) 
    return 0.0; 
  else 
    return E;
}

//Methods Dealing With Current State
//As mentioned, when the algorithm finds a solution state as it goes from one 
//converged solution to the next. As it attempts to find these solutions it goes 
//through a number of trial steps (each setTrialStrain() is invoked in each of these 
//steps). Once it finds a trial step that is on the solution path it will stop and 
//invoke commitState() on the material. Any state variables that the material uses 
//needs to be updated at this time. Should the algorithm fail to find a solution it 
//may return to the last converged step or indeed the start. You the developer must 
//provide code so that your material can indeed go back to these states and report 
//correct getTangent() and getStress() values for subsequent analysis atte,pts.

int 
GNGMaterial::commitState(void)
{
	//update state variables for next step
	
// LOADING
if (trialStrain > epsP) {
    
    // BEYOND YIELD (3)
    if (trialStrain >= epsY) {
        epsE = trialStrain - trialStress/E; // UPDATE X AXIS CROSSING
		
		if (epsP > epsY) { //UPDATE CUMULATIVE PLASTIC DEMAND
			pdemand = pdemand + trialStrain - epsP;
		}
		else {
			pdemand = pdemand + trialStrain - epsY;
		}	
		
    }
	
// UNLOADING
    
} else {
    
    // BELOW ENGAGEMENT THRESHOLD (5)
    if (trialStrain <= epsE) {
        if (trialStrain < (epsE - P)) { // CHECK FOR RATCHETING
		
			epsE = epsE - P; //*****LIMITED TO SINGLE RATCHET*****// max dy/dt appears to be < 5e-4
			
            epsY = epsE + sigY/E; // NEW YIELD STRAIN
			
			nratchet = nratchet + 1;
        }
        
    // ELASTIC RECOVERY (4)
    } else { 
        if (sigP > sigY) {
            sigY = sigP; // NEW YIELD STRESS
            epsY = epsE + sigY/E; // NEW YIELD STRAIN
        }
    }
}
	
	epsP = trialStrain;
	sigP = trialStress;

    commitStrain = trialStrain;

    return 0;
	
}


int 
GNGMaterial::revertToLastCommit(void)
{
	
    trialStrain = commitStrain;

    return 0;
}


int 
GNGMaterial::revertToStart(void)
{	
    commitStrain = 0.0;
    trialStrain = 0.0;
	
	pdemand = 0.0;
	nratchet = 0;
	
	epsP = 0.0;
	sigP = 0.0;
	epsE = 0.0;
	
	epsY = epsE + sigY/E;

    return 0;
}

//getCopy() Method
//This is the method called by each element or section to get unique copies of a material.

UniaxialMaterial *
GNGMaterial::getCopy(void)
{
    GNGMaterial *theCopy = new GNGMaterial(this->getTag(),E,sigY,P,eta);

    theCopy->commitStrain = commitStrain;
    theCopy->trialStrain = trialStrain;	
	
	theCopy->epsP = epsP;
	theCopy->sigP = sigP;
    theCopy-> epsE = epsE;
	theCopy->epsY = epsY;
	theCopy->pdemand = pdemand;
	theCopy->nratchet = nratchet;

	theCopy->trialStress = trialStress;
	theCopy->trialTangent = trialTangent;	
	
    return theCopy;
}

//Methods Dealing With Databases/Parallel Processing
//There are two methods provided which are required if the user wishes to use the
//database or parallel processing features of the OpenSees applications. If neither
//are to be used, the developer need simply return a negative value in both methods.
//The idea is that the material must pack up it's information using Vector and ID
//objects and send it off to a Channel object. On the flip side, the receiving blank
//element must receive the same Vector and ID data, unpack it and set the variables.

int 
GNGMaterial::sendSelf(int cTag, Channel &theChannel)
{
	//we place all the data needed to define the material and its state
	//into a vector object
  static Vector data(12);
  data(0) = this->getTag();
  data(1) = commitStrain;
  data(2) = E;
  data(3) = sigY;
  data(4) = P;
  data(5) = eta;
  data(6) = epsY;
  data(7) = epsE;
  data(8) = epsP;
  data(9) = sigP;
  data(10) = pdemand;
  data(11) = nratchet;

  //send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "GNGMaterial::sendSelf() - failed to send data" << endln;
    return -1;
  }
  
  return 0;
}

int 
GNGMaterial::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
	//receive the vector object from the channel which defines material
	//parameters and state
  static Vector data(12);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "GNGMaterial::recvSelf() - failed to recv data" << endln;
    return -1;
  }
  else {
    this->setTag((int)data(0));
    commitStrain = data(1);
    trialStrain = commitStrain;
    E = data(2);
    sigY = data(3);
    P = data(4);
    eta = data(5);
    epsY = data(6);
    epsE = data(7);
	epsP = data(8);
	sigP = data(9);
	pdemand = data(10);
	nratchet = (int)data(11);
	
  }

  return 0;
}

//Methods Dealing With Output
//Information is obtained by the user when the print command is invoked by the user 
//and also when the user issues the recorder command. When the print command is invoked 
//the Print method is invoked. This method simply prints information about the element, 
//and then asks the material to do likewise.

void 
GNGMaterial::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "GNG tag: " << this->getTag() << endln;
		s << "  E: " << E << ", kinematic hardening ratio: " << eta << endln;
		s << "  sigY: " << sigY << endln;
		s << "  P: " << P << endln;
		s << " plastic demand: " << pdemand << endln;
		s << " ratchet count: " << nratchet << endln;
	}
	
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"GNG\", ";
		s << "\"E\": " << E << ", ";
		s << "\"eta\": " << eta << ", ";
		s << "\"sigY\": " << sigY << ", ";
		s << "\"P\": " << P << ", ";
		s << "\"plastic demand\": " << pdemand << ", ";
		s << "\"ratchet count\": " << nratchet << ", ";
	}
}

//Responses available to recorders

Response* 
GNGMaterial::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{
  //by default, See if the response is one of the defaults
  Response *res =  UniaxialMaterial::setResponse(argv, argc, theOutput);
  if (res != 0)
    return res;

  if (strcmp(argv[0],"demand") == 0) {
    return new MaterialResponse(this, 11, 0.0);
  }
  else if (strcmp(argv[0],"ratchetCount") == 0) {
    return new MaterialResponse(this, 12, 0.0);
  } 	

  return 0;
}

int 
GNGMaterial::getResponse(int responseID, Information &matInfo)
{
  if (responseID==11) {
    return matInfo.setDouble(pdemand);
  }
  else if (responseID==12) {
    return matInfo.setDouble(nratchet);
  }
  
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return UniaxialMaterial::getResponse(responseID, matInfo);
}

int
GNGMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"E") == 0 || strcmp(argv[0],"k") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"Fy") == 0 || strcmp(argv[0],"fy") == 0 || strcmp(argv[0],"sigY") == 0) {
    param.setValue(sigY);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"eta") == 0) {
    param.setValue(eta);
    return param.addObject(4, this);
  }
  
  return 0;
}

int
GNGMaterial::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    E = info.theDouble;
    break;
  case 2:
    sigY = info.theDouble;
    break;
  case 4:
    eta = info.theDouble;
    break;    
  default:
    return -1;
  }

  return 0;
}


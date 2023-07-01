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

// $Date: 2018/11/25 
                                                                        
// Written: Yi Xiao Tongji University/University of Washington 
// Modified from GNGMaterial by Jook
//
// Ratchet Material 
//
// What: "@(#) Ratchet.C, revA"

#include <stdlib.h>   //add

#include <elementAPI.h>
#include <Ratchet.h>

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

static int numRatchet = 0;

void *
OPS_Ratchet(void)
{
  // print out some KUDO's
  if (numRatchet == 0) {
    opserr << "Ratchet device installed in this structure!\n";
    numRatchet =1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 3) {
	opserr << "Invalid #args,  want: uniaxialMaterial Ratchet tag E toothSize";
	return 0;
    }
  
    int tag;
    double dData[2];

    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
	opserr << "WARNING invalid tag for uniaxialMaterial Ratchet" << endln;
	return 0;
    }

    numData = OPS_GetNumRemainingInputArgs();
    if(numData > 3) numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
	opserr << "Invalid E or toothSize data for uniaxial Ratchet \n";
	return 0;	
    }

    // Parsing was successful, allocate the material
    theMaterial = new Ratchet(tag, dData[0], dData[1]);
    if (theMaterial == 0) {
	opserr << "WARNING could not create uniaxialMaterial of type Ratchet\n";
	return 0;
    }

    return theMaterial;
}


//full constructor
Ratchet::Ratchet(int tag, double e, double tSize)//, int accum)
:UniaxialMaterial(tag, MAT_TAG_Ratchet),
 commitStrain(0.0),commitStress(0.0),commitTangent(e),commitEngageStrain(0.0),
 trialStrain(0.0),trialStress(0.0),trialTangent(e),
 E(e),toothSize(tSize),engageStrain(0.0),currentStrain(0.0),
 nratchet(0),commitNratchet(0)
{
  if (E <= 0.0) {
    opserr << "Ratchet::Ratchet -- E is zero or negative\n";
	  exit(-1);
  } 
}

//null constructor
Ratchet::Ratchet()
:UniaxialMaterial(0, MAT_TAG_Ratchet),
 commitStrain(0.0),commitStress(0.0),commitTangent(0.0),commitEngageStrain(0.0),
 trialStrain(0.0),trialStress(0.0),trialTangent(0.0),
 E(0.0),toothSize(0.0),engageStrain(0.0),currentStrain(0.0),
 nratchet(0),commitNratchet(0)
{

}

Ratchet::~Ratchet()
{
  // does nothing
}

int 
Ratchet::setTrialStrain(double strain, double strainRate)
{
 //   if (fabs(trialStrain - strain) < DBL_EPSILON)
 //     return 0;
 // set the trial strain
   trialStrain = strain;
  
   double preciValue = toothSize/10000000;
   
      // BELOW ENGAGEMENT THRESHOLD
    if (trialStrain < engageStrain-preciValue) {
    	
        trialStress = 0; // NO STRESS
		    trialTangent = 0;   
   
    } else {
        trialStress = E*(trialStrain-engageStrain); // ELASTIC
		    trialTangent = E;
    }
  
  
  //if (trialStress < 0) {
	  
	 // trialStress = 0;
	  
  //}
  
  return 0;
}

double 
Ratchet::getStrain(void)
{
  return trialStrain;
}

double 
Ratchet::getStress(void)
{
  return trialStress;
}

double 
Ratchet::getTangent(void)
{
  return trialTangent;
}

double 
Ratchet::getInitialTangent(void)
{
	double preciValue = toothSize/10000000;
  if (currentStrain < engageStrain-preciValue)
    return 0.0; 
  else 
    return E;
}

int 
Ratchet::commitState(void)
{
   	//update state variables for next step
	
	double preciValue = toothSize/10000000;
	//UNLOADING
	if (trialStrain <= currentStrain){
		
		if (trialStrain <= engageStrain){
			
			while (trialStrain <= (engageStrain - toothSize + preciValue)){
				
				engageStrain = engageStrain - toothSize;
				nratchet = nratchet +1;
				
			}
			
		} 
		
	}
	
	  currentStrain = trialStrain;

    commitStrain = trialStrain;
    commitStress = trialStress;
    commitTangent = trialTangent;
    commitEngageStrain = engageStrain;
    commitNratchet = nratchet;
    return 0;
}	


int 
Ratchet::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialStress = commitStress;
  trialTangent = commitTangent;
  engageStrain = commitEngageStrain;
  nratchet = commitNratchet;
  
  currentStrain = commitStrain;
  
  return 0;
}


int 
Ratchet::revertToStart(void)
{
  commitStrain = 0.0;
  commitStress = 0.0;
  commitTangent = 0.0;
  commitEngageStrain = 0.0;
  trialStrain = 0.0;
  trialStress = 0.0;
  trialTangent = E;
  
  engageStrain = 0.0;
  currentStrain = 0.0;
  
  nratchet = 0;
  commitNratchet = 0; 

    return 0;
}


UniaxialMaterial *
Ratchet::getCopy(void)
{
  Ratchet *theCopy =
    new Ratchet(this->getTag(),E,toothSize);
  theCopy->commitStrain = this->commitStrain;
  theCopy->commitStress = this->commitStress;
  theCopy->commitTangent = this->commitTangent;
  theCopy->commitEngageStrain = this->commitEngageStrain;
  theCopy->trialStrain = this->trialStrain;
  theCopy->trialStress = this->trialStress;
  theCopy->trialTangent = this->trialTangent;
  
  theCopy->engageStrain = this->engageStrain;  
	theCopy->currentStrain = this->currentStrain;
  theCopy->nratchet = this->nratchet;
  theCopy->commitNratchet = this->commitNratchet;
	 
  return theCopy;
}


int 
Ratchet::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
Ratchet::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
   return -1;
}

void 
Ratchet::Print(OPS_Stream &s, int flag)
{
    s << "  Ratchet tag: " << this->getTag() << endln;
    s << "  currentStrain: " << currentStrain << endln;
    s << "  getStrain: " << getStrain() << endln;   
    s << "  engageStrain: " << engageStrain << endln;  
		s << "  getTangent: " << getTangent() << endln;
		s << "  E: " << E << endln;
		s << "  toothSize: " << toothSize << endln;

		s << "  ratchet count: " << nratchet << endln;
		
		s << "  trialStrain -(engageStrain - toothSize + DBL_EPSILON): " << trialStrain -(engageStrain - toothSize + DBL_EPSILON) << endln;
}


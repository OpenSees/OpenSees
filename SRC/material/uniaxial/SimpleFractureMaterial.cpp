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
                                                                        
// Written: fmk
// Created: Feb 2015

#include <stdlib.h>
#include <string.h>

#include <SimpleFractureMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>
#include <cmath>
#include <elementAPI.h>

void *
OPS_SimpleFractureMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  double maxStrain = 1.0e16;
  int    iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 3) {
    opserr << "WARNING invalid uniaxialMaterial SimpleFracture $tag $otherTag $maxStrain>" << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial SimpleFracture $tag $otherTag $maxStrain" << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag:  uniaxialMaterial SimpleFracture $tag $otherTag $max: " << iData[0] << endln;
    return 0;	
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &maxStrain) != 0) {
    opserr << "WARNING invalid maxStrain: uniaxialMaterial  SimpleFracture $tag $otherTag $maxStrain" << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new SimpleFractureMaterial(iData[0], *theOtherMaterial, maxStrain);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type SimpleFractureMaterial\n";
    return 0;
  }

  return theMaterial;
}

SimpleFractureMaterial::SimpleFractureMaterial(int tag, UniaxialMaterial &material, double max)
  :UniaxialMaterial(tag,MAT_TAG_SimpleFractureMaterial), theMaterial(0),
   maxStrain(max), TstartCompStrain(0), CstartCompStrain(0), Tfailed(false), Cfailed(false)
{
  theMaterial = material.getCopy();

  Cstress = theMaterial->getStress();
  Ctangent = theMaterial->getTangent();
  Cstrain = theMaterial->getStrain();
  Tstress = Cstress;
  Ttangent = Ctangent;
  Tstrain = Cstrain;

  if (theMaterial == 0) {
    opserr <<  "SimpleFractureMaterial::SimpleFractureMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

SimpleFractureMaterial::SimpleFractureMaterial()
  :UniaxialMaterial(0,MAT_TAG_SimpleFractureMaterial), theMaterial(0),
   maxStrain(0), TstartCompStrain(0), CstartCompStrain(0), Tfailed(false), Cfailed(false)
{
  Cstress = 0;
  Ctangent = 0;
  Cstrain = 0;
  Tstress = Cstress;
  Ttangent = Ctangent;
  Tstrain = Cstrain;
}

SimpleFractureMaterial::~SimpleFractureMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int 
SimpleFractureMaterial::setTrialStrain(double strain, double strainRate)
{  
  // for now do this until finialize method below
  return this->setTrialStrain(strain, 0, strainRate);
}


int 
SimpleFractureMaterial::setTrialStrain(double strain, double temp, double strainRate)
{
  // determine based on last converged step
  Tfailed = Cfailed;
  TstartCompStrain = CstartCompStrain;

  Tstress = Cstress;
  Tstrain = strain;
  theMaterial->revertToLastCommit();

  if (Tfailed == true && strain >= TstartCompStrain) {

    // material has failed and strain is greater than strain where goes into compression

    Ttangent = 0.;
    Tstress = 0.;
    return 0;

  } else if (Tfailed == false && strain > maxStrain) {
    //    opserr << "(Tfailed == false && strain > maxStrain)\n";
    // material fails for first time, strain is greater than strain where goes into compression

    Tfailed = true;
    Ttangent = 0.;
    Tstress = 0.;
      
    // use divide & conquer to find strain where crosses into compression (0 stress by unloaading material in small increments)
    theMaterial->setTrialStrain(maxStrain);
    double tStress = theMaterial->getStress();
    double dStrain = fabs(strain/1e4);
    while (tStress > 0.) {
      strain -= dStrain;
      theMaterial->setTrialStrain(strain, temp, strainRate);
      tStress = theMaterial->getStress();
    }
    TstartCompStrain = strain;
    return 0;

  } else if (Tfailed == true && strain < TstartCompStrain) {

    //    opserr << "(Tfailed == true && strain < TstartCompStrain): " << TstartCompStrain << endln;
    // set strain in material, if tensile, revert to last & increment ro where crosses, set TsratCompStrain
    theMaterial->setTrialStrain(strain, temp, strainRate);
    Tstress = theMaterial->getStress();
    Ttangent = theMaterial->getTangent();
    double dStrain = fabs(strain/1e4);
    if (Tstress > 0) {
      double tStress = Tstress;
      Ttangent = 0.;
      Tstress = 0.;
      while (tStress > 0.) {
	strain -= dStrain;
	theMaterial->setTrialStrain(strain, temp, strainRate);
	tStress = theMaterial->getStress();
      }
      TstartCompStrain = strain;
    } 

  } else {
    //    opserr << "ELSE\n";
    theMaterial->setTrialStrain(strain, temp, strainRate);
    Ttangent = theMaterial->getTangent();
    Tstress = theMaterial->getStress();
    Tfailed = false;
    
  }
  return 0;
}


double 
SimpleFractureMaterial::getStress(void)
{
  return Tstress;
}

double 
SimpleFractureMaterial::getTangent(void)
{
  return Ttangent;
}

double 
SimpleFractureMaterial::getDampTangent(void)
{
  return theMaterial->getDampTangent();
}


double 
SimpleFractureMaterial::getStrain(void)
{
  return Tstrain;
}

double 
SimpleFractureMaterial::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
SimpleFractureMaterial::commitState(void)
{	
  Cfailed = Tfailed;
  Cstress = Tstress;
  Ctangent = Ttangent;
  Cstrain = Tstrain;
  CstartCompStrain = TstartCompStrain;
  
  return theMaterial->commitState();
}

int 
SimpleFractureMaterial::revertToLastCommit(void)
{
  Tfailed = Cfailed;
  Tstress = Cstress;
  Ttangent = Ctangent;
  Tstrain = Cstrain;
  TstartCompStrain = CstartCompStrain;

  return theMaterial->revertToLastCommit();
}

int 
SimpleFractureMaterial::revertToStart(void)
{
  Tfailed = false;
  Cstrain = 0;
  theMaterial->revertToStart();
  Ctangent = theMaterial->getTangent();
  Cstress = theMaterial->getStress();
  return 0;
}

UniaxialMaterial *
SimpleFractureMaterial::getCopy(void)
{
  SimpleFractureMaterial *theCopy = 
    new SimpleFractureMaterial(this->getTag(), *theMaterial, maxStrain);
        
  return theCopy;
}

int 
SimpleFractureMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();

  int matDbTag = theMaterial->getDbTag();
  if ( matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "SimpleFractureMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(6);
  dataVec(0) = maxStrain;
  if (Cfailed == true)
    dataVec(1) = 1.0;
  else
    dataVec(1) = 0.0;

  dataVec(2) = Cstress;
  dataVec(3) = Cstrain;
  dataVec(4) = Ctangent;
  dataVec(5) = CstartCompStrain;


  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "SimpleFractureMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "SimpleFractureMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
SimpleFractureMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "SimpleFractureMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "SimpleFractureMaterial::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));


  static Vector dataVec(6);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "SimpleFractureMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  maxStrain = dataVec(0);
  
  if (dataVec(2) == 1.0)
    Cfailed = true;
  else
    Cfailed = false;

  Cstress =   dataVec(2);
  Cstrain = dataVec(3);
  Ctangent = dataVec(4);
  CstartCompStrain = dataVec(5);

  this->revertToLastCommit();

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "SimpleFractureMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
SimpleFractureMaterial::Print(OPS_Stream &s, int flag)
{
  s << "SimpleFractureMaterial tag: " << this->getTag() << endln;
  s << "\tMaterial: " << theMaterial->getTag() << endln;
  s << "\tMax strain: " << maxStrain << endln;
}

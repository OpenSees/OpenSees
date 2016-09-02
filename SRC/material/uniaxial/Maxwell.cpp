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
                                                                        
// Written: Dimitrios G. Lignos, PhD, Assistant Professor, McGill University 
// Created: February, 2011
// Revision: A
//
// Description: This file contains the class implementation for Maxwell Model

#include <math.h>

#include <elementAPI.h>
#include <Maxwell.h>
#include <Vector.h>
#include <Channel.h>
#include <string.h>

#include <OPS_Globals.h>

static int numMaxwellMaterials = 0;

void *
OPS_Maxwell()
{
  if (numMaxwellMaterials == 0) {
    numMaxwellMaterials++;
    opserr << "Maxwell Model - D.Lignos, McGill University\n";
  }
  
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  
  int    iData[2];
  double dData[4];
  iData[1] = 0;

  int numRData = OPS_GetNumRemainingInputArgs();
  if (numRData != 5 && numRData != 6) {
    opserr << "Invalid #args for command uniaxialMaterial Maxwell\n";
    return 0;
  }

  // Check tag
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  Maxwell tag" << endln;
    return 0;
  }
  // Check if we have 4 input variables for K, C, Alpha, L
  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial Maxwell tag? K? C? Alpha? Length L?"<< endln;
    return 0;	
  }

  if (numRData == 6) {
    const char *cArray = OPS_GetString();
    // OPS_GetStringCopy(&cArray);
    if ((strcmp(cArray, "-returnD") == 0) || (strcmp(cArray, "-D") == 0)) 
      iData[1] = 1;
    delete [] cArray;
  }      
  
  // Parsing was successful, allocate the material with zero index
  theMaterial = new Maxwell(iData[0], 
			    dData[0], dData[1], dData[2], dData[3], 
			    iData[1]);
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Maxwell Material\n";
    return 0;
  }
  
  return theMaterial;
}


Maxwell::Maxwell(int tag, double k, double c, double a, double l, int retD)
  :UniaxialMaterial(tag,MAT_TAG_Maxwell), K(k), C(c), Alpha(a), L(l), returnD(retD)
{
  if (Alpha < 0.0) {
    opserr << "Maxwell::Maxwell -- Alpha < 0.0, setting to 1.0\n";
    Alpha = 1.0;
  }
  // Sets all history and state variables to initial values
  
  // initialize variables
  
  Tstrain = 0.0;
  Tstress = 0.0;
  Ttangent = K;
  
  Cstrain = 0.0;
  Cstress = 0.0;
  Ctangent = K;
  
}

Maxwell::Maxwell()
 :UniaxialMaterial(0,MAT_TAG_Maxwell),
  K(0.0), C(0.0), Alpha(0.0), L(0.0)
{
  // Do nothing
}

Maxwell::~Maxwell()
{
  // does nothing
}

int 
Maxwell::setTrialStrain(double strain, double strainRate)
{
  // Set Total Strain
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;
  
  // Determine change in strain from last converged state
  double dStrain = strain - Cstrain;
  
  // Constant based on Viscous Damper Viscocity C and Axial Stiffness K
  double tR = (C/pow(L,Alpha))/K;
  
  // Integration Constants as a function of incremental time dt
  // Note: The integration time dt = ops_Dt is coming in from the OPS_Globals
  double A = Tstress * ( exp(-ops_Dt/tR) - 1.0);
  double B = (K/2.0) * (1.0 + exp(-ops_Dt/tR) );
  
  // Incremental Stress Calculation (considers the history of strain in A and B
  double Dstress = A + B * dStrain;
  
  // Total Stress Calculation (Tstress) is increment Stress Plus Previous Total Stress
  Tstress = Dstress + Tstress;
  Tstrain = strain;

  //  double DTangent = Tstress/Tstrain;
  
  return 0;
}

double 
Maxwell::getStress(void)
{
  
  return  Tstress;
}

double 
Maxwell::getTangent(void)
{
  return K;
}

double 
Maxwell::getInitialTangent(void)
{
  return K;
}

double
Maxwell::getDampTangent(void)
{
  if (returnD == 1) {
    double DTangent = Tstress/Tstrain;
    return DTangent;
  }

  return 0.0;
}


double 
Maxwell::getStrain(void)
{
  return Tstrain;
}

double 
Maxwell::getStrainRate(void)
{
  return 0;
}

int 
Maxwell::commitState(void)
{
  //commit trial  variables
  Cstrain = Tstrain;
  Cstress = Tstress;
  Ctangent = Ttangent;
  
  return 0;
}

int 
Maxwell::revertToLastCommit(void)
{
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;
  
  return 0;
}

int 
Maxwell::revertToStart(void)
{
  // Initialize state variables
  Tstrain=0.0;
  Tstress=0.0;
  Ttangent = K;
  
  Cstrain=0.0;
  Cstress = 0.0;
  Ctangent = K;
  
  return 0;
}

UniaxialMaterial *
Maxwell::getCopy(void)
{
  Maxwell *theCopy = new Maxwell(this->getTag(), K, C, Alpha, L, returnD);
  // Converged state variables
  theCopy->Cstrain = Cstrain;
  theCopy->Cstress = Cstress;
  theCopy->Ctangent = Ctangent;
  
  // Trial state variables
  theCopy->Tstrain = Tstrain;
  theCopy->Tstress = Tstress;	
  theCopy->Ttangent = Ttangent;
  
  return theCopy;
}

int 
Maxwell::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(9);
  data(0) = this->getTag();
  data(8) = returnD;

  // Material properties
  data(1) = K;
  data(2) = C;
  data(3) = Alpha;
  data(4) = L;
  
  // State variables from last converged state
  data(5) = Cstrain;
  data(6) = Cstress;
  data(7) = Ctangent;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Maxwell::sendSelf() - failed to send data\n";

  return res;
}

int 
Maxwell::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(9);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);

  
  if (res < 0) {
      opserr << "Maxwell::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {

    this->setTag((int)data(0));
    returnD = (int)data(8);

    // Material properties
    K = data(1);
    C = data(2);
    Alpha = data(3);
    L = data(4);
    
    // State variables from last converged state 
    Cstrain = data(5);
    Cstress = data(6);
    Ctangent = data(7);
    
    //Copy converged state values into trial values
    Tstrain = Cstrain;
    Tstress = Cstress;
    Ttangent = Ctangent;
  }
    
  return res;
}

void 
Maxwell::Print(OPS_Stream &s, int flag)
{
  s << "Maxwell tag: " << this->getTag() << endln;
  s << "  K: " << K << endln;	
  s << "  C: " << C << endln;
  s << "  Alpha: " << Alpha << endln;
  s << "  Length: " << L << endln;
}



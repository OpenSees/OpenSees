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

// $Revision: 1.3 $
// $Date: 2006-03-16 19:28:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FatigueMaterial.cpp,v $

// Written: Patxi 
// Created: Aug 2003
//
// Description: This file contains the class definition for 
// FatigueMaterial.  FatigueMaterial wraps a UniaxialMaterial
// and imposes fatigue limits. More information about this material can
// be found in the doctoral dissertation of Patxi Uriz:
//
//   Uriz, Patxi, "Towards Earthquake Resistant Design of 
//      Concentrically Braced Steel Frames," Ph.D. Dissertation, 
//      Structural Engineering, Mechanics, and Materials, Civil 
//      and Envrironmental Engineering, University of California, 
//      Berkeley, December 2005
//
// Modifications to the code have been added by Kevin Mackie, Ph.D., as
// shown below.

#include <stdlib.h>
#include <MaterialResponse.h>
#include <Information.h>


#include <FatigueMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

FatigueMaterial::FatigueMaterial(int tag, UniaxialMaterial &material,
				 double dmax, double E_0, double slope_m, 
				 double epsmin, double epsmax )
  :UniaxialMaterial(tag,MAT_TAG_Fatigue), theMaterial(0), 
   Cfailed(false), trialStrain(0)
{
  DI  = 0; //Damage index
  X   = 0; //Range in consideration
  Y   = 0; //Previous Adjacent Range
  A   = 0; //Peak or valley 1
  B   = 0; //Peak or valley 2
  C   = 0; //Peak or valley 2
  D   = 0; //Peak or valley 4
  PCC = 0; /*Previous Cycle counter flag if >1 then previous 'n' 
	     cycles did not flag a complete cycle */
  R1F = 0; //Flag for first  peak count
  R2F = 0; //Flag for second peak count
  CS  = 0; //Current Slope
  PS  = 0; //Previous slope
  EP  = 0; //Previous Strain
  SF  = 0; /*Start Flag = 0 if very first strain, (i.e. when initializing)
	     = 1 otherwise */
  DL  = 0; //Damage if current strain was last peak.

  if ( dmax > 1.0 || dmax < 0.0 ) {
    opserr << "FatigueMaterial::FatigueMaterial - Dmax must be between 0 and 1, assuming Dmax = 1\n";
    Dmax = 1;
  } else 
    Dmax      = dmax;
  
  E0        = E_0;
  m         = slope_m;
  minStrain = epsmin;
  maxStrain = epsmax;
  
  theMaterial = material.getCopy();
  
  if (theMaterial == 0) {
    opserr <<  "FatigueMaterial::FatigueMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

FatigueMaterial::FatigueMaterial()
  :UniaxialMaterial(0,MAT_TAG_Fatigue), theMaterial(0), 
   Cfailed(false), trialStrain(0)
{
  DI  = 0; //Damage index
  X   = 0; //Range in consideration
  Y   = 0; //Previous Adjacent Range
  A   = 0; //Peak or valley 1
  B   = 0; //Peak or valley 2
  C   = 0; //Peak or valley 2
  D   = 0; //Peak or valley 4
  PCC = 0; /*Previous Cycle counter flag if >1 then previous 'n' 
	     cycles did not flag a complete cycle */
  R1F = 0; //Flag for first  peak count
  R2F = 0; //Flag for second peak count
  CS  = 0; //Current Slope
  PS  = 0; //Previous slope
  EP  = 0; //Previous Strain
  SF  = 0; /*Start Flag = 0 if very first strain, (i.e. when initializing)
	     = 1 otherwise */
  DL  = 0; //Damage if current strain was last peak.

  Dmax    = 0;
  E0      = 0; 
  m       = 0;
  minStrain    = 0;
  maxStrain    = 0;
}

FatigueMaterial::~FatigueMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

static int sign(double a) {
  if (a < 0)
    return -1;
  else if (a == 0)
    return 0;
  else
    return 1;
}

int 
FatigueMaterial::setTrialStrain(double strain, double strainRate)
{
  if (Cfailed) {
    trialStrain = strain;
    // return 0;
    return theMaterial->setTrialStrain(strain, strainRate);
  } else {
    Cfailed = false;
    trialStrain = strain;
    return theMaterial->setTrialStrain(strain, strainRate);
  }
}

double 
FatigueMaterial::getStress(void)
{
  double modifier = 1.0;
  double damageloc = 1.0-Dmax+DL;
  if (Cfailed)
    // Reduce stress to 0.0 
    return theMaterial->getStress()*1.0e-8;
  
  /*
  // This portion of the code was added by Kevin Mackie, Ph.D.
  //  This is appropriate for steel material.  Uncomment to use.
  else if ( damageloc <= 0.9 )
    modifier = 1.0-725.0/2937.0*pow(damageloc,2);
  else 
    modifier = 8.0*(1.0-damageloc);
  
  if (modifier <= 0)
    modifier = 1.0e-8;
  
  else
    return theMaterial->getStress()*modifier;
  */

  else
    return theMaterial -> getStress();
}

double 
FatigueMaterial::getTangent(void)
{
  double modifier = 1.0;
  double damageloc = 1.0-Dmax+DL;
  if (Cfailed)
    // Reduce tangent to 0.0 
    return 1.0e-8*theMaterial->getInitialTangent();

    /* 
    // This portion of the code was added by Kevin Mackie, Ph.D.
    //  This is appropriate for steel material.  Uncomment to use.
    else if ( damageloc <= 0.9 )
      modifier = 1.0-725.0/2937.0*pow(damageloc,2);
    else
      modifier = 8.0*(1.0-damageloc);
  
    if (modifier <= 0)
      // modifier = 1.0e-8;
      modifier = 1.0e-3;
    */

  else
    return theMaterial->getTangent()*modifier;
}

double 
FatigueMaterial::getDampTangent(void)
{
  if (Cfailed)
    return 0.0;
  else
    return theMaterial->getDampTangent();
}



double 
FatigueMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

double 
FatigueMaterial::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
FatigueMaterial::commitState(void)
{	

  // NOTE: Do not accumulate damage if peaks are too small (e.g. if X < 1e-10)
  // may get floating point errors.  This is essentially a filter for 
  // strain cycles smaller than 1e-10 strain.  
  //  THIS FATIGE MATERIAL CODE WAS NOT INTENDED FOR HIGH CYCLE FATIGUE, 
  //  THE LINEAR ACCUMULATION OF DAMAGE IS NOT AS APPROPRIATE FOR HIGH
  //  CYCLE FATIGUE. 


  // No need to continue if the uniaxial material copy 
  // has already failed.
  if (Cfailed) {
    return 0;
  }

  //Simple check to see if we reached max strain capacities
  if (trialStrain >= maxStrain || trialStrain <= minStrain) { 
      Cfailed = true;
      opserr << "FatigueMaterial: material tag " << this->getTag() << " failed from excessive strain\n";
      DI = Dmax;
      DL = Dmax;
      return 0;
  }

  //Initialize the fatigue parameters if they have 
  // not been initialized yet
  if (SF == 0) {

    A   = trialStrain;
    SF  = 1  ;
    EP  = trialStrain;
    // Initialize other params if not done so already
    PCC = 0;
    B   = 0;
    C   = 0;
    D   = 0;

  }

  /* Now we need to determine if we are at a peak or not 
     If we are, then we need to do some calcs to determine
     the amount of damage suffered. If we are not at a peak, we need to
     pretend like we are at a peak, so that we can calculate the damage
     as if it WERE a peak.   
  */

  // Determine the slope of the strain hysteresis
  if ( EP == trialStrain ) {
    CS = PS;         // No real slope here....
  } else {
    CS = trialStrain - EP;   // Determine Current Slope
  }


  // If we are at a peak or a valley, then check for damage
  if (sign(PS) != sign(CS) && sign(PS) != 0 ) {

    if ( R1F == 0 )  {    // mark second peak

      B = EP; 
      Y = fabs(B-A);
      R1F = 1;

    } else {   // start at least the third peak

      // begin modified Rainflow cycle counting
      if (PCC == 1) { 
	
	D = EP;
	X = fabs(D-C);
	
      } else {
	
	C = EP;
	X = fabs(C-B);
	
      }

      if (X < Y) {

	PCC = PCC + 1;

	if (PCC == 1) {
	  Y = fabs(C-B);
	} else if (PCC == 2 ) {
	  // Count X = |D-C| as a 1.0 cycle
	  DI = DI + 1.0 / fabs(pow( (X/E0) , 1/m )) ;
	  // Reset parameters
	  D = 0;
	  C = 0;
	  Y = fabs(B-A);
	  PCC = 0;
	}

      } else {
	
	if (PCC == 1 ) {
	  
	  // Count Y = |C-B| as a 1.0 cycle
	  DI = DI + 1.0 / fabs(pow( (Y/E0) , 1/m ));
	  // Reser parameters
	  B = D;
	  C = 0;
	  D = 0;
	  Y = fabs(B-A);
	  PCC = 0;
	  
	} else {
	  
	  // Count Y = |A-B| as a 0.5 cycle
	  DI = DI + 0.5 / fabs(pow( (Y/E0), 1/m ));
	  // Reset parameters
	  A = B;
	  B = C;
	  C = 0;
	  D = 0;
	  Y = X;
	  PCC = 0;
	  
	}
	  
      }
	
    }

    // Flag failure if we have reached that point
    if (DI >= Dmax )  {
      // Most likely will not fail at this point, more 
      // likely at the psuedo peak. But this step is
      // is important for accumulating damage
      Cfailed = true;
      opserr << "FatigueMaterial: material tag " << this->getTag() << " failed at peak\n";
      DL=DI;
    } else {
      Cfailed = false;
      DL=DI;
    }

    // Modified by Patxi 3/5/2006
  // } else {
  }
  if (Cfailed == false) {
    
    // Now check for damage, although we may not be at a peak at all.
    // Store temporary damage only as if it were the last peak: DL
    // Commit to DI only if failure occurs.
    if (B == 0 && C == 0 &&  D == 0) {
      
      // If we have not yet found the second peak
      X = fabs(trialStrain - A);
      
      if (fabs(X) < 1e-10) {
	DL = DI ;
      } else {
	DL = DI +  0.5 / fabs(pow( (X/E0), 1/m ));
      }
      
    } else if (B != 0 && C == 0 &&  D == 0) {
      
      // On our way to find point C. Range Y defined, no X yet
      X = fabs(trialStrain - B);
      
      if (fabs(X) < 1e-10) {
	DL = DI;
      } else {
	DL = DI +  0.5 / fabs(pow( (X/E0) , 1/m ));
      }	
      
      if (fabs(Y) < 1e-10) {
	DL = DL;
      } else {
	DL = DL +  0.5 / fabs(pow( (Y/E0) , 1/m ));
      }
      
    } else if (B != 0 && C != 0 &&  D == 0) {
      
      // Two ranges stored, but no cycles for either stored
      //   Make sure we get the potential |D-A| range.
      X = fabs(trialStrain-A);
      
      if (fabs(Y) < 1e-10) {
	DL = DI;
      } else {
	DL = DI +  1.0 / fabs(pow( (Y/E0) , 1/m ));
      } 
      
      if (fabs(X) < 1e-10) {
	DL = DL;
      } else {
	DL = DL +  0.5 / fabs(pow( (X/E0) , 1/m ));
      }
      
    }
    
    // Did we fail before a peak?
    double mStress = theMaterial->getStress();
    if (DL > Dmax && mStress > 0.0 ) {
      DI = DL;
      Cfailed = true;
      opserr << "FatigueMaterial: material tag " << this->getTag() << " failed at pseudo peak\n";
    } else {
      Cfailed = false;
    }
    
  }

  PS = CS;            // Previous Slope
  EP = trialStrain;   // Keep track of previous strain
  
  // Check if failed at current step
  if (Cfailed) {
    return 0;
  }
  else 
    return theMaterial->commitState();

}


int 
FatigueMaterial::revertToLastCommit(void)
{
  // Check if failed at last step
  if (Cfailed)
    return 0;
  else
    return theMaterial->revertToLastCommit();
}

// Kevin Mackie, Ph.D. is responsible for the addition of 
// revertToStart, sendSelf, recvSelf.  Thanks to Kevin for making this
// code more user friendly!

int 
FatigueMaterial::revertToStart(void)
{

  Cfailed = false;
  DI  = 0; //Damage index
  X   = 0; //Range in consideration
  Y   = 0; //Previous Adjacent Range
  A   = 0; //Peak or valley 1
  B   = 0; //Peak or valley 2
  C   = 0; //Peak or valley 2
  D   = 0; //Peak or valley 4
  PCC = 0; /*Previous Cycle counter flag if >1 then previous 'n' 
	     cycles did not flag a complete cycle */
  R1F = 0; //Flag for first  peak count
  R2F = 0; //Flag for second peak count
  CS  = 0; //Current Slope
  PS  = 0; //Previous slope
  EP  = 0; //Previous Strain
  SF  = 0; /*Start Flag = 0 if very first strain, (i.e. when initializing)
	     = 1 otherwise */
  DL  = 0; //Damage if current strain was last peak.

  Dmax    = 0;
  E0      = 0; 
  m       = 0;
  minStrain    = 0;
  maxStrain    = 0;

  return theMaterial->revertToStart();
}

UniaxialMaterial *
FatigueMaterial::getCopy(void)
{
  FatigueMaterial *theCopy = 
    new FatigueMaterial(this->getTag(), *theMaterial, Dmax, E0, m ,minStrain, maxStrain);

  theCopy->Cfailed = Cfailed;
  theCopy->trialStrain = trialStrain;

  return theCopy;
}

int 
FatigueMaterial::sendSelf(int cTag, Channel &theChannel)
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
    opserr << "FatigueMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(21);
  dataVec(0)  = DI;
  dataVec(1)  = X;
  dataVec(2)  = Y;
  dataVec(3)  = A;
  dataVec(4)  = B;
  dataVec(5)  = C;
  dataVec(6)  = D;
  dataVec(7)  = PCC;
  dataVec(8)  = R1F;
  dataVec(9)  = R2F;
  dataVec(10) = CS;
  dataVec(11) = PS;
  dataVec(12) = EP;
  dataVec(13) = SF;
  dataVec(14) = DL;
  dataVec(15) = Dmax;
  dataVec(16) = E0;
  dataVec(17) = m;
  dataVec(18) = minStrain;
  dataVec(19) = maxStrain;

  if (Cfailed == true)
    dataVec(20) = 1.0;
  else
    dataVec(20) = 0.0;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "FatigueMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "FatigueMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
FatigueMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "FatigueMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "FatigueMaterial::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(21);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "FatigueMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  DI   = dataVec(0);
  X    = dataVec(1);
  Y    = dataVec(2);
  A    = dataVec(3);
  B    = dataVec(4);
  C    = dataVec(5);
  D    = dataVec(6);
  PCC  = int(dataVec(7));
  R1F  = int(dataVec(8));
  R2F  = int(dataVec(9));
  CS   = dataVec(10);
  PS   = dataVec(11);
  EP   = dataVec(12);
  SF   = int(dataVec(13));
  DL   = dataVec(14);
  Dmax = dataVec(15);
  E0   = dataVec(16);
  m    = dataVec(17);
  minStrain = dataVec(18);
  maxStrain = dataVec(19);

  if (dataVec(20) == 1.0)
    Cfailed = true;
  else
    Cfailed = false;

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "FatigueMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

//Printing of current damage.  NOTE:  The damage that is returned
//  is damage at the psuedo peak, DL, ( if not at a peak when the 
//  print function  is called). The generic print will print both 
//  the damage recorded at the last peak, along with current damage

void 
FatigueMaterial::Print(OPS_Stream &s, int flag)
{
  if (flag == 100) {
    s << DL << endln;
  } else {
    s << "FatigueMaterial tag: " << this->getTag() << endln;
    s << "\tMaterial: " << theMaterial->getTag() << endln;
    s << "\tDI: " << DI << " Dmax: " << Dmax << endln;
    s << "\tE0: " << E0 <<  " m: " << m  << endln;
    s << "\tDL: " << DL << endln;

  }
}

Response* 
FatigueMaterial::setResponse(const char **argv, int argc, Information &matInfo)
{
  if (argc == 0) 
    return 0;

  // stress
  if (strcmp(argv[0],"stress") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  
  // tangent
  else if (strcmp(argv[0],"tangent") == 0)
    return new MaterialResponse(this, 2, this->getTangent());

  // strain
  else if (strcmp(argv[0],"strain") == 0)
    return new MaterialResponse(this, 3, this->getStrain());

  // strain
  else if ((strcmp(argv[0],"stressStrain") == 0) || 
	   (strcmp(argv[0],"stressANDstrain") == 0)) {
    return new MaterialResponse(this, 4, Vector(2));

  }

  else if (strcmp(argv[0],"damage") == 0)
    return new MaterialResponse(this, 5, DL);

  // otherwise unknown
  else {
    return 0;
  }
}

int 
FatigueMaterial::getResponse(int responseID, Information &matInfo)
{
  static Vector stressStrain(2);
  // each subclass must implement its own stuff    
  switch (responseID) {
    case 1:
      matInfo.setDouble(this->getStress());
      return 0;
      
    case 2:
      matInfo.setDouble(this->getTangent());
      return 0;      

    case 3:
      matInfo.setDouble(this->getStrain());
      return 0;      
    
    case 4:
      stressStrain(0) = this->getStress();
      stressStrain(1) = this->getStrain();
      matInfo.setVector(stressStrain);
      return 0;

    case 5:
      matInfo.setDouble(DL);
      return 0;      
      
  default:      
    return -1;
  }
}



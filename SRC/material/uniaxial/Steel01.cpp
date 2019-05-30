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
                                                                        
// $Revision: 1.20 $
// $Date: 2008-08-26 16:35:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel01.cpp,v $
                                                                        
// Written: MHS 
// Created: 06/99
// Revision: A
//
// Description: This file contains the class implementation for 
// Steel01. 
//
// What: "@(#) Steel01.C, revA"


#include <Steel01.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <string.h>

#include <math.h>
#include <float.h>


#include <elementAPI.h>
#include <OPS_Globals.h>


void *
OPS_Steel01()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[7];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel01 tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel01 " << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel01 " << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 3) {
    dData[3] = STEEL_01_DEFAULT_A1;
    dData[4] = STEEL_01_DEFAULT_A2;
    dData[5] = STEEL_01_DEFAULT_A3;
    dData[6] = STEEL_01_DEFAULT_A4;
  }

  // Parsing was successful, allocate the material
  theMaterial = new Steel01(iData[0], dData[0], dData[1], 
			    dData[2], dData[3], dData[4], 
			    dData[5], dData[6]);

  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel01 Material\n";
    return 0;
  }

  return theMaterial;
}



Steel01::Steel01
(int tag, double FY, double E, double B,
 double A1, double A2, double A3, double A4):
   UniaxialMaterial(tag,MAT_TAG_Steel01),
   fy(FY), E0(E), b(B), a1(A1), a2(A2), a3(A3), a4(A4)
{
   // Sets all history and state variables to initial values
   // History variables
	Energy = 0;	//by SAJalali
	
	CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0;

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////
}

Steel01::Steel01():UniaxialMaterial(0,MAT_TAG_Steel01),
 fy(0.0), E0(0.0), b(0.0), a1(0.0), a2(0.0), a3(0.0), a4(0.0)
{
	Energy = 0;	//by SAJalali

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

}

Steel01::~Steel01 ()
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	if (SHVs != 0) 
		delete SHVs;
// AddingSensitivity:END //////////////////////////////////////
}

int Steel01::setTrialStrain (double strain, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }

   return 0;
}

int Steel01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

void Steel01::determineTrialState (double dStrain)
{
      double fyOneMinusB = fy * (1.0 - b);

      double Esh = b*E0;
      double epsy = fy/E0;
      
      double c1 = Esh*Tstrain;
      
      double c2 = TshiftN*fyOneMinusB;

      double c3 = TshiftP*fyOneMinusB;

      double c = Cstress + E0*dStrain;

      /**********************************************************
         removal of the following lines due to problems with
	 optimization may be required (e.g. on gnucc compiler
         with optimization turned on & -ffloat-store option not
         used) .. replace them with line that follows but which 
         now requires 2 function calls to achieve same result !!
      ************************************************************/

      double c1c3 = c1 + c3;

      if (c1c3 < c)
	Tstress = c1c3;
      else
	Tstress = c;

      double c1c2 = c1-c2;

      if (c1c2 > Tstress)
	Tstress = c1c2;

      /* ***********************************************************
      and replace them with:

      Tstress = fmax((c1-c2), fmin((c1+c3),c));
      **************************************************************/

      if (fabs(Tstress-c) < DBL_EPSILON)
	  Ttangent = E0;
      else
	Ttangent = Esh;

      //
      // Determine if a load reversal has occurred due to the trial strain
      //

      // Determine initial loading condition
      if (Tloading == 0 && dStrain != 0.0) {
	  if (dStrain > 0.0)
	    Tloading = 1;
	  else
	    Tloading = -1;
      }

      // Transition from loading to unloading, i.e. positive strain increment
      // to negative strain increment
      if (Tloading == 1 && dStrain < 0.0) {
	  Tloading = -1;
	  if (Cstrain > TmaxStrain)
	    TmaxStrain = Cstrain;
	  TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
      }

      // Transition from unloading to loading, i.e. negative strain increment
      // to positive strain increment
      if (Tloading == -1 && dStrain > 0.0) {
	  Tloading = 1;
	  if (Cstrain < TminStrain)
	    TminStrain = Cstrain;
	  TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
      }
}

void Steel01::detectLoadReversal (double dStrain)
{
   // Determine initial loading condition
   if (Tloading == 0 && dStrain != 0.0)
   {
      if (dStrain > 0.0)
         Tloading = 1;
      else
         Tloading = -1;
   }

   double epsy = fy/E0;

   // Transition from loading to unloading, i.e. positive strain increment
   // to negative strain increment
   if (Tloading == 1 && dStrain < 0.0)
   {
      Tloading = -1;
      if (Cstrain > TmaxStrain)
         TmaxStrain = Cstrain;
      TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
   }

   // Transition from unloading to loading, i.e. negative strain increment
   // to positive strain increment
   if (Tloading == -1 && dStrain > 0.0)
   {
      Tloading = 1;
      if (Cstrain < TminStrain)
         TminStrain = Cstrain;
      TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
   }
}

double Steel01::getStrain ()
{
   return Tstrain;
}

double Steel01::getStress ()
{
   return Tstress;
}

double Steel01::getTangent ()
{
   return Ttangent;
}

int Steel01::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CmaxStrain = TmaxStrain;
   CshiftP = TshiftP;
   CshiftN = TshiftN;
   Cloading = Tloading;

   // State variables
   //by SAJalali
   Energy += 0.5*(Tstress + Cstress)*(Tstrain - Cstrain);

   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int Steel01::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;

   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int Steel01::revertToStart ()
{
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0;

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0;

// AddingSensitivity:BEGIN /////////////////////////////////
	if (SHVs != 0) 
		SHVs->Zero();
// AddingSensitivity:END //////////////////////////////////

   return 0;
}

UniaxialMaterial* Steel01::getCopy ()
{
   Steel01* theCopy = new Steel01(this->getTag(), fy, E0, b,
				  a1, a2, a3, a4);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CmaxStrain = CmaxStrain;
   theCopy->CshiftP = CshiftP;
   theCopy->CshiftN = CshiftN;
   theCopy->Cloading = Cloading;

   // Trial history variables
   theCopy->TminStrain = TminStrain;
   theCopy->TmaxStrain = TmaxStrain;
   theCopy->TshiftP = TshiftP;
   theCopy->TshiftN = TshiftN;
   theCopy->Tloading = Tloading;

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

int Steel01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(16);
   data(0) = this->getTag();

   // Material properties
   data(1) = fy;
   data(2) = E0;
   data(3) = b;
   data(4) = a1;
   data(5) = a2;
   data(6) = a3;
   data(7) = a4;

   // History variables from last converged state
   data(8) = CminStrain;
   data(9) = CmaxStrain;
   data(10) = CshiftP;
   data(11) = CshiftN;
   data(12) = Cloading;

   // State variables from last converged state
   data(13) = Cstrain;
   data(14) = Cstress;
   data(15) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "Steel01::sendSelf() - failed to send data\n";

   return res;
}

int Steel01::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(16);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "Steel01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      fy = data(1);
      E0 = data(2);
      b = data(3);
      a1 = data(4);
      a2 = data(5);
      a3 = data(6);
      a4 = data(7);

      // History variables from last converged state
      CminStrain = data(8);
      CmaxStrain = data(9);
      CshiftP = data(10);
      CshiftN = data(11);
      Cloading = int(data(12));

      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
      TminStrain = CminStrain;
      TmaxStrain = CmaxStrain;
      TshiftP = CshiftP;
      TshiftN = CshiftN;
      Tloading = Cloading;

      // State variables from last converged state
      Cstrain = data(13);
      Cstress = data(14);
      Ctangent = data(15);      

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }
    
   return res;
}

void Steel01::Print (OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {    
    s << "Steel01 tag: " << this->getTag() << endln;
    s << "  fy: " << fy << " ";
    s << "  E0: " << E0 << " ";
    s << "   b: " << b << " ";
    s << "  a1: " << a1 << " ";
    s << "  a2: " << a2 << " ";
    s << "  a3: " << a3 << " ";
    s << "  a4: " << a4 << " ";
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": \"" << this->getTag() << "\", ";
	s << "\"type\": \"Steel01\", ";
	s << "\"E\": " << E0 << ", ";
	s << "\"fy\": " << fy << ", ";
    s << "\"b\": " << b << ", ";
    s << "\"a1\": " << a1 << ", ";
    s << "\"a2\": " << a2 << ", ";
    s << "\"a3\": " << a3 << ", ";
    s << "\"a4\": " << a4 << "}";
  }
  
}




// AddingSensitivity:BEGIN ///////////////////////////////////
int
Steel01::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0 || strcmp(argv[0],"Fy") == 0) {
    param.setValue(fy);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E0);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"b") == 0) {
    param.setValue(b);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"a1") == 0) {
    param.setValue(a1);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"a2") == 0) {
    param.setValue(a2);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"a3") == 0) {
    param.setValue(a3);
    return param.addObject(6, this);
  }
  if (strcmp(argv[0],"a4") == 0) {
    param.setValue(a4);
    return param.addObject(7, this);
  }

  return -1;
}



int
Steel01::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->fy = info.theDouble;
		break;
	case 2:
		this->E0 = info.theDouble;
		break;
	case 3:
		this->b = info.theDouble;
		break;
	case 4:
		this->a1 = info.theDouble;
		break;
	case 5:
		this->a2 = info.theDouble;
		break;
	case 6:
		this->a3 = info.theDouble;
		break;
	case 7:
		this->a4 = info.theDouble;
		break;
	default:
		return -1;
	}

	Ttangent = E0;          // Initial stiffness

	return 0;
}




int
Steel01::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}



double
Steel01::getStressSensitivity(int gradIndex, bool conditional)
{
	// Initialize return value
	double gradient = 0.0;


	// Pick up sensitivity history variables
	double CstrainSensitivity = 0.0;
	double CstressSensitivity = 0.0;
	if (SHVs != 0) {
		CstrainSensitivity = (*SHVs)(0,gradIndex);
		CstressSensitivity = (*SHVs)(1,gradIndex);
	}


	// Assign values to parameter derivatives (depending on what's random)
	double fySensitivity = 0.0;
	double E0Sensitivity = 0.0;
	double bSensitivity = 0.0;
	if (parameterID == 1) {
		fySensitivity = 1.0;
	}
	else if (parameterID == 2) {
		E0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		bSensitivity = 1.0;
	}


	// Compute min and max stress
	double Tstress;
	double dStrain = Tstrain-Cstrain;
	double sigmaElastic = Cstress + E0*dStrain;
	double fyOneMinusB = fy * (1.0 - b);
	double Esh = b*E0;
	double c1 = Esh*Tstrain;
	double c2 = TshiftN*fyOneMinusB;
	double c3 = TshiftP*fyOneMinusB;
	double sigmaMax = c1+c3;
	double sigmaMin = c1-c2;


	// Evaluate stress sensitivity 
	if ( (sigmaMax < sigmaElastic) && (fabs(sigmaMax-sigmaElastic)>1e-5) ) {
		Tstress = sigmaMax;
		gradient = E0Sensitivity*b*Tstrain 
				 + E0*bSensitivity*Tstrain
				 + TshiftP*(fySensitivity*(1-b)-fy*bSensitivity);
	}
	else {
		Tstress = sigmaElastic;
		gradient = CstressSensitivity 
			     + E0Sensitivity*(Tstrain-Cstrain)
				 - E0*CstrainSensitivity;
	}
	if (sigmaMin > Tstress) {
		gradient = E0Sensitivity*b*Tstrain
			     + E0*bSensitivity*Tstrain
				 - TshiftN*(fySensitivity*(1-b)-fy*bSensitivity);
	}

	return gradient;
}




double
Steel01::getInitialTangentSensitivity(int gradIndex)
{
	// For now, assume that this is only called for initial stiffness 
	if (parameterID == 2) {
		return 1.0; 
	}
	else {
		return 0.0;
	}
}


int
Steel01::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(2,numGrads);
	}


	// Initialize unconditaional stress sensitivity
	double gradient = 0.0;


	// Pick up sensitivity history variables
	double CstrainSensitivity = 0.0;
	double CstressSensitivity	 = 0.0;
	if (SHVs != 0) {
		CstrainSensitivity = (*SHVs)(0,gradIndex);
		CstressSensitivity = (*SHVs)(1,gradIndex);
	}


	// Assign values to parameter derivatives (depending on what's random)
	double fySensitivity = 0.0;
	double E0Sensitivity = 0.0;
	double bSensitivity = 0.0;
	if (parameterID == 1) {
		fySensitivity = 1.0;
	}
	else if (parameterID == 2) {
		E0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		bSensitivity = 1.0;
	}


	// Compute min and max stress
	double Tstress;
	double dStrain = Tstrain-Cstrain;
	double sigmaElastic = Cstress + E0*dStrain;
	double fyOneMinusB = fy * (1.0 - b);
	double Esh = b*E0;
	double c1 = Esh*Tstrain;
	double c2 = TshiftN*fyOneMinusB;
	double c3 = TshiftP*fyOneMinusB;
	double sigmaMax = c1+c3;
	double sigmaMin = c1-c2;


	// Evaluate stress sensitivity ('gradient')
	if ( (sigmaMax < sigmaElastic) && (fabs(sigmaMax-sigmaElastic)>1e-5) ) {
		Tstress = sigmaMax;
		gradient = E0Sensitivity*b*Tstrain 
				 + E0*bSensitivity*Tstrain
				 + E0*b*TstrainSensitivity
				 + TshiftP*(fySensitivity*(1-b)-fy*bSensitivity);
	}
	else {
		Tstress = sigmaElastic;
		gradient = CstressSensitivity 
			     + E0Sensitivity*(Tstrain-Cstrain)
				 + E0*(TstrainSensitivity-CstrainSensitivity);
	}
	if (sigmaMin > Tstress) {
		gradient = E0Sensitivity*b*Tstrain
			     + E0*bSensitivity*Tstrain
			     + E0*b*TstrainSensitivity
				 - TshiftN*(fySensitivity*(1-b)-fy*bSensitivity);
	}


	// Commit history variables
	(*SHVs)(0,gradIndex) = TstrainSensitivity;
	(*SHVs)(1,gradIndex) = gradient;

	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////


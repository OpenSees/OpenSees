// TendonL01.cpp
// Written: ALaskar
// Created: 2007.8  
//
// Description: This file contains the class implementation for
// TendonL01
// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.

#include "TendonL01.h"
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include <elementAPI.h>
#define OPS_Export 
#include <MaterialResponse.h>

OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_TendonL01Material)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs < 6) {
      opserr << "Invalid Args want: uniaxialMaterial TendonL01 tag? fpy? Eps? fpu? rou? epsp? <ac?> <rc?>" << endln;
      return 0;	
  }

  int    iData[1];
  double dData[7];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial TendonL01 tag" << endln;
    return 0;
  }

  numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 5) {
    if (OPS_GetDoubleInput(&numRemainingArgs, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial TendonL01 tag? fpy? Eps? fpu? rou? epsp? <ac?> <rc?>" << endln;
      return 0;	
    } else
      theMaterial = new TendonL01(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4]);
  } else if (numRemainingArgs == 7) {
    if (OPS_GetDoubleInput(&numRemainingArgs, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial TendonL01 tag? fpy? Eps? fpu? rou? epsp? <ac?> <rc?>" << endln;
      return 0;	
    } else
      theMaterial = new TendonL01(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6]);
  } else {
    opserr << "Invalid Args want: uniaxialMaterial TendonL01 tag? fpy? Eps? fpu? rou? epsp? <ac?> <rc?>" << endln;
    return 0;
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type TendonL01\n";
    return 0;
  }

  return theMaterial;
}


TendonL01::TendonL01
(int tag, double FPY, double E, double FPU, double ROU, double EPSP, double AC, double RC):
UniaxialMaterial(tag, MAT_TAG_TendonL01),
fpy(FPY), Eps(E), fpu(FPU), rou(ROU), epsp(EPSP), ac(AC), rc(RC)
{
    tt1 = 0.0; // temp add
	tt2 = 0.0; // temp add by zhong
	ttStrain = 0.0;

	//if ( fpc < 0.0 )
	//{
	//	fpc = - fpc; // Make fpc positive to get sqrt(fpc)
	//}
	// Sets all history and state variables to initial values
	this->revertToStart();
    
}

TendonL01::TendonL01(): UniaxialMaterial(0, MAT_TAG_TendonL01),
fpy(0.0), Eps(0.0), fpu(0.0), rou(0.0), epsp(0), ac(0.0), rc(0.0)
{
	this->revertToStart();
}

TendonL01::~TendonL01()
{
	// does nothing;
}

int TendonL01::revertToStart()
{
    // History variables
    CminStrain = 0.0;  
    CmaxStrain = 0.0;  
    CloadingState = 0; 
	CloopPathState = 0;
	
    reverseFromTenEnvelopeStrain = 0.0;  
    reverseFromTenEnvelopeStress = 0.0;  
	approachToComEnvelopeStrain = 0.0;  
	approachToComEnvelopeStress = 0.0;   
    
    reverseFromComEnvelopeStrain = 0.0;  
    reverseFromComEnvelopeStress = 0.0;  
	approachToTenEnvelopeStrain = 0.0;   
	approachToTenEnvelopeStress = 0.0;  
	
	int i;
	for (i=0; i<SIZE1; i++)
	{
		CreverseTopStrain[i] = 0.0;    
		CreverseTopStress[i] = 0.0;
		CreverseBottomStrain[i] = 0.0;
		CreverseBottomStress[i] = 0.0;
	}
	CreverseTopNum = 0;             
	CreverseBottomNum = 0;   

	for (i=0; i<SIZE1; i++)
	{
		TreverseTopStrain[i] = 0.0;    
		TreverseTopStress[i] = 0.0;
		TreverseBottomStrain[i] = 0.0;
		TreverseBottomStress[i] = 0.0;
	}
	TreverseTopNum = 0;             
	TreverseBottomNum = 0;  
	

	downPathPointOneStrain = 0.0;    
	downPathPointTwoStrain = 0.0;
	downPathPointTwoStress = 0.0;
    upPathPointOneStrain = 0.0;      
	upPathPointTwoStrain = 0.0;
	upPathPointTwoStress = 0.0;

    // State variables
	TminStrain = 0.0;  
    TmaxStrain = 0.0;  
    TloadingState = 0;
	TloopPathState = 0;

    Cstrain = 0.0;
    Cstress = 0.0;
    Ctangent = Eps;

    Tstrain = 0.0;
    Tstress = 0.0;
    Ttangent = Eps;

    return 0;
}


int TendonL01::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;  
   TmaxStrain = CmaxStrain;  
   TloadingState = CloadingState;
   TloopPathState = CloopPathState;

   for (int i=0; i<SIZE1; i++)
   {
    	TreverseTopStrain[i] = CreverseTopStrain[i];    
		TreverseTopStress[i] = CreverseTopStress[i];
		TreverseBottomStrain[i] = CreverseBottomStrain[i];
		TreverseBottomStress[i] = CreverseBottomStress[i];
   }
   TreverseTopNum = CreverseTopNum;             
   TreverseBottomNum = CreverseBottomNum; 

   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}




int TendonL01::commitState ()
{
   // Reset trial history variables to last committed state
    CminStrain = TminStrain;  
    CmaxStrain = TmaxStrain;  
    CloadingState = TloadingState;
	CloopPathState = TloopPathState;

   for (int i=0; i<SIZE1; i++)
   {
    	CreverseTopStrain[i] = TreverseTopStrain[i];    
		CreverseTopStress[i] = TreverseTopStress[i];
		CreverseBottomStrain[i] = TreverseBottomStrain[i];
		CreverseBottomStress[i] = TreverseBottomStress[i];
   }
   CreverseTopNum = TreverseTopNum;             
   CreverseBottomNum = TreverseBottomNum; 

   // Reset trial state variables to last committed state
    Cstrain = Tstrain;
    Cstress = Tstress;
    Ctangent = Ttangent;

   return 0;
}


double TendonL01::getCommittedStrain()
{
   return Cstrain;
}


double TendonL01::getStrain ()
{
   return Tstrain;
}

double TendonL01::getStress ()
{
   return Tstress;
}

double TendonL01::getTangent ()
{
   // temp add
   if ( Ttangent == 0.0 )
   {
   	   opserr << " TendonL01:getTangent() -- Ttangent = 0.0\n";
	   opserr<< " Tstrain = " << this->getStrain() << endln;
       opserr<< " Tstress = " << this->getStress() << endln; 
	   opserr<< " CloadingState = " <<CloadingState << endln;
       opserr<< " CloopPathState = " << CloopPathState << endln;
       opserr<< " TloadingState = " <<TloadingState << endln;
       opserr<< " TloopPathState = " << TloopPathState << endln;
	   opserr<< " Cstrain = " << Cstrain << endln;
	   opserr<< " Cstress = " << Cstress << endln;
	   opserr<< " dStrain = " << ttStrain << endln;
	   opserr<< " TreverseTopStrain[TreverseTopNum] = " << TreverseTopStrain[TreverseTopNum] << endln;
       opserr<< " TreverseBottomStrain[TreverseBottomNum] = " << TreverseBottomStrain[TreverseBottomNum] << endln;
	   opserr<< " TreverseBottomNum = " << TreverseBottomNum << endln;
	   opserr<< " approachToComEnvelopeStrain = " << approachToComEnvelopeStrain << endln;

   }

   return Ttangent;
}

double TendonL01::getSecant ()
{
	double TSecant;
	if ( Tstrain == 0.0 )
	{
		TSecant = Eps;
	}
    else
	{
		TSecant = Tstress/Tstrain;
	}
	return TSecant;
}


int TendonL01::setTrialStrain (double strain, double strainRate)
{
   // Reset history variables to last converged state	
    TminStrain = CminStrain;  
    TmaxStrain = CmaxStrain;  
    TloadingState = CloadingState;
	TloopPathState = CloopPathState;

   for (int i=0; i<SIZE1; i++)
   {
    	TreverseTopStrain[i] = CreverseTopStrain[i];    
		TreverseTopStress[i] = CreverseTopStress[i];
		TreverseBottomStrain[i] = CreverseBottomStrain[i];
		TreverseBottomStress[i] = CreverseBottomStress[i];
   }
   TreverseTopNum = CreverseTopNum;             
   TreverseBottomNum = CreverseBottomNum; 

	// Set trial strain
   Tstrain = strain + epsp;

   // Determine change in strain from last converged state
   double dStrain = Tstrain - Cstrain;

   // Calculate the trial state given the trial strain
   //   if (fabs(dStrain) > DBL_EPSILON)
   if (fabs(dStrain) > 1e-10)
   determineTrialState (dStrain);

   ttStrain = dStrain;

   return 0;
}



int TendonL01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
	TminStrain = CminStrain;  
    TmaxStrain = CmaxStrain;  
    TloadingState = CloadingState;
	TloopPathState = CloopPathState;

   for (int i=0; i<SIZE1; i++)
   {
    	TreverseTopStrain[i] = CreverseTopStrain[i];    
		TreverseTopStress[i] = CreverseTopStress[i];
		TreverseBottomStrain[i] = CreverseBottomStrain[i];
		TreverseBottomStress[i] = CreverseBottomStress[i];
   }
   TreverseTopNum = CreverseTopNum;             
   TreverseBottomNum = CreverseBottomNum; 

   // Set trial strain
   Tstrain = strain + epsp;

   // Determine change in strain from last converged state
   double dStrain;
   dStrain = Tstrain - Cstrain;

   // Calculate the trial state given the trial strain
   // if (fabs(dStrain) > DBL_EPSILON) 
   determineTrialState (dStrain);

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}


void TendonL01::determineTrialState (double dStrain)
{
    tt1 = Tstrain;
	tt2 = 0.0;
	ttStrain = dStrain;

	//Set parameter values
	//double epsy = fy/E0;
	//double fcr = 0.31*sqrt(fpc);
	//if ( rou < 0.0025 ) rou = 0.0025;
	//double B = pow((fcr/fy),1.5)/rou;
	//double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);
	double epsy = fpy/Eps;
	double epsn = 0.7 * fpu / Eps;

	if ( TloadingState == 0 ) //Start from initial point
	{
		TloadingState = 1;
		initialEnvelope();
	}

	else if ( TloadingState == 1 ) // at initial envelope
	{
		if ( TmaxStrain > epsn && dStrain < 0 ) // reverse from tensile part after yield
		{
			reverseFromTenEnvelope();
            if ( Tstrain < approachToComEnvelopeStrain ) // directly reverse to compressive envelope
			{
				TloadingState = 3;
				compressionEnvelope();
			}
			else // unloading path, pathdown
			{
				TloadingState = 4;
				
				TreverseTopNum = 0;
				TreverseBottomNum = 0;

				TreverseTopStrain[TreverseTopNum] = reverseFromTenEnvelopeStrain;
				TreverseTopStress[TreverseTopNum] = reverseFromTenEnvelopeStress;
				TreverseBottomStrain[TreverseBottomNum] = approachToComEnvelopeStrain;
				TreverseBottomStress[TreverseBottomNum] = approachToComEnvelopeStress;

				determineDownPathPoint(); 
				downPath();
			}
		}
		else if ( TminStrain < -epsy && dStrain > 0 ) // reverse from compressive part after yield
		{   
			reverseFromComEnvelope();
            if ( Tstrain > approachToTenEnvelopeStrain ) // direct reverse to tensile envelope
			{
				TloadingState = 2;
				tensionEnvelope();
			}
			else // reloading path, pathup
			{
				TloadingState = 5;

				TreverseTopNum = 0;
				TreverseBottomNum = 0;
				
				TreverseTopStrain[TreverseTopNum] = approachToTenEnvelopeStrain;
				TreverseTopStress[TreverseTopNum] = approachToTenEnvelopeStress;
				TreverseBottomStrain[TreverseBottomNum] = reverseFromComEnvelopeStrain;
				TreverseBottomStress[TreverseBottomNum] = reverseFromComEnvelopeStress;                

				determineUpPathPoint();
				upPath();
			}
		}
		else // continue on initial envelope
		{
			initialEnvelope();
		}
	} // else if TloadingState == 1

	else if ( TloadingState == 2 )
	{
		if ( dStrain >= 0 ) // Continue on tension envelope
		{
			tensionEnvelope();
		}
		else // Reverse from tension envelope
		{
			reverseFromTenEnvelope();
			if ( Tstrain < approachToComEnvelopeStrain ) // direct reverse to compression envelope
			{
				TloadingState = 3;
				compressionEnvelope();
			}
			else // unloading path, pathdown
			{
				TloadingState = 4;

				TreverseTopNum = 0;
				TreverseBottomNum = 0;
				
				TreverseTopStrain[TreverseTopNum] = reverseFromTenEnvelopeStrain;
				TreverseTopStress[TreverseTopNum] = reverseFromTenEnvelopeStress;
				TreverseBottomStrain[TreverseBottomNum] = approachToComEnvelopeStrain;
				TreverseBottomStress[TreverseBottomNum] = approachToComEnvelopeStress;

				determineDownPathPoint(); 
				downPath();
			}
		}
	} // end of TloadingState = 2

	else if ( TloadingState == 3 )
	{
		if ( dStrain <= 0 ) // Continue on compression envelope
		{
			compressionEnvelope();
		}
		else // Reverse from compression envelope
		{
			reverseFromComEnvelope();
			if ( Tstrain > approachToTenEnvelopeStrain ) // direct reverse to tensile envelope
			{
				TloadingState = 2;
				tensionEnvelope();
			}
			else // reloading path, pathup
			{
				TloadingState = 5;
				
				TreverseTopNum = 0;
				TreverseBottomNum = 0;
				
				TreverseTopStrain[TreverseTopNum] = approachToTenEnvelopeStrain;
				TreverseTopStress[TreverseTopNum] = approachToTenEnvelopeStress;
				TreverseBottomStrain[TreverseBottomNum] = reverseFromComEnvelopeStrain;
				TreverseBottomStress[TreverseBottomNum] = reverseFromComEnvelopeStress;           

				determineUpPathPoint();
				upPath();
			}
		}
	} // end of TloadingState = 3


	else if ( TloadingState == 4 ) // In the loop reversed from tension envelope
	{
		if ( (Tstrain>reverseFromTenEnvelopeStrain)||(fabs(Tstrain-reverseFromTenEnvelopeStrain)<1e-6) ) // Back to tension envelope
		{
			reverseLoopSetZero( ); 
			TloadingState = 2;
			tensionEnvelope();
		}
		else if (( Tstrain<approachToComEnvelopeStrain)||(fabs(Tstrain-approachToComEnvelopeStrain)<1e-6) ) // Approach to compression envelope
		{
			reverseLoopSetZero( );
			TloadingState = 3;
			compressionEnvelope();
		}
		else
		{
            //opserr << " reverseTopStrain[reverseTopNum] = " << reverseTopStrain[reverseTopNum] << endln;
			// temp add
			//if ( Tstrain > TreverseTopStrain[TreverseTopNum] )
			//{
				//opserr << " TendonL01 : determinTrialState \n";
				//opserr << " Tstrain = " << Tstrain << endln;
				//opserr << " reverseTopNum = " << reverseTopNum << endln;
				//opserr << " reverseTopStrain[0] = " << reverseTopStrain[0] << endln;
				//opserr << " reverseTopStrain[reverseTopNum] = " << reverseTopStrain[reverseTopNum] << endln;
				//opserr << " reverseTopStress[reverseTopNum] = " << reverseTopStress[reverseTopNum] << endln;
				//opserr << " reverseFromTenEnvelopeStrain = " << reverseFromTenEnvelopeStrain << endln;
			//}
            // temp add finish

			determineTrialLoop(dStrain); // Still in TloadingState = 4
		}
	}

	else if ( TloadingState == 5 ) // In loop reverse from compression envelope
	{
		if ( Tstrain < reverseFromComEnvelopeStrain ) // Back to compression envelope
		{
			reverseLoopSetZero();
			TloadingState = 3;
			compressionEnvelope();
		}
		else if ( Tstrain > approachToTenEnvelopeStrain ) // Approach to tension envelope
		{
			reverseLoopSetZero();
			TloadingState = 2;
			tensionEnvelope();
		}
		else
		{
			determineTrialLoop(dStrain); // Still in TloadingState = 5
		}
	}

	else
	{
		opserr << "TendonL01::determineTrialState -- Improper TloadingState : " 
			 << TloadingState << " for TendonL01\n" ;
	}

	if ( TmaxStrain < Tstrain )		TmaxStrain = Tstrain;
	if ( TminStrain > Tstrain )		TminStrain = Tstrain;
} //end of determineTrialState( )



void TendonL01::initialEnvelope( )
{
	//Set parameter values
	//double epsy = fy/E0;
	//double fcr = 0.31*sqrt(fpc);
	//if ( rou < 0.0025 ) rou = 0.0025;
	//double B = pow((fcr/fy),1.5)/rou;
	//double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);
	double epsn = 0.7 * fpu / Eps;
	double epsy = fpy / Eps;
	double Epsr = 1.046 * Eps;
	double fpur = 0.963 * fpu;
	
	if ( Tstrain > epsn )
	{
		Tstress = Epsr/pow((1+pow(Epsr*Tstrain/fpur,5)),0.2)*Tstrain; 
		/*Ttangent = Tstress/Tstrain;
		if (Tstrain == 0)*/
		Ttangent = Epsr/pow((1+pow(Epsr*Tstrain/fpur,5)),1.2); 
	}
	else if ( Tstrain < - epsy )
	{
		Tstress = 0.001*Eps*(Tstrain+epsy) -fpy;
		/*Ttangent = Tstress/Tstrain;
		if (Tstrain == 0)*/
		Ttangent = 0.001*Eps;         // set Ttangent !=0 to avoid numerical problem
	}
	else
	{
		Tstress = Eps*Tstrain;
		/*Ttangent = Tstress/Tstrain;
		if (Tstrain == 0)*/
		Ttangent = Eps;
	}	
} 

void TendonL01::tensionEnvelope( )
{
	//double epsy = fy/E0;
	//double fcr = 0.31*sqrt(fpc);
	//if ( rou < 0.0025 ) rou = 0.0025;
	//double B = pow((fcr/fy),1.5)/rou;
	//double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);
	double epsn = 0.7 * fpu / Eps;
	double epsy = fpy / Eps;
	double Epsr = 1.046 * Eps;
	double fpur = 0.963 * fpu;

	if ( Tstrain > epsn )
	{
		Tstress = Epsr/pow((1+pow(Epsr*Tstrain/fpur,5)),0.2)*Tstrain; 
		/*Ttangent = Tstress/Tstrain;
		if (Tstrain == 0)*/
		Ttangent = Epsr/pow((1+pow(Epsr*Tstrain/fpur,5)),1.2);	
	}
	else
	{
		Tstress = 0.001*Eps*(Tstrain-epsn) +Eps*epsn;
		/*Ttangent = Tstress/Tstrain;
		if (Tstrain == 0)*/
		Ttangent = 0.001*Eps;          // set Ttangent !=0 to avoid numerical problem
	}
}

void TendonL01::compressionEnvelope( )
{
	double epsy = fpy/Eps;
	Tstress = 0.001*Eps*(Tstrain+epsy) -fpy;
	/*Ttangent = Tstress/Tstrain;
	if (Tstrain == 0)*/
	Ttangent = 0.001*Eps; 
}
	
void TendonL01::determineTrialLoop(double dStrain)
{

	//if ( Tstrain > TreverseTopStrain[TreverseTopNum] )
	//{		
		    //opserr << " reverseTopNum = " << reverseTopNum << endln;
			//opserr << " reverseBottomNum = " << reverseBottomNum << endln;
			//opserr << "\n TendonL01 : determinTrialLoop \n";
			//opserr << " Tstrain = " << Tstrain << endln;
			//opserr << " reverseTopStrain[reverseTopNum] = " << reverseTopStrain[reverseTopNum] << endln;
			//opserr << " reverseTopStress[reverseTopNum] = " << reverseTopStress[reverseTopNum] << endln;
			//opserr << " reverseFromTenEnvelopeStrain = " << reverseFromTenEnvelopeStrain << endln;
	//}
    // temp add finish

	if ( Tstrain > TreverseTopStrain[TreverseTopNum] ) // go upper, expand to outerloop
	{ 
        while ( Tstrain > TreverseTopStrain[TreverseTopNum] ) 
		{
			if ( TreverseTopNum > 0 )
			{
				TreverseTopStrain[TreverseTopNum] = 0.0;
				TreverseTopStress[TreverseTopNum] = 0.0;
				TreverseTopNum = TreverseTopNum -1;
			}
			if ( TreverseBottomNum > 0 )
			{
				TreverseBottomStrain[TreverseBottomNum] = 0.0;
				TreverseBottomStress[TreverseBottomNum] = 0.0;
				TreverseBottomNum = TreverseBottomNum -1;
			}
		}
		   
    	determineUpPathPoint();
		upPath();
	}

	else if ( Tstrain < TreverseBottomStrain[TreverseBottomNum] ) // go lower, expand to outerloop
	{		
        while ( Tstrain < TreverseBottomStrain[TreverseBottomNum] ) 
		{
          if ( TreverseTopNum > 0 )
		  {
			TreverseTopStrain[TreverseTopNum] = 0.0;
			TreverseTopStress[TreverseTopNum] = 0.0;
			TreverseTopNum = TreverseTopNum -1;
		  }
          if ( TreverseBottomNum > 0 )
		  {
			TreverseBottomStrain[TreverseBottomNum] = 0.0;
			TreverseBottomStress[TreverseBottomNum] = 0.0;
			TreverseBottomNum = TreverseBottomNum -1;
		  }
		}

		determineDownPathPoint();
		downPath();
	}
    
	else // create a new inner loop or continue on previous loop
	{

		 if (( TloopPathState == 2 || TloopPathState == 3 ) && ( dStrain > 0)) // reverse up at down path
		{
			TreverseBottomNum = TreverseBottomNum + 1;
			if ( TreverseBottomNum > (SIZE1-1) )
			{
				opserr << " TendonL01::determineTrialLoop -- overflowed the size of the array storing the loop!\n"
					 << " Size of the array : " << SIZE1 << "\n";
			}		
			TreverseBottomStrain[TreverseBottomNum] = Cstrain;
			TreverseBottomStress[TreverseBottomNum] = Cstress;
			determineUpPathPoint();
			upPath();
		}

		else if (( TloopPathState == 5 || TloopPathState == 6 ) && ( dStrain < 0 )) // reverse down at up path
		{
			TreverseTopNum = TreverseTopNum + 1;
			if ( TreverseTopNum > (SIZE1-1) )
			{
				opserr << " TendonL01::determineTrialLoop -- overflowed the size of the array storing the loop!\n"
					 << " Size of the array : " << SIZE1 << "\n";			
			}
			TreverseTopStrain[TreverseTopNum] = Cstrain;
			TreverseTopStress[TreverseTopNum] = Cstress;
			determineDownPathPoint();
			downPath();

			//opserr << " TloopPathState = " << TloopPathState << endln;
			//opserr << " Tstrain = " << Tstrain << endln;
			//opserr << " dStrain = " << dStrain << endln;
			//opserr << " TreverseTopStrain[TreverseTopNum] = " << TreverseTopStrain[TreverseTopNum] << endln;

		}

		else // continue on previous path of the loop
		{
			 if ( TloopPathState == 1 || TloopPathState == 2 || TloopPathState == 3 ) // downpath
			{
				determineDownPathPoint();
				downPath();
			}

			else if ( TloopPathState == 4 || TloopPathState == 5 || TloopPathState == 6 ) //up path
			{
				determineUpPathPoint();
				upPath();
			}
			else
			{
				  opserr << " TendonL01::determineTrialLoop -- improper TloopPathState : " 
					   << TloopPathState << "\n" ;
			}
		}

	}  // Inside loop


}



void TendonL01::determineDownPathPoint( )
{
		
	//double epsy = fy/E0;
	//double fcr = 0.31*sqrt(fpc);
	//if ( rou < 0.0025 ) rou = 0.0025;
	//double B = pow((fcr/fy),1.5)/rou;
	//double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);
	double epsn = 0.7 * fpu / Eps;
	double epsy = fpy / Eps;
	double Epsr = 1.046 * Eps;
	double fpur = 0.963 * fpu;

    double topStrain = TreverseTopStrain[TreverseTopNum];
	double topStress = TreverseTopStress[TreverseTopNum];

    double bottomStrain = TreverseBottomStrain[TreverseBottomNum];
	double bottomStress = TreverseBottomStress[TreverseBottomNum];

	double maxStrain;
	if ( fabs(topStrain) > fabs (bottomStrain) ) {
		maxStrain = topStrain;
	} else {
		maxStrain = bottomStrain;
	}

	double Kp;
	if ( ( maxStrain > epsn ) || ( maxStrain < 0 ))
	{
		Kp = fabs(( maxStrain - epsn )/ epsn);
	}
	else
	{
		Kp = fabs(( -maxStrain - epsn )/epsn);
	}
	double A = ac * pow(Kp, -0.1);
	double R = rc * pow(Kp, -0.2);

//c1

	// get downPathPointOneStress & downPathPointOneStrain
	downPathPointOneStress = 0.0;
	double tempOne1 = pow( fabs((downPathPointOneStress-topStress)/fpy), R-1);
	double tempOne2 = pow(A, -R)*tempOne1;
	downPathPointOneStrain = topStrain + (downPathPointOneStress-topStress)*(1+tempOne2)/Eps;	


	// get downPathPointTwoStress & downPathPointTwoStrain 
 
	if ( TreverseBottomStress[TreverseBottomNum] > -0.65*fpy )
	{
        downPathPointTwoStrain = TreverseBottomStrain[TreverseBottomNum];
		downPathPointTwoStress = TreverseBottomStress[TreverseBottomNum];
	}
	else
	{
		downPathPointTwoStress = -0.65*fpy;
		double tempTwo1 = pow( fabs((downPathPointTwoStress-topStress)/fpy), R-1);
		double tempTwo2 = pow(A, -R)*tempTwo1;
		downPathPointTwoStrain = topStrain + (downPathPointTwoStress-topStress)*(1+tempTwo2)/Eps;		

	}
}

void TendonL01::determineUpPathPoint( )
{
	//double epsy = fy/E0;
	//double fcr = 0.31*sqrt(fpc);
	//if ( rou < 0.0025 ) rou = 0.0025;
	//double B = pow((fcr/fy),1.5)/rou;
	//double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);
	double epsn = 0.7 * fpu / Eps;
	double epsy = fpy / Eps;
	double Epsr = 1.046 * Eps;
	double fpur = 0.963 * fpu;

	double topStrain = TreverseTopStrain[TreverseTopNum];
	double topStress = TreverseTopStress[TreverseTopNum];
    
    double bottomStrain = TreverseBottomStrain[TreverseBottomNum];
	double bottomStress = TreverseBottomStress[TreverseBottomNum];

	double maxStrain;
	if ( fabs(topStrain) > fabs (bottomStrain) ) {
		maxStrain = topStrain;
	} else {
		maxStrain = bottomStrain;
	}

	double Kp;
	if ( ( maxStrain > epsn ) || ( maxStrain < 0 ))
	{
		Kp = fabs(( maxStrain - epsn )/ epsn);
	}
	else
	{
		Kp = fabs(( -maxStrain - epsn )/epsn);
	}

	double A = ac * pow(Kp, -0.1);
	double R = rc * pow(Kp, -0.2);


    // get upPathPointOneStrain & upPathPointOneStress
   
	upPathPointOneStress = 0.0;
	double tempOne1 = pow( fabs((upPathPointOneStress-bottomStress)/fpy), R-1); 
	double tempOne2 = pow(A, -R)*tempOne1;
	upPathPointOneStrain = bottomStrain + (upPathPointOneStress-bottomStress)*(1+tempOne2)/Eps;

    // get upPathPointTwoStrain & upPathPointTwoStress	

	if ( TreverseTopStress[TreverseTopNum] < 0.65*fpy )
	{
        upPathPointTwoStrain = TreverseTopStrain[TreverseTopNum];
		upPathPointTwoStress = TreverseTopStress[TreverseTopNum];
	}
	else
	{
		upPathPointTwoStress = 0.65*fpy;
		double tempTwo1 = pow( fabs((upPathPointTwoStress-bottomStress)/fpy), R-1); 
		double tempTwo2 = pow(A, -R)*tempTwo1;
		upPathPointTwoStrain = bottomStrain + (upPathPointTwoStress-bottomStress)*(1+tempTwo2)/Eps;
	}
}

void TendonL01::downPath( )
{
	double topStrain = TreverseTopStrain[TreverseTopNum];
    double topStress = TreverseTopStress[TreverseTopNum];
	double bottomStrain = TreverseBottomStrain[TreverseBottomNum];
	double bottomStress = TreverseBottomStress[TreverseBottomNum];
	if ( Tstrain >= downPathPointOneStrain ) 
	{
// c3
		TloopPathState = 1;
		double slope1 = (topStress - downPathPointOneStress)/(topStrain-downPathPointOneStrain);
		Ttangent = slope1;
		Tstress = slope1*( Tstrain- topStrain) + topStress;		
	}
	else if ( Tstrain < downPathPointOneStrain && Tstrain >= downPathPointTwoStrain )
	{
		TloopPathState = 2;
		double slope2 = (downPathPointTwoStress - downPathPointOneStress)/(downPathPointTwoStrain-downPathPointOneStrain);
		Tstress = slope2*(Tstrain - downPathPointOneStrain) + downPathPointOneStress;
		Ttangent = slope2;
	}
	else
	{
		TloopPathState = 3;
		double slope3 = (bottomStress-downPathPointTwoStress)/(bottomStrain-downPathPointTwoStrain);
		Tstress = slope3*(Tstrain - downPathPointTwoStrain) +downPathPointTwoStress;
		Ttangent = slope3;
	}

}

void TendonL01::upPath( )
{
	double topStrain = TreverseTopStrain[TreverseTopNum];
    double topStress = TreverseTopStress[TreverseTopNum];
	double bottomStrain = TreverseBottomStrain[TreverseBottomNum];
	double bottomStress = TreverseBottomStress[TreverseBottomNum];
	double slope;
	if ( Tstrain <= upPathPointOneStrain ) 
	{
// c4
		TloopPathState = 4;
		slope = (bottomStress-upPathPointOneStress)/(bottomStrain-upPathPointOneStrain);
		Ttangent = slope;
		Tstress = slope*( Tstrain- bottomStrain) + bottomStress;
	}
	else if ( Tstrain > upPathPointOneStrain && Tstrain <= upPathPointTwoStrain )
	{
		TloopPathState = 5;
		slope = (upPathPointTwoStress-upPathPointOneStress)/(upPathPointTwoStrain-upPathPointOneStrain);
		Tstress = slope*(Tstrain - upPathPointOneStrain) + upPathPointOneStress;
		Ttangent = slope;
	}
	else
	{
		TloopPathState = 6;
		slope = (topStress-upPathPointTwoStress)/(topStrain-upPathPointTwoStrain);
		Tstress = slope*(Tstrain - upPathPointTwoStrain) + upPathPointTwoStress;
		Ttangent = slope;
	}
    
	tt2 = upPathPointTwoStrain; // temp add by zhong

}

void TendonL01::reverseFromTenEnvelope( )
{   
	// Get strain and stress of point reverse from tension envelope
	reverseFromTenEnvelopeStrain = Cstrain;
	reverseFromTenEnvelopeStress = Cstress;
 
	// Get strain and stress of point if approach to compression envelope
	//double epsy = fy/E0;
	//double fcr = 0.31*sqrt(fpc);
	//if ( rou < 0.0025 ) rou = 0.0025;
	//double B = pow((fcr/fy),1.5)/rou;
	//double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);
	double epsn = 0.7 * fpu / Eps;
	double epsy = fpy / Eps;
	double Epsr = 1.046 * Eps;
	double fpur = 0.963 * fpu;
     
	double Kp;
	if ( ( reverseFromTenEnvelopeStrain > epsn ) || ( reverseFromTenEnvelopeStrain < 0 ))
	{
		Kp = fabs(( reverseFromTenEnvelopeStrain - epsn )/ epsn);
	}
	else
	{
		Kp = fabs(( -reverseFromTenEnvelopeStrain - epsn )/epsn);
	}
   
   
	//double Kp = fabs(( TreverseFromTenEnvelopeStrain - epsn )/ epsn);
	double A = ac * pow(Kp, -0.1);
	double R = rc * pow(Kp, -0.2);
	double temp1 = pow( fabs((fpy+reverseFromTenEnvelopeStress)/fpy), R-1);
	double temp2 = pow(A, -R)*temp1;

	approachToComEnvelopeStrain = reverseFromTenEnvelopeStrain 
		                           + (-fpy-reverseFromTenEnvelopeStress)*(1+temp2)/Eps;
	// get approximate point of approachToComEnvelopeStress
	approachToComEnvelopeStress = 0.001*Eps*(approachToComEnvelopeStrain+epsy)-fpy;

}

void TendonL01::reverseFromComEnvelope( )
{
	// Get strain and stress of point reverse from compression envelope
	reverseFromComEnvelopeStrain = Cstrain;
	reverseFromComEnvelopeStress = Cstress;
 
	// Get strain and stress of point
	// Assume approaching to the part of tensile envelope when strain > epsn
	//double epsy = fy/E0;
	//double fcr = 0.31*sqrt(fpc);
	//if ( rou < 0.0025 ) rou = 0.0025;
	//double B = pow((fcr/fy),1.5)/rou;
	//double epsn = epsy*(0.91-2.0*B)/(0.98-0.25*B);
	//double fn = E0*epsn;
	//double Kp = fabs(( TreverseFromComEnvelopeStrain - epsn )/ epsn);
	double epsn = 0.7 * fpu / Eps;
	double epsy = fpy / Eps;
	double Epsr = 1.046 * Eps;
	double fpur = 0.963 * fpu;
	double fn = Eps * epsn;
	double Kp;
	if ( ( reverseFromComEnvelopeStrain > epsn ) || ( reverseFromComEnvelopeStrain < 0 ))
	{
		Kp = fabs(( reverseFromComEnvelopeStrain - epsn )/ epsn);
	}
	else
	{
		Kp = fabs(( -reverseFromComEnvelopeStrain - epsn )/epsn);
	}

	double A = ac * pow(Kp, -0.1);
	double R = rc * pow(Kp, -0.2);

	double tempStrainOne = reverseFromComEnvelopeStrain - reverseFromComEnvelopeStress/Eps;
	double temp1 = pow( fabs(( 0.65*fpy -reverseFromComEnvelopeStress)/fpy), R-1);
	double temp2 = pow(A, -R)*temp1;
	double tempStrainTwo = reverseFromComEnvelopeStrain 
		                           + (0.65*fpy-reverseFromComEnvelopeStress)*(1+temp2)/Eps;
	double slopeOneTwo = -0.65*fpy/(tempStrainOne-tempStrainTwo);
	
    //If intersection occurs at pathOneTwo and tensionEnvelope, calculate intersection stress
    double x, fx;
	x = tempStrainOne;
	fx = slopeOneTwo*(x-tempStrainOne) - Epsr/pow((1+pow(Epsr*x/fpur,5)),0.2)*x;

	while ( fabs(fx) > 1e-2 ) {
 	        x += 0.0001;
	        fx = slopeOneTwo*(x-tempStrainOne) - Epsr/pow((1+pow(Epsr*x/fpur,5)),0.2)*x;
		}
	double tempInterOneEnvelopeStrain = x;
    double tempInterOneEnvelopeStress = slopeOneTwo*(tempInterOneEnvelopeStrain-tempStrainOne);
	
   
   	if ( tempInterOneEnvelopeStress < 0.65*fpy ) // Intersection occurs at LoopPath5 and tensionEnvelope
	{
		approachToTenEnvelopeStrain = tempInterOneEnvelopeStrain;
        approachToTenEnvelopeStress = tempInterOneEnvelopeStress;
	}
	else // Intersection occurs at LoopPath6 and tensionEnvelope
	{
		double slope = 0.25*slopeOneTwo;
		x = tempStrainTwo;
        fx = slope*(x-tempStrainTwo) + 0.65*fpy - Epsr/pow((1+pow(Epsr*x/fpur,5)),0.2)*x;
		while ( fabs(fx) > 1e-2 ) {
 	        x += 0.0001;
	        fx = slope*(x-tempStrainTwo) + 0.65*fpy - Epsr/pow((1+pow(Epsr*x/fpur,5)),0.2)*x;
		}

		approachToTenEnvelopeStrain = x;
		approachToTenEnvelopeStress = slope*(approachToTenEnvelopeStrain-tempStrainTwo) + 0.65*fpy;
	}
    
    // Check if approach to the tension envelope when strain < epsn
	if ( approachToTenEnvelopeStrain < epsn )
	{
		temp1 = pow( fabs((fn-reverseFromComEnvelopeStress)/fpy), R-1);
		temp2 = pow(A, -R)*temp1;
		approachToTenEnvelopeStrain = reverseFromComEnvelopeStrain 
		                           + (fn-reverseFromComEnvelopeStress)*(1+temp2)/Eps;
		// get approximate point of approachToTenEnvelopeStress
		approachToTenEnvelopeStress = 0.001*Eps*(approachToTenEnvelopeStrain-epsn) + fn;
	}
}

void TendonL01::reverseLoopSetZero( )
{
	TloopPathState = 0;
	TreverseTopNum = 0;
	TreverseBottomNum = 0;
	for (int i=0; i<SIZE1; i++)
	{
		TreverseTopStrain[i] = 0.0;    
		TreverseTopStress[i] = 0.0;
		TreverseBottomStrain[i] = 0.0;
		TreverseBottomStress[i] = 0.0;
	}

}




UniaxialMaterial* TendonL01::getCopy ()
{
   TendonL01* theCopy = new TendonL01(this->getTag(), fpy, Eps, fpu, rou, epsp, ac, rc);     
   
   // History variables
   theCopy->CminStrain = CminStrain;  
   theCopy->CmaxStrain = CmaxStrain;   
   theCopy->CloadingState = CloadingState; 
   theCopy->CloopPathState = CloopPathState;
	
   theCopy->reverseFromTenEnvelopeStrain = reverseFromTenEnvelopeStrain;  
   theCopy->reverseFromTenEnvelopeStress = reverseFromTenEnvelopeStress;  
   theCopy->approachToComEnvelopeStrain = approachToComEnvelopeStrain;  
   theCopy->approachToComEnvelopeStress = approachToComEnvelopeStress;   
    
   theCopy->reverseFromComEnvelopeStrain = reverseFromComEnvelopeStrain;  
   theCopy->reverseFromComEnvelopeStress = reverseFromComEnvelopeStress;  
   theCopy->approachToTenEnvelopeStrain = approachToTenEnvelopeStrain;   
   theCopy->approachToTenEnvelopeStress = approachToTenEnvelopeStress;  
	
   int i;
   for (i=0; i<SIZE1; i++)
   {
		theCopy->CreverseTopStrain[i] = CreverseTopStrain[i];    
		theCopy->CreverseTopStress[i] = CreverseTopStress[i];
		theCopy->CreverseBottomStrain[i] = CreverseBottomStrain[i];
		theCopy->CreverseBottomStress[i] = CreverseBottomStress[i];
   }
   theCopy->CreverseTopNum = CreverseTopNum;             
   theCopy->CreverseBottomNum = CreverseBottomNum;          

   for ( i=0; i<SIZE1; i++)
   {
		theCopy->TreverseTopStrain[i] = TreverseTopStrain[i];    
		theCopy->TreverseTopStress[i] = TreverseTopStress[i];
		theCopy->TreverseBottomStrain[i] = TreverseBottomStrain[i];
		theCopy->TreverseBottomStress[i] = TreverseBottomStress[i];
   }
   theCopy->TreverseTopNum = TreverseTopNum;             
   theCopy->TreverseBottomNum = TreverseBottomNum;  


   theCopy->downPathPointOneStrain = downPathPointOneStrain;    
   theCopy->downPathPointTwoStrain = downPathPointTwoStrain;
   theCopy->downPathPointTwoStress = downPathPointTwoStress;
   theCopy->upPathPointOneStrain = upPathPointOneStrain;      
   theCopy->upPathPointTwoStrain = upPathPointTwoStrain;
   theCopy->upPathPointTwoStress = upPathPointTwoStress;

   // State variables
   theCopy->TminStrain = TminStrain;  
   theCopy->TmaxStrain = TmaxStrain;  
   theCopy->TloadingState = TloadingState;
   theCopy->TloopPathState = TloopPathState;

   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;

   return theCopy;
}


int TendonL01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(151);
   data(0) = this->getTag();

   // Material properties
   data(1) = fpy;
   data(2) = Eps;
   data(3) = fpu;
   data(4) = rou;
   data(5) = epsp;
   data(6) = ac;
   data(7) = rc;


   // History variables from last converged state
   data(8) = CminStrain;
   data(9) = CmaxStrain;
   data(10) = CloadingState;
   data(11) = CloopPathState;

   data(12) = reverseFromTenEnvelopeStrain;
   data(13) = reverseFromTenEnvelopeStress;
   data(14) = approachToComEnvelopeStrain;
   data(15) = approachToComEnvelopeStress;
   
   data(16) = reverseFromComEnvelopeStrain;
   data(17) = reverseFromComEnvelopeStress;
   data(18) = approachToTenEnvelopeStrain;
   data(19) = approachToTenEnvelopeStress; 

   // Here SIZE1 = 30;
   for (int i=0; i<SIZE1; i++)
	{
		data(20+i) = CreverseTopStrain[i];    
		data(50+i) = CreverseTopStress[i];
		data(80+i) = CreverseBottomStrain[i];
		data(110+i)= CreverseBottomStress[i];
	}

	data(140) = CreverseTopNum;
	data(141) = CreverseBottomNum;
	data(142) = downPathPointOneStrain;
	data(143) = downPathPointTwoStrain;
	data(144) = downPathPointTwoStress;
	data(145) = upPathPointOneStrain;
	data(146) = upPathPointTwoStrain;
	data(147) = upPathPointTwoStress;

    // State variables from last converged state
    data(148) = Cstrain;
    data(149) = Cstress;
    data(150) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "TendonL01::sendSelf() - failed to send data\n";

   return res;

} //sendSelf



int TendonL01::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(.65);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "TendonL01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      fpy = data(1);
      Eps = data(2);
      fpu = data(3);
      rou = data(4);
	  epsp = data(5);
      ac = data(6);
      rc = data(7);
      

      // History variables from last converged state

      CminStrain = data(8);
      CmaxStrain = data(9);
      CloadingState = int(data(10));
      CloopPathState = int(data(11));

      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
      TminStrain = CminStrain;  
      TmaxStrain = CmaxStrain;  
      TloadingState = CloadingState;
	  TloopPathState = CloopPathState;	

      reverseFromTenEnvelopeStrain = data(12);
      reverseFromTenEnvelopeStress = data(13);
      approachToComEnvelopeStrain = data(14);
      approachToComEnvelopeStress = data(15);
   
      reverseFromComEnvelopeStrain = data(16);
      reverseFromComEnvelopeStress = data(17);
      approachToTenEnvelopeStrain = data(18);
      approachToTenEnvelopeStress = data(19);
 

    	// Here SIZE1 = 30;
      for (int i=0; i<SIZE1; i++)
	  {
	  	CreverseTopStrain[i] = data(20+i);    
		CreverseTopStress[i] = data(50+i);
		CreverseBottomStrain[i] = data(80+i);
		CreverseBottomStress[i] = data(110+i);
	  }  

      CreverseTopNum = int(data(140));
	  CreverseBottomNum = int(data(141));
	  downPathPointOneStrain = data(142);
	  downPathPointTwoStrain = data(143);
	  downPathPointTwoStress = data(144);
	  upPathPointOneStrain = data(145);
	  upPathPointTwoStrain = data(146);
	  upPathPointTwoStress = data(147);

      // State variables from last converged state
      Cstrain = data(148);
      Cstress = data(149);
      Ctangent = data(150);      

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }
    
   return res;
}


Response* 
TendonL01::setResponse(const char **argv, int argc,
			 OPS_Stream &theOutput)
{
  Response *theResponse = 0;

  if (strcmp(argv[0],"getCommittedStrain") == 0) {
    double data = 0.0;
    theResponse = new MaterialResponse(this, 100, data);
  } else
    return this->UniaxialMaterial::setResponse(argv, argc, theOutput);

  return theResponse;
}
 
int 
TendonL01::getResponse(int responseID, Information &matInfo)
{
  if (responseID == 100) {
    matInfo.theDouble = this->getCommittedStrain();
  } else
    return this->UniaxialMaterial::getResponse(responseID, matInfo);

  return 0;
}

void TendonL01::Print (OPS_Stream& s, int flag)
{
   s << "TendonL01 tag: " << this->getTag() << endln;
  // comment for temporary
  // s << "  fy: " << fy << endln;
  // s << "  E0: " << E0 << endln;
  // s << "  fpc:  " << fpc << endln;
  // s << "  rou: " << rou << endln;
  //  s << "  ac: " << ac << endln;
  // s << "  rc: " << rc << endln; 

   // add for check purpose
   s<< " Strain = " << this->getStrain() << endln;
   s<< " Stress = " << this->getStress() << endln;
   s<< " Tangent = " << this->getTangent() << endln;
   s<< " LoadingState = " <<TloadingState << endln;
   s<< " LoopPathState = " << TloopPathState << endln;
   //s<< " Tstrain = " << tt1 << endln;
   //s<< " dstrain = " << ttStrain << endln;
   //s<< " approachToComEnvelopeStrain = "<<approachToComEnvelopeStrain << endln;
   //s<< " reverseFromTenEnvelopeStrain = "<<reverseFromTenEnvelopeStrain << endln;
   //s<< " tt2 = " << tt2 << endln;
   //s<< " ttStrain = " << ttStrain << endln;
   //s<< " TreverseTopNum = " << TreverseTopNum << endln;
   //s<< " TreverseTopStrain[TreverseTopNum] = " << TreverseTopStrain[TreverseTopNum] << endln;
   //s<< " TreverseTopStrain[1] = " << TreverseTopStrain[1] << endln;
}





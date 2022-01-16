// SteelZ01.cpp
// Written: JZhong
// Created: 2003.7  Revised 2005.4
//
// Description: This file contains the class implementation for
// SteelZ01
// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.

#include "SteelZ01.h"
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include <MaterialResponse.h>
#include <elementAPI.h>
#define OPS_Export 

OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_SteelZ01Material)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs < 5) {
    opserr << "Invalid Args want: uniaxialMaterial SteelZ01 tag? fy? E0? fpc? rou? <ac?> <rc?>" << endln;
    return 0;	
  }

  int    iData[1];
  double dData[6];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial SteelZ01 tag" << endln;
    return 0;
  }

  numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 4) {
    if (OPS_GetDoubleInput(&numRemainingArgs, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial SteelZ01 tag? fy? E0? fpc? rou? <ac?> <rc?>" << endln;
      return 0;	
    } else
      theMaterial = new SteelZ01(iData[0], dData[0], dData[1], dData[2], dData[3]);
  } else if (numRemainingArgs == 6) {
    if (OPS_GetDoubleInput(&numRemainingArgs, dData) != 0) {
      opserr << "Invalid Args want: uniaxialMaterial SteelZ01 tag? fy? E0? fpc? rou? <ac?> <rc?>" << endln;
      return 0;	
    } else
      theMaterial = new SteelZ01(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);
  } else {

    return 0;
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type SteelZ01\n";
    return 0;
  }

  return theMaterial;
}


SteelZ01::SteelZ01
(int tag, double FY, double E, double FPC, double ROU, double AC, double RC):
UniaxialMaterial(tag, MAT_TAG_SteelZ01),
fy(FY), E0(E), fpc(FPC), rou(ROU), ac(AC), rc(RC)
{
  tt1 = 0.0; // temp add
  tt2 = 0.0; // temp add by zhong
  ttStrain = 0.0;
  
  if ( fpc < 0.0 )
    {
      fpc = - fpc; // Make fpc positive to get sqrt(fpc)
    }
  // Sets all history and state variables to initial values
  this->revertToStart();
  
}

SteelZ01::SteelZ01(): UniaxialMaterial(0, MAT_TAG_SteelZ01),
fy(0.0), E0(0.0), fpc(0.0), rou(0.0), ac(0.0), rc(0.0)
{
	this->revertToStart();
}

SteelZ01::~SteelZ01()
{
	// does nothing;
}

int SteelZ01::revertToStart()
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
	for (i=0; i<SIZE; i++)
	{
		CreverseTopStrain[i] = 0.0;    
		CreverseTopStress[i] = 0.0;
		CreverseBottomStrain[i] = 0.0;
		CreverseBottomStress[i] = 0.0;
	}
	CreverseTopNum = 0;             
	CreverseBottomNum = 0;   

	for (i=0; i<SIZE; i++)
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
    Ctangent = E0;

    Tstrain = 0.0;
    Tstress = 0.0;
	Ttangent = E0;

    return 0;
}


int SteelZ01::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;  
   TmaxStrain = CmaxStrain;  
   TloadingState = CloadingState;
   TloopPathState = CloopPathState;

   for (int i=0; i<SIZE; i++)
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




int SteelZ01::commitState ()
{
   // Reset trial history variables to last committed state
    CminStrain = TminStrain;  
    CmaxStrain = TmaxStrain;  
    CloadingState = TloadingState;
	CloopPathState = TloopPathState;

   for (int i=0; i<SIZE; i++)
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


double SteelZ01::getCommittedStrain()
{
   return Cstrain;
}


double SteelZ01::getStrain ()
{
   return Tstrain;
}

double SteelZ01::getStress ()
{
   return Tstress;
}

double SteelZ01::getTangent ()
{
   // temp add
   if ( Ttangent == 0.0 )
   {
   	   opserr << " SteelZ01:getTangent() -- Ttangent = 0.0\n";
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

double SteelZ01::getSecant ()
{
	double TSecant;
	if ( Tstrain == 0.0 )
	{
		TSecant = E0;
	}
    else
	{
		TSecant = Tstress/Tstrain;
	}
	return TSecant;
}


int SteelZ01::setTrialStrain (double strain, double strainRate)
{
   // Reset history variables to last converged state	
    TminStrain = CminStrain;  
    TmaxStrain = CmaxStrain;  
    TloadingState = CloadingState;
	TloopPathState = CloopPathState;

   for (int i=0; i<SIZE; i++)
   {
    	TreverseTopStrain[i] = CreverseTopStrain[i];    
		TreverseTopStress[i] = CreverseTopStress[i];
		TreverseBottomStrain[i] = CreverseBottomStrain[i];
		TreverseBottomStress[i] = CreverseBottomStress[i];
   }
   TreverseTopNum = CreverseTopNum;             
   TreverseBottomNum = CreverseBottomNum; 

	// Set trial strain
   Tstrain = strain;

   // Determine change in strain from last converged state
   double dStrain = Tstrain - Cstrain;

   // Calculate the trial state given the trial strain
   //   if (fabs(dStrain) > DBL_EPSILON)
   if (fabs(dStrain) > 1e-10)
   determineTrialState (dStrain);

   ttStrain = dStrain;

   return 0;
}



int SteelZ01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
	TminStrain = CminStrain;  
    TmaxStrain = CmaxStrain;  
    TloadingState = CloadingState;
	TloopPathState = CloopPathState;

   for (int i=0; i<SIZE; i++)
   {
    	TreverseTopStrain[i] = CreverseTopStrain[i];    
		TreverseTopStress[i] = CreverseTopStress[i];
		TreverseBottomStrain[i] = CreverseBottomStrain[i];
		TreverseBottomStress[i] = CreverseBottomStress[i];
   }
   TreverseTopNum = CreverseTopNum;             
   TreverseBottomNum = CreverseBottomNum; 

   // Set trial strain
   Tstrain = strain;

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


void SteelZ01::determineTrialState (double dStrain)
{
    tt1 = Tstrain;
	tt2 = 0.0;
	ttStrain = dStrain;

	//Set parameter values
	double epsy = fy/E0;
	double fcr = 0.31*sqrt(fpc);
	if ( rou < 0.0025 ) rou = 0.0025;
	double B = pow((fcr/fy),1.5)/rou;
	double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);

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
				//opserr << " steelZ01 : determinTrialState \n";
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
		opserr << "SteelZ01::determineTrialState -- Improper TloadingState : " 
			 << TloadingState << " for SteelZ01\n" ;
	}

	if ( TmaxStrain < Tstrain )		TmaxStrain = Tstrain;
	if ( TminStrain > Tstrain )		TminStrain = Tstrain;
} //end of determineTrialState( )



void SteelZ01::initialEnvelope( )
{
	//Set parameter values
	double epsy = fy/E0;
	double fcr = 0.31*sqrt(fpc);
	if ( rou < 0.0025 ) rou = 0.0025;
	double B = pow((fcr/fy),1.5)/rou;
	double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);
	
	if ( Tstrain > epsn )
	{
		Tstress = fy * (0.91-2.0*B) + E0*(0.25*B+0.02)*Tstrain; //corrected
		/*Ttangent = Tstress/Tstrain;
		if (Tstrain == 0) */
		Ttangent = E0*(0.25*B+0.02); //corrected
	}
	else if ( Tstrain < - epsy )
	{
		Tstress = 0.001*E0*(Tstrain+epsy) -fy;
		/*Ttangent = Tstress/Tstrain;
		if (Tstrain == 0)*/ 
		Ttangent = 0.001*E0;         // set Ttangent !=0 to avoid numerical problem
	}
	else
	{
		Tstress = E0*Tstrain;
		/*Ttangent = Tstress/Tstrain;
		if (Tstrain == 0) */
		Ttangent = E0;
	}	
} 

void SteelZ01::tensionEnvelope( )
{
	double epsy = fy/E0;
	double fcr = 0.31*sqrt(fpc);
	if ( rou < 0.0025 ) rou = 0.0025;
	double B = pow((fcr/fy),1.5)/rou;
	double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);

	if ( Tstrain > epsn )
	{
		Tstress = fy * (0.91-2.0*B) + E0*(0.25*B+0.02)*Tstrain; //corrected
		/*Ttangent = Tstress/Tstrain;
		if (Tstrain == 0) */
		Ttangent = (0.25*B+0.02)*E0; //corrected	
	}
	else
	{
		Tstress = 0.001*E0*(Tstrain-epsn) +E0*epsn;
		/*Ttangent = Tstress/Tstrain;
		if (Tstrain == 0)*/
		Ttangent = 0.001*E0;          // set Ttangent !=0 to avoid numerical problem
	}
}

void SteelZ01::compressionEnvelope( )
{
	double epsy = fy/E0;
	Tstress = 0.001*E0*(Tstrain+epsy) -fy;
	/*Ttangent = Tstress/Tstrain;
	if (Tstrain == 0)*/
	Ttangent = 0.001*E0; 
}
	
void SteelZ01::determineTrialLoop(double dStrain)
{

	//if ( Tstrain > TreverseTopStrain[TreverseTopNum] )
	//{		
		    //opserr << " reverseTopNum = " << reverseTopNum << endln;
			//opserr << " reverseBottomNum = " << reverseBottomNum << endln;
			//opserr << "\n steelZ01 : determinTrialLoop \n";
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
			if ( TreverseBottomNum > (SIZE-1) )
			{
				opserr << " SteelZ01::determineTrialLoop -- overflowed the size of the array storing the loop!\n"
					 << " Size of the array : " << SIZE << "\n";
			}		
			TreverseBottomStrain[TreverseBottomNum] = Cstrain;
			TreverseBottomStress[TreverseBottomNum] = Cstress;
			determineUpPathPoint();
			upPath();
		}

		else if (( TloopPathState == 5 || TloopPathState == 6 ) && ( dStrain < 0 )) // reverse down at up path
		{
			TreverseTopNum = TreverseTopNum + 1;
			if ( TreverseTopNum > (SIZE-1) )
			{
				opserr << " SteelZ01::determineTrialLoop -- overflowed the size of the array storing the loop!\n"
					 << " Size of the array : " << SIZE << "\n";			
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
				  opserr << " SteelZ01::determineTrialLoop -- improper TloopPathState : " 
					   << TloopPathState << "\n" ;
			}
		}

	}  // Inside loop


}



void SteelZ01::determineDownPathPoint( )
{
		
	double epsy = fy/E0;
	double fcr = 0.31*sqrt(fpc);
	if ( rou < 0.0025 ) rou = 0.0025;
	double B = pow((fcr/fy),1.5)/rou;
	double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);
    
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
	double tempOne1 = pow( fabs((downPathPointOneStress-topStress)/fy), R-1);
	double tempOne2 = pow(A, -R)*tempOne1;
	downPathPointOneStrain = topStrain + (downPathPointOneStress-topStress)*(1+tempOne2)/E0;	


	// get downPathPointTwoStress & downPathPointTwoStrain 
 
	if ( TreverseBottomStress[TreverseBottomNum] > -0.65*fy )
	{
        downPathPointTwoStrain = TreverseBottomStrain[TreverseBottomNum];
		downPathPointTwoStress = TreverseBottomStress[TreverseBottomNum];
	}
	else
	{
		downPathPointTwoStress = -0.65*fy;
		double tempTwo1 = pow( fabs((downPathPointTwoStress-topStress)/fy), R-1);
		double tempTwo2 = pow(A, -R)*tempTwo1;
		downPathPointTwoStrain = topStrain + (downPathPointTwoStress-topStress)*(1+tempTwo2)/E0;		

	}
}

void SteelZ01::determineUpPathPoint( )
{
	double epsy = fy/E0;
	double fcr = 0.31*sqrt(fpc);
	if ( rou < 0.0025 ) rou = 0.0025;
	double B = pow((fcr/fy),1.5)/rou;
	double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);

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
	double tempOne1 = pow( fabs((upPathPointOneStress-bottomStress)/fy), R-1); 
	double tempOne2 = pow(A, -R)*tempOne1;
	upPathPointOneStrain = bottomStrain + (upPathPointOneStress-bottomStress)*(1+tempOne2)/E0;

    // get upPathPointTwoStrain & upPathPointTwoStress	

	if ( TreverseTopStress[TreverseTopNum] < 0.65*fy )
	{
        upPathPointTwoStrain = TreverseTopStrain[TreverseTopNum];
		upPathPointTwoStress = TreverseTopStress[TreverseTopNum];
	}
	else
	{
		upPathPointTwoStress = 0.65*fy;
		double tempTwo1 = pow( fabs((upPathPointTwoStress-bottomStress)/fy), R-1); 
		double tempTwo2 = pow(A, -R)*tempTwo1;
		upPathPointTwoStrain = bottomStrain + (upPathPointTwoStress-bottomStress)*(1+tempTwo2)/E0;
	}
}

void SteelZ01::downPath( )
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

void SteelZ01::upPath( )
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

void SteelZ01::reverseFromTenEnvelope( )
{   
	// Get strain and stress of point reverse from tension envelope
	reverseFromTenEnvelopeStrain = Cstrain;
	reverseFromTenEnvelopeStress = Cstress;
 
	// Get strain and stress of point if approach to compression envelope
	double epsy = fy/E0;
	double fcr = 0.31*sqrt(fpc);
	if ( rou < 0.0025 ) rou = 0.0025;
	double B = pow((fcr/fy),1.5)/rou;
	double epsn = epsy * (0.91-2.0*B) / (0.98-0.25*B);
     
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
	double temp1 = pow( fabs((fy+reverseFromTenEnvelopeStress)/fy), R-1);
	double temp2 = pow(A, -R)*temp1;

	approachToComEnvelopeStrain = reverseFromTenEnvelopeStrain 
		                           + (-fy-reverseFromTenEnvelopeStress)*(1+temp2)/E0;
	// get approximate point of approachToComEnvelopeStress
	approachToComEnvelopeStress = 0.001*E0*(approachToComEnvelopeStrain+epsy)-fy;

}

void SteelZ01::reverseFromComEnvelope( )
{
	// Get strain and stress of point reverse from compression envelope
	reverseFromComEnvelopeStrain = Cstrain;
	reverseFromComEnvelopeStress = Cstress;
 
	// Get strain and stress of point
	// Assume approaching to the part of tensile envelope when strain > epsn
	double epsy = fy/E0;
	double fcr = 0.31*sqrt(fpc);
	if ( rou < 0.0025 ) rou = 0.0025;
	double B = pow((fcr/fy),1.5)/rou;
	double epsn = epsy*(0.91-2.0*B)/(0.98-0.25*B);
	double fn = E0*epsn;
	//double Kp = fabs(( TreverseFromComEnvelopeStrain - epsn )/ epsn);
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

	double tempStrainOne = reverseFromComEnvelopeStrain - reverseFromComEnvelopeStress/E0;
	double temp1 = pow( fabs(( 0.65*fy -reverseFromComEnvelopeStress)/fy), R-1);
	double temp2 = pow(A, -R)*temp1;
	double tempStrainTwo = reverseFromComEnvelopeStrain 
		                           + (0.65*fy-reverseFromComEnvelopeStress)*(1+temp2)/E0;
	double slopeOneTwo = -0.65*fy/(tempStrainOne-tempStrainTwo);
	//If intersection occurs at pathOneTwo and tensionEnvelope, calculate intersection stress
	double tempInterOneEnvelopeStrain = (slopeOneTwo*tempStrainOne + fy*(0.91-2.0*B))/(slopeOneTwo-0.25*B*E0-0.02*E0);//corrected
	double tempInterOneEnvelopeStress = slopeOneTwo*(tempInterOneEnvelopeStrain-tempStrainOne);
   
   	if ( tempInterOneEnvelopeStress < 0.65*fy ) // Intersection occurs at LoopPath5 and tensionEnvelope
	{
		approachToTenEnvelopeStrain = tempInterOneEnvelopeStrain;
        approachToTenEnvelopeStress = tempInterOneEnvelopeStress;
	}
	else // Intersection occurs at LoopPath6 and tensionEnvelope
	{
		double slope = 0.25*slopeOneTwo;
		approachToTenEnvelopeStrain = (slope*tempStrainTwo + fy*(0.26-2.0*B))/(slope-0.25*B*E0-0.02*E0);//corrected
		approachToTenEnvelopeStress = slope*(approachToTenEnvelopeStrain-tempStrainTwo) + 0.65*fy;
	}
    
    // Check if approach to the tension envelope when strain < epsn
	if ( approachToTenEnvelopeStrain < epsn )
	{
		temp1 = pow( fabs((fn-reverseFromComEnvelopeStress)/fy), R-1);
		temp2 = pow(A, -R)*temp1;
		approachToTenEnvelopeStrain = reverseFromComEnvelopeStrain 
		                           + (fn-reverseFromComEnvelopeStress)*(1+temp2)/E0;
		// get approximate point of approachToTenEnvelopeStress
		approachToTenEnvelopeStress = 0.001*E0*(approachToTenEnvelopeStrain-epsn) + fn;
	}
}

void SteelZ01::reverseLoopSetZero( )
{
	TloopPathState = 0;
	TreverseTopNum = 0;
	TreverseBottomNum = 0;
	for (int i=0; i<SIZE; i++)
	{
		TreverseTopStrain[i] = 0.0;    
		TreverseTopStress[i] = 0.0;
		TreverseBottomStrain[i] = 0.0;
		TreverseBottomStress[i] = 0.0;
	}

}




UniaxialMaterial* SteelZ01::getCopy ()
{
   SteelZ01* theCopy = new SteelZ01(this->getTag(), fy, E0, fpc, rou, ac,rc);     
   
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
   for (i=0; i<SIZE; i++)
   {
		theCopy->CreverseTopStrain[i] = CreverseTopStrain[i];    
		theCopy->CreverseTopStress[i] = CreverseTopStress[i];
		theCopy->CreverseBottomStrain[i] = CreverseBottomStrain[i];
		theCopy->CreverseBottomStress[i] = CreverseBottomStress[i];
   }
   theCopy->CreverseTopNum = CreverseTopNum;             
   theCopy->CreverseBottomNum = CreverseBottomNum;          

   for ( i=0; i<SIZE; i++)
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


int SteelZ01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(149);
   data(0) = this->getTag();

   // Material properties
   data(1) = fy;
   data(2) = E0;
   data(3) = fpc;
   data(4) = rou;
   data(5) = ac;
   data(6) = rc;


   // History variables from last converged state
   data(7) = CminStrain;
   data(8) = CmaxStrain;
   data(9) = CloadingState;
   data(10) = CloopPathState;

   data(11) = reverseFromTenEnvelopeStrain;
   data(12) = reverseFromTenEnvelopeStress;
   data(13) = approachToComEnvelopeStrain;
   data(14) = approachToComEnvelopeStress;
   
   data(15) = reverseFromComEnvelopeStrain;
   data(16) = reverseFromComEnvelopeStress;
   data(17) = approachToTenEnvelopeStrain;
   data(18) = approachToTenEnvelopeStress; 

   // Here SIZE = 30;
   for (int i=0; i<SIZE; i++)
	{
		data(19+i) = CreverseTopStrain[i];    
		data(49+i) = CreverseTopStress[i];
		data(79+i) = CreverseBottomStrain[i];
		data(109+i)= CreverseBottomStress[i];
	}

	data(139) = CreverseTopNum;
	data(140) = CreverseBottomNum;
	data(141) = downPathPointOneStrain;
	data(142) = downPathPointTwoStrain;
	data(143) = downPathPointTwoStress;
	data(144) = upPathPointOneStrain;
	data(145) = upPathPointTwoStrain;
	data(146) = upPathPointTwoStress;

    // State variables from last converged state
    data(147) = Cstrain;
    data(148) = Cstress;
    data(149) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "SteelZ01::sendSelf() - failed to send data\n";

   return res;

} //sendSelf



int SteelZ01::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(149);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "SteelZ01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      fy = data(1);
      E0 = data(2);
      fpc = data(3);
      rou = data(4);
      ac = data(5);
      rc = data(6);
      

      // History variables from last converged state

      CminStrain = data(7);
      CmaxStrain = data(8);
      CloadingState = int(data(9));
      CloopPathState = int(data(10));

      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
      TminStrain = CminStrain;  
      TmaxStrain = CmaxStrain;  
      TloadingState = CloadingState;
	  TloopPathState = CloopPathState;	

      reverseFromTenEnvelopeStrain = data(11);
      reverseFromTenEnvelopeStress = data(12);
      approachToComEnvelopeStrain = data(13);
      approachToComEnvelopeStress = data(14);
   
      reverseFromComEnvelopeStrain = data(15);
      reverseFromComEnvelopeStress = data(16);
      approachToTenEnvelopeStrain = data(17);
      approachToTenEnvelopeStress = data(18);
 

    	// Here SIZE = 30;
      for (int i=0; i<SIZE; i++)
	  {
	  	CreverseTopStrain[i] = data(19+i);    
		CreverseTopStress[i] = data(49+i);
		CreverseBottomStrain[i] = data(79+i);
		CreverseBottomStress[i] = data(109+i);
	  }  

      CreverseTopNum = int(data(139));
	  CreverseBottomNum = int(data(140));
	  downPathPointOneStrain = data(141);
	  downPathPointTwoStrain = data(142);
	  downPathPointTwoStress = data(143);
	  upPathPointOneStrain = data(144);
	  upPathPointTwoStrain = data(145);
	  upPathPointTwoStress = data(146);

      // State variables from last converged state
      Cstrain = data(147);
      Cstress = data(148);
      Ctangent = data(149);      

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }
    
   return res;
}


Response* 
SteelZ01::setResponse(const char **argv, int argc,
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
SteelZ01::getResponse(int responseID, Information &matInfo)
{
  if (responseID == 100) {
    matInfo.theDouble = this->getCommittedStrain();
  } else
    return this->UniaxialMaterial::getResponse(responseID, matInfo);

  return 0;
}

void SteelZ01::Print (OPS_Stream& s, int flag)
{
   s << "SteelZ01 tag: " << this->getTag() << endln;
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





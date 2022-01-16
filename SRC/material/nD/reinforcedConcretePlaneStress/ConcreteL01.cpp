// File: ConcreteL01.cpp
// Written: ALaskar
// Created: 2007.8
// Description: This file contains the class implementation for
// ConcreteL01
// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.


#include "ConcreteL01.h"
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

OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_ConcreteL01Material)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[2];
  int numData = 1;


  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs < 3) {
  opserr << "Want: uniaxialMaterial ConcreteL01 tag? fpc? epsc0?" << endln;
    return 0;	
  }


  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConcreteL01 tag" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial ConcreteL01 tag? fpc? epsc0?" << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  theMaterial = new ConcreteL01(iData[0], dData[0], dData[1]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ConcreteL01\n";
    return 0;
  }

  return theMaterial;
}

ConcreteL01::ConcreteL01
(int tag, double FPC, double EPSC0): 
UniaxialMaterial(tag, MAT_TAG_ConcreteL01), fpc(FPC), epsc0(EPSC0)
{
	// Make all concrete parameters negative
	if (fpc > 0.0)
		fpc = -fpc;

	if (epsc0 > 0.0)
		epsc0 = -epsc0;
	
	// Set trial values
	this->revertToStart();
 }

ConcreteL01::ConcreteL01():
UniaxialMaterial(0, MAT_TAG_ConcreteL01), fpc(0.0), epsc0(0.0)
{
	// Set trial values
	this->revertToStart();
 }

ConcreteL01::~ConcreteL01 ()
{
   // Does nothing
}

//int ConcreteL01::setTrialStrain (double x, double k, double Dfactor, double BETA, double EPSLONTP, double strain, double strainRate)
int ConcreteL01::setTrialStrain (double strain, double strainRate)
{
  /* FMK
   X = x;
   K = k;
   D = Dfactor;
   beta = BETA;
   epslonTP = EPSLONTP;
  */
   
   /* Make all prestressing stress and strain negative
	if (sigmaCI > 0.0)
		sigmaCI = -sigmaCI;

	if (epslonCI > 0.0)
		epslonCI = -epslonCI;
	*/
   
   fbeta = 1 - fabs(beta)/24;
   Wp = 1.15 + fabs(beta) * (0.09 * fabs(beta) - 1) / 6;
   
   //Calculate softerning effect: zeta if epslonTP >0 
   if (epslonTP > 0.0)
	{
	    // add K into zeta, K is delta
		zeta =  (K) * 5.8 / sqrt( -fpc * ( 1.0 + 400.0 * epslonTP) ) * fbeta * Wp;
		if ( zeta >= 0.9 )
		{ 
			zeta = 0.9;
		}
		if ( zeta <= 0.25 ) //min zeta
		{
			zeta = 0.25;
		}
	}
   else
	{
		zeta = 1.0;
    }
 
   // Reset history variables to last converged state
   TloadingState = CloadingState;  

   // Set trial strain
   Tstrain = strain;

   // Determine change in strain from last converged state
   double dStrain = Tstrain - Cstrain;

   // Calculate the trial state given the trial strain
   //   if (fabs(dStrain) > DBL_EPSILON)   
   determineTrialState (dStrain); 
   
   return 0;
}


//int ConcreteL01::setTrial (double x, double k, double Dfactor, double BETA, double EPSLONTP, double strain, double &stress, double &tangent, double strainRate)
int ConcreteL01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
  /* FMK
   X = x;
   K = k;
   D = Dfactor;
   beta = BETA;
   epslonTP = EPSLONTP;
  */


   /* Make all prestressing stress and strain negative
	if (sigmaCI > 0.0)
		sigmaCI = -sigmaCI;

	if (epslonCI > 0.0)
		epslonCI = -epslonCI;
		*/

   //Calculate softerning effect: zeta if epslonTP >0 
   fbeta = 1 - fabs(beta)/24;
   Wp = 1.15 + fabs(beta) * (0.09 * fabs(beta) - 1) / 6;

   if (epslonTP > 0.0)
	{
	    // add K into zeta
		zeta =  (K) * 5.8 / sqrt( -fpc * ( 1.0 + 400.0 * epslonTP) ) * fbeta * Wp;
		if ( zeta >= 0.9 ) // max zeta
		{ 
			zeta = 0.9;
		}
		if ( zeta <= 0.25 )
		{
			zeta = 0.25; // min zeta
		}
	}
   else
	{
		zeta = 1.0;
    }

   // Reset history variables to last converged state
   TloadingState = CloadingState;

   // Set trial strain
   Tstrain = strain;

   // Determine change in strain from last converged state
   double dStrain = Tstrain - Cstrain;


   // Calculate the trial state given the trial strain
   // if (fabs(dStrain) > DBL_EPSILON) 
   determineTrialState (dStrain); 

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}


double ConcreteL01::getStrain ()
{
   return Tstrain;
}

double ConcreteL01::getStress ()
{
   return Tstress;
}

double ConcreteL01::getTangent ()
{
   return Ttangent;
}

double ConcreteL01::getSecant ()
{
	double TSecant;
	if ( Tstrain == 0.0 )
	{
		TSecant = 1.7*fpc/epsc0;
	}
    else
	{
		TSecant = Tstress/Tstrain;
	}
	return TSecant;
}

double ConcreteL01::getZeta ()
{
	return zeta;
}

int ConcreteL01::commitState ()
{
   // History variables
   CloadingState = TloadingState;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int ConcreteL01::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TloadingState = CloadingState;

   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}


int ConcreteL01::revertToStart ()
{
   // History variables
   
   TloadingState = 0;
   CloadingState = 0;
   
   reloadPath = 0;

   reverseFromOneStrain = 0.0;
   reverseFromOneStress = 0.0;
   reverseFromTwoStrain = 0.0; 
   reverseFromTwoStress = 0.0;
   reverseFromFourStrain = 0.0; 
   reverseFromFourStress = 0.0;   

   interFiveSevenStrain = 0.0;

   approachFiveToComStrain = 0.0;
   approachSixToComStrain = 0.0;
 
   // State variables
   zeta = 1.0;
   beta = 1.0;
   epslonTP = 0.0;
   double Ec0 = 1.4*fpc/epsc0;
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Ec0;

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = Ec0;

   return 0;
}

void ConcreteL01::determineTrialState (double dStrain)
{
	double Ec0 = 1.4*fpc/epsc0;
	if ( TloadingState == 0 )
	{
		envelope( );		
	}
	
	if ( TloadingState == 1 )
	{
		if ( dStrain < 0 ) //Continues on envelope
			envelope( );
		else //Reverse occurs
		{
			reverseFromOneStrain = Cstrain;
			reverseFromOneStress = Cstress;

			TloadingState = 5;
			reloadPath = 1;
        
        	interFiveSevenStrain = reverseFromOneStrain - reverseFromOneStress/Ec0;
			getApproachFiveToComStrain ( );
			pathFive();

			if ( Tstress > 0 ) //Reach path 7
			{
				TloadingState = 7;
				pathSeven();
			}
		}
	}


	else if ( TloadingState == 2 )
	{
		if ( dStrain < 0 ) 
			envelope( );
		else
		{
			reverseFromTwoStrain = Cstrain;
			reverseFromTwoStress = Cstress;

			TloadingState = 5;
			reloadPath = 2;

			interFiveSevenStrain = reverseFromTwoStrain - reverseFromTwoStress/(0.8*Ec0);
            getApproachFiveToComStrain ( );
			pathFive();

			if ( Tstress > 0 ) //Reach path 7
			{
				TloadingState = 7;
				pathSeven();
			}
		}
	}

	else if ( TloadingState == 3 )
	{
		envelope ( );
	}

	else if ( TloadingState == 4 )
	{
		if ( dStrain > 0)
		{
			envelope ( );
		}
		else
		{
			reverseFromFourStrain = Cstrain;
			reverseFromFourStress = Cstress;
			TloadingState = 6;

			if ( reloadPath != 0) {
				getApproachSixToComStrain ( );
			}
			pathSix();
		}
	}

	else if ( TloadingState == 5 )
	{
		getApproachFiveToComStrain( );
		pathFive();
		if (dStrain<0)
		{
			if ( Tstrain < approachFiveToComStrain )
			{
				envelope( );
			}
		}
		else //dStrain > 0
		{
			if ( Tstress > 0 ) //Reach path 7
			{
				TloadingState = 7;
				pathSeven();
			}
		}
	}	

	else if ( TloadingState == 6 )
	{  
    	if ( reloadPath != 0) {
				getApproachSixToComStrain ( );
		}
		pathSix();
	}

	else if ( TloadingState == 7 )
	{
		pathSeven();
		if (dStrain < 0) {
			//Check if go back to path 5
			if ( Tstrain < interFiveSevenStrain ) 
			{
				TloadingState = 5;
				pathFive();
				getApproachFiveToComStrain ( );
				if ( Tstrain < approachFiveToComStrain )
				{
			    	envelope( );
				}                    
			}
		}
	} // if TloadingState ==7
	
	else
		opserr << " ConcreteL01::determineTrialState -- impropter TloadingState: " 
		     << TloadingState << "\n";

} // end trialState



void ConcreteL01::pathFive()
{
	double Ec0 = 1.4*fpc/epsc0;

    if ( reloadPath == 1 )
	{
		Tstress = Ec0 * ( Tstrain - reverseFromOneStrain ) + reverseFromOneStress;
		Ttangent = Ec0;
	} // if ( reloadPath == 1 )

	else if ( reloadPath == 2 )
	{
		double slope =  0.8*Ec0;
	    Tstress = slope * ( Tstrain - reverseFromTwoStrain ) + reverseFromTwoStress;
		Ttangent = slope;
	} // if ( reloadPath == 2 )

	else
	{
		opserr << " ConcreteL01::pathFive -- improper reloadPath : " 
			 << reloadPath << "\n";
	}
}



void ConcreteL01::pathSeven()
{
	double ecr = 0.00008; //crack strain
	double fcr = 0.31 * sqrt( -fpc ); //crack strength

    // Approach to reverse point at descending branch of tensile envelope
    if ( reverseFromFourStrain > ecr ) 
	{
		if ( Tstrain > reverseFromFourStrain ) {
			envelope();
		}
		else {
			double slope1 = reverseFromFourStress / ( reverseFromFourStrain - interFiveSevenStrain );
		    Tstress = slope1 * ( Tstrain - interFiveSevenStrain );
		    Ttangent = slope1;
		}
	}

	// Approach to peak point of the tensile envelope
	else
	{
		if ( Tstrain > ecr ) {
			envelope();
		}
		else {
			double slope2 =  fcr / ( ecr - interFiveSevenStrain );
			Tstress = slope2 * ( Tstrain - interFiveSevenStrain );
			Ttangent = slope2;
		}
		//if ( Tstrain > 0 ) {
		//	envelope();
		//}
		//else {
		//	Tstress = 0.0;
		//	Ttangent = 0.0;
		//}
	}
}


void ConcreteL01::pathSix()
{  
   double Ec = 3875.0 * sqrt( -fpc );
   double fcr = 0.31 * sqrt( -fpc ); //crack strength
   double ecr = 0.00008; //crack strain
   
   // point A
   double epslonA = reverseFromFourStrain - reverseFromFourStress/Ec; 
   // point B
   double epslonB = reverseFromFourStrain/3.0;
   // point C
   double epslonC;
   double stressC;
   stressC = -1.5*fcr + 0.8*reverseFromFourStress;
   //stressC = -fcr;
   double temp = 1 - stressC/zeta/fpc;
   if ( temp<0 ) {
	   opserr << " ConcreteL01::pathSix -- can not get epslonC \n";
	   epslonC = 0.0;
   }
   else {
	   epslonC = zeta*epsc0*(1-sqrt(temp));
   }
   
   double slope;

   if ( Tstrain > reverseFromFourStrain ) {
	   envelope();
   }
   else if ((Tstrain <= reverseFromFourStrain) && (Tstrain > epslonB)) {
	   slope = (reverseFromFourStress+0.2*fcr)/(reverseFromFourStrain-epslonB);
	   Tstress = reverseFromFourStress + slope*(Tstrain-reverseFromFourStrain);
	   Ttangent = slope;
   }
   else if ((Tstrain <= epslonB) && (Tstrain > epslonC)) {
	   slope = (stressC+0.2*fcr)/(epslonC-epslonB);
	   Tstress = slope*(Tstrain-epslonB)-0.2*fcr;
	   Ttangent = slope;
   }
   else {
		if ( reloadPath == 0 ) // no history stress in compressive zone
		{
			envelope();
		}
		else // has history stress in compressive zone, approach to history point
		{
            if ( reloadPath == 1 ) 
			{
				slope = (reverseFromOneStress-stressC) / (reverseFromOneStrain-epslonC);
				Tstress = slope * (Tstrain-epslonC) + stressC;
				Ttangent = slope;
			}
			else // ( reloadPath == 2 )
			{
				slope = 0.93* (reverseFromTwoStress-stressC) / (reverseFromTwoStrain-epslonC);
				Tstress = slope * (Tstrain-epslonC) + stressC; 
				Ttangent = slope; 
			}  	
			
			// Check if reach compression envelope
			if ( Tstrain < approachSixToComStrain  )
			{				
				envelope();			
			}
		}
   } // else

}


void ConcreteL01::envelope( )
{
	//double fcr = 0.5167 * sqrt( -fpc ); //crack strength for hament's test
	double fcr = 0.31 * sqrt( -fpc ); //crack strength
	double ecr = 0.00008; //crack strain
	double Ec0 = 1.4*fpc/epsc0; //decompression modulus
	//double ecx = epslonCI - sigmaCI/Ec0; //decompression strain
	double Ec = fcr / ecr;  //(ecr - ecx); //modulus of concrete post decompression
	
	if ( Tstrain >= 0 ) // Tension
	{
		if ( Tstrain <= ecr ) //decompression
		{
			//double Ec = 6458.0 * sqrt( -fpc ); // for hament's test
			Tstress = Ec * Tstrain; //(Tstrain - epslonCI) + sigmaCI;
			/*Ttangent = Tstress/Tstrain;
		    if (Tstrain == 0) */
			Ttangent = Ec;
			TloadingState = 3;
		}
		/*else if ( Tstrain > ecx && Tstrain <= ecr ) //decompression
		{
			//double Ec = 6458.0 * sqrt( -fpc ); // for hament's test
			//double Ec = 3875.0 * sqrt( -fpc );
			Tstress = Ec * (Tstrain - ecx);
			Ttangent = Ec;
			TloadingState = 3;
		}*/
		else //Descending branch
		{
			Tstress = fcr * pow((ecr/Tstrain), 0.4);
			/*Ttangent = Tstress/Tstrain;
		    if (Tstrain == 0)*/
			Ttangent = - fcr * 0.4 * pow(ecr, 0.4) * pow(Tstrain, -1.4);
			TloadingState = 4;
		}
	}
	else //Compression
	{
		if ( Tstrain >= zeta*epsc0 ) //Ascending branch
        {	
			TloadingState = 1;
			double eta = Tstrain/(zeta*epsc0);
			Tstress = D * zeta * fpc * (2*eta-eta*eta);
			double Ec0 = 1.4*fpc/epsc0;
			/*Ttangent = Tstress/Tstrain;
		    if (Tstrain == 0)*/ 
			Ttangent = D*Ec0*2/1.4*(1.0-eta);
			if (Tstress >= 0.84 * D * zeta * fpc)
			{
				Tstress = Ec0 * Tstrain;
				/*Ttangent = Tstress/Tstrain;
		        if (Tstrain == 0) */
				Ttangent = Ec0;
			}
		}
		else //Descending branch
		{
			TloadingState = 2;
			double temp = (Tstrain/(zeta*epsc0)-1.0)/(4.0/zeta-1.0);
			
			// initially descending branch
			// Tstress = D * zeta * fpc * ( 1.0 - temp*temp );
			// Ttangent = - D*2.0*fpc*temp/epsc0/(4.0/zeta-1.0);

			// change by using X as powder, 2004/08, X is n
            Tstress = D * zeta * fpc * ( 1.0 - pow(temp, X) );
			/*Ttangent = Tstress/Tstrain;
		    if (Tstrain == 0) */
			Ttangent = - D*fpc*(X)*pow(temp,-1.0+X)/epsc0/(4.0/zeta-1.0);

			//check if Tstress > 0.2*zeta*fpc, change to platum part
			if ( Tstress > D*0.2*zeta*fpc )
			{
				Tstress = D*0.2*zeta*fpc;
				/*Ttangent = Tstress/Tstrain;
		        if (Tstrain == 0)*/
				Ttangent = 0.0;
			}

		}
	}

}


void ConcreteL01::getApproachFiveToComStrain( )
{
    approachFiveToComStrain = 0.0; 

    double Ec0 = 1.4*fpc/epsc0;
    double tempB;  // temporary values for solving equations
    double tempC;  // temporary values for solving equations
    double tempK = 0.0;  // slope of the line y=kx+b
    double b = 0.0;      // distance of the line

    double fiveToOneStrain = 0.0;
    double fiveToTwoStrain = 0.0;
    double fiveToTwoStress = 0.0;

	if ( reloadPath == 1 )
	{
		tempK = Ec0;
		b = - tempK*reverseFromOneStrain + reverseFromOneStress ;
	}
	else if ( reloadPath == 2 )
	{
		tempK = 0.8*Ec0;
		b = - tempK*reverseFromTwoStrain + reverseFromTwoStress ;
	}
	else // error reloadPath
	{
	    opserr << " ConcreteL01::getApproachFiveToComStrain -- improper reloadPath! \n";
	}

    // First, get fiveToOneStrain: intersection of pathFive to ascending branch
	tempB = (tempK-D*Ec0)*zeta*epsc0*epsc0/(D*fpc);
	tempC = b*zeta*epsc0*epsc0/(D*fpc);

	if ( tempB*tempB-4.0*tempC < 0.0 ) 
	{
        opserr << " ConcreteL01::getApproachFiveToComStrain -- can not get root of equation: sqrt(x) x<0! \n";
	}

	fiveToOneStrain = -0.5*tempB - 0.5*sqrt(tempB*tempB-4.0*tempC);

	if ( reverseFromOneStress > D*zeta*fpc) // intersection occurs at ascending branch
	{
       approachFiveToComStrain = reverseFromOneStrain;	   
	}
	else if ( fiveToOneStrain > zeta*epsc0) // intersection occurs at ascending branch
	{
       approachFiveToComStrain = fiveToOneStrain;	   
	}
    else  // intersection occurs at descending branch
	{        
		// Second, get intersection of pathFive to descending branch		
		
		// Initial solution for interection when no X
		//tempB = -2.0*zeta*epsc0 + tempK*pow((4.0*epsc0-zeta*epsc0),2.0)/(D*zeta*fpc);
		//tempC = pow(zeta*epsc0,2.0) + pow((4.0*epsc0-zeta*epsc0),2.0)*( b/(D*zeta*fpc) - 1.0);
       
		//fiveToTwoStrain = -0.5*tempB - 0.5*sqrt(tempB*tempB-4.0*tempC);
		//fiveToTwoStress = tempK * fiveToTwoStrain + b ;


        // Iteration is used to get intersection when there is X, 2004/08

 	    double xn, xnn; // x(i), x(i+1)
	    double fxnn;    // f(xi) = stress(descending)- (kx+b)
	    double fxnp;    // slope of f(xi)

	    xnn = 1.5*zeta*epsc0;  // initial values 
	    fxnn = zeta*D*fpc - zeta*D*fpc*pow(xnn/(zeta*epsc0)-1.0,X)/pow(4.0/zeta-1.0,X) -tempK*xnn -b;
 
	    int counter = 0;

	    if ( ( tempK*zeta*epsc0+b ) < zeta*D*fpc ) {
	        opserr << " ConcreteL01::getApproachFiveToComStrain -- No intersection of reloading path with descending branch! \n";	
	        counter = 50;
		}

	    while ( ( fabs(fxnn) > 1e-4 ) && ( counter < 50 )) {
 	        xn = xnn;
	        fxnp = - (X)*D*fpc*pow(xnn/(zeta*epsc0)-1,-1+X)/pow(4/zeta-1.0,X)/epsc0 - tempK;
	        xnn = xn - fxnn/fxnp;
            fxnn = zeta*D*fpc - zeta*D*fpc*pow((xnn/(zeta*epsc0)-1.0),X)/pow(4.0/zeta-1.0,X) -tempK*xnn -b;
  	        counter++;
		}
 
        if ( counter == 50 ) {
            opserr << " ConcreteL01::getApproachFiveToComStrain -- overflow the iteration limit! \n";
        } else {
			fiveToTwoStrain = xnn;
            fiveToTwoStress = tempK*xnn+b;
        } 

		if ( fiveToTwoStress > D*0.2*zeta*fpc ) // intersection occurs at platum, 0.2*zeta*fpc		  
		{
			approachFiveToComStrain = ( D*0.2*zeta*fpc -b ) / tempK;		
		}
		else // itersection occurs at descending branch
		{
		    approachFiveToComStrain = fiveToTwoStrain;
		}

	}  // itersection occurs at descending branch or platum

	if ( approachFiveToComStrain == 0.0 )
	{
		opserr << " ConcreteL01::getApproachFiveToComStrain -- can not get approachFiveToComStrain! \n";
		opserr << " approachFiveToComStrain = " << approachFiveToComStrain << endln;
		opserr << " reloadPath = " << reloadPath << endln;
		opserr << " zeta = " << zeta << endln;
		opserr << " reverseFromOneStrain = " << reverseFromOneStrain << endln;
		opserr << " reverseFromOneStress = " << reverseFromOneStress << endln;
		opserr << " reverseFromTwoStrain = " << reverseFromTwoStrain << endln;
		opserr << " reverseFromTwoStress = " << reverseFromTwoStress << endln;
		opserr << " fiveToOneStrain = " << fiveToOneStrain << endln;
		opserr << " fiveToTwoStrain = " << fiveToTwoStrain << endln;
		
	}

}






void ConcreteL01::getApproachSixToComStrain( )
{
    approachSixToComStrain = 0.0;

	double Ec0 = 1.4*fpc/epsc0;
    //double tempB;  // temporary values for solving equations
    //double tempC;
    double tempK = 0.0;

	double sixToOneStrain = 0.0;
	double sixToTwoStrain = 0.0;
	double sixToTwoStress = 0.0;

	if ( reloadPath == 1 )
	{
		tempK = reverseFromOneStress / reverseFromOneStrain;
	}
	else if ( reloadPath == 2 )
	{
		tempK = 0.93 * reverseFromTwoStress / reverseFromTwoStrain;
	}
	else // error reloadPath
	{
	    opserr << " ConcreteL01::getApproachSixToComStrain -- improper reloadPath! \n";
	}

    // First, get sixToOneStrain: intersection of pathSix to ascending branch
	sixToOneStrain = (D*Ec0-tempK)*zeta*epsc0*epsc0/(D*fpc); 

	if ( reverseFromOneStress > zeta*D*fpc) // intersection occurs at ascending branch
	{
       approachSixToComStrain = reverseFromOneStrain;
	}
	else if ( sixToOneStrain > zeta*epsc0) // intersection occurs at ascending branch
	{
       approachSixToComStrain = sixToOneStrain;
	}
    else  // intersection occurs at descending branch
	{        
		// Second, get intersection of pathSix to descending branch		
		
		// Initial solution to get itersection
		//tempB = -2.0*zeta*epsc0 + tempK*pow((4.0*epsc0-zeta*epsc0),2.0)/(D*zeta*fpc);
		//tempC = pow(zeta*epsc0,2.0)- pow((4.0*epsc0-zeta*epsc0),2.0);
       
		//sixToTwoStrain = -0.5*tempB - 0.5*sqrt(tempB*tempB-4.0*tempC);
		//sixToTwoStress = tempK * sixToTwoStrain;

		// iteration is used when there is X, 2004/08

        double xn, xnn; // x(i), x(i+1)
	    double fxnn;    // f(xi) = stress(descending)- (kx+b)
	    double fxnp;    // slope of f(xi)

	    xnn = 1.5*zeta*epsc0;  // initial values 
	    fxnn = zeta*D*fpc - zeta*D*fpc*pow(xnn/(zeta*epsc0)-1.0, X)/pow(4.0/zeta-1.0, X) -tempK*xnn;
 
	    int counter = 0;

	    if ( ( tempK*zeta*epsc0 ) < zeta*D*fpc ) {
	        opserr << " ConcreteL01::getApproachFiveToComStrain -- No intersection of reloading path with descending branch! \n";	
	        counter = 50;
		}

	    while ( ( fabs(fxnn) > 1e-4 ) && ( counter < 50 )) {
 	        xn = xnn;
	        fxnp = - (X)*D*fpc*pow(xnn/(zeta*epsc0)-1,-1+X)/pow(4/zeta-1.0,X)/epsc0 - tempK;
	        xnn = xn - fxnn/fxnp;
            fxnn = zeta*D*fpc - zeta*D*fpc*pow((xnn/(zeta*epsc0)-1.0), X)/pow(4.0/zeta-1.0, X) -tempK*xnn;
  	        counter++;
		}
 
        if ( counter == 50 ) {
            opserr << " ConcreteL01::getApproachSixToComStrain -- overflow the iteration limit! \n";
        } else {
            sixToTwoStrain = xnn;
            sixToTwoStress = tempK*xnn;
        } 

		if ( sixToTwoStress > D*0.2*zeta*fpc ) // intersection occurs at platum, 0.2*zeta*fpc		  
		{
			approachSixToComStrain = D*0.2*zeta*fpc / tempK;
		} 
		else // intersection occurs at descending branch
		{
			approachSixToComStrain = sixToTwoStrain;
		}
	
	}

	if ( approachSixToComStrain == 0.0 )
	{
		opserr << " ConcreteL01::getApproachSixToComStrain -- can not get approachSixToComStrain! \n";
	}
}




double ConcreteL01::getPD ()
{
  double PD; 
  //PD = 0.0;
  
  double tempRatio;	
  if ( epslonTP <= 0.0 )
    {
      PD = 0.0;
    }
  else
    {
      if     ( TloadingState == 1 )
	{
	  tempRatio = Tstrain/(zeta*epsc0);
	  PD = - D * fbeta * Wp * 1160.0 * sqrt(-fpc)* pow((1+400.0*epslonTP), -1.5)
	    * pow(tempRatio,2.0);
	}
      else if ( TloadingState == 2 )
	{
	  if ( Ttangent == 0.0 ) // at the end platum part of descending branch FMK CHANGED FROM = 0.0
	    {
	      PD = 0.0;
	    }
	  else
	    {
	      tempRatio = Tstrain/(zeta*epsc0);
	      PD = - D * fbeta * Wp *1160.0 * sqrt(-fpc)* pow((1+400.0*epslonTP), -1.5) 
		* (1.0 -(tempRatio-1)/pow(4.0/zeta-1.0,3.0)*(1.0-12.0/zeta+(4.0/zeta+1.0)*tempRatio)); 
	      // 1160, 400
	    }
	}		
      else
	{
	  PD = 0.0;
	}
      if ( (zeta == 0.9) || (zeta == 0.25)) // zeta = max or min value
	{
	  PD = 0.0;
	}
    }
  
  return PD;
}



UniaxialMaterial* ConcreteL01::getCopy ()
{
   ConcreteL01* theCopy = new ConcreteL01(this->getTag(), fpc, epsc0);  
   
   // History variables
   theCopy->TloadingState = TloadingState;
   theCopy->CloadingState = CloadingState;
   
   theCopy->reloadPath = reloadPath;

   theCopy->reverseFromOneStrain = reverseFromOneStrain;
   theCopy->reverseFromOneStress = reverseFromOneStress;
   theCopy->reverseFromTwoStrain = reverseFromTwoStrain; 
   theCopy->reverseFromTwoStress = reverseFromTwoStress;
   theCopy->reverseFromFourStrain = reverseFromFourStrain; 
   theCopy->reverseFromFourStress = reverseFromFourStress;   

   theCopy->interFiveSevenStrain = interFiveSevenStrain;

   theCopy->approachFiveToComStrain = approachFiveToComStrain;
   theCopy->approachSixToComStrain = approachSixToComStrain;
 
   // State variables
   theCopy->zeta = zeta;
   theCopy->beta = beta;
   theCopy->epslonTP = epslonTP;
   
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;

   theCopy->D = D;

   return theCopy;
}


int ConcreteL01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(21);
   data(0) = this->getTag();

   // Material properties
   data(1) = fpc;
   data(2) = epsc0; 
   data(3) = zeta;
   data(4) = beta;
   data(5) = epslonTP;

   // History variables from last converged state
   data(6) = CloadingState;
   data(7) = reloadPath;
   data(8) = reverseFromOneStrain;
   data(9) = reverseFromOneStress;
   data(10) = reverseFromTwoStrain; 
   data(11) = reverseFromTwoStress;
   data(12) = reverseFromFourStrain; 
   data(13) = reverseFromFourStress;
   data(14) = interFiveSevenStrain;
   data(15) = approachFiveToComStrain;
   data(16) = approachSixToComStrain;

   // State variables from last converged state
   data(17) = Cstrain;
   data(18) = Cstress;
   data(19) = Ctangent;

   data(20) = D;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "ConcreteL01::sendSelf() - failed to send data\n";

   return res;
}


int ConcreteL01::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(21);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "ConcreteL01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties 
      fpc = data(1);
      epsc0 = data(2); 
      zeta = data(3);
      beta = data(4);
      epslonTP = data(5);

      // History variables from last converged state
	  CloadingState = int(data(6));
      reloadPath = int(data(7));
      reverseFromOneStrain = data(8);
      reverseFromOneStress = data(9);
      reverseFromTwoStrain = data(10); 
      reverseFromTwoStress = data(11);
      reverseFromFourStrain = data(12); 
      reverseFromFourStress = data(13);
      interFiveSevenStrain = data(14);
      approachFiveToComStrain = data(15);
      approachSixToComStrain = data(16);

	  // Copy converged history values into trial values since data is only
      // sent (received) after convergence	  
      TloadingState = CloadingState;

      // State variables from last converged state
      Cstrain = data(17);
      Cstress = data(18);
      Ctangent = data(19);

	  D = data(20);

      // Set trial state variables
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }

   return res;
}


Response* 
ConcreteL01::setResponse(const char **argv, int argc,
			 OPS_Stream &theOutput)
{
  Response *theResponse = 0;

  if (strcmp(argv[0],"getPD") == 0) {
    double data = 0.0;
    theResponse = new MaterialResponse(this, 100, data);
  } else if (strcmp(argv[0],"setWallVar") == 0) {
    theResponse = new MaterialResponse(this, 101, Vector(5));
  } else
    return this->UniaxialMaterial::setResponse(argv, argc, theOutput);

  return theResponse;
}
 
int 
ConcreteL01::getResponse(int responseID, Information &matInfo)
{
  if (responseID == 100) {
    matInfo.theDouble = this->getPD();
  } else if (responseID == 101){
    Vector *theVector = matInfo.theVector;
    X = (*theVector)(0);
    K = (*theVector)(1);
    D = (*theVector)(2);
    beta = (*theVector)(3);
    epslonTP = (*theVector)(4);   
  } else
    return this->UniaxialMaterial::getResponse(responseID, matInfo);

  return 0;
}


void ConcreteL01::Print (OPS_Stream& s, int flag)
{
   s << "ConcreteL01, tag: " << this->getTag() << endln;
   //s << "  fpc: " << fpc << endln;
   //s << "  epsc0: " << epsc0 << endln;

   s << " strain: " << this->getStrain() << endln;
   s << " stress: " << this->getStress() << endln;
   s << " tangent: " << this->getTangent() << endln;
   //s << " PD: " << this->getPD() << endln;
   s << " zeta: " << zeta << endln;
   s << " D: " << D << endln;
   s << " TloadingState: " << TloadingState << endln;
   s << " reverseFromFourStrain: " << reverseFromFourStrain << endln;   

}

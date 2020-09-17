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
                                                                        
// $Revision: 1.2 $
// $Date: 2010-01-20 20:34:28 $
                                                                        
// File: ~/Bond_SP01.cpp
//
// Written: 		Jian Zhao, Iowa State University 		04/2004
// Revision:		Jian Zhao, University of Wisconsin, Milwaukee 		04/2006
// Revision:		Jian Zhao, University of Wisconsin, Milwaukee. Fixed the jump in the negative direction
//
// Description: This file contains the class definition for Uniaxial material Bond_SP01. 
// Bond_SP01: Strain penetration of fully anchored rebars w/o bond damage.


#include <Bond_SP01.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_Bond_SP01(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numInput = OPS_GetNumRemainingInputArgs();
  if (numInput != 7 && numInput != 11) {
    opserr << "Invalid #args,  uniaxialMaterial Bond_SP01 tag? fy? sy? fu? su? b? R?";
    opserr << " <Cd? db? fc? la?>" << endln;	
    return 0;
  }
  
  int iData[1];
  double dData[10];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
    return 0;
  }

  numData = numInput-1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
    return 0;
  }

  if (numInput == 7)
    theMaterial = new Bond_SP01 (iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], 
				 dData[5]);
  else
    theMaterial = new Bond_SP01 (iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], 
				 dData[5], dData[6], dData[7], dData[8], dData[9]);

  return theMaterial;
}


Bond_SP01::Bond_SP01
(int tag, double FY, double SY, double FU, double SU, double KZ, double r, double CD, double DB, double FC, double LA):
  UniaxialMaterial(tag,MAT_TAG_Bond_SP01),
  fy(FY), sy(SY),  fu(FU), su(SU), Kz(KZ), R(r), Cd(CD), db(DB), fc(FC), lba(LA)
{
  // Check units.  Need parameters in ksi and in.
  if ( fy >= 1000 || sy >= 1)
    opserr << "WARNING: For the Strain-Penetration Model: input values in ksi and in." << endln;
  
  
  // Set parameters for bar stress-slip envelope 
  Cr = 1.01;						//% need to be large than 1 pretty arbitrary
  Ks = pow(R,Kz/2.5);		//% pretty arbitrary
  slvrg = pow(12.0/30.0,1/0.4)*0.04;			//%the slip corresponding to the virgin friction in local bond-slip model
  
  // Assume symmetric envelope. This needs to be changed later because end-bearing participate when under compression
  E0 = fy/sy;
  
  // bond condition (for future use) 
  la = fy*db*1000.0/40.0/pow(fc*1000,0.5);	// effective anchorage length (not being used inversions 1.x)
  
  // Set all history and state variables to initial values
  this->revertToStart ();

}

Bond_SP01::Bond_SP01
(int tag, double FY, double SY, double FU, double SU, double KZ, double r): 
UniaxialMaterial(tag,MAT_TAG_Bond_SP01),
 fy(FY), sy(SY), fu(FU), su(SU), Kz(KZ), R(r), Cd(0.0), db(1.0), fc(4.35), lba(0.0)
{
	// Check units.  Need parameters in ksi and in.
	if (fy >= 1000 || sy >= 1)
		opserr << "WARNING: WARNING: For the Strain-Penetration Model: input values in ksi and in." << endln;

	
	// Set bar stress-slip envelope  (yield point)
	Cr = 1.01;				//% need to be large than 1 pretty arbitrary
	Ks = pow(R,Kz/2.5);		//% pretty arbitrary
	slvrg = pow(12.0/30.0,1/0.4)*0.04;			//%the slip corresponding to the virgin friction in local bond-slip model

	// Assume symmetric envelope 
	E0 = fy/sy;

	// bond condition (for calculation) 
    la = fy*db*1000.0/40.0/pow(fc*1000,0.5);	// effective anchorage length

	// Set all history and state variables to initial values
	this->revertToStart ();
}

Bond_SP01::Bond_SP01 ():
UniaxialMaterial(0,MAT_TAG_Bond_SP01),
fy(64.0), sy(0.021), fu(102.2), su(1.0), Kz(0.3), R(.8), Cd(0.0), db(1.0), fc(4.35), lba(0.0)
{
	// constructor w/o parameter.  does nothing

}

Bond_SP01::~Bond_SP01 ()
{
}

int Bond_SP01::setTrialStrain (double strain, double strainRate)
{
	// Reset history variables to last converged state
	TRSlip = CRSlip;			// Return slip
	TRLoad = CRLoad;			// Return load
	TRSlope = CRSlope;			// Return slope
	TmaxHSlip = CmaxHSlip;		// Maximum slip in tension
	TminHSlip = CminHSlip;		// Maximum slip in compression
	Tloading = Cloading;		// Loading flag
	TYieldFlag = CYieldFlag;	// Yield flag

	// Set trial slip
	Tslip = strain;

	// Determine change in slip from last converged state
	double dslip = Tslip - Cslip;

	// Calculate the trial state given the trial Strain
	determineTrialState (Tslip, dslip);

	return 0;
}

int Bond_SP01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
	// Reset history variables to last converged state
	TRSlip = CRSlip;			// Return slip
	TRLoad = CRLoad;			// Return load
	TRSlope = CRSlope;			// Return slope
	TmaxHSlip = CmaxHSlip;		// Maximum slip in tension
	TminHSlip = CminHSlip;		// Maximum slip in compression
	Tloading = Cloading;		// Loading flag
	TYieldFlag = CYieldFlag;	// Yield flag

	// Set trial Strain
	Tslip = strain;

	// Determine change in Strain from last converged state
	double dslip;
	dslip = Tslip - Cslip;

	// Calculate the trial state given the trial Strain
	determineTrialState (Tslip, dslip);

	stress = Tload;
	tangent = Ttangent;

	return 0;
}

void Bond_SP01::determineTrialState (double ts, double dslip)
{

	double maxrs;					// maximum return slip in tension in history
	double maxrl;					// maximum return load in tension in history
	double minrs;					// maximum return slip in compression in history
	double minrl;					// maximum return load in compression in history
	double rslip;					// last return slip 
	double rload;					// last return load
	double rsvg;					// last return vergin slip
	double trs;						// this return slip
	double trl;						// this return load
	double trsvg;					// this return vergin slip
	double Eun;						// unloading slope (tension --> compression)
	double Ere;						// reloading slope (compression --> tension)
	double Er;						// return slope 
	double kkk = 0.38;				// unloading slope (damage)
	double R1;						// temporary virables

	double ss;						// normalized slip
	double Sy;						// dummy yield slip to keep reloadig slope E0
	double ssy;						// ratio of slip to Sy
	double suy;						// ratio of dummy su to Sy
	double ft;						// normalized load
	double st;						// normalized slope
	
	if (fabs(dslip) <= DBL_EPSILON)             //ignore trivial slip change
    {
		Tload = Cload;
		Ttangent = Ctangent;
		return;
	}

	if (Tloading == 0)						//the initial step when loading/unloading is not determined
	{
		Tload = getEnvelopeStress(ts);
		if (dslip > 0.0)					//positive change in slip --> loading (Tloading = 1)
		{
			Tloading = 1;
			CminHSlip = -slvrg;				//avoid divided by zero later in determine reloading path
		}
		else								//negative change in slip --> unloading (Tloading = -1)
		{
			Tloading = -1;
			CmaxHSlip = slvrg;				//avoid divided by zero later in determine reloading path
		}
		return;
	}

	if (TYieldFlag == 0)					//no yielding in previous steps 
	{
		Tload = getEnvelopeStress(ts);
		if (Tloading > 0.0)									// (Cloading>0): was loading (tension)
		{
			if (dslip < 0.0)											//% (dslip<0): turn unloading (compression)
			{
				Tloading = -1;
				TRSlip = Cslip;
				TRLoad = Cload;
				TRSlope = E0;
				if (Cslip > TmaxHSlip)
					TmaxHSlip = Cslip;
			}
		}
		else								//(Tloading<0): was unloading (compression)
		{
			if (dslip > 0.0)					//%(dslip<0): turn reloading (tension)
			{
				Tloading = 1;
				TRSlip = Cslip;
				TRLoad = Cload;
				TRSlope = E0;
				if (Cslip < TminHSlip)
					TminHSlip = Cslip;
			}
		}
		return;
	}

	//set limits for elastic unloading/reloading: before the curve hits the x-axis
	maxrs = TmaxHSlip;
	maxrl = getEnvelopeStress(maxrs);
	Eun = E0;

	minrs = TminHSlip;
	minrl = getEnvelopeStress(minrs);
	Ere = E0;

	rslip = TRSlip;
	rload = TRLoad;
	Er = TRSlope;
	rsvg = rslip-rload/Er;

	//get load for a slip
	if (Tloading > 0.0)									// (Cloading>0): was loading (tension)
	{
		if (dslip > 0.0)											//% (dslip>0): keep loading (tension)
		{
			if (ts >= maxrs)										//%% larger than the max. slip -> envelope stress
			{
				Tload = getEnvelopeStress(ts);
			}
			else													//%% ts < maxrs
			{
				if (rsvg > rslip)										//%%  turn from compression to tension
				{
					if (ts >= rsvg)										//%%% curve reloading (tension)
					{
						Sy = maxrl/Ere;
						suy = (maxrs-rsvg)/Sy;
						ssy = (ts-rsvg)/Sy; 
						ss = ssy/(suy-ssy);
						R1 = R+(1.01-R)*pow(((ts-rsvg)/(maxrs-rsvg)),(1/R/R));			//%%%% to avoid infinite slope 
						ft = ss/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1));
						Tload = ft*maxrl;
						st = pow(suy,(1-R1))/pow(suy-ssy,2)/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1+1));
						Ttangent = st*Ere;
					}
					else							//%% ts < rsvg
					{
						Tload = rload*(ts-rsvg)/(rslip-rsvg);
						Ttangent = Ere;
					}
				}
				else										//%%  rsvg < rslip turn from linear tension to tension
				{
					Sy = maxrl/Er;
					suy = (maxrs-rslip)/Sy;
					ssy = (ts-rslip)/Sy; 
					ss = ssy/(suy-ssy);
					R1 = R+(1.01-R)*pow(((ts-rslip)/(maxrs-rslip)),(1/R/R));	//%%%% to avoid infinite slope 
					ft = ss/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1));
					Tload = ft*(maxrl-rload)+rload;
					st = pow(suy,(1-R1))/pow(suy-ssy,2)/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1+1));
					Ttangent = st*Er;
				}						
			}
		}
		else								//%(dslip<0): turn unloading (compression)
		{
			Tloading = -1;
			TRSlip = Cslip;
			TRLoad = Cload;
			TRSlope = E0;
			if (Cslip > TmaxHSlip)
			{
				TmaxHSlip = Cslip;
				maxrs = TmaxHSlip;
				maxrl = getEnvelopeStress(maxrs);
				Eun = E0;
			}

			trs = TRSlip;
			trl = TRLoad;
			Er = TRSlope;
			trsvg = trs-trl/Er;
        
			if (ts <= minrs)
			{
				Tload = getEnvelopeStress(ts);
			}
			else
			{
                if (trsvg < trs)						//%% from tension to compression
				{
					if (ts <= trsvg)					//%%linear unloading (compression)
					{
						Sy = minrl/Eun;
						suy = (minrs-trsvg)/Sy;
						ssy = (ts-trsvg)/Sy; 
						ss = ssy/(suy-ssy);
						R1 = R+(1.01-R)*pow(((ts-trsvg)/(minrs-trsvg)),(1/R/R));			//%%%% to avoid infinite slope 
						ft = ss/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1));
						Tload = ft*minrl;
						st = pow(suy,(1-R1))/pow(suy-ssy,2)/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1+1));
						Ttangent = st*Eun;
					}
					else
					{
						Tload = trl*(ts-trsvg)/(trs-trsvg);
						Ttangent = Eun;
					}
				}
				else							//%% from compression to compression
				{
					Sy = minrl/Er;
					suy = (minrs-trs)/Sy;
					ssy = (ts-trs)/Sy; 
					ss = ssy/(suy-ssy);
					R1 = R+(1.01-R)*pow(((ts-trs)/(minrs-trs)),(1/R/R));	//%%%% to avoid infinite slope 
					ft = ss/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1));
					Tload = ft*(minrl-trl)+trl;
					st = pow(suy,(1-R1))/pow(suy-ssy,2)/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1+1));
					Ttangent = st*Er;
				}
			}
		}
	}
	else								//(Tloading<0): was unloading (compression)
	{
		if (dslip < 0.0)					//%(dslip<0): keep unloading (compression)
		{
			if (ts <= minrs)					//%%envelope stress
			{
				Tload = getEnvelopeStress(ts);
			}
			else								//%% ts > minrs								
			{
				if (rsvg < rslip)					//%%tension to compression and keep compression
				{
					if (ts <= rsvg)					//%% ts < rsvg curve unloading (compression)
					{
						Sy = minrl/Eun;
						suy = (minrs-rsvg)/Sy;
						ssy = (ts-rsvg)/Sy; 
						ss = ssy/(suy-ssy);
						R1 = R+(1.01-R)*pow(((ts-rsvg)/(minrs-rsvg)),(1/R/R));			//%%%% to avoid infinite slope 
						ft = ss/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1));
						Tload = ft*minrl;
						st = pow(suy,(1-R1))/pow(suy-ssy,2)/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1+1));
						Ttangent = st*Eun;
					}
					else				//%%linear unloading (compression)
					{
						Tload = rload*(ts-rsvg)/(rslip-rsvg);
						Ttangent = Eun;
					}

				}
				else							//%%%  unloading from linear reloading and keep unloading
				{
					Sy = minrl/Er;
					suy = (minrs-rslip)/Sy;
					ssy = (ts-rslip)/Sy; 
					ss = ssy/(suy-ssy);
					R1 = R+(1.01-R)*pow(((ts-rslip)/(minrs-rslip)),(1/R/R));	//%%%% to avoid infinite slope 
					ft = ss/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1));
					Tload = ft*(minrl-rload)+rload;
					st = pow(suy,(1-R1))/pow(suy-ssy,2)/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1+1));
					Ttangent = st*Er;
				}
			}
		}
		else								//(dslip>0): turn loading (tension)
		{
			Tloading = 1;
			TRSlip = Cslip;
			TRLoad = Cload;
			TRSlope = E0;
			if (Cslip < TminHSlip)
			{
				TminHSlip = Cslip;
				minrs = TminHSlip;
				minrl = getEnvelopeStress(minrs);
				Ere = E0;
			}

			trs = TRSlip;
			trl = TRLoad;
			Er = TRSlope;
			trsvg = trs-trl/Er;
        
			if (ts >= maxrs)						//%%envelope stress
			{
				Tload = getEnvelopeStress(ts);
			}
			else								//%% ts > trs
			{
				if (trsvg > trs)			//%%compression to tension
				{
					if (ts >= trsvg)					//%%curve reloading (tension)
					{
						Sy = maxrl/Ere;
						suy = (maxrs-trsvg)/Sy;
						ssy = (ts-trsvg)/Sy; 
						ss = ssy/(suy-ssy);
						R1 = R+(1.01-R)*pow(((ts-trsvg)/(maxrs-trsvg)),(1/R/R));			//%%%% to avoid infinite slope 
						ft = ss/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1));
						Tload = ft*maxrl;
						st = pow(suy,(1-R1))/pow(suy-ssy,2)/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1+1));
						Ttangent = st*Ere;
					}
					else							//%% ts < trsvg linear reloading (conpression)
					{
						Tload = trl*(ts-trsvg)/(trs-trsvg);
						Ttangent = Ere;
					}
				}
				else						
				{
					Sy = maxrl/Er;
					suy = (maxrs-trs)/Sy;
					ssy = (ts-trs)/Sy; 
					ss = ssy/(suy-ssy);
					R1 = R+(1.01-R)*pow(((ts-trs)/(maxrs-trs)),(1/R/R));	//%%%% to avoid infinite slope 
					ft = ss/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1));
					Tload = ft*(maxrl-trl)+trl;
					st = pow(suy,(1-R1))/pow(suy-ssy,2)/pow((pow(1/suy,R1)+pow(ss,R1)),(1/R1+1));
					Ttangent = st*Er;
				}
			}
        }
    }
}

void Bond_SP01::detectStressReversal (double dslip)
{
	if (Tloading > 0.0)									// (Cloading>0): was loading (tension)
	{
		if (dslip < 0.0)											//% (dslip<0): turn unloading (compression)
		{
			Tloading = -1;
			TRSlip = Cslip;
			TRLoad = Cload;
			TRSlope = E0;
			if (Cslip > TmaxHSlip)
				TmaxHSlip = Cslip;
		}
	}
	else								//(Tloading<0): was unloading (compression)
	{
		if (dslip > 0.0)					//%(dslip<0): turn reloading (tension)
		{
			Tloading = 1;
			TRSlip = Cslip;
			TRLoad = Cload;
			TRSlope = E0;
			if (Cslip < TminHSlip)
				TminHSlip = Cslip;
		}
	}
}

double Bond_SP01::getEnvelopeStress (double s)
{
	double ss;					// normalized slip
	double ft;					// normalized load
	double st;					// normalized slope
	double ssy;					// the ratio of s-Sy to Sy
	double suy;					// the ratio of su-Sy to Sy
	double load;
	
	if (fabs(s) < DBL_EPSILON)
	{
		load = 0.0;
		Ttangent = E0;
		return load;
	}

	if (s > 0.0)				//tension side
	{
		if (s <= sy)
		{
			load = s*E0;                     //%stress
			Ttangent = E0;                   //%slope
		}
		else 
		{
			TYieldFlag = 1;
			if (s < su)
			{
				ssy = (s-sy)/sy;
				suy = (su-sy)/sy;
				ss = ssy/(suy-ssy);
				ft = ss/pow((pow(1/suy/Kz,Cr)+pow(ss,Cr)),1/Cr);
				load = fy+ft*(fu-fy);
				st = (pow(suy,(1-Cr))/pow(Kz,Cr)/pow((suy-ssy),2))/pow((pow(1/suy/Kz,Cr)+pow(ss,Cr)),(1/Cr+1));
				Ttangent = st*E0;
			}
			else
			{
				load = fu;
				Ttangent = 0.0;
			}
		}
	}
	else						//compression side
	{
		if (s >= -sy)
		{
			load = s*E0;
			Ttangent = E0;
		}
		else 
		{
			TYieldFlag = 1;
			if (s > -su)
			{
				ssy = (s-(-sy))/(-sy);
				suy = (su-sy)/sy;
				ss = ssy/(suy-ssy);
				ft = ss/pow((pow(1/suy/Kz,Cr)+pow(ss,Cr)),1/Cr);
				load = -fy+ft*(-fu+fy);
				st = (pow(suy,(1-Cr))/pow(Kz,Cr)/pow((suy-ssy),2))/pow((pow(1/suy/Kz,Cr)+pow(ss,Cr)),(1/Cr+1));
				Ttangent = st*E0;
			}
			else
			{
				load = -fu;
				Ttangent = 0.0;
			}
		}
	}

	return load;
}

double Bond_SP01::getStrain ()
{
   return Tslip;
}

double Bond_SP01::getStress ()
{
   return Tload;
}

double Bond_SP01::getTangent ()
{
   return Ttangent;
}

double Bond_SP01::getInitialTangent()
{
	return fy/sy;
}

int Bond_SP01::commitState ()
{
	// History variables
	CRSlip = TRSlip;			// Return slip
	CRLoad = TRLoad;			// Return load
	CRSlope = TRSlope;			// Return slope
	CmaxHSlip = TmaxHSlip;		// Maximum slip in tension
	CminHSlip = TminHSlip;		// Maximum slip in compression
	Cloading = Tloading;		// Loading flag
	CYieldFlag = TYieldFlag;	// Yield flag

	// State variables
	Cslip = Tslip;
	Cload = Tload;
	Ctangent = Ttangent;

	return 0;
}

int Bond_SP01::revertToLastCommit ()
{
	// Reset trial history variables to last committed state
	TRSlip = CRSlip;			// Return slip
	TRLoad = CRLoad;			// Return load
	TRSlope = CRSlope;			// Return slope
	TmaxHSlip = CmaxHSlip;		// Maximum slip in tension
	TminHSlip = CminHSlip;		// Maximum slip in compression
	Tloading = Cloading;		// Loading flag
	TYieldFlag = CYieldFlag;	// Yield flag

	// Reset trial state variables to last committed state
	Tslip = Cslip;
	Tload = Cload;
	Ttangent = Ctangent;

	return 0;
}

int Bond_SP01::revertToStart ()
{
	// History variables
	CRSlip = 0.0;
	CRLoad = 0.0;
	CRSlope = E0;			
	CmaxHSlip = 0.0;
	CminHSlip = 0.0;
	Cloading = 0;
	CYieldFlag = 0;

	TRSlip = 0.0;
	TRLoad = 0.0;
	TRSlope = E0;			
	TmaxHSlip = 0.0;
	TminHSlip = 0.0;
	Tloading = 0;
	TYieldFlag = 0;

	// State variables
	Cslip = 0.0;
	Cload = 0.0;
	Ctangent = E0;		 

	Tslip = 0.0;
	Tload = 0.0;
	Ttangent = E0; 

	return 0;
}

UniaxialMaterial* Bond_SP01::getCopy ()
{
   Bond_SP01* theCopy = new Bond_SP01(this->getTag(), fy, sy, fu, su, Kz, R, Cd, db, fc, lba);

   // Converged history variables
   theCopy->CRSlip = CRSlip;
   theCopy->CRLoad = CRLoad;
   theCopy->CRSlope = CRSlope;
   theCopy->CmaxHSlip = CmaxHSlip;
   theCopy->CminHSlip = CminHSlip;
   theCopy->Cloading = Cloading;
   theCopy->CYieldFlag = CYieldFlag;

   // Trial history variables
   theCopy->TRSlip = TRSlip;
   theCopy->TRLoad = TRLoad;
   theCopy->TRSlope = TRSlope;
   theCopy->TmaxHSlip = TmaxHSlip;
   theCopy->TminHSlip = TminHSlip;
   theCopy->Tloading = Tloading;
   theCopy->TYieldFlag = TYieldFlag;

   // Converged state variables
   theCopy->Cslip = Cslip;
   theCopy->Cload = Cload;
   theCopy->Ctangent = Ctangent;

   // Trial state variables
   theCopy->Tslip = Tslip;
   theCopy->Tload = Tload;
   theCopy->Ttangent = Ttangent;

   // Other global virables
   theCopy->Cr = Cr;
   theCopy->Ks = Ks;
   theCopy->E0 = E0;
   theCopy->la = la;

   return theCopy;
}

int Bond_SP01::sendSelf (int cTag, Channel& theChannel)
{
  int res = 0;
  static Vector dData(26);

  dData(0) = this->getTag();
  dData(1) =db;
  dData(2) =fc;
  dData(3) =lba;
  dData(4) =la;
  dData(5) =sy;
  dData(6) =su;
  dData(7) =fy;
  dData(8) =fu;
  dData(9) =E0;
  dData(10) =Kz;
  dData(11) =Cr;
  dData(12) =Ks;
  dData(13) =slvrg;
  dData(14) =R;		
  dData(15) =Cd;	
  dData(16) =CRSlip;	
  dData(17) =CRLoad;	
  dData(18) =CRSlope;	
  dData(19) =CmaxHSlip;	
  dData(20) =CminHSlip;	
  dData(21) =Cloading;	
  dData(22) =  CYieldFlag;
  dData(23) =Cslip;
  dData(24) =Cload;
  dData(25) =Ctangent;

  res = theChannel.sendVector(this->getDbTag(), cTag, dData);
  if (res < 0) 
    opserr << "Bond_SP01::sendSelf() - failed to send data\n";

  return res;
}

int Bond_SP01::recvSelf (int cTag, Channel& theChannel,
			 FEM_ObjectBroker& theBroker)
{
  int res = 0;
  static Vector dData(26);

  res = theChannel.recvVector(this->getDbTag(), cTag, dData);
  if (res < 0) 
    opserr << "Bond_SP01::sendSelf() - failed to send data\n";

  this->setTag(int(dData(0)));
  db = dData(1);
  fc= dData(2);
  lba= dData(3);
  la= dData(4);
  sy= dData(5);
  su= dData(6);
  fy= dData(7);
  fu= dData(8);
  E0= dData(9);
  Kz= dData(10);
  Cr= dData(11);
  Ks= dData(12);
  slvrg= dData(13);
  R= dData(14);		
  Cd= dData(15);	
  CRSlip= dData(16);	
  CRLoad= dData(17);	
  CRSlope= dData(18);	
  CmaxHSlip= dData(19);	
  CminHSlip= dData(20);	
  Cloading= int(dData(21));	
  CYieldFlag= int(dData(22));
  Cslip= dData(23);
  Cload= dData(24);
  Ctangent = dData(25);
  
  this->revertToLastCommit();
    
  return res;
}

void Bond_SP01::Print (OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "Bond_SP01 tag: " << this->getTag() << endln;
        s << "  sy: " << sy << ", ";
        s << "  fy: " << fy << ", ";
        s << "  su: " << su << ", ";
        s << "  fu: " << fu << ", ";
        s << "  Kz: " << Kz << ", ";
        s << "  R: " << R << ", ";
        s << "  Cd: " << Cd << ", ";
        s << "  db: " << db << ", ";
        s << "  fc: " << fc << ", ";
        s << "  lba:" << lba;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"Bond_SP01\", ";
        s << "\"sy\": " << sy << ", ";
        s << "\"fy\": " << fy << ", ";
        s << "\"su\": " << su << ", ";
        s << "\"fu\": " << fu << ", ";
        s << "\"Kz\": " << Kz << ", ";
        s << "\"R\": " << R << ", ";
        s << "\"Cd\": " << Cd << ", ";
        s << "\"db\": " << db << ", ";
        s << "\"fc\": " << fc << ", ";
        s << "\"lba\": " << lba << "}";
    }
}



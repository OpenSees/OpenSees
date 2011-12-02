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
                                                                        
// $Revision: 1.1 $
// $Date: 2006-05-24 21:10:28 $

// Written: 		Jian Zhao, Iowa State University 		04/2004
// Revision:		Jian Zhao, University of Wisconsin, Milwaukee 		04/2006
//
// Description: This file contains the class defination for Uniaxial material Bond_SP01. 
//		Bond_SP01: Strain penetration of rebars at footings w/o bond damage.


#include <Bond_SP01.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>

Bond_SP01::Bond_SP01
(int tag, double FY, double SY, double FU, double SU, double KZ, double r, double CD, double DB, double FC, double LA):
UniaxialMaterial(tag,MAT_TAG_Bond_SP01),
fy(FY), sy(SY),  fu(FU), su(SU), Kz(KZ), R(r), Cd(CD), db(DB), fc(FC), lba(LA)
{
	// Check units.  Need parameters in ksi and in.
	if ( fy >= 1000 || sy >= 1)
		opserr << "WARNING: For the Strain-Penetration Model: input values in ksi and in." << endln;

	
	// Set bar stress-slip envelope 
	Cr = 1.01;				//% need to be large than 1 pretty arbitary
	Ks = pow(R,Kz/2.5);		//% pretty arbitary

	// Assume symmetric envelope 
	E0 = fy/sy;

	// bond condition (for future use) 
    la = fy*db*1000.0/40.0/pow(fc*1000,0.5);	// effective anchorage length

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
	Cr = 1.01;				//% need to be large than 1 pretty arbitary
	Ks = pow(R,Kz/2.5);		//% pretty arbitary

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

	double maxrs;				// maximum return slip in tension in history
	double maxrl;				// maximum return load in tension in history
	double maxvgs;				// maximum return vergin slip in tension in history
	double minrs;				// maximum return slip in compression in history
	double minrl;				// maximum return load in compression in history
	double minvgs;				// maximum return vergin slip in compression in history
	double rslip;				// last return slip 
	double rload;				// last return load
	double rsvg;				// last return vergin slip
	double trs;					// this return slip
	double trl;					// this return load
	double trsvg;				// this return vergin slip
	double Eun;					// unloading slope (tension --> compression)
	double Ere;					// reloading slope (compression --> tension)
	double Er;					// return slope 
	double kkk = 0.38;			// unloading slope (damage)
	double templ, temps;		// temporary virables

	double ss;					// normalized slip
	double Sy;					// dummy yield slip to keep reloadig slope E0
	double ssy;					// ratio of slip to Sy
	double suy;					// ratio of dummy su to Sy
	double ft;					// normalized load
	double st;					// normalized slope
	
	if (fabs(dslip) <= DBL_EPSILON)             //ignore trivial slip change
    {
		Tload = Cload;
		return;
	}

	if (Tloading == 0)
	{
		Tload = getEnvelopeStress(ts);
		if (dslip > 0.0) 
		{
			Tloading = 1;
			CminHSlip = -ts;			//avoid divided by zero
		}
		else
		{
			Tloading = -1;
			CmaxHSlip = -ts;			//avoid divided by zero
		}
		return;
	}

	if (TYieldFlag == 0)
	{
		if (fabs(ts) <= sy)
		{
			Tload = ts*E0;
			Ttangent = E0;
		}
		else
		{
			TYieldFlag = 1;
			Tload = getEnvelopeStress(ts);
		}

		if (Tloading > 0 && dslip < 0.0 && Cslip > TmaxHSlip)
			TmaxHSlip = Cslip;
		if (Tloading < 0 && dslip > 0.0 && Cslip < TminHSlip)
			TminHSlip = Cslip;

		return;
	}

	//set limits
	maxrs = TmaxHSlip;
	maxrl = getEnvelopeStress(maxrs);
	if (maxrs > sy)
	{
		Eun = E0/pow(maxrs/sy,kkk);
	}
	else
	{
		Eun = E0;
	}
	maxvgs = maxrs-maxrl/Eun;

	minrs = TminHSlip;
	minrl = getEnvelopeStress(minrs);
	if (minrs < -sy)
	{
		Ere = E0/pow(-minrs/sy,kkk);
	}
	else
	{
		Ere = E0;
	}
	minvgs = minrs-minrl/Ere;

	rslip = TRSlip;
	rload = TRLoad;
	Er = TRSlope;
	rsvg = rslip-rload/Er;

	//get load for a slip
	if (Tloading > 0.0)					//(Cloading>0): was loading (tension)
	{
		if (dslip > 0.0)					//%(dslip>0): keep loading (tension)
		{
			if (ts >= maxrs)					//%%envelope stress
			{
				Tload = getEnvelopeStress(ts);
			}
			else								//%% ts < maxrs
			{
				if (ts >= rsvg)					//%%curve reloading (tension)
				{
					Sy = maxrl/Ere;
					suy = (maxrs-rsvg)/Sy;
					if (ts-rsvg <= Ks*(maxrs-rsvg))
					{
						ssy = (ts-rsvg)/Sy; 
						ss = ssy/(suy-ssy);
						ft = ss/pow((pow(1/suy,R)+pow(ss,R)),(1/R));
						Tload = ft*maxrl;
						st = pow(suy,(1-R))/pow(suy-ssy,2)/pow((pow(1/suy,R)+pow(ss,R)),(1/R+1));
						Ttangent = st*Ere;
					}
					else
					{
						ssy = Ks*(maxrs-rsvg)/Sy; 
						ss = ssy/(suy-ssy);
						ft = ss/pow((pow(1/suy,R)+pow(ss,R)),(1/R));
						templ = ft*maxrl;
						temps = Ks*(maxrs-rsvg)+rsvg;
						Tload = templ+(maxrl-templ)/(maxrs-temps)*(ts-temps);
						Ttangent = (maxrl-templ)/(maxrs-temps);
					}
				}
				else							//%% ts < rsvg
				{
					if (ts >= rslip)			//%%linear reloading (tension)
					{
						Tload = rload*(ts-rsvg)/(rslip-rsvg);
						Ttangent = Er;
					}
					else						//%%ts > rsvg not possible for safety
					{
						Tload = rload;
					}
                }
			}
		}
		else								//%(dslip<0): turn unloading (compression)
		{
			Tloading = -1;
			TRSlip = Cslip;
			TRLoad = Cload;
			if (Cslip > TmaxHSlip)
			{
				TmaxHSlip = Cslip;
	            if (TmaxHSlip > sy)
				{
					TRSlope = E0/pow(TmaxHSlip/sy, kkk);
				}
				else
				{
					TRSlope = E0;
				}
			}
			Er = TRSlope;

			trs = TRSlip;
			trl = TRLoad;
			trsvg = trs-trl/Er;
        
			if (trsvg > 0.0 && CminHSlip == 0.0)	//%%guess a minimum return slip
				CminHSlip = -sy*trsvg/su;

			if (ts > trs)						//%%not possible for safety
			{
				Tload = trl;
			}
			else
			{
				if (ts > trsvg)					//%%linear unloading (compression)
				{
					Tload = trl*(ts-trsvg)/(trs-trsvg);
					Ttangent = Er;
				}
				else
				{
					if (ts > minrs)				//%%curve unloading (compression)
					{
						Sy = minrl/Eun;
						suy = (minrs-trsvg)/Sy;
						if (ts-trsvg >= Ks*(minrs-trsvg))
						{
							ssy = (ts-trsvg)/Sy; 
							ss = ssy/(suy-ssy);
							ft = ss/pow((pow(1/suy,R)+pow(ss,R)),(1/R));
							Tload = ft*minrl;
							st = pow(suy,(1-R))/pow(suy-ssy,2)/pow((pow(1/suy,R)+pow(ss,R)),(1/R+1));
							Ttangent = st*Eun;
						}
					else
						{
							ssy = Ks*(minrs-trsvg)/Sy; 
							ss = ssy/(suy-ssy);
							ft = ss/pow((pow(1/suy,R)+pow(ss,R)),(1/R));
							templ = ft*minrl;
							temps = Ks*(minrs-trsvg)+trsvg;
							Tload = templ+(minrl-templ)/(minrs-temps)*(ts-temps);
							Ttangent = (minrl-templ)/(minrs-temps);
						}
					}
					else						//%%envelope stress
					{
						Tload = getEnvelopeStress(ts);
					}
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
				if (ts <= rsvg)					//%%curve unloading (compression)
				{
					Sy = minrl/Eun;
					suy = (minrs-rsvg)/Sy;
					if (ts-rsvg >= Ks*(minrs-rsvg))
					{
						ssy = (ts-rsvg)/Sy; 
						ss = ssy/(suy-ssy);
						ft = ss/pow((pow(1/suy,R)+pow(ss,R)),(1/R));
						Tload = ft*minrl;
						st = pow(suy,(1-R))/pow(suy-ssy,2)/pow((pow(1/suy,R)+pow(ss,R)),(1/R+1));
						Ttangent = st*Eun;
					}
					else
					{
						ssy = Ks*(minrs-rsvg)/Sy; 
						ss = ssy/(suy-ssy);
						ft = ss/pow((pow(1/suy,R)+pow(ss,R)),(1/R));
						templ = ft*minrl;
						temps = Ks*(minrs-rsvg)+rsvg;
						Tload = templ+(minrl-templ)/(minrs-temps)*(ts-temps);
						Ttangent = (minrl-templ)/(minrs-temps);
					}
				}
				else							//%% ts > rsvg
				{
					if (ts <= rslip)			//%%linear unloading (compression)
					{
						Tload = rload*(ts-rsvg)/(rslip-rsvg);
						Ttangent = E0;
					}
					else						//%%not possible for safety
					{
                    Tload = rload;
					}
				}
			}
		}
		else								//(dslip>0): turn loading (tension)
		{
			Tloading = 1;
			TRSlip = Cslip;
			TRLoad = Cload;
			if (Cslip < TminHSlip)
			{
				TminHSlip = Cslip;
	            if (TminHSlip < -sy)
				{
					TRSlope = E0/pow(-TminHSlip/sy, kkk);
				}
				else
				{
					TRSlope = E0;
				}
			}
			Er = TRSlope;

			trs = TRSlip;
			trl = TRLoad;
			trsvg = trs-trl/Er;
        
			if (trsvg < 0.0 && CmaxHSlip == 0.0)
				CmaxHSlip = sy*trsvg/(-su);

			if (ts < trs)						//%%not possible for safety
			{
				Tload = trl;
			}
			else								//%% ts > trs
			{
				if (ts < trsvg)					//%%linear unloading (compression)
				{
					Tload = trl*(ts-trsvg)/(trs-trsvg);
					Ttangent = Er;
				}
				else							//%% ts > trsvg
				{
					if (ts < maxrs)				//%%curve reloading (compression)
					{
						Sy = maxrl/Ere;
						suy = (maxrs-rsvg)/Sy;
						if (ts-trsvg <= Ks*(maxrs-trsvg))
						{
							ssy = (ts-trsvg)/Sy; 
							ss = ssy/(suy-ssy);
							ft = ss/pow((pow(1/suy,R)+pow(ss,R)),(1/R));
							Tload = ft*maxrl;
							st = pow(suy,(1-R))/pow(suy-ssy,2)/pow((pow(1/suy,R)+pow(ss,R)),(1/R+1));
							Ttangent = st*Ere;
						}
						else
						{
							ssy = Ks*(maxrs-trsvg)/Sy; 
							ss = ssy/(suy-ssy);
							ft = ss/pow((pow(1/suy,R)+pow(ss,R)),(1/R));
							templ = ft*maxrl;
							temps = Ks*(maxrs-trsvg)+trsvg;
							Tload = templ+(maxrl-templ)/(maxrs-temps)*(ts-temps);
							Ttangent = (maxrl-templ)/(maxrs-temps);
						}
					}
					else						//%%envelope stress
					{
						Tload = getEnvelopeStress(ts);
					}
				}
			}
			if (Cslip < TminHSlip)
			{
				TminHSlip = Cslip;
			}
        }
    }

	// Detect load reversals and determine damage indexes
	//this->detectloadReversal (dslip);
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

int Bond_SP01::sendSelf (int commitTag, Channel& theChannel)
{
   return -1;
}

int Bond_SP01::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   return -1;
}

void Bond_SP01::Print (OPS_Stream& s, int flag)
{
   s << "Bond_SP01 tag: " << this->getTag() << endln;
   s << "  sy: " << sy << " ";
   s << "  fy: " << fy << " ";
   s << "  su: " << su << " ";
   s << "  fu: " << fu << " ";
   s << "  Kz: " << Kz << " ";
   s << "  R: " << R << " ";
   s << "  Cd: " << Cd << " ";
   s << "  db: " << db << " ";
   s << "  fc: " << fc << " ";
   s << "  lba:" << lba << " ";
}



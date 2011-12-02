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
// $Date: 2006-08-03 23:44:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/CloughHenry.cpp,v $
//
//
// CloughHenry.cpp: implementation of the CloughHenry class from Fortran version.
// Originally from SNAP PROGRAM by Prof H.K. Krawinkler
//
// Written: Arash Altoontash, Gregory Deierlein, 12/01
// Revised: 03/02
//
// Purpose: This file contains the implementation for the CloughHenry class.
//
//////////////////////////////////////////////////////////////////////

#include <CloughHenry.h>
#include <stdlib.h>
#include <Channel.h>
#include <math.h>

#include <string.h>

#define DEBG 0

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CloughHenry::CloughHenry(int tag, Vector inputParam )
:UniaxialMaterial(tag,0)
{
  if( (inputParam.Size()) < 16) 
    opserr << "Error: CloughHenry(): inputParam, size <16\n" << "\a";
  
  /*
    Input parameters
    
    elstk      =  initial elastic stiffness
    fyieldPos  =  positive yield strength	fyieldNeg  =  yield strength in compression
    alpha      =  strain hardening ratio (fraction of elstk)
    elstk      =  initial elastic stiffness
    Resfac	   =  residual stress after collapse
    dyieldPos  =  positive yield displacement
    dyieldNeg  =  negative yield displacement
    
    ecaps, ecapd, ecapk, ecapa = parameter expressing the hystetic
    energy dissipacion capacity.
    Enrgts, Enrgtd, Enrgtk, Enrgta = hysteretic energy dissipation capacity
  */
  
  elstk		= inputParam[0];
  fyieldPos	= inputParam[1];
  fyieldNeg	= inputParam[2];
  alpha		= inputParam[3];
  Resfac		= inputParam[4];
  capSlope	= inputParam[5];
  capDispPos	= inputParam[6];
  capDispNeg	= inputParam[7];
  ecaps		= inputParam[8];
  ecapk		= inputParam[9];
  ecapa		= inputParam[10];
  ecapd		= inputParam[11];
  cs			= inputParam[12];
  ck			= inputParam[13];
  ca			= inputParam[14];
  cd			= inputParam[15];
  
  // Error check
  
  if ( ecaps < 0.0 || ecapk < 0.0 || ecapa < 0.0 || ecapd < 0.0)
    opserr << "Error: CloughHenry::CloughHenry  : All gamma values must be >= 0\n" << "\a";		
  
  if ( cs < 0.0 || ck < 0.0 || ca < 0.0 || cd < 0.0 )
    opserr << "Error: CloughHenry::CloughHenry  : All 'c' values must be >= 0\n" << "\a";
  
  if ( capSlope > 0.0 )
    opserr << "Error: CloughHenry::CloughHenry  : CapSlope must be < 0\n" << "\a";
  
  if ( Resfac <  0.0  || Resfac > 1.0)
    opserr << "Error: CloughHenry::CloughHenry  : Residual must be > 0 and <= 1\n" << "\a";
  
  if ( alpha > 0.8 || alpha < -0.8 )
    opserr << "Error: CloughHenry::CloughHenry  : alpha must be < 0.8 and > -0.8\n" << "\a";	
  
  if ( alpha == capSlope )
    opserr << "Error: CloughHenry::CloughHenry  : Error: alpha Hard. can not be equal to alphaCap\n" << "\a";	
  
  // Initialize history data
  this->revertToStart();
  
}


CloughHenry::CloughHenry()
  :UniaxialMaterial(0,0) 
{
}



CloughHenry::~CloughHenry()
{

}


int CloughHenry::revertToStart()
{

	dyieldPos = fyieldPos/elstk;
	dyieldNeg = fyieldNeg/elstk;
	Enrgts = fyieldPos*dyieldPos*ecaps;
	Enrgtk = fyieldPos*dyieldPos*ecapk;
	Enrgta = fyieldPos*dyieldPos*ecapa;
	Enrgtd = fyieldPos*dyieldPos*ecapd;
	
	double ekhard = elstk*alpha;
	double fPeakPos = fyieldPos+ekhard*(capDispPos-dyieldPos);
	double fPeakNeg = fyieldNeg+ekhard*(capDispNeg-dyieldNeg);	

	hsTrial[0] = 0.0;									// d
	hsTrial[1] = 0.0;									// f
	hsTrial[2]  = elstk;								// ek
	hsTrial[3]  = elstk;								// ekunload
	hsTrial[4]  = elstk;								// ekexcurs
	hsTrial[5]  = 0.0;									// Enrgtot
	hsTrial[6]  = 0.0;									// Enrgc
	hsTrial[7]  = 0.0;									// sp
	hsTrial[8]  = 0.0;									// sn
	hsTrial[9]  = 0.0;									// kon
	hsTrial[10] = dyieldPos;							// dmax
	hsTrial[11] = dyieldNeg;							// dmin
	hsTrial[12] = fyieldPos;							// fyPos
	hsTrial[13] = fyieldNeg;							// fyNeg
	hsTrial[14] = capDispPos;							// cpPos
	hsTrial[15] = capDispNeg;							// cpNeg
	hsTrial[16] = 0.0;									// dlstPos
	hsTrial[17] = 0.0;									// flstPos
	hsTrial[18] = 0.0;									// dlstNeg
	hsTrial[19] = 0.0;									// flstNeg
	hsTrial[20] = alpha;								// alphaPos
	hsTrial[21] = alpha;								// alphaNeg
	hsTrial[22] = -capSlope*elstk*capDispPos+fPeakPos;	// fCapRefPos
	hsTrial[23] = -capSlope*elstk*capDispNeg+fPeakNeg;	// fCapRefNeg : indicates cap reference point
		
	for( int i=0 ; i<24; i++) {
		hsCommit[i] = hsTrial[i];
		hsLastCommit[i] = hsTrial[i];
	}

	return 0;
}


void CloughHenry::Print(OPS_Stream &s, int flag)
{
	s << "BondSlipMaterial Tag: " << this->getTag() << endln;
	s << "D : " << hsTrial[0] << endln;
	s << "F : " << hsTrial[1] << endln;
	s << "EK: " << hsTrial[2]  << endln;
	
	s << endln;
}


int CloughHenry::revertToLastCommit()
{
  
  for(int i=0; i<24; i++) hsTrial[i] = hsLastCommit[i];
  
  return 0;
}


double CloughHenry::getTangent()
{
  if (hsTrial[1] < 0)
    return hsTrial[2];
  else
    return 0;
}

double CloughHenry::getInitialTangent (void)
{
  return elstk;
}

double CloughHenry::getStress()
{
  if (hsTrial[1] < 0)
    return hsTrial[1];
  else
    return 0;
}


double CloughHenry::getStrain (void)
{

	return hsTrial[0];
}


int CloughHenry::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
	return 0;
}


int CloughHenry::sendSelf(int cTag, Channel &theChannel)
{
	return 0;
}


UniaxialMaterial *CloughHenry::getCopy(void)
{
	Vector inp(16);
	
	inp[0]  = 	elstk;
	inp[1]  =  	fyieldPos;
	inp[2]  =  	fyieldNeg;
	inp[3]  = 	alpha;
	inp[4]  = 	Resfac;
	inp[5]  = 	capSlope;
	inp[6]  = 	capDispPos;
	inp[7]  = 	capDispNeg;
	inp[8]  = 	ecaps;
	inp[9]  = 	ecapk;
	inp[10] = 	ecapa;
	inp[11] = 	ecapd;
	inp[12] =	cs;
	inp[13] = 	ck;
	inp[14] = 	ca;
	inp[15] = 	cd;

  
  CloughHenry *theCopy = new CloughHenry(this->getTag(), inp );
  
  for (int i=0; i<24; i++) {
    theCopy->hsTrial[i] = hsTrial[i];
    theCopy->hsLastCommit[i] = hsLastCommit[i];
  }
  
  return theCopy;
}


int CloughHenry::setTrialStrain( double d, double strainRate)
{

	int kon,idiv,Unl,flagDeg,flgstop;

	double ekP,ek,ekunload,sp,sn,dP,fP,f;
	double deltaD,err,f1,f2,dmax,dmin,fmax,fmin,ekt;
	double betas,betak,betaa;
	double Enrgtot,Enrgc,Enrgi,dyPos,dyNeg;
	double fyPos,fyNeg;
	double betad,cpPos,cpNeg,a2;
	double dlstPos,flstPos,dlstNeg,flstNeg,ekc,RSE,tst,ekexcurs;
	double dCap1Pos,dCap2Pos,dCap1Neg,dCap2Neg;
	double ekhardNeg,ekhardPos,alphaPos,alphaNeg;
	double fCapRefPos,fCapRefNeg;
	
	//	Relationship between basic variables and hsTrial array
	idiv = 1;
	flagDeg = 0;
	flgstop = 0;
	Unl = 1;
	err = 0.0;

	dP			= hsLastCommit[0];
	fP			= hsLastCommit[1];
	ekP			= hsLastCommit[2];
	ekunload	= hsLastCommit[3];
	ekexcurs	= hsLastCommit[4];
	Enrgtot		= hsLastCommit[5];
	Enrgc		= hsLastCommit[6];
	sp			= hsLastCommit[7];
	sn			= hsLastCommit[8];
	kon			= (int) hsLastCommit[9];
	dmax		= hsLastCommit[10];
	dmin		= hsLastCommit[11];
	fyPos		= hsLastCommit[12];
	fyNeg		= hsLastCommit[13];
	cpPos		= hsLastCommit[14];
	cpNeg		= hsLastCommit[15];
	dlstPos		= hsLastCommit[16];
	flstPos		= hsLastCommit[17];
	dlstNeg		= hsLastCommit[18];
	flstNeg		= hsLastCommit[19];
	alphaPos	= hsLastCommit[20];
	alphaNeg	= hsLastCommit[21];
	fCapRefPos	= hsLastCommit[22];
	fCapRefNeg	= hsLastCommit[23];

	ekhardPos = alphaPos * elstk;
	ekhardNeg = alphaNeg * elstk;
	deltaD = d-dP;

		betas = betak = betaa = betad = 0.0;
	//	Loop to initialize parameters 
	
	if ( kon == 0 ) {
		if( deltaD >= 0.0) {
			kon = 1;
		}
		else {
			kon = 2;
		}
	}
	
	//	STARTS BIG LOOP
	//	Positive Delta indicating loading
	
	if ( deltaD >= 0.0 ) {
		
		if ( kon==2 ) {
			kon = 1;
			Unl = 0;
			RSE = 0.5 * fP * fP / ekunload;
			if ( (Enrgc-RSE) <= 0.0 || (Enrgtk-(Enrgtot-RSE)) <0.0) RSE = 0.0;
			a2 = Enrgtk - ( Enrgtot - RSE );			
			if( a2 <= 0.0 && Enrgtk != 0.0)
				opserr << "Warning: CloughHenry::SetTrial  : Maximum energy capacity has been reached for stiffness degradation\n" << "\a";	

			if( ecapk != 0.0) {
				betak = pow ( ((Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE))) , ck );
				ekunload = ekexcurs * ( 1 - betak );
				if( ekunload <= ekhardNeg ) 
					opserr << "Warning: CloughHenry::SetTrial  : Maximum energy capacity has been reached for stiffness degradation\n" << "\a";	
			}
			
			//	Determination of sn according to the hysteresis status
			if( ekunload <= 1.e-7 )
				opserr << "Warning: CloughHenry::SetTrial  : Total stiffness loss\n" << "\a";	
			
			if( fP < 0.0) {
				tst = dP - fP / ekunload;
				if( fabs(dmax-dyieldPos) >= 1.e-10 && fabs(tst) <= 1.e-10) {
					sn = 1.e-9;
				}
				else {
					sn = dP - fP / ekunload;
				}
			}
	    if( fabs(dmin-dP) <= 1.e-10 ) sp = sn + 1.0e-10;
		}
		
		//	LOADING
		//	Push envelope
		
		if ( d >= dmax ) {
			this->envelPosCap (fyPos,alphaPos,capSlope,cpPos,d,&f,&ek);
			dmax = d;
			fmax = f;
			dlstPos = dmax+ 1.0e-10;
			flstPos = f;
		}
		else if ( fabs(sn) > 1.0e-10) {
			this->envelPosCap (fyPos,alphaPos,capSlope,cpPos,dmax,&fmax,&ekt);
			if ( d<=sn) {
				ek = ekunload;
				f  = fP+ek*deltaD;
				if ( Unl == 0 && fabs(ek-ekP) > 1.0e-10 && dP != dmin) {
					dlstNeg = dP;
					flstNeg = fP;
				}
			}
			else {
				ek = fmax/(dmax-sn);
				
				if (ek>=ekunload) {
					flgstop = 1;
					opserr << "Unloading stiffness < reloading stiff";
				}
				
				f2 = ek *( d - sn );
				if( dlstPos > sn && dlstPos < dmax ) {
					ekc = flstPos / ( dlstPos - sn );
					if( ekc > ek && flstPos < fmax ) {
						
						if( d < dlstPos ) {
							ek = flstPos/(dlstPos-sn);
							f2 = ek*(d-sn);
						}
						else {
							ek = (fmax-flstPos)/(dmax-dlstPos);
							f2 = flstPos+ek*(d-dlstPos);
						}
					}
				}
				
				f1 = fP + ekunload * deltaD;
				f = (f1 < f2) ? f1 : f2;
				
				if ( fabs(f-f1) < 1.0e-10 ) ek=ekunload;
			}
		}
		else {
			if ( d > 0.0 ) {
				this->envelPosCap (fyPos,alphaPos,capSlope,cpPos,d,&f,&ek);
			}
			else {
				this->envelNegCap (fyNeg,alphaNeg,capSlope,cpNeg,d,&f,&ek);
			}
		}
	}
	else {
		if (kon==1) {
			kon = 2;
			Unl = 0;
			RSE = 0.5 * fP * fP / ekunload;
			if ( (Enrgc-RSE) <= 0.0 || (Enrgtk-(Enrgtot-RSE)) <0.0 ) RSE = 0.0;

			a2 = Enrgtk - ( Enrgtot - RSE );

			if( a2 <= 0.0 && Enrgtk != 0.0) flgstop=1;
			
			if(ecapk != 0.0) {
				betak = pow ( ((Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE))) , ck);
				ekunload = ekexcurs*(1-betak);
				if(ekunload<=ekhardPos) flgstop=1;
			}
			
			//	Determination of sn according to the hysteresis status
			
			if( fP > 0.0 ) {
				tst = dP-fP/ekunload;
				if( fabs(dmin-dyieldNeg) >= 1.e-10 && fabs(tst) <= 1.e-10 ) {
					sp = 1.e-9;
				}
				else {
					sp = dP-fP/ekunload;
				}
			}
			
			if( fabs(dmax-dP) <= 1.e-10 ) sn=sp-1.0e-10;
		}
		
		//	UNLOADING
		//	Push envelope
		
		if (d <= dmin) {
			this->envelNegCap (fyNeg,alphaNeg,capSlope,cpNeg,d,&f,&ek);
			dmin = d;
			fmin = f;
			dlstNeg = dmin - 1.0e-10;
			flstNeg = f;
			
		}
		else if ( fabs(sp) > 1.0e-10) {
			this->envelNegCap (fyNeg,alphaNeg,capSlope,cpNeg,dmin,&fmin,&ekt);
			if ( d>=sp) {
				ek = ekunload;
				f = fP+ek*deltaD;
				if( Unl==0 && fabs(ek-ekP) > 1.0e-10 && dP != dmax ) {
					dlstPos = dP;
					flstPos = fP;
				}
			}
			else {
				ek = fmin/(dmin-sp);
				
				if ( ek >= ekunload ) {
					flgstop = 1;
					opserr << "Unloading stiffness < reloading stiff\n";
				}
				
				f2 = ek * ( d - sp );
				if( dlstNeg < sp && dlstNeg > dmin ) {
					ekc = flstNeg / ( dlstNeg - sp );
					
					if( ekc > ek && flstNeg > fmin ) {
						if( d > dlstNeg ) {
							ek = flstNeg / ( dlstNeg - sp );
							f2 = ek * ( d - sp );
						}
						else {
							ek = ( fmin - flstNeg ) / ( dmin - dlstNeg );
							f2 = flstNeg + ek * ( d - dlstNeg );
						}
					}
				}
				
				f1 = fP + ekunload * deltaD;
				f = (f1 > f2) ? f1 : f2;
				if ( fabs(f-f1) < 1.e-10 ) ek=ekunload;
			}
		}
		else {
			if ( d > 0.0 ) {
				this->envelPosCap (fyPos,alphaPos,capSlope,cpPos,d,&f,&ek);
			}
			else {
				this->envelNegCap (fyNeg,alphaNeg,capSlope,cpNeg,d,&f,&ek);
			}
		}
	}
	
	//	ENDS BIG LOOP ---------------------------------------------------
	
	//	Flag to deteriorate parameters on the opposite side of the loop --	

	if( f * fP < 0.0) {
		if( fP >0.0 && dmax >dyieldPos) flagDeg = 1;	// positive half cycle
		if( fP < 0.0 && dmin < dyieldNeg) flagDeg = 2;	// negative half cycle 
	}	

	//	ENERGY CALCULATIONS ---------------------------------------------

		Enrgi = 0.5 * ( f + fP ) * deltaD;
		Enrgc = Enrgc + Enrgi;
		Enrgtot = Enrgtot + Enrgi;
	
	//	UPDATING OF DETERIORATION PARAMETERS ----------------------------
	//	Check the remaining capacity of the system
	
	if( flagDeg == 1 || flagDeg == 2 ) {
		if( ( Enrgtot >= Enrgts && Enrgts != 0.0) || ( Enrgtot >= Enrgtk && Enrgtk != 0.0) || 
			( Enrgtot >= Enrgta && Enrgta != 0.0) || ( Enrgtot >= Enrgtd && Enrgtd != 0.0))
			opserr << "Total Energy greater than capacity\n";
		
		//	Update beta values for strength, acc. stiff. and capping

		if( ecaps != 0.0 ) betas = pow ((Enrgc/(Enrgts-Enrgtot)) , cs );
		if( betas>=1.0 ){
			opserr << "Warning: CloughHenry::SetTrial  : Total Strength loss\n" << "\a";	
			betas = 1.0;
		}

		if( ecapa != 0.0 ) betaa = pow ((Enrgc/(Enrgta-Enrgtot)) , ca );
		if( betaa>=1.0 ){
			opserr << "Warning: CloughHenry::SetTrial  : Total accelerated stiffness loss\n" << "\a";
			betaa = 1.0;
		}

		if( ecapd != 0.0 ) betad = pow ((Enrgc/(Enrgtd-Enrgtot)) , cd );
		if( betad>=1.0){
			opserr << "Warning: CloughHenry::SetTrial  : Total capping loss\n" << "\a";	
			betad = 1.0;
		}
	
		//	Update values for the next half cycle
		ekexcurs = ekunload;
		Enrgc = 0.0;
		
		//	Deteriorate parameters for the next half cycle
		
		if( deltaD < 0.0 ) {
			
			fyNeg = fyNeg * ( 1 - betas );
			alphaNeg = alphaNeg * ( 1 - betas );
			fCapRefNeg = fCapRefNeg * ( 1 - betad );
			dmin = dmin * ( 1 + betaa );
			
			dyNeg = fyNeg / elstk;
			ekhardNeg = alphaNeg * elstk;
			
			dCap1Neg = fCapRefNeg/(elstk-capSlope*elstk);
			dCap2Neg = (fCapRefNeg+ekhardNeg*dyNeg-fyNeg)/(ekhardNeg-capSlope*elstk);
			cpNeg = ( dCap1Neg < dCap2Neg ) ? dCap1Neg : dCap2Neg;
		}
		else {
			fyPos = fyPos * ( 1 - betas );
			alphaPos = alphaPos * ( 1 - betas );
			fCapRefPos = fCapRefPos * ( 1 - betad );
			dmax = dmax * ( 1 + betaa );

			dyPos = fyPos / elstk;
			ekhardPos = alphaPos * elstk;

			dCap1Pos = fCapRefPos / ( elstk - capSlope * elstk );
			dCap2Pos = (fCapRefPos+ekhardPos*dyPos-fyPos)/(ekhardPos-capSlope*elstk);
			cpPos= (dCap1Pos > dCap2Pos ) ? dCap1Pos : dCap2Pos;
		}
		flagDeg = 0;
	}
	
	// Relationship between basic variables and hsTrial array	for next cycle
	
	hsTrial[0] = d;
	hsTrial[1] = f;
	hsTrial[2] = ek;
	hsTrial[3] = ekunload;
	hsTrial[4] = ekexcurs;
	hsTrial[5] = Enrgtot;
	hsTrial[6] = Enrgc;
	hsTrial[7] = sp;
	hsTrial[8] = sn;
	hsTrial[9] = (double) kon;
	hsTrial[10] = dmax;
	hsTrial[11] = dmin;
	hsTrial[12] = fyPos;
	hsTrial[13] = fyNeg;
	hsTrial[14] = cpPos;
	hsTrial[15] = cpNeg;
	hsTrial[16] = dlstPos;
	hsTrial[17] = flstPos;
	hsTrial[18] = dlstNeg;
	hsTrial[19] = flstNeg;
	hsTrial[20] = alphaPos;
	hsTrial[21] = alphaNeg;
	hsTrial[22] = fCapRefPos;
	hsTrial[23] = fCapRefNeg;
	
	return 0;
}


int CloughHenry::commitState()
{
  int i;
  for( i=0; i<24; i++ ) hsLastCommit[i] = hsTrial[i];
  
  this->recordInfo();
  
  return 0;
}


void CloughHenry::recordInfo(int cond )
{

}


void CloughHenry::envelPosCap( double fy, double alphaPos, double alphaCap,
			    double cpDsp, double d, double *f, double *ek )
{
  
  double dy,Res,rcap,dres;
  
  dy = fy / elstk;
  
  if (dy < cpDsp)
    {
      Res = Resfac * fyieldPos;
      rcap = fy+alphaPos * elstk * ( cpDsp - dy );
      dres = cpDsp + ( Res - rcap ) / ( alphaCap * elstk );
      
      if ( d < 0.0 )  
	{
	  *f = 0.0;
	  *ek = 0.0;
	}
      else
	{
	  if ( d <= dy )
	    {
	      *ek = elstk;
	      *f = (*ek) * d;
	    }
	  else
	    {
	      if( d <= cpDsp )
		{
		  *ek = elstk * alphaPos;
		  *f = fy + (*ek) * ( d - dy );
		}
	      else
		{
		  if( d <=  dres )
		    {
		      *ek = alphaCap * elstk;
		      *f = rcap + (*ek) * ( d - cpDsp );
		    }
		  else
		    {
		      *ek = 0;
		      // ek = 1.e-10;
		      *f = Res + d * (*ek);
		    }
		}
	    }
	}
    }
  else
    {
      rcap = elstk * cpDsp;
      Res = Resfac * rcap;
      dres = cpDsp + ( Res - rcap ) / ( alphaCap * elstk );
      
      if ( d < 0.0 ) 
	{
	  *f = 0.0;
	  *ek = 0.0;
	}
      else
	{
	  if( d <= cpDsp )
	    {
	      *ek = elstk;
	      *f = (*ek) * d;
	    }
	  else
	    {
	      if( d <= dres )
		{
		  *ek = alphaCap * elstk;
		  *f = rcap + (*ek) * ( d - cpDsp );
		}
	      else
		{
		  *ek = 0;
		  // ek = 1.e-10;
		  *f = Res + d * (*ek);
		}
	    }
	}
    }
  
  return;
}


void CloughHenry::envelNegCap( double fy, double alphaNeg, double alphaCap,
			    double cpDsp, double d, double *f, double *ek)
{
  
  double dy,Res,rcap,dres;
  
  dy = fy / elstk;
  
  if( dy > cpDsp )
    {
      
      Res = Resfac * fyieldNeg;
      rcap = fy + alphaNeg * elstk * ( cpDsp - dy );
      dres = cpDsp + ( Res - rcap ) / ( alphaCap * elstk );
      
      if (d > 0.0) 
	{
	  *f = 0.0;
	  *ek = 0.0;
	}
      else
	{
	  if ( d >= dy )
	    {
	      *ek = elstk;
	      *f = (*ek) * d;
	    }
	  else
	    {
	      if ( d >= cpDsp )
		{
		  *ek = elstk * alphaNeg;
		  *f = fy + (*ek) * ( d - dy );
		}
	      else
		{
		  if ( d >= dres )
		    {
		      *ek = elstk * alphaCap;
		      *f = rcap + (*ek) * ( d - cpDsp );
		    }
		  else
		    {
		      *ek = 0;
		      // *ek = 1.d-10
		      *f = Res + (*ek) * d;
		    }
		}
	    }
	}
    }
  else
    {
      rcap = elstk * cpDsp;
      Res = Resfac * rcap;
      dres = cpDsp + ( Res - rcap ) / ( alphaCap * elstk );
      
      if (d > 0.0) 
	{
	  *f = 0.0;
	  *ek = 0.0;
	}
      else
	{
	  if ( d >= cpDsp )
	    {
	      *ek = elstk;
	      *f = (*ek) * d;
	    }
	  else
	    {
	      if ( d >= dres )
		{
		  *ek = elstk * alphaCap;
		  *f = rcap + (*ek) * ( d - cpDsp );
		}
	      else
		{
		  *ek = 0;
		  // *ek = 1.e-10;
		  *f = Res + (*ek) * d;
		}
	    }
	}
    }
  return;
}

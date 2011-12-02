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
                                                                        
// $Revision: 1.4 $
// $Date: 2004-09-01 03:53:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/Pinching.cpp,v $
//
//
// Pinching.cpp: implementation of the Pinching class from Fortran version.
// Originally from SNAP PROGRAM by Prof H.K. Krawinkler
//
// Written: Arash Altoontash, Gregory Deierlein, 12/01
// Revised: 05/03
//
// Purpose: This file contains the implementation for the Pinching class.
//
//////////////////////////////////////////////////////////////////////

#include <Pinching.h>
#include <stdlib.h>
#include <Channel.h>
#include <math.h>

#include <string.h>

#define DEBG 0

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Pinching::Pinching(int tag, Vector inputParam )
  :UniaxialMaterial(tag,MAT_TAG_Pinching)
{
	if( (inputParam.Size()) < 19) 
		opserr << "Error: Pinching(): inputParam, size <19\n" << "\a";
  
	/*
 	Input parameters

	elstk      =  initial elastic stiffness
	fyieldPos  =  positive yield strength
	fyieldNeg  =  yield strength in compression
	alpha      =  strain hardening ratio (fraction of elstk)
	elstk      =  initial elastic stiffness
	Resfac	   =  residual stress after collapse
	dyieldPos  =  positive yield displacement
	dyieldNeg  =  negative yield displacement

	ecaps, ecapd, ecapk, ecapa = parameter expressing the hystetic
		    energy dissipacion capacity.
	Enrgts, Enrgtd, Enrgtk, Enrgta = hysteretic energy dissipation
			capacity
	*/
	elstk			= inputParam[0];
	fyieldPos		= inputParam[1];
	fyieldNeg		= inputParam[2];
	alpha			= inputParam[3];
	Resfac			= inputParam[4];
	capSlope		= inputParam[5];
	capDispPos		= inputParam[6];
	capDispNeg		= inputParam[7];
	fpPos			= inputParam[8];
	fpNeg			= inputParam[9];
	a_pinch			= inputParam[10];
	ecaps			= inputParam[11];
	ecapk			= inputParam[12];
	ecapa			= inputParam[13];
	ecapd			= inputParam[14];
	cs				= inputParam[15];
	ck				= inputParam[16];
	ca				= inputParam[17];
	cd				= inputParam[18];
  
	// Error check
	if ( ecaps < 0.0 || ecapk < 0.0 || ecapa < 0.0 || ecapd < 0.0)
    {
		opserr << "Error: Pinching::Pinching  : All gamma values must be >= 0\n" << "\a";
    }		
	
	if ( cs < 0.0 || ck < 0.0 || ca < 0.0 || cd < 0.0 )
    {
		opserr << "Error: Pinching::Pinching  : All 'c' values must be >= 0\n" << "\a";
    }
	
	if ( capSlope > 0.0 )
    {
		opserr << "Error: Pinching::Pinching  : CapSlope must be < 0\n" << "\a";
    }		
	
	if ( Resfac <  0.0  || Resfac > 1.0)
    {
		opserr << "Error: Pinching::Pinching  : Residual must be > 0 and <= 1\n" << "\a";
    }		
	
	if ( a_pinch < 0.0 || a_pinch > 1.0)
    {
		opserr << "Error: Pinching::Pinching  : kappad (dev. point)must be > 0 and <= 1\n" << "\a";
    }
	
	if ( alpha > 0.8 || alpha < -0.8 )
    {
		opserr << "Error: Pinching::Pinching  : alpha must be < 0.8 and > -0.8\n" << "\a";	
    }
	
	if ( alpha == capSlope )
    {
		opserr << "Error: Pinching::Pinching  : Error: alpha Hard. can not be equal to alphaCap\n" << "\a";	
    }
	
	if ( DEBG ==1 )
    {
		// Open an output file for debugging
		char FileName[20];						// debugging
		sprintf(FileName, "Pinching%d.out", tag);
		OutputFile = fopen ( FileName , "w" );	// debugging
		fprintf( OutputFile , "Constructor\n" );	// debugging
    }
	
	// Initialice history data
	this->revertToStart();
}


Pinching::Pinching()
  :UniaxialMaterial(0,MAT_TAG_Pinching) 
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Empty constructor\n" );	// debugging
}



Pinching::~Pinching()
{
  if ( DEBG ==1 )
    {
      fprintf( OutputFile , "Distructor\n" );		// debugging
      fclose( OutputFile );						// debugging
    }
}


int Pinching::revertToStart()
{
	dyieldPos = fyieldPos / elstk;
	dyieldNeg = fyieldNeg / elstk;

	Enrgts = fyieldPos * dyieldPos * ecaps;
	Enrgta = fyieldPos * dyieldPos * ecapa;
	Enrgtk = fyieldPos * dyieldPos * ecapk;
	Enrgtd = fyieldPos * dyieldPos * ecapd;
	
	double ekhard = elstk * alpha;

	double fPeakPos = fyieldPos + ekhard * ( capDispPos - dyieldPos );
	double fPeakNeg = fyieldNeg + ekhard * ( capDispNeg - dyieldNeg );
	
	hsTrial[0] = 0.0;			// d
	hsTrial[1] = 0.0;			// f
	hsTrial[2] = elstk;		// ek
	hsTrial[3] = elstk;		// ekunload
	hsTrial[4] = elstk;		// ekexcurs
	hsTrial[5] = 0.0;			// Enrgtot
	hsTrial[6] = 0.0;			// Enrgc
	hsTrial[7] = 0.0;			// sp
	hsTrial[8] = 0.0;			// sn
	hsTrial[9] = 0.0;			// kon
	hsTrial[10] = dyieldPos;	// dmax
	hsTrial[11] = dyieldNeg;	// dmin
	hsTrial[12] = fyieldPos;	// fyPos
	hsTrial[13] = fyieldNeg;	// fyNeg
	hsTrial[14] = capDispPos;	// cpPos
	hsTrial[15] = capDispNeg;	// cpNeg
	hsTrial[16] = fyieldPos;	// fmax
	hsTrial[17] = fyieldNeg;	// fmin
	hsTrial[18] = alpha;		// alphaPos
	hsTrial[19] = alpha;		// alphaNeg
	hsTrial[20] = -capSlope * elstk * capDispPos + fPeakPos;	// fCapRefPos
	hsTrial[21] = -capSlope * elstk * capDispNeg + fPeakNeg;	// fCapRefNeg
	
	for( int i=0 ; i<22; i++) {
		hsCommit[i] = hsTrial[i];
		hsLastCommit[i] = hsTrial[i];
	}
	
	return 0;
}


void Pinching::Print(OPS_Stream &s, int flag)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Print\n" );	// debugging
	s << "BondSlipMaterial Tag: " << this->getTag() << endln;
	s << "D : " << hsTrial[0] << endln;
	s << "F : " << hsTrial[1] << endln;
	s << "EK: " << hsTrial[2]  << endln;
	s << "Input Parameters:\n";
	s << endln;
}


int Pinching::revertToLastCommit()
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Revert to last commit\n" );	// debugging
	
	for(int i=0; i<22; i++) {
		hsTrial[i] = hsCommit[i];
		hsCommit[i] = hsLastCommit[i];
	}
	
	return 0;
}


double Pinching::getTangent()
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Get tangent\n" );	// debugging
	return hsTrial[2];
}

double Pinching::getInitialTangent (void)
{
	return elstk;
}

double Pinching::getStress()
{
	return hsTrial[1];
}


double Pinching::getStrain (void)
{
	return hsTrial[0];
}


int Pinching::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
	return 0;
}


int Pinching::sendSelf(int cTag, Channel &theChannel)
{
	return 0;
}


UniaxialMaterial *Pinching::getCopy(void)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Get copy\n" );	// debugging
	Vector inp(19);
	
	inp[0] = elstk;
	inp[1] = fyieldPos;
	inp[2] = fyieldNeg;
	inp[3] = alpha;
	inp[4] = Resfac;
	inp[5] = capSlope;
	inp[6] = capDispPos;
	inp[7] = capDispNeg;
	inp[8] = fpPos;
	inp[9] = fpNeg;
	inp[10] = a_pinch;
	inp[11] = ecaps;
	inp[12] = ecapk;
	inp[13] = ecapa;
	inp[14] = ecapd;
	inp[15] = cs;
	inp[16] = ck;
	inp[17] = ca;
	inp[18] = cd;
	
	Pinching *theCopy = new Pinching(this->getTag(), inp );

	for (int i=0; i<22; i++) {
		theCopy->hsTrial[i] = hsTrial[i];
		theCopy->hsCommit[i] = hsCommit[i];
		theCopy->hsLastCommit[i] = hsLastCommit[i];
	}
	
	return theCopy;
}


int Pinching::setTrialStrain( double d, double strainRate)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Set trial displacement\n" );	// debugging
	int flagDeg , flgstop , Unl, kon;
	double deltaD;
	double ekP,ekunload,ekexcurs,Enrgtot,Enrgc,sp,sn,dmax,dmin,fyPos,fyNeg,fP,dP;
	double ek,f,RSE,a2,betaa,betad,betak,betas,ekLim;
	double ekt,f1,f2,fpinch,dCap1Pos,dCap2Pos,dCap1Neg,dCap2Neg;
	double dyPos,dyNeg,dch,ekpinch, cpPos, cpNeg;
	double Enrgi,ekhardNeg,ekhardPos,fmax,fmin;
	double alphaPos, alphaNeg, fCapRefPos, fCapRefNeg;
	
	//	Updating of information in each new call to this Model ----------

	flagDeg = 0;
	flgstop = 0;
	Unl = 1;
	
  	dP			= hsCommit[0];
	fP			= hsCommit[1];
	ekP			= hsCommit[2];
 	ekunload	= hsCommit[3];
	ekexcurs	= hsCommit[4];
	Enrgtot		= hsCommit[5];
	Enrgc		= hsCommit[6];
	sp			= hsCommit[7];
	sn			= hsCommit[8];
	kon			= (int) hsCommit[9];
	dmax		= hsCommit[10];
	dmin		= hsCommit[11];
	fyPos		= hsCommit[12];
	fyNeg		= hsCommit[13];			   
	cpPos		= hsCommit[14];
	cpNeg		= hsCommit[15];
	fmax		= hsCommit[16];
	fmin		= hsCommit[17];
	alphaPos	= hsCommit[18];
	alphaNeg	= hsCommit[19];
	fCapRefPos	= hsCommit[20];
	fCapRefNeg	= hsCommit[21];
	
	ekhardPos	= alphaPos * elstk;
	ekhardNeg	= alphaNeg * elstk;
	deltaD		= d-dP;

	betas = betak = betaa = betad = 0.0;

	//	Loop to initialize parameters 
	if ( kon == 0 )
    {
		if( deltaD >= 0.0)
		{
			kon = 1;
		}
		else
		{
			kon = 2;
		}
    }

	//	Modification to the origin of the deviation point
	// a_pinch = 1 - a_pinch;


	//	STARTS BIG LOOP
	//	Positive Delta indicating loading

	if ( deltaD >= 0.0 ) {
		if ( kon == 2 ) {
			kon = 1;
			Unl = 0;
			RSE = 0.5 * fP * fP / ekunload;
			if ( (Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE)) <0.0) RSE = 0.0;
			a2	= Enrgtk - ( Enrgtot - RSE );
			
			if( a2 <= 0.0 && Enrgtk != 0.0) flgstop = 1;
			
			if( ecapk != 0.0) {
				betak = pow( ( (Enrgc-RSE) / (Enrgtk - (Enrgtot-RSE) ) ) , ck );
				ekunload = ekexcurs * ( 1 - betak );
				ekLim = ( fmax - fmin ) / ( dmax - dmin );
				if( ekunload <= ekLim ) ekunload = ekLim;
				//	if( ekunload <= ekhardNeg ) flgstop = 1;
			}
			
			//	Determination of sn according to the hysteresis status
			
			if( ekunload <= 1.e-7 ) flgstop = 1;
			
			if(fP < 0.0) {
				if( (fabs(dmax-dyieldPos) >= 1.e-10) && (fabs(dP - fP / ekunload) <= 1.e-10) ) {
					sn = 1.e-9;
				}
				else {
					sn = dP - fP / ekunload;
				}
			}

	    if( fabs( dmin - dP ) <= 1.e-10) sp = sn + 1.e-10;
		}
		
		//	LOADING
		//	Push envelope
		
		if ( d >= dmax ) {
			this->envelPosCap ( fyPos, alphaPos, capSlope, cpPos, d, &f, &ek );
			dmax = d;
			fmax = f;
		} else if ( fabs(sn) > 1.e-10) {
			this->envelPosCap ( fyPos, alphaPos, capSlope, cpPos, dmax, &fmax, &ekt );
			dch = dmax - fmax / ekunload;
			double fpdegPos = fmax * fpPos;
			ekpinch = fpdegPos / ( dmax - sn );
			fpinch = ekpinch * ( a_pinch * dch - sn );
			
			if(sn <= a_pinch * dch ) {
				if ( d < sn ) {
					ek = ekunload;
					f  = fP + ek * deltaD;
				}
				else if ( d >= sn  &&  d < a_pinch * dch ) {
					ek = ekpinch;
					f2 = ek * ( d - sn );
					f1 = fP + ekunload * deltaD;
					f= (f1<f2) ? f1 : f2;
					if( ekunload < ek ) flgstop = 1;
					if ( fabs(f-f1) < 1.e-10) ek=ekunload;
				}
				else {
					ek = (fmax-fpinch) / (dmax - a_pinch * dch);
					f2 = fpinch + ek * ( d - a_pinch * dch );
					f1 = fP + ekunload * deltaD;
					f = (f1<f2) ? f1 : f2;
					if( ekunload < ek) flgstop = 1;
					if ( fabs( f - f1 )  <  1.e-10) ek = ekunload;
				}
			}
			else if( sn  >  a_pinch * dch) {
				//	If d is larger than sn
				
				if( d  <  sn ) {
					ek = ekunload;
					f = fP + ek * deltaD;
				}
				else {
					ek = fmax / ( dmax - sn );
					f2 = ek * ( d - sn );
					f1 = fP + ekunload * deltaD;
					f = (f1 < f2) ? f1 : f2;
					if ( ekunload  <  ek) flgstop = 1;
					if ( fabs(f-f1)  <  1.e-10) ek = ekunload;
				}
			}
		}
		else {
			if (d  >  0.0) { 
				this->envelPosCap ( fyPos, alphaPos, capSlope, cpPos, d, &f, &ek );
			}
			else {
				this->envelNegCap ( fyNeg, alphaNeg, capSlope, cpNeg, d, &f, &ek );
			}
		}
		
		//	UNLOADING (deltaD<0) ---------------------------------------------
	} 
	else {

		if (kon  == 1) {
			kon = 2;
			Unl = 0;
			RSE = 0.5 * fP * fP / ekunload;
			if ( (Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE)) <0.0) RSE = 0.0;
			a2 = Enrgtk - ( Enrgtot - RSE );
			
			if( a2 <= 0.0  &&  Enrgtk!= 0.0) flgstop = 1;
			if( ecapk != 0.0 ) {
				betak = pow( ((Enrgc-RSE) / (Enrgtk-(Enrgtot-RSE))) ,ck );
				ekunload = ekexcurs * ( 1 - betak );
				ekLim= ( fmax - fmin ) / ( dmax - dmin );
				if( ekunload <= ekLim) ekunload = ekLim;
				//	if( ekunload <= ekhardPos) flgstop = 1;
			}
			
			// Determination of sn according to the hysteresis status

			if( ekunload <= 1.e-7) return 0;
			
			if( fP > 0.0) {
				if( fabs( dmin-dyieldNeg ) >= 1.e-10  && fabs(dP - fP / ekunload) <= 1.e-10) {
					sp = 1.e-9;
				}
				else
				{
					sp = dP - fP / ekunload;
				}
			}
			
			if( fabs( dmax - dP ) <= 1.e-10 ) sn = sp - 1.e-10;
			
		}
		
		//	UNLOADING (deltaD<0) ---------------------------------------------
		
		//	Push envelope
		
		if ( d < dmin ) {
			this->envelNegCap ( fyNeg, alphaNeg, capSlope, cpNeg, d, &f, &ek );
			dmin = d;
			fmin = f;
		}
		else if ( fabs(sp) > 1.e-10 ) {
			this->envelNegCap ( fyNeg, alphaNeg, capSlope, cpNeg, dmin, &fmin, &ekt );
			dch = dmin - fmin / ekunload;
			double fpdegNeg = fmin * fpNeg;
			ekpinch = fpdegNeg / ( dmin - sp );
			fpinch = ekpinch * ( a_pinch * dch - sp );
			
			if( sp >= a_pinch * dch ) {
				if ( d > sp ) {
					ek = ekunload;
					f  = fP + ek * deltaD;
				}
				else if ( d <= sp  &&  d > a_pinch * dch) {
					ek = ekpinch;
					f2 = ek * ( d - sp );
					f1 = fP + ekunload * deltaD;
					
					f = (f1>f2) ? f1 : f2;
					if ( ekunload < ek ) flgstop = 1;
					if ( fabs ( f - f1 ) < 1.e-10 ) ek = ekunload;
				}
				else {
					ek = ( fmin - fpinch ) / ( dmin - a_pinch * dch );
					f2 = fpinch + ek * ( d - a_pinch * dch );
					f1 = fP + ekunload * deltaD ;
					f = (f1>f2) ? f1 : f2;
					if( ekunload < ek ) flgstop = 1;
					if ( fabs( f - f1 ) < 1.e-10 ) ek = ekunload;
				}
			}
			else if( sp < a_pinch * dch ) {
				if( d > sp ) {
					ek = ekunload;
					f = fP + ek * deltaD;
				}
				else {
					ek = fmin / ( dmin - sp );
					f2 = ek * ( d - sp );
					f1 = fP + ekunload * deltaD;
					f  = (f1 > f2) ? f1 : f2;
					if ( ekunload < ek) flgstop = 1;
					if ( fabs(f-f1) < 1.e-10) ek = ekunload;
					
				}
			}
		}
		
		else {
			if (d > 0.0) {
				this->envelPosCap ( fyPos, alphaPos, capSlope, cpPos, d, &f, &ek );
			}
			else {
				this->envelNegCap ( fyNeg, alphaNeg, capSlope, cpNeg, d, &f, &ek );
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
	//	For the static analysis
	
	Enrgi = 0.5 * ( f + fP ) *deltaD;
	Enrgc = Enrgc + Enrgi;
	Enrgtot = Enrgtot + Enrgi;
		
	//	UPDATING OF DETERIORATION PARAMETERS ----------------------------
	//	Check the remaining capacity of the system
	
	if( flagDeg == 1 || flagDeg == 2 ) {
		
		if( ( Enrgtot >= Enrgts && Enrgts != 0.0) || ( Enrgtot >= Enrgtk && Enrgtk != 0.0) || 
			( Enrgtot >= Enrgta && Enrgta != 0.0) || ( Enrgtot >= Enrgtd && Enrgtd != 0.0)) {
			
			opserr << "Total Energy greater than capacity\n";
			flgstop=1;
		}
		//	Update beta values for strength, acc. stiff. and capping

		if( ecaps != 0.0 ) betas = pow ((Enrgc/(Enrgts-Enrgtot)) , cs );
		if( ecapa != 0.0 ) betaa = pow ((Enrgc/(Enrgta-Enrgtot)) , ca );
		if( ecapd != 0.0 ) betad = pow ((Enrgc/(Enrgtd-Enrgtot)) , cd );
	
		if( betas>=1.0 || betak>=1.0 || betaa>=1.0) {
			opserr << "Beta greater than one\n";
			flgstop=1;
		}
			
		// Initialize energy of the cycle and Kstif for next loop -------
		ekexcurs = ekunload;
		Enrgc = 0.0;
		
		//Deteriorate parameters for the next half cycle ---------------

		if( deltaD < 0.0) {
			fyNeg = fyNeg * ( 1 - betas);
			alphaNeg = alphaNeg * ( 1 - betas );
			fCapRefNeg = fCapRefNeg * ( 1 - betad );
			dmin = dmin * ( 1 + betaa );

			dyNeg = fyNeg / elstk;
			ekhardNeg = alphaNeg * elstk;

			dCap1Neg = fCapRefNeg / (elstk-capSlope*elstk);
			dCap2Neg = (fCapRefNeg+ekhardNeg*dyNeg-fyNeg) / (ekhardNeg-capSlope*elstk);
			cpNeg = (dCap1Neg < dCap2Neg) ? dCap1Neg : dCap2Neg;
			
		}
		
		else {
			fyPos = fyPos * ( 1 - betas );
	  		alphaPos = alphaPos * ( 1 - betas );
			fCapRefPos = fCapRefPos * ( 1 - betad );
			dmax = dmax * ( 1 + betaa );

			dyPos = fyPos / elstk;
			ekhardPos = alphaPos * elstk;

			dCap1Pos = fCapRefPos / (elstk-capSlope*elstk);
			dCap2Pos = (fCapRefPos+ekhardPos*dyPos-fyPos) / (ekhardPos-capSlope*elstk);
			cpPos= (dCap1Pos>dCap2Pos) ? dCap1Pos:dCap2Pos;
			
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
	hsTrial[16] = fmax;
	hsTrial[17] = fmin;
	hsTrial[18] = alphaPos;
	hsTrial[19] = alphaNeg;
	hsTrial[20] = fCapRefPos;
	hsTrial[21] = fCapRefNeg;

	return 0;
}


int Pinching::commitState()
{
  if ( DEBG ==1 ) fprintf( OutputFile , "Commit State\n" );	// debugging
  int i;
  for( i=0; i<22 ; i++ ) {
	  hsLastCommit[i] = hsCommit[i];
	  hsCommit[i] = hsTrial[i];
  }
  
  this->recordInfo();
  
  return 0;
}


void Pinching::recordInfo(int cond )
{

}


void Pinching::envelPosCap( double fy, double alphaPos, double alphaCap,
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


void Pinching::envelNegCap( double fy, double alphaNeg, double alphaCap,
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

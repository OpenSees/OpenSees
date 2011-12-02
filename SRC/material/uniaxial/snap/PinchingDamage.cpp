/* ****************************************************************** 
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
// $Date: 2004-09-01 03:53:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/PinchingDamage.cpp,v $
//
//
// PinchingDamage.cpp: implementation of the PinchingDamage class from Fortran version.
// Originally from SNAP PROGRAM by Prof H.K. Krawinkler
//
// Written: A. Altoontash & Prof. G. Deierlein 12/01
// Revised: 05/03
//
// Purpose: This file contains the implementation for the PinchingDamage class.
//
//////////////////////////////////////////////////////////////////////

#include <PinchingDamage.h>
#include <DamageModel.h>
#include <stdlib.h>
#include <Channel.h>
#include <math.h>

#include <string.h>

#define DEBG 0

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PinchingDamage::PinchingDamage(int tag, Vector inputParam , DamageModel *strength, DamageModel *stiffness,DamageModel *accelerated,DamageModel *capping)
  :UniaxialMaterial(tag,MAT_TAG_SnapPinchingDamage)
{
	if( (inputParam.Size()) < 11) 
		opserr << "Error: PinchingDamage(): inputParam, size <19\n" << "\a";
  
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
  
	// Error check
	
	if ( capSlope > 0.0 )
    {
		opserr << "Error: PinchingDamage::PinchingDamage  : CapSlope must be < 0\n" << "\a";
    }		
	
	if ( Resfac <  0.0  || Resfac > 1.0)
    {
		opserr << "Error: PinchingDamage::PinchingDamage  : Residual must be > 0 and <= 1\n" << "\a";
    }		
	
	if ( a_pinch < 0.0 || a_pinch > 1.0)
    {
		opserr << "Error: PinchingDamage::PinchingDamage  : kappad (dev. point)must be > 0 and <= 1\n" << "\a";
    }
	
	if ( alpha > 0.8 || alpha < -0.8 )
    {
		opserr << "Error: PinchingDamage::PinchingDamage  : alpha must be < 0.8 and > -0.8\n" << "\a";	
    }
	
	if ( alpha == capSlope )
    {
		opserr << "Error: PinchingDamage::PinchingDamage  : Error: alpha Hard. can not be equal to alphaCap\n" << "\a";	
    }
  
  	StrDamage = StfDamage = AccDamage = CapDamage = NULL;
	if ( strength != NULL )
	{
		StrDamage = strength->getCopy();
		if ( StrDamage == NULL ) {
			opserr << "Error: CloughDamage::CloughDamage  : Can not make a copy of strength damage model\n" << "\a";	
			exit(-1);
		}
	}

	if ( stiffness != NULL )
	{
		StfDamage = stiffness->getCopy();
		if ( StfDamage == NULL ) {
			opserr << "Error: CloughDamage::CloughDamage  : Can not make a copy of stiffness damage model\n" << "\a";	
			exit(-1);
		}
	}
	
	if ( accelerated != NULL )
	{
		AccDamage = accelerated->getCopy();
		if ( AccDamage == NULL ) {
			opserr << "Error: CloughDamage::CloughDamage  : Can not make a copy of accelerated stiffness degradation damage model\n" << "\a";	
			exit(-1);
		}
	}	

	if ( capping != NULL )
	{
		CapDamage = capping->getCopy();
		if ( CapDamage == NULL ) {
			opserr << "Error: CloughDamage::CloughDamage  : Can not make a copy of capping damage model\n" << "\a";	
			exit(-1);
		}
	}

	
	if ( DEBG ==1 )
    {
		// Open an output file for debugging
		char FileName[20];						// debugging
		sprintf(FileName, "PinchingDamage%d.out", tag);
		OutputFile = fopen ( FileName , "w" );	// debugging
		fprintf( OutputFile , "Constructor\n" );	// debugging
    }
	
	// Initialice history data
	this->revertToStart();
}


PinchingDamage::PinchingDamage()
  :UniaxialMaterial(0,MAT_TAG_SnapPinchingDamage) 
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Empty constructor\n" );	// debugging
}



PinchingDamage::~PinchingDamage()
{
  if ( DEBG ==1 )
    {
      fprintf( OutputFile , "Distructor\n" );		// debugging
      fclose( OutputFile );						// debugging
    }

   	if ( StrDamage != 0 ) delete StrDamage;
	if ( StfDamage != 0 ) delete StfDamage;
	if ( AccDamage != 0 ) delete AccDamage;
	if ( CapDamage != 0 ) delete CapDamage;
}


int PinchingDamage::revertToStart()
{
	dyieldPos = fyieldPos / elstk;
	dyieldNeg = fyieldNeg / elstk;
	
	double ekhard = elstk * alpha;

	double fPeakPos = fyieldPos + ekhard * ( capDispPos - dyieldPos );
	double fPeakNeg = fyieldNeg + ekhard * ( capDispNeg - dyieldNeg );
	
	hsCommit[0] = 0.0;			// d
	hsCommit[1] = 0.0;			// f
	hsCommit[2] = elstk;		// ek
	hsCommit[3] = elstk;		// ekunload
	hsCommit[4] = elstk;		// ekexcurs
	hsCommit[5] = 0.0;			// Enrgtot
	hsCommit[6] = 0.0;			// Enrgc
	hsCommit[7] = 0.0;			// sp
	hsCommit[8] = 0.0;			// sn
	hsCommit[9] = 0.0;			// kon
	hsCommit[10] = dyieldPos;	// dmax
	hsCommit[11] = dyieldNeg;	// dmin
	hsCommit[12] = fyieldPos;	// fyPos
	hsCommit[13] = fyieldNeg;	// fyNeg
	hsCommit[14] = capDispPos;	// cpPos
	hsCommit[15] = capDispNeg;	// cpNeg
	hsCommit[16] = fyieldPos;	// fmax
	hsCommit[17] = fyieldNeg;	// fmin
	hsCommit[18] = alpha;		// alphaPos
	hsCommit[19] = alpha;		// alphaNeg
	hsCommit[20] = -capSlope * elstk * capDispPos + fPeakPos;	// fCapRefPos
	hsCommit[21] = -capSlope * elstk * capDispNeg + fPeakNeg;	// fCapRefNeg
	hsCommit[22] = dyieldPos;	// dmaxDeg
	hsCommit[23] = dyieldNeg;	// dminDeg
	for( int i=0 ; i<24; i++) {
		hsTrial[i]		= hsCommit[i];
		hsLastCommit[i] = hsCommit[i];
	}
	if ( StrDamage != NULL ) StrDamage->revertToStart();
	if ( StfDamage != NULL ) StfDamage->revertToStart();
	if ( AccDamage != NULL ) AccDamage->revertToStart();
	if ( CapDamage != NULL ) CapDamage->revertToStart();
	
	return 0;
}


void PinchingDamage::Print(OPS_Stream &s, int flag)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Print\n" );	// debugging
	s << "BondSlipMaterial Tag: " << this->getTag() << endln;
	s << "D : " << hsTrial[0] << endln;
	s << "F : " << hsTrial[1] << endln;
	s << "EK: " << hsTrial[2]  << endln;
	s << "Input Parameters:\n";
	s << endln;
}


int PinchingDamage::revertToLastCommit()
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Revert to last commit\n" );	// debugging
	
	for(int i=0; i<24; i++) {
		hsTrial[i]	= hsCommit[i];
		hsCommit[i] = hsLastCommit[i];
	}
	if ( StrDamage != NULL ) StrDamage->revertToLastCommit();
	if ( StfDamage != NULL ) StfDamage->revertToLastCommit();
	if ( AccDamage != NULL ) AccDamage->revertToLastCommit();
	if ( CapDamage != NULL ) CapDamage->revertToLastCommit();
	
	return 0;
}


double PinchingDamage::getTangent()
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Get tangent\n" );	// debugging
	return hsTrial[2];
}

double PinchingDamage::getInitialTangent (void)
{
	return elstk;
}

double PinchingDamage::getStress()
{
	return hsTrial[1];
}


double PinchingDamage::getStrain (void)
{
	return hsTrial[0];
}


int PinchingDamage::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
	return 0;
}


int PinchingDamage::sendSelf(int cTag, Channel &theChannel)
{
	return 0;
}


UniaxialMaterial *PinchingDamage::getCopy(void)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Get copy\n" );	// debugging
	Vector inp(11);
	
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

	
	PinchingDamage *theCopy = new PinchingDamage(this->getTag(), inp, StrDamage, StfDamage, AccDamage, CapDamage  );

	for (int i=0; i<24; i++) {
		theCopy->hsTrial[i] = hsTrial[i];
		theCopy->hsCommit[i] = hsCommit[i];
		theCopy->hsLastCommit[i] = hsLastCommit[i];
	}
	
	return theCopy;
}


int PinchingDamage::setTrialStrain( double d, double strainRate)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Set trial displacement\n" );	// debugging
	int flagDeg , flgstop , Unl, kon;
	double deltaD;
	double ekP,ekunload,ekexcurs,Enrgtot,Enrgc,sp,sn,dmax,dmin,fyPos,fyNeg,fP,dP;
	double ek,f,betaa,betad,betak,betas,ekLim;
	double ekt,f1,f2,fpinch,dCap1Pos,dCap2Pos,dCap1Neg,dCap2Neg;
	double dyPos,dyNeg,dch,ekpinch, cpPos, cpNeg;
	double ekhardNeg,ekhardPos,fmax,fmin,dmaxDeg,dminDeg;
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
	dmaxDeg		= hsCommit[22];
	dminDeg		= hsCommit[23];
	
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
			
			if( StfDamage != NULL ) {				
	
				betak = StfDamage->getDamage();
				if( betak>=1.0 ) {
					opserr << "Total loss for stiffness degradation\n";
					betak = 1.0;
				}
				ekunload = elstk * ( 1 - betak );
				ekLim = ( fmax - fmin ) / ( dmaxDeg - dminDeg );
				if( ekunload <= ekLim ) ekunload = ekLim;
				//	if( ekunload <= ekhardNeg ) flgstop = 1;
			}
			
			//	Determination of sn according to the hysteresis status
			
			if( ekunload <= 1.e-7 ) flgstop = 1;
			
			if(fP < 0.0) {
				if( (fabs(dmaxDeg-dyieldPos) >= 1.e-10) && (fabs(dP - fP / ekunload) <= 1.e-10) ) {
					sn = 1.e-9;
				}
				else {
					sn = dP - fP / ekunload;
				}
			}

	    if( fabs( dminDeg - dP ) <= 1.e-10) sp = sn + 1.e-10;
		}
		
		//	LOADING
		//	Push envelope
		
		if ( d >= dmaxDeg ) {
			this->envelPosCap ( fyPos, alphaPos, capSlope, cpPos, d, &f, &ek );
			dmax = d;
			dmaxDeg = d;
			fmax = f;
		} else if ( fabs(sn) > 1.e-10) {
			this->envelPosCap ( fyPos, alphaPos, capSlope, cpPos, dmaxDeg, &fmax, &ekt );
			dch = dmaxDeg - fmax / ekunload;
			double fpdegPos = fmax * fpPos;
			ekpinch = fpdegPos / ( dmaxDeg - sn );
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
					ek = (fmax-fpinch) / (dmaxDeg - a_pinch * dch);
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
					ek = fmax / ( dmaxDeg - sn );
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

		if ( kon  == 1) {
			kon = 2;
			Unl = 0;

			if( StfDamage != NULL ) {

				betak = StfDamage->getDamage();
				if( betak>=1.0 ) {
					opserr << "Total loss for stiffness degradation\n";
					betak = 1.0;
				}
				ekunload = elstk * ( 1 - betak );
				ekLim= ( fmax - fmin ) / ( dmaxDeg - dminDeg );
				if( ekunload <= ekLim) ekunload = ekLim;
				//	if( ekunload <= ekhardPos) flgstop = 1;
			}
			
			// Determination of sn according to the hysteresis status

			if( ekunload <= 1.e-7) return 0;
			
			if( fP > 0.0) {
				if( fabs( dminDeg-dyieldNeg ) >= 1.e-10  && fabs(dP - fP / ekunload) <= 1.e-10) {
					sp = 1.e-9;
				}
				else
				{
					sp = dP - fP / ekunload;
				}
			}
			
			if( fabs( dmaxDeg - dP ) <= 1.e-10 ) sn = sp - 1.e-10;
			
		}
		
		//	UNLOADING (deltaD<0) ---------------------------------------------
		
		//	Push envelope
		
		if ( d < dminDeg ) {
			this->envelNegCap ( fyNeg, alphaNeg, capSlope, cpNeg, d, &f, &ek );
			dmin = d;
			dminDeg = d;
			fmin = f;
		}
		else if ( fabs(sp) > 1.e-10 ) {
			this->envelNegCap ( fyNeg, alphaNeg, capSlope, cpNeg, dminDeg, &fmin, &ekt );
			dch = dminDeg - fmin / ekunload;
			double fpdegNeg = fmin * fpNeg;
			ekpinch = fpdegNeg / ( dminDeg - sp );
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
					ek = ( fmin - fpinch ) / ( dminDeg - a_pinch * dch );
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
					ek = fmin / ( dminDeg - sp );
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

	//	DAMAGE CALCULATIONS ---------------------------------------------
	// Calling the damage object
	if ( StfDamage != NULL ) { 
		betak = StfDamage->getDamage();
		if( fabs(betak) >= 1.0) betak = 1.0;
	}

	if ( StrDamage != NULL ) { 
		betas = StrDamage->getDamage();
		if( fabs(betas) >= 1.0) betas = 1.0;
	}

	if ( AccDamage != NULL ) {
		betaa = AccDamage->getDamage();
		if( fabs(betaa) >= 1.0) betaa = 1.0;
	}
	if ( CapDamage != NULL ) {
		betad = CapDamage->getDamage();
		if( fabs(betad) >= 1.0) betad = 1.0;
	}

	//	Flag to deteriorate parameters on the opposite side of the loop --	

	if( f * fP < 0.0) {
		if( fP >0.0 && dmaxDeg >dyieldPos) flagDeg = 1;	// positive half cycle
		if( fP < 0.0 && dminDeg < dyieldNeg) flagDeg = 2;	// negative half cycle 
	}
		
	//	UPDATING OF DETERIORATION PARAMETERS ----------------------------
	//	Check the remaining capacity of the system
	
	if( flagDeg == 1 || flagDeg == 2 ) {
		// Initialize energy of the cycle and Kstif for next loop -------
		ekexcurs = ekunload;
		Enrgc = 0.0;
		
		//Deteriorate parameters for the next half cycle ---------------

		if( deltaD < 0.0) {
			
			//Update beta values for strength, acc. stiff. and capping
			if( StrDamage != NULL ) betas = StrDamage->getNegDamage();
			if( betas>=1.0 ) {
				opserr << "Total loss for strength degradation\n";
				betas = 1.0;
			}
			if( AccDamage != NULL ) betaa = AccDamage->getNegDamage();
			if( betaa>=1.0 ) {
				opserr << "Total loss for accelerated stiffness degradation\n";
				betaa = 1.0;
			}
			if( CapDamage != NULL ) betad = CapDamage->getNegDamage();
			if( betad>=1.0 ) {
				opserr << "Total loss for capping degradation\n";
				betad = 1.0;
			}
		

			fyNeg = fyieldNeg * ( 1 - betas );
			alphaNeg = alpha * ( 1 - betas );
			fCapRefNeg = (-capSlope * elstk * capDispNeg + fyieldNeg + elstk * alpha * ( capDispNeg - dyieldNeg) ) * ( 1 - betad );
			dminDeg = dmin * ( 1 + betaa );

			dyNeg = fyNeg / elstk;
			ekhardNeg = alphaNeg * elstk;

			dCap1Neg = fCapRefNeg / (elstk-capSlope*elstk);
			dCap2Neg = (fCapRefNeg+ekhardNeg*dyNeg-fyNeg) / (ekhardNeg-capSlope*elstk);
			cpNeg = (dCap1Neg < dCap2Neg) ? dCap1Neg : dCap2Neg;
		}
		
		else {
			
			//Update beta values for strength, acc. stiff. and capping

			if( StrDamage != NULL ) betas = StrDamage->getPosDamage();
			if( betas>=1.0 ) {
				opserr << "Total loss for strength degradation\n";
				betas = 1.0;
			}
			if( AccDamage != NULL ) betaa = AccDamage->getPosDamage();
			if( betaa>=1.0 ) {
				opserr << "Total loss for accelerated stiffness degradation\n";
				betaa = 1.0;
			}
			if( CapDamage != NULL ) betad = CapDamage->getPosDamage();
			if( betad>=1.0 ) {
				opserr << "Total loss for capping degradation\n";
				betad = 1.0;
			}
		
			fyPos = fyieldPos * ( 1 - betas );
	  		alphaPos = alpha * ( 1 - betas );
			fCapRefPos = (-capSlope * elstk * capDispPos + fyieldPos + elstk * alpha * ( capDispPos - dyieldPos ) ) * ( 1 - betad );
			dmaxDeg = dmax * ( 1 + betaa );

			dyPos = fyPos / elstk;
			ekhardPos = alphaPos * elstk;

			dCap1Pos = fCapRefPos / (elstk-capSlope*elstk);
			dCap2Pos = (fCapRefPos+ekhardPos*dyPos-fyPos) / (ekhardPos-capSlope*elstk);
			cpPos= (dCap1Pos>dCap2Pos) ? dCap1Pos:dCap2Pos;
		}
		flagDeg = 0;
	}
	
	// Relationship between basic variables and hstory array for next cycle
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
	hsTrial[22] = dmaxDeg;
	hsTrial[23] = dminDeg;

	return 0;
}


int PinchingDamage::commitState()
{
  if ( DEBG ==1 ) fprintf( OutputFile , "Commit State\n" );	// debugging
  int i;
  for( i=0; i<24 ; i++ ) {
	  hsLastCommit[i] = hsCommit[i];
	  hsCommit[i] = hsTrial[i];
  }
  
  	// Calling the damage object
	Vector InforForDamage(3);
	InforForDamage(0) = hsCommit[0];
	InforForDamage(1) = hsCommit[1];
	InforForDamage(2) = hsCommit[3];

	if ( StfDamage != NULL ) { 
		StfDamage->setTrial(InforForDamage);
		StfDamage->commitState();
	}

	InforForDamage(2) = 0.0;

	if ( StrDamage != NULL ) { 
		StrDamage->setTrial(InforForDamage);
		StrDamage->commitState();
	}

	if ( AccDamage != NULL ) {
		AccDamage->setTrial(InforForDamage);
		AccDamage->commitState();
	}
	if ( CapDamage != NULL ) {
		CapDamage->setTrial(InforForDamage);
		CapDamage->commitState();
	}

  this->recordInfo();
  
  return 0;
}


void PinchingDamage::recordInfo(int cond )
{

}


void PinchingDamage::envelPosCap( double fy, double alphaPos, double alphaCap,
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


void PinchingDamage::envelNegCap( double fy, double alphaNeg, double alphaCap,
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

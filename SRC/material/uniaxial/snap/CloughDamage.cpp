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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/CloughDamage.cpp,v $
//
//
// CloughDamage.cpp: implementation of the CloughDamage class from Fortran version.
// Originally from SNAP PROGRAM by Prof H.K. Krawinkler
//
// Written: A. Altoontash & Prof. G. Deierlein 12/01
// Revised: 03/02
//
// Purpose: This file contains the implementation for the CloughDamage class.
//
//////////////////////////////////////////////////////////////////////

#include <CloughDamage.h>
#include <DamageModel.h>
#include <stdlib.h>
#include <Channel.h>
#include <math.h>

#include <string.h>

#define DEBG 0

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CloughDamage::CloughDamage(int tag, Vector inputParam, DamageModel *strength, DamageModel *stiffness,DamageModel *accelerated,DamageModel *capping )
:UniaxialMaterial(tag,MAT_TAG_SnapCloughDamage)
{
	if( (inputParam.Size()) < 8) 
		opserr << "Error: CloughDamage(): inputParam, size <16\n" << "\a";
	
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

	
	// Error check
	
	if ( capSlope > 0.0 )
		opserr << "Error: CloughDamage::CloughDamage  : CapSlope must be < 0\n" << "\a";
	
	if ( Resfac <  0.0  || Resfac > 1.0)
		opserr << "Error: CloughDamage::CloughDamage  : Residual must be > 0 and <= 1\n" << "\a";
	
	if ( alpha > 0.8 || alpha < -0.8 )
		opserr << "Error: CloughDamage::CloughDamage  : alpha must be < 0.8 and > -0.8\n" << "\a";	
	
	if ( alpha == capSlope )
		opserr << "Error: CloughDamage::CloughDamage  : Error: alpha Hard. can not be equal to alphaCap\n" << "\a";	
	
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
			opserr << "Error: CloughDamage::CloughDamage  : Can not make a copy of accelerated stiffness damage model\n" << "\a";	
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

	// Initialize history data
	this->revertToStart();

}


CloughDamage::CloughDamage()
  :UniaxialMaterial(0,MAT_TAG_SnapCloughDamage) 
{

}



CloughDamage::~CloughDamage()
{
	
  	if ( StrDamage != 0 ) delete StrDamage;
	if ( StfDamage != 0 ) delete StfDamage;
	if ( AccDamage != 0 ) delete AccDamage;
	if ( CapDamage != 0 ) delete CapDamage;
}



int CloughDamage::revertToStart()
{

	dyieldPos = fyieldPos/elstk;
	dyieldNeg = fyieldNeg/elstk;
	
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
		hsCommit[i]		= hsTrial[i];
		hsLastCommit[i] = hsTrial[i];
	}
	if ( StrDamage != NULL ) StrDamage->revertToStart();
	if ( StfDamage != NULL ) StfDamage->revertToStart();
	if ( AccDamage != NULL ) AccDamage->revertToStart();
	if ( CapDamage != NULL ) CapDamage->revertToStart();

	return 0;
}

void CloughDamage::Print(OPS_Stream &s, int flag)
{

	s << "BondSlipMaterial Tag: " << this->getTag() << endln;
	s << "D : " << hsTrial[0] << endln;
	s << "F : " << hsTrial[1] << endln;
	s << "EK: " << hsTrial[2]  << endln;
	
	s << endln;
}


int CloughDamage::revertToLastCommit()
{
  
  for(int i=0; i<24; i++) {
	  hsTrial[i] = hsCommit[i];
	  hsCommit[i] = hsLastCommit[i];
  }
  if ( StrDamage != NULL ) StrDamage->revertToLastCommit();
  if ( StfDamage != NULL ) StfDamage->revertToLastCommit();
  if ( AccDamage != NULL ) AccDamage->revertToLastCommit();
  if ( CapDamage != NULL ) CapDamage->revertToLastCommit();

  return 0;
}


double CloughDamage::getTangent()
{

	return hsTrial[2];
}

double CloughDamage::getInitialTangent (void)
{

	return elstk;
}

double CloughDamage::getStress()
{

	return hsTrial[1];
}


double CloughDamage::getStrain (void)
{

	return hsTrial[0];
}


int CloughDamage::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{

	return 0;
}


int CloughDamage::sendSelf(int cTag, Channel &theChannel)
{

	return 0;
}


UniaxialMaterial *CloughDamage::getCopy(void)
{

	Vector inp(8);
	
	inp[0]  = 	elstk;
	inp[1]  =  	fyieldPos;
	inp[2]  =  	fyieldNeg;
	inp[3]  = 	alpha;
	inp[4]  = 	Resfac;
	inp[5]  = 	capSlope;
	inp[6]  = 	capDispPos;
	inp[7]  = 	capDispNeg;
  
  CloughDamage *theCopy = new CloughDamage(this->getTag(), inp ,StrDamage, StfDamage, AccDamage, CapDamage);
  
  for (int i=0; i<24; i++) {
    theCopy->hsTrial[i] = hsTrial[i];
	theCopy->hsCommit[i] = hsCommit[i];
    theCopy->hsLastCommit[i] = hsLastCommit[i];
  }
  
  return theCopy;
}


int CloughDamage::setTrialStrain( double d, double strainRate)
{


	int kon,idiv,Unl,flagDeg,flgstop;

	double ekP,ek,ekunload,sp,sn,dP,fP,f;
	double deltaD,err,f1,f2,dmax,dmin,fmax,fmin,ekt;
	double betas,betak,betaa;
	double Enrgtot,Enrgc,dyPos,dyNeg;
	double fyPos,fyNeg;
	double betad,cpPos,cpNeg;
	double dlstPos,flstPos,dlstNeg,flstNeg,ekc,tst,ekexcurs;
	double dCap1Pos,dCap2Pos,dCap1Neg,dCap2Neg;
	double ekhardNeg,ekhardPos,alphaPos,alphaNeg;
	double fCapRefPos,fCapRefNeg;
	
	//	Relationship between basic variables and hsTrial array
	idiv = 1;
	flagDeg = 0;
	flgstop = 0;
	Unl = 1;
	err = 0.0;

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
	dlstPos		= hsCommit[16];
	flstPos		= hsCommit[17];
	dlstNeg		= hsCommit[18];
	flstNeg		= hsCommit[19];
	alphaPos	= hsCommit[20];
	alphaNeg	= hsCommit[21];
	fCapRefPos	= hsCommit[22];
	fCapRefNeg	= hsCommit[23];

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

			if( StfDamage != NULL ) {
		
				betak = StfDamage->getDamage();
				if( betak >= 1.0 ) {
					opserr << "Total loss for stiffness degradation\n";
					betak = 1.0;
				}
				ekunload = ekexcurs * ( 1 - betak );
				if( ekunload <= ekhardNeg ) flgstop=1;
			}
			
			//	Determination of sn according to the hysteresis status
			if( ekunload <= 1.e-7 ) flgstop=1;
			
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
			
			if( StfDamage != NULL ) {
	
				betak = StfDamage->getDamage();
				if( betak >= 1.0 ) {
					opserr << "Total loss for stiffness degradation\n";
					betak = 1.0;
				}
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
	
	//	DAMAGE CALCULATIONS ---------------------------------------------
	// Calling the damage object
	if ( StfDamage != NULL ) { 
		betak = StrDamage->getDamage();
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
		if( fP >0.0 && dmax >dyieldPos) flagDeg = 1;	// positive half cycle
		if( fP < 0.0 && dmin < dyieldNeg) flagDeg = 2;	// negative half cycle 
	}	

	//	UPDATING OF DETERIORATION PARAMETERS ----------------------------
	//	Check the remaining capacity of the system
	
	if( flagDeg == 1 || flagDeg == 2 ) {

		//	Update beta values for strength, acc. stiff. and capping

		if( StrDamage != NULL ) betas = StrDamage->getDamage();
		if( betas>=1.0 ) {
			opserr << "Total loss for strength degradation\n";
			betas = 1.0;
		}
		if( AccDamage != NULL ) betaa = AccDamage->getDamage();
		if( betaa>=1.0 ) {
			opserr << "Total loss for accelerated stiffness degradation\n";
			betaa = 1.0;
		}
		if( CapDamage != NULL ) betad = CapDamage->getDamage();
		if( betad>=1.0 ) {
			opserr << "Total loss for capping degradation\n";
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


int CloughDamage::commitState()
{

  int i;
  for( i=0; i<24; i++ ) {
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


void CloughDamage::recordInfo(int cond )
{

}


void CloughDamage::envelPosCap( double fy, double alphaPos, double alphaCap,
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


void CloughDamage::envelNegCap( double fy, double alphaNeg, double alphaCap,
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

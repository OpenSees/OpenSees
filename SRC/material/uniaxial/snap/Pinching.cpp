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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:01:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/Pinching.cpp,v $
//
//
// Pinching.cpp: implementation of the Pinching class from Fortran version.
// Originally from SNAP PROGRAM by Prof H.K. Krawinkler
//
// Written: A. Altoontash & Prof. G. Deierlein 12/01
// Revised: 03/02
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
	alpha      =  strain hardening ratio (fraction of elstk)

	ecaps, ecapd, ecapk, ecapa = parameter expressing the hystetic
		    energy dissipacion capacity.
	Enrgts, Enrgtd, Enrgtk, Enrgta = hysteretic energy dissipation
			capacity
	*/
  elstk		= inputParam[0];
  fyieldPos	= inputParam[1];
  fyieldNeg	= inputParam[2];
  alpha		= inputParam[3];
  Resfac		= inputParam[4];
  capSlope	= inputParam[5];
  capDispPos	= inputParam[6];
  capDispNeg	= inputParam[7];
  fpPos		= inputParam[8];
  fpNeg		= inputParam[9];
  a_pinch		= inputParam[10];
  ecaps		= inputParam[11];
  ecapk		= inputParam[12];
  ecapa		= inputParam[13];
  ecapd		= inputParam[14];
  cs			= inputParam[15];
  ck			= inputParam[16];
  ca			= inputParam[17];
  cd			= inputParam[18];
  
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
      opserr << "Error: Pinching::Pinching  : Error: alfa Hard. can not be equal to alfaCap\n" << "\a";	
    }
  
  
  if ( DEBG ==1 )
    {
      // Open an output file for debugging
      //char buffer[10];						// debugging
      char FileName[20];						// debugging
      //strcpy( FileName, "Pinching");			// debugging
      //_itoa(tag,buffer,10);					// debugging
      //strcat( FileName , buffer);				// debugging
      //strcat( FileName , ".out");				// debugging
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
  if ( DEBG ==1 ) fprintf( OutputFile , "Revert to start\n" );		// debugging
  ekhard = elstk * alpha;
  alfaNeg = alpha;
  alfaPos = alpha;
  dyieldPos = fyieldPos / elstk;
  dyieldNeg = fyieldNeg / elstk;
  Enrgts = fyieldPos * dyieldPos * ecaps;
  Enrgta = fyieldPos * dyieldPos * ecapa;
  Enrgtk = fyieldPos * dyieldPos * ecapk;
  Enrgtd = fyieldPos * dyieldPos * ecapd;
  
  fPeakPos = fyieldPos + ekhard * ( capDispPos - dyieldPos );
  fPeakNeg = fyieldNeg + ekhard * ( capDispNeg - dyieldNeg );
  
  fpdegPos = fpPos;
  fpdegNeg = fpNeg;
  cpPos = capDispPos;
  cpNeg = capDispNeg;
  
  fCapRefPos = -capSlope * elstk * capDispPos + fPeakPos;
  fCapRefNeg = -capSlope * elstk * capDispNeg + fPeakNeg;
  
  flagDir = 0;
  
  hstv[0] = elstk;		// ekP
  hstv[1] = elstk;		// ekunload
  hstv[2] = elstk;		// ekexcurs
  hstv[3] = 0.0;			// Enrgtot
  hstv[4] = 0.0;			// Enrgc
  hstv[5] = 0.0;			// sp
  hstv[6] = 0.0;			// sn
  hstv[7] = 0.0;			// kon
  hstv[8] = dyieldPos;	// dmax
  hstv[9] = dyieldNeg;	// dmin
  hstv[10] = fyieldPos;	// fyPos
  hstv[11] = fyieldNeg;	// fyNeg
  hstv[12] = capDispPos;	// cpPos
  hstv[13] = capDispNeg;	// cpNeg
  hstv[14] = 0.0;			// f
  hstv[15] = 0.0;			// d
  hstv[16] = fyieldPos;	// fmax
  hstv[17] = fyieldNeg;	// fmin
  
  for( int i=0 ; i<18 ; i++) hsLastCommit[i] = hstv[i];
  
  return 0;
}


void Pinching::Print(OPS_Stream &s, int flag)
{
  if ( DEBG ==1 ) fprintf( OutputFile , "Print\n" );	// debugging
  s << "BondSlipMaterial Tag: " << this->getTag() << endln;
  s << "D : " << hstv[15] << endln;
  s << "EK: " << hstv[0]  << endln;
  s << "F : " << hstv[14] << endln;
  s << "Input Parameters:\n";
  s << endln;
}


int Pinching::revertToLastCommit()
{
  if ( DEBG ==1 ) fprintf( OutputFile , "Revert to last commit\n" );	// debugging
  
  for(int i=0; i<18; i++) hstv[i] = hsLastCommit[i];
  
  return 0;
}


double Pinching::getTangent()
{
  if ( DEBG ==1 ) fprintf( OutputFile , "Get tangent\n" );	// debugging
  return hstv[0];
}

double Pinching::getInitialTangent (void)
{
  if ( DEBG ==1 ) fprintf( OutputFile , "Get initial tangent\n" );	// debugging
  return elstk;
}

double Pinching::getStress()
{
  if ( DEBG ==1 ) fprintf( OutputFile , "Get stress\n" );	// debugging
  return hstv[14];
}


double Pinching::getStrain (void)
{
  if ( DEBG ==1 ) fprintf( OutputFile , "Get strain\n" );	// debugging
  return hstv[15];
}


int Pinching::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  if ( DEBG ==1 ) fprintf( OutputFile , "Receive self\n" );	// debugging
  return 0;
}


int Pinching::sendSelf(int cTag, Channel &theChannel)
{
  if ( DEBG ==1 ) fprintf( OutputFile , "Send self\n" );	// debugging
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
  
  for (int i=0; i<18; i++) {
    theCopy->hstv[i] = hstv[i];
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
  double dyPos,dyNeg,testValue,dch,ekpinch;
  double Enrgi,ekhardNeg,ekhardPos,fmax,fmin;
  
  flagDeg = 0;
  flgstop = 0;
  Unl = 1;
  
  
  ekP			= hsLastCommit[0];
  ekunload	= hsLastCommit[1];
  ekexcurs	= hsLastCommit[2];
  Enrgtot		= hsLastCommit[3];
  Enrgc		= hsLastCommit[4];
  sp			= hsLastCommit[5];
  sn			= hsLastCommit[6];
  kon			= (int) hsLastCommit[7];
  dmax		= hsLastCommit[8];
  dmin		= hsLastCommit[9];
  fyPos		= hsLastCommit[10];
  fyNeg		= hsLastCommit[11];			   
  cpPos		= hsLastCommit[12];
  cpNeg		= hsLastCommit[13];
  fP			= hsLastCommit[14];
  dP			= hsLastCommit[15];
  fmax		= hsLastCommit[16];
  fmin		= hsLastCommit[17];
  
  
  deltaD = d - dP;
  
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
  
  
  int flag_anal = 1;
  double tol = 1e-7;
  
  //  Modification to the origin of the deviation point
  a_pinch = 1 - a_pinch;
  
  // STARTS BIG LOOP
  // Positive Delta indicating loading
  
  if ( deltaD >= 0.0 )
    {
      if ( kon == 2 )
	{
	  kon = 1;
	  Unl = 0;
	  RSE = 0.5 * fP * fP / ekunload;
	  a2	= Enrgtk - ( Enrgtot - RSE );
	  
	  if( a2 <= 0.0 && Enrgtk != 0.0) flgstop = 1;
	  if( ecapk != 0.0 ) 
	    {
	      betak = pow( ( ( Enrgc - RSE ) / ( Enrgtk - ( Enrgtot - RSE ) ) ) , ck );
	      ekunload = ekexcurs * ( 1 - betak );
	      ekLim = ( fmax - fmin ) / ( dmax - dmin );
	      if ( ekunload <= ekLim ) ekunload = ekLim;
	      //if ( ekunload <= ekhardNeg ) flgstop = 1;
	    }
	  
	  // Determination of sn according to the hysteresis status
	  if( ekunload <= 1.e-7) flgstop = 1;
	  if( fP < 0.0 ) 
	    {
	      testValue = dP - fP / ekunload;
	      if( fabs ( dmax - dyieldPos ) >= 1.e-10 && fabs( testValue ) <= 1.e-10) 
		{
		  sn = 1.e-9;
		}
	      else
		{
		  sn = dP - fP / ekunload;
		}
	    }
	  if( fabs( dmin - dP ) <= 1.e-10) sp = sn + 1.e-10;
	}
      
      
      // LOADING
      
      // Push envelope
      
      if ( d >= dmax ) 
	{
	  
	  this->envelPosCap ( fyPos, alfaPos, capSlope, cpPos, d, &f, &ek );
	  dmax = d;
	  fmax = f;
	}
      else
	{
	  if ( fabs( sn ) > 1.e-10) 
	    {
	      this->envelPosCap (fyPos,alfaPos,capSlope,cpPos,dmax,&fmax,&ekt);
	      dch = dmax - fmax / ekunload;
	      fpdegPos = fmax * fpPos;
	      ekpinch = fpdegPos / ( dmax - sn );
	      fpinch = ekpinch * ( a_pinch * dch - sn );
	      
	      if( sn <= a_pinch * dch ) 
		{
		  if ( d < sn ) 
		    {
		      ek = ekunload;
		      f  = fP + ek * deltaD;
		    }
		  else
		    {
		      if ( d >= sn && d < a_pinch * dch )
			{
			  ek = ekpinch;
			  f2 = ek * ( d - sn );
			  f1 = fP + ekunload * deltaD;
			  f= ( f1 < f2 ) ? f1 : f2;
			  if( ekunload < ek ) flgstop = 1;
			  if ( fabs( f - f1 ) < 1.e-10 ) ek = ekunload;
			}
		      else
			{
			  ek = ( fmax - fpinch ) / ( dmax - a_pinch * dch );
			  f2 = fpinch + ek * ( d - a_pinch * dch);
			  f1 = fP + ekunload * deltaD;
			  f = ( f1 < f2 ) ? f1 : f2;
			  if( ekunload < ek ) flgstop = 1;
			  if ( fabs( f - f1 ) < 1.e-10 ) ek = ekunload;
			}
		    }
		}
	      else
		{
		  if( sn > a_pinch * dch )
		    {
		      // If d is larger than sn
		      
		      if( d < sn ) 
			{
			  ek = ekunload;
			  f = fP + ek * deltaD;
			}
		      else
			{
			  ek = fmax / ( dmax - sn );
			  f2 = ek * ( d - sn );
			  f1 = fP + ekunload * deltaD;
			  f = ( f1 < f2 ) ? f1: f2;
			  if( ekunload < ek ) flgstop = 1 ;
			  if( fabs( f - f1 ) < 1.e-10 ) ek = ekunload;
			}    
		    }
		}
	    }
	  else
	    {
	      if ( d > 0.0) 
		{
		  this->envelPosCap (fyPos,alfaPos,capSlope,cpPos,d,&f,&ek);
		}
	      else
		{
		  this->envelNegCap (fyNeg,alfaNeg,capSlope,cpNeg,d,&f,&ek);
		}
	    }
	}
    }
  //UNLOADING (deltaD<0) ---------------------------------------------
  else
    {
      if (kon == 1) 
	{
	  kon = 2;
	  Unl = 0;
	  RSE = 0.5 * fP * fP / ekunload;
	  a2 = Enrgtk - ( Enrgtot - RSE );
	  
	  if( a2 <= 0.0 && Enrgtk != 0.0) flgstop=1;
	  if( ecapk != 0.0 ) 
	    {
	      betak = pow( ( ( Enrgc - RSE ) / ( Enrgtk - ( Enrgtot - RSE ) ) ) , ck );
	      ekunload = ekexcurs * ( 1 - betak );
	      ekLim= ( fmax - fmin ) / ( dmax - dmin );
	      if( ekunload <= ekLim ) ekunload = ekLim;
	      //if( ekunload <= ekhardPos ) flgstop = 1;
	    }
	  
	  // Determination of sn according to the hysteresis status
	  
	  if( ekunload <= 1.e-7 ) return 0;
	  if( fP > 0.0 ) 
	    {
	      testValue = dP - fP / ekunload;
	      if( fabs( dmin - dyieldNeg ) >= 1.e-10 && fabs( testValue ) <= 1.e-10 )
		{
		  sp = 1.e-9;
		}
	      else
		{
		  sp = dP - fP / ekunload;
		}
	    }
	  if( fabs( dmax - dP ) <= 1.e-10 ) sn = sp - 1.e-10;
	}
      
      // UNLOADING (deltaD<0) ---------------------------------------------
      
      // Push envelope
      if ( d < dmin ) 
	{
	  this-> envelNegCap (fyNeg,alfaNeg,capSlope,cpNeg,d,&f,&ek);
	  dmin = d;
	  fmin = f;
	}
      else
	{
	  if ( fabs( sp ) > 1.e-10 )
	    {
	      this->envelNegCap (fyNeg,alfaNeg,capSlope,cpNeg,dmin,&fmin,&ekt);
	      dch = dmin - fmin / ekunload;
	      fpdegNeg = fmin * fpNeg;
	      ekpinch = fpdegNeg / ( dmin - sp );
	      fpinch = ekpinch * ( a_pinch * dch - sp );
	      if( sp >= a_pinch * dch )
		{
		  if ( d > sp ) 
		    {
		      ek = ekunload;
		      f  = fP + ek * deltaD;
		    }
		  else
		    {
		      if ( d <= sp && d > a_pinch * dch ) 
			{
			  ek = ekpinch;
			  f2 = ek * ( d - sp );
			  f1 = fP + ekunload * deltaD;
			  f= (f1>f2) ? f1 :f2;
			  if( ekunload < ek ) flgstop = 1;
			  if ( fabs( f - f1 ) < 1.e-10 ) ek = ekunload;
			}
		      else
			{
			  ek = ( fmin - fpinch ) / ( dmin - a_pinch * dch );
			  f2 = fpinch + ek * ( d - a_pinch * dch );
			  f1 = fP + ekunload * deltaD  ;
			  f = (f1>f2) ? f1 : f2;
			  if( ekunload < ek ) flgstop = 1;
			  if ( fabs( f - f1 ) < 1.e-10 ) ek = ekunload;
			}
		    }
		}
	      else
		{
		  if( sp < a_pinch * dch )
		    {
		      if( d > sp ) 
			{
			  ek = ekunload;
			  f = fP + ek * deltaD;
			}
		      else   
			{
			  ek = fmin / ( dmin - sp );
			  f2 = ek * ( d - sp );
			  f1 = fP + ekunload * deltaD;
			  f  = (f1>f2) ? f1 : f2;
			  if( ekunload < ek) flgstop = 1;
			  if( fabs( f - f1 ) < 1.e-10 ) ek = ekunload ;
			}
		    }
		}
	    }
	  else
	    {
	      if ( d > 0.0) 
		{
		  this->envelPosCap ( fyPos,alfaPos,capSlope,cpPos,d,&f,&ek);
		}
	      else
		{
		  this->envelNegCap (fyNeg,alfaNeg,capSlope,cpNeg,d,&f,&ek);
		}
	    }
	}
    }
  
  // ENDS BIG LOOP ---------------------------------------------------		 	     
  
  // ENERGY CALCULATIONS ---------------------------------------------
  
  Enrgi = 0.5 * ( f + fP ) * deltaD;
  Enrgc = Enrgc + Enrgi;
  Enrgtot = Enrgtot + Enrgi;
  
  
  // UPDATING OF DETERIORATION PARAMETERS ----------------------------
  // Check the remaining capacity of the system
  if ( ecaps != 0.0 || ecapa != 0.0 || ecapd != 0.0 )
    {
      if( ( Enrgtot >= Enrgts && Enrgts != 0.0)||
	  ( Enrgtot >= Enrgtk && Enrgtk != 0.0)||
	  ( Enrgtot >= Enrgta && Enrgta != 0.0)||
	  ( Enrgtot >= Enrgtd && Enrgtd != 0.0))  
	{
	  opserr << "Error: Pinching::setTrialStrain  : Total Energy greater than capacity" << "\a";
	  flgstop = 1;
	}
      
      // Update beta values for strength, acc. stiff. and capping
      
      if(ecaps != 0.0) betas = pow( ( Enrgc / ( Enrgts - Enrgtot ) ) , cs);
      if(ecapa != 0.0) betaa = pow( ( Enrgc / ( Enrgta - Enrgtot ) ) , ca);
      if(ecapd != 0.0) betad = pow( ( Enrgc / ( Enrgtd - Enrgtot ) ) , cd);
      
      if( betas >= 1.0 || betak >= 1.0 || betaa >= 1.0 )   
	{
	  opserr << "Error: Pinching::setTrialStrain  : Beta greater than one" << "\a";
	  flgstop=1;
	}
      
      // Update values for the next half cycle
      ekexcurs = ekunload;
      Enrgc = 0.0;
      
      // Deteriorate parameters for the next half cycle
      
      if( deltaD < 0.0) 
	{
	  if(ecaps != 0.0)
	    {
	      fyNeg = fyNeg * ( 1.0 - betas );
	      fpdegNeg = fpdegNeg * ( 1.0 - betas );
	      alfaNeg = alfaNeg * ( 1.0 - betas );
	    }
	  if(ecapd != 0.0) fCapRefNeg = fCapRefNeg * ( 1.0 - betad );
	  if(ecapa != 0.0) dmin = dmin * ( 1.0 + betaa );
	  
	  dyNeg = fyNeg / elstk;
	  ekhardNeg = alfaNeg * elstk;
	  
	  dCap1Neg = fCapRefNeg / ( elstk - capSlope * elstk );
	  dCap2Neg = ( fCapRefNeg + ekhardNeg * dyNeg - fyNeg ) / ( ekhardNeg - capSlope * elstk );
	  cpNeg = ( dCap1Neg < dCap2Neg ) ? dCap1Neg : dCap2Neg;
	  
	  // else if(deltaD.ge.0.d0)
	}
      else
	{
	  if(ecaps != 0.0)
	    {
	      fyPos = fyPos * ( 1.0 - betas );
	      fpdegPos = fpdegPos * ( 1.0 - betas );
	      alfaPos = alfaPos * ( 1.0 - betas );
	    }
	  if(ecapd != 0.0) fCapRefPos = fCapRefPos * ( 1.0 - betad );
	  if(ecapa != 0.0) dmax = dmax * ( 1.0 + betaa );
	  
	  dyPos = fyPos / elstk;
	  ekhardPos = alfaPos * elstk;
	  
	  dCap1Pos = fCapRefPos / ( elstk - capSlope * elstk );
	  dCap2Pos = ( fCapRefPos + ekhardPos * dyPos - fyPos ) / ( ekhardPos - capSlope * elstk );
	  cpPos= ( dCap1Pos > dCap2Pos ) ? dCap1Pos : dCap2Pos;
	}
    }
  
  // Relationship between basic variables and hstv array	for next cycle
  
  hstv[0] = ek;
  hstv[1] = ekunload;
  hstv[2] = ekexcurs;
  hstv[3] = Enrgtot;
  hstv[4] = Enrgc;
  hstv[5] = sp;
  hstv[6] = sn;
  hstv[7] = (double) kon;
  hstv[8] = dmax;
  hstv[9] = dmin;
  hstv[10] = fyPos;
  hstv[11] = fyNeg;			   
  hstv[12] = cpPos;
  hstv[13] = cpNeg;
  hstv[14] = f;
  hstv[15] = d;
  hstv[16] = fmax;
  hstv[17] = fmin;
  
  return 0;
}


int Pinching::commitState()
{
  if ( DEBG ==1 ) fprintf( OutputFile , "Commit State\n" );	// debugging
  int i;
  for( i=0; i<18; i++ ) hsLastCommit[i] = hstv[i];
  
  this->recordInfo();
  
  return 0;
}


void Pinching::recordInfo(int cond )
{

}


void Pinching::envelPosCap( double fy, double alfaPos, double alfaCap,
			    double cpDsp, double d, double *f, double *ek )
{
  
  double dy,Res,rcap,dres;
  
  dy = fy / elstk;
  
  if (dy < cpDsp)
    {
      Res = Resfac * fyieldPos;
      rcap = fy+alfaPos * elstk * ( cpDsp - dy );
      dres = cpDsp + ( Res - rcap ) / ( alfaCap * elstk );
      
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
		  *ek = elstk * alfaPos;
		  *f = fy + (*ek) * ( d - dy );
		}
	      else
		{
		  if( d <=  dres )
		    {
		      *ek = alfaCap * elstk;
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
      dres = cpDsp + ( Res - rcap ) / ( alfaCap * elstk );
      
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
		  *ek = alfaCap * elstk;
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


void Pinching::envelNegCap( double fy, double alfaNeg, double alfaCap,
			    double cpDsp, double d, double *f, double *ek)
{
  
  double dy,Res,rcap,dres;
  
  dy = fy / elstk;
  
  if( dy > cpDsp )
    {
      
      Res = Resfac * fyieldNeg;
      rcap = fy + alfaNeg * elstk * ( cpDsp - dy );
      dres = cpDsp + ( Res - rcap ) / ( alfaCap * elstk );
      
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
		  *ek = elstk * alfaNeg;
		  *f = fy + (*ek) * ( d - dy );
		}
	      else
		{
		  if ( d >= dres )
		    {
		      *ek = elstk * alfaCap;
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
      dres = cpDsp + ( Res - rcap ) / ( alfaCap * elstk );
      
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
		  *ek = elstk * alfaCap;
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

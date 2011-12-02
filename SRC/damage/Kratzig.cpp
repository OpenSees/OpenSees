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
// $Date: 2008-04-14 22:38:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/damage/Kratzig.cpp,v $
                                                                        
// Written: AA,GGD 
// Created: 08/02
// Revision: AA
//
// Description: This file contains the class implementation for Damage model  .
//

#include <Kratzig.h>
#include <DamageResponse.h>
#include <math.h>

Kratzig::Kratzig( int tag, double ultimatePosVal ,  double ultimateNegVal )
:DamageModel(tag,DMG_TAG_Kratzig),
UltimatePosValue(ultimatePosVal) , UltimateNegValue(ultimateNegVal)
{
	if ( UltimatePosValue <= 0.0 ) 	{
		opserr << "Kratzig::Kratzig : Incorrect arguments for the damage model"<<endln;
		exit (-1);
	}
	if ( UltimateNegValue < 0.0 ) UltimateNegValue = fabs( UltimateNegValue);
	if ( UltimateNegValue == 0.0 ) UltimateNegValue = UltimatePosValue;
	
	this->revertToStart();
}


Kratzig::Kratzig()
:DamageModel(0,DMG_TAG_Kratzig)
{
	// Does nothing
}


Kratzig::~Kratzig()
{
	// Does nothing
}


int
Kratzig::setTrial ( const Vector & trialVector )
{
	// Trial step variables
	double TDefo, TForce, TKunload, TSumPosFHC, TPosPHC, TSumNegFHC, TNegPHC, TMaxPosDefo, TMinNegDefo, TDamage;
	
	// retrieve history data
	double	CDefo			= CommitInfo[0];
	double	CForce			= CommitInfo[1];
	double	CKunload		= CommitInfo[2];
	double	CSumPosFHC		= CommitInfo[3];
	double	CPosPHC			= CommitInfo[4];
	double	CSumNegFHC		= CommitInfo[5];
	double	CNegPHC			= CommitInfo[6];
	double	CMaxPosDefo		= CommitInfo[7];
	double	CMinNegDefo		= CommitInfo[8];
	double	CDamage			= CommitInfo[9];

	if ( trialVector.Size() < 3 ) {
	  opserr << "WARNING: Kratzig::setTrial Wrong vector size for trial data" << endln;
	  return -1;
	}
		
	TDefo		= trialVector(0);
	TForce		= trialVector(1);
	TKunload	= trialVector(2);

	if ( TDefo > 0.0 ) {
	  if ( CDefo < 0.0 ) {
	    
	    double ZeroDispForce;
	    if ( fabs ( TDefo - CDefo ) > 1.0e-6 ) {
	      ZeroDispForce = CForce - (TForce - CForce) * CDefo / (TDefo - CDefo );
	    } else {
	      ZeroDispForce = 0.5*( CForce + TForce );
	    }
	    
	    TSumNegFHC = CSumNegFHC + 0.5*( ZeroDispForce + CForce )*( 0.0 - CDefo );
	    
	    if ( TDefo <= CMaxPosDefo ) {
	      TSumPosFHC = CSumPosFHC + 0.5*( TForce + ZeroDispForce )*( TDefo - 0.0);
	      TPosPHC = CPosPHC;
	    } else
	      {
		TMaxPosDefo = TDefo;
		double MaxDispForce = CForce +(TForce - CForce ) * ( TMaxPosDefo - CDefo)/( TDefo - CDefo );
		TPosPHC = CPosPHC + 0.5*( TForce + MaxDispForce )*( TDefo - TMaxPosDefo);
		TSumPosFHC = CSumPosFHC + 0.5*( MaxDispForce + CForce )*( TMaxPosDefo - CDefo);
	      }
	    
	  } else {
	    
	    if ( TDefo <= CMaxPosDefo ) {
	      TSumPosFHC = CSumPosFHC + 0.5*( TForce + CDefo )*( TDefo - CDefo );
	      TPosPHC = CPosPHC;
	    } else
	      {
		double MaxDispForce = CForce +(TForce - CForce ) * ( CMaxPosDefo - CDefo)/( TDefo - CDefo );
		TPosPHC = CPosPHC + 0.5*( TForce + MaxDispForce )*( TDefo - CMaxPosDefo);
		TSumPosFHC = CSumPosFHC + 0.5*( MaxDispForce + CForce )*( CMaxPosDefo - CDefo);
		TMaxPosDefo = TDefo;
	      }			
	  }
	  
	} else if (TDefo < 0.0 ) {
	  
	  if ( CDefo > 0.0 ) {
	    
	    double ZeroDispForce;
	    if ( fabs( TDefo - CDefo ) > 1.0e-6 ) {
	      ZeroDispForce = CForce - (TForce - CForce) * CDefo / (TDefo - CDefo );
	    } else {
	      ZeroDispForce = 0.5*( CForce + TForce );
	    }
	    
	    TSumPosFHC = CSumPosFHC + 0.5*( ZeroDispForce + CForce )*( 0.0 - CDefo );
	    
	    if ( TDefo >= CMinNegDefo ) {
	      TSumNegFHC = CSumNegFHC + 0.5*( TForce + ZeroDispForce )*( TDefo - 0.0);
	      TNegPHC = CNegPHC;
	    } else
	      {
		TMinNegDefo = TDefo;
		double MinDispForce = CForce +(TForce - CForce ) * ( TMinNegDefo - CDefo)/( TDefo - CDefo );
		TNegPHC = CNegPHC + 0.5*( TForce + MinDispForce )*( TDefo - TMinNegDefo);
		TSumNegFHC = CSumNegFHC + 0.5*( MinDispForce + CForce )*( TMinNegDefo - CDefo);
	      }
	    
	  } else {
	    
	    if ( TDefo >= CMinNegDefo ) {
	      TSumNegFHC = CSumNegFHC + 0.5*( TForce + CDefo )*( TDefo - CDefo );
	      TNegPHC = CNegPHC;
	    } else
	      {
		double MinDispForce = CForce +(TForce - CForce ) * ( CMinNegDefo - CDefo)/( TDefo - CDefo );
		TNegPHC = CNegPHC + 0.5*( TForce + MinDispForce )*( TDefo - CMinNegDefo);
		TSumNegFHC = CSumNegFHC + 0.5*( MinDispForce + CForce )*( CMinNegDefo - CDefo);
		TMinNegDefo = TDefo;
	      }
	  }
	  
	} else { 
	  // TDefo = 0.0
	  if ( CDefo < 0.0 ) {
	    
	    TSumNegFHC = CSumNegFHC + 0.5*( TForce + CForce )*( TDefo - CDefo );
	  } else if ( CDefo > 0.0 ) {
	    TSumPosFHC = CSumPosFHC + 0.5*( TForce + CForce )*( TDefo - CDefo );
	  } else TSumPosFHC = CSumPosFHC;
	}
	
	double PosDamage , NegDamage;
	PosDamage = ( TPosPHC + TSumPosFHC ) / ( UltimatePosValue + TSumPosFHC ) ;
	
	NegDamage = ( TNegPHC + TSumNegFHC ) / ( UltimateNegValue + TSumNegFHC ) ;
	
	TDamage = PosDamage + NegDamage - PosDamage * NegDamage;

	if ( TDamage < CDamage ) TDamage = CDamage;

	
	TrialInfo[0] = TDefo;
	TrialInfo[1] = TForce;
	TrialInfo[2] = TKunload;
	TrialInfo[3] = TSumPosFHC;
	TrialInfo[4] = TPosPHC;
	TrialInfo[5] = TSumNegFHC;
	TrialInfo[6] = TNegPHC;
	TrialInfo[7] = TMaxPosDefo;
	TrialInfo[8] = TMinNegDefo;
	TrialInfo[9] = TDamage;

	return 0;
}


double
Kratzig::getDamage (void)
{
	return TrialInfo[9];
}
    

double Kratzig::getPosDamage (void)
{
	double PosDamage = ( TrialInfo[4] + TrialInfo[3] ) / ( UltimatePosValue + TrialInfo[3] ) ;
	return PosDamage;
}


double Kratzig::getNegDamage (void)
{
	double 	NegDamage = ( fabs(TrialInfo[6]) + fabs(TrialInfo[5]) ) / ( fabs(UltimateNegValue) + fabs(TrialInfo[5]) ) ;
	return NegDamage;
}


int
Kratzig::commitState (void)
{
	for ( int i=0 ; i<10 ; i++ )
	{
		LastCommitInfo[i] = CommitInfo[i];
		CommitInfo[i] = TrialInfo[i];
	}

	return 0;
}

int
Kratzig::revertToLastCommit (void)
{
	for ( int i=0 ; i<10 ; i++ )
	{
		CommitInfo[i] = LastCommitInfo[i];
	}

	return 0;
}

int
Kratzig::revertToStart (void)
{
	for ( int i = 0 ; i< 10 ; i++ ){
		TrialInfo[i] = 0.0;
		CommitInfo[i] = 0.0;
		LastCommitInfo[i] = 0.0;
	}

	return 0;
}


DamageModel*
Kratzig::getCopy (void)
{
	Kratzig *theCopy = new Kratzig(this->getTag() , UltimatePosValue , UltimateNegValue );
    
	for ( int i=0 ; i<10 ; i++ )
	{
		theCopy->TrialInfo[i] = TrialInfo[i];
		theCopy->CommitInfo[i] = CommitInfo[i];
		theCopy->LastCommitInfo[i] = LastCommitInfo[i];
	}

	return theCopy;
}

Response*
Kratzig::setResponse(const char **argv, int argc, OPS_Stream &info)
{
//
// we compare argv[0] for known response types for the Truss
//

	if ( strcmp(argv[0],"damage") == 0 || strcmp(argv[0],"damageindex") == 0 )
    return new DamageResponse( this , 1 , 0.0 );

	else if ( strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformation") == 0 )
    return new DamageResponse( this , 2 , 0.0 );

	else if (strcmp(argv[0],"trial") == 0 || strcmp(argv[0],"trialinfo") == 0 )
    return new DamageResponse( this , 3 , Vector(4) );

	else 
		return 0;

}

int 
Kratzig::getResponse(int responseID, Information &info)
{
	switch (responseID) {
	case -1:
		return -1;
	
	case 1:
		return info.setDouble( this->getDamage() );

	case 2:
		return info.setDouble( TrialInfo[0] );

	case 3:
		if(info.theVector!=0)
		{
			(*(info.theVector))(0) = TrialInfo[4];
			(*(info.theVector))(1) = TrialInfo[3];
			(*(info.theVector))(2) = TrialInfo[6];
			(*(info.theVector))(3) = TrialInfo[5];
		}
		return 0;

	default:
		return -1;
	}
}


int
Kratzig::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}


int
Kratzig::recvSelf(int commitTag, Channel &theChannel,
								FEM_ObjectBroker &theBroker)
{
	return 0;
}


void
Kratzig::Print(OPS_Stream &s, int flag )
{
    s << "CumulativePeak tag: " << this->getTag() << endln;
	s << " UltimatePosValue: " << UltimatePosValue << " UltimateNegValue: " << UltimateNegValue << endln;
}


int Kratzig::setInputResponse ( Element *elem , const char **argv , int argc, int ndof )
{
	return -1;
}

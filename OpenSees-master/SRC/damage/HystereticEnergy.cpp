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
// $Source: /usr/local/cvs/OpenSees/SRC/damage/HystereticEnergy.cpp,v $

// Written: Arash Altoontash, Gregory Deierlein 
// Created: 06/03
// Revision: AA
//
// Description: This file contains the class implementation for HystereticEnergy 
// damage model. HystereticEnergy damage model calculates the damage index based on
// last half cycle energy and the cumulative energy. If an unloading stiffness is
// introduced, the calculation would be based on plastic energy
//

#include <HystereticEnergy.h>
#include <DamageResponse.h>
#include <math.h>

HystereticEnergy::HystereticEnergy (int tag, double Etot , double Cpow)
:DamageModel(tag,DMG_TAG_HystereticEnergy),
Etotal( Etot ) , Cpower( Cpow )
{
  if ( Etot <= 0 || Cpow <= 0 )
    opserr << "DamageModel::DamageModel : Incorrect arguments for the damage model";
  
  this->revertToStart();
}

HystereticEnergy::HystereticEnergy()
:DamageModel(0,DMG_TAG_HystereticEnergy)
{
	// Does nothing
}


HystereticEnergy::~HystereticEnergy()
{
  // Does nothing
}


int
HystereticEnergy::setTrial (const Vector & trialVector )
{	
  double TDisp, TForce, TKunload, TEnrgTot, TEnrgc, TExcurDmg, TCyclicDmg;
  double CDisp, CForce, CKunload, CEnrgTot, CEnrgc, CExcurDmg, CCyclicDmg;
  
  // retrieve history data
  CDisp		= CommitInfo[0];
  CForce		= CommitInfo[1];
  CKunload	= CommitInfo[2];
  CEnrgTot	= CommitInfo[3];
  CEnrgc		= CommitInfo[4];
  CExcurDmg	= CommitInfo[5];
  CCyclicDmg	= CommitInfo[6];
  
  if ( trialVector.Size() < 3 ) {
    opserr << "WARNING: HystereticEnergy::setTrial Wrong vector size for trial data" << endln;
    return -1;
  }
  
  TDisp		= trialVector[0];
  TForce		= trialVector[1];
  TKunload	= trialVector[2];
  
  if ( TKunload < 0.0 ) {
    opserr << "WARNING: HystereticEnergy::setTrial negative unloading stiffness specified" << endln;
    return -1;
  }
  
  if ( TForce == 0.0 )
    {
      // submitt the cyclic damage
      TCyclicDmg = CCyclicDmg + CExcurDmg - CExcurDmg * CCyclicDmg;
      // a new excursion just started
      TEnrgc = 0.0;
      TEnrgTot = CEnrgTot;
    } 
  else 
    if ( CForce * TForce < 0.0 )
      {
	double ZeroForceDisp;
	if ( fabs( CForce + TForce ) < 1.0e-6 )
	  {
	    ZeroForceDisp = 0.5 * ( TDisp + CDisp );
	  }
	else
	  {
	    ZeroForceDisp = ( CForce * TDisp + TForce * CDisp ) / (CForce + TForce );
	  }
	
	// calculate the energies to the end of last cycle
	TEnrgc = CEnrgc + 0.5 * CForce * ( ZeroForceDisp - CDisp );
	TEnrgTot = CEnrgTot + 0.5 * CForce * ( ZeroForceDisp - CDisp );
	TExcurDmg = pow( TEnrgc / ( Etotal - TEnrgTot) , Cpower );
	TCyclicDmg = CCyclicDmg + TExcurDmg - TExcurDmg * CCyclicDmg;
	
	// then add the new cycle
	TEnrgc = 0.5 * TForce * ( TDisp - ZeroForceDisp );
	TEnrgTot = CEnrgTot + 0.5 * TForce * ( TDisp - ZeroForceDisp );
	
      } else {
	TEnrgc = CEnrgc + 0.5 * ( TForce + CForce ) * ( TDisp - CDisp );
	TEnrgTot = CEnrgTot + 0.5 * ( TForce + CForce ) * ( TDisp - CDisp );
	TCyclicDmg = CCyclicDmg;
      }
  
  double RSE = 0.0 ;
  if ( TKunload != 0.0 )
    {
      // Calculate and deduct the elastic energy from the total energy
      RSE =  0.5 * TForce * TForce / TKunload;	
      if ( (TEnrgc - RSE) < 0.0 ) RSE = 0.0;
      if ( (TEnrgTot -RSE) < 0.0 ) RSE = 0.0;
    }
  
  TExcurDmg = pow( (TEnrgc -RSE) / ( (Etotal -RSE) - (TEnrgTot -RSE) ) , Cpower );
  
  
  TrialInfo[0] = TDisp;
  TrialInfo[1] = TForce;
  TrialInfo[2] = TKunload;
  TrialInfo[3] = TEnrgTot;
  TrialInfo[4] = TEnrgc;
  TrialInfo[5] = TExcurDmg;
  TrialInfo[6] = TCyclicDmg;
  return 0;
}


int
HystereticEnergy::setTrial ()
{
  opserr << "WARNING: HystereticEnergy::setTrial Wrong Method called" << endln;
  opserr << "HystereticEnergy Model uses vector based setTrial method" << endln;
  return -1;
}


double
HystereticEnergy::getDamage (void)
{
  TrialInfo[7] =  CommitInfo[6] + TrialInfo[5] - TrialInfo[5] * CommitInfo[6];
  if ( TrialInfo[7] < CommitInfo[7] ) TrialInfo[7] = CommitInfo[7];
  return TrialInfo[7];
}


double HystereticEnergy::getPosDamage (void)
{
  return this->getDamage();
}


double HystereticEnergy::getNegDamage (void)
{
	return this->getDamage();
}


int
HystereticEnergy::commitState (void)
{
	for ( int i=0 ; i<8 ; i++ )
	{
		LastCommitInfo[i] = CommitInfo[i];
		CommitInfo[i] = TrialInfo[i];
	}

	return 0;
}


int
HystereticEnergy::revertToLastCommit (void)
{
	for ( int i=0 ; i<8 ; i++ )
	{
		CommitInfo[i] = LastCommitInfo[i];
	}

	return 0;
}


int
HystereticEnergy::revertToStart (void)
{
	for ( int i = 0 ; i< 8 ; i++ ){
		TrialInfo[i] = 0.0;
		CommitInfo[i] = 0.0;
		LastCommitInfo[i] = 0.0;
	}

	return 0;
}


DamageModel*
HystereticEnergy::getCopy (void)
{
	HystereticEnergy *theCopy = new HystereticEnergy(this->getTag(), Etotal , Cpower);
 
	for ( int i=0 ; i<8 ; i++ )
	{
		theCopy->TrialInfo[i] = TrialInfo[i];
		theCopy->CommitInfo[i] = CommitInfo[i];
		theCopy->LastCommitInfo[i] = LastCommitInfo[i];
	}

	return theCopy;
}



Response*
HystereticEnergy::setResponse(const char **argv, int argc, OPS_Stream  &info)
{
//
// we compare argv[0] for known response types for the Truss
//

  if ( strcmp(argv[0],"damage") == 0 || strcmp(argv[0],"damageindex") == 0 )
    return new DamageResponse( this , 1 , 0.0 );
  
  else if (strcmp(argv[0],"trial") == 0 || strcmp(argv[0],"trialinfo") == 0 )
    return new DamageResponse( this , 2 , Vector(7) );
  
  else 
    return 0;

}


int 
HystereticEnergy::getResponse(int responseID, Information &info)
{
	switch (responseID) {
	case -1:
		return -1;
	
	case 1:
		return info.setDouble( this->getDamage() );

	case 2:
		if(info.theVector!=0)
		{
			for (int i = 0 ; i < 8 ; i++ ) (*(info.theVector))(i) = TrialInfo[i];
		}
		return 0;

	default:
		return -1;
	}
}


int
HystereticEnergy::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}


int
HystereticEnergy::recvSelf(int commitTag, Channel &theChannel,
								FEM_ObjectBroker &theBroker)
{
  return -1;
}


void
HystereticEnergy::Print(OPS_Stream &s, int flag )
{
    s << "HystereticEnergy tag: " << this->getTag() << endln;
    s << "  Etotal: " << Etotal << " Cpower: " << Cpower << endln;
}



int HystereticEnergy::setInputResponse ( Element *elem , const char **argv , int argc, int ndof )
{
	return -1;
}

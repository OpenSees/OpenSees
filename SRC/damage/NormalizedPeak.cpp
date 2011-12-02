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
// $Source: /usr/local/cvs/OpenSees/SRC/damage/NormalizedPeak.cpp,v $

// Written: Arash Altoontash, Gregory Deierlein
// Created: 08/02
// Revision: AA
//
// Description: This file contains the class implementation for Damage model  .
//

#include <NormalizedPeak.h>
#include <DamageResponse.h>
#include <Element.h>
#include <math.h>

NormalizedPeak::NormalizedPeak(int tag, double maxVal, double minVal , const char *argv)
:DamageModel(tag,DMG_TAG_NormalizedPeak), damagetype( NotSpecified ), 
MaxValue(maxVal), MinValue(minVal), TrialDmg(0.0), CommitDmg(0.0), LCommitDmg(0.0),
TrialVector(3), CommitVector(3), LCommitVector(3)
{
  if ( MaxValue < 0.0 || MinValue > 0.0 || argv == NULL )
    {
      opserr << "NormalizedPeak::NormalizedPeak : Incorrect arguments for the damage model";
      exit (-1);
    }
  
  strcpy ( damagename, argv);
  
  
  if ( strcmp( damagename , "force" ) == 0 || strcmp( damagename , "Force" ) == 0 )
    {
      damagetype = Force;
    }
  else if ( strcmp( damagename , "strain" ) == 0 || strcmp( damagename , "Strain" ) == 0 ||
	    strcmp( damagename , "defo" ) == 0 || strcmp( damagename , "deformation" ) == 0 ||
	    strcmp( damagename , "Deformation" ) == 0 )
    {
      damagetype = Deformation;
    }
  else if ( strcmp( damagename , "plasticDefo" ) == 0 || strcmp( damagename , "PlasticDefo" ) == 0 ||
	    strcmp( damagename , "plasticStrain" ) == 0 || strcmp( damagename , "PlasticStrain" ) == 0 ||
	    strcmp( damagename , "plasticDeformation" ) == 0 || strcmp( damagename , "PlasticDeformation" ) == 0 )
    {
      damagetype = PlasticDefo;
    }
  else if ( strcmp( damagename , "energy" ) == 0 || strcmp( damagename , "Energy" ) == 0 ||
	    strcmp( damagename , "totalEnergy" ) == 0 || strcmp( damagename , "TotalEnergy" ) == 0 )
    {
      damagetype = TotalEnergy;
    }
  else if ( strcmp( damagename , "plasticDefo" ) == 0 || strcmp( damagename , "PlasticDefo" ) == 0 ||
	    strcmp( damagename , "plasticStrain" ) == 0 || strcmp( damagename , "PlasticStrain" ) == 0 ||
	    strcmp( damagename , "plasticDeformation" ) == 0 || strcmp( damagename , "PlasticDeformation" ) == 0 ){
    damagetype = 	PlasticEnergy;
  }
  else
    {
      opserr << "NormalizedPeak::NormalizedPeak : The damage type specified is not supported";
      exit (-1);	
    }
  
  this->revertToStart();
}


NormalizedPeak::NormalizedPeak()
  :DamageModel(0,DMG_TAG_NormalizedPeak)
{
  // Does nothing
}

NormalizedPeak::~NormalizedPeak()
{
  // Does nothing
}


DamageModel*
NormalizedPeak::getCopy (void)
{
  NormalizedPeak *theCopy = new NormalizedPeak( this->getTag(), MaxValue, MinValue , damagename);
  
  theCopy->TrialScalar = TrialScalar;
  theCopy->TrialDmg = TrialDmg;
  theCopy->CommitScalar = CommitScalar;
  theCopy->CommitDmg = CommitDmg;
  theCopy->LCommitScalar = LCommitScalar;
  theCopy->LCommitDmg = LCommitDmg;
  
  for ( int i=0 ; i < 3 ; i++ )
    {
      (theCopy->TrialVector)(i) = TrialVector(i);
      (theCopy->CommitVector)(i) = CommitVector(i);
      (theCopy->LCommitVector)(i) = LCommitVector(i);
    } 
  
  return theCopy;
}


    
int
NormalizedPeak::commitState (void)
{
  LCommitScalar = CommitScalar;
  LCommitDmg = CommitDmg;
  LCommitVector = CommitVector;
  
  CommitScalar = TrialScalar;
  CommitDmg = TrialDmg;
  CommitVector = TrialVector;
  
  return 0;
}

int
NormalizedPeak::revertToLastCommit (void)
{
  CommitScalar = LCommitScalar;
  CommitDmg = LCommitDmg;
  
  CommitVector = LCommitVector;
  
  return 0;
}

int
NormalizedPeak::revertToStart (void)
{
  TrialScalar = CommitScalar = LCommitScalar = 0.0;
  TrialDmg = CommitDmg = LCommitDmg = 0.0;
  TrialVector.Zero();
  CommitVector.Zero();
  LCommitVector.Zero();
  
  return 0;
}


int
NormalizedPeak::setTrial (const Vector &trialVector )
{
  if ( trialVector.Size() < 3 ) {
    opserr << "WARNING: NormalizedPeak::setTrial Wrong vector size for trial data" << endln;
    return -1;
  }
  
  TrialVector = trialVector;
  
  double TrialDefo = trialVector(0);
  double TrialForce = trialVector(1);
  double TrialKU = trialVector(2);
  
  TrialScalar = 0.0;
  
  switch( damagetype )
    {
    case Force:
      TrialScalar = TrialVector(1);
      break;
    case Deformation:
      TrialScalar = TrialVector(0);
      break;
    case PlasticDefo:
      if ( TrialVector(2) != 0.0 )
	{ TrialScalar = TrialVector(0) - TrialVector(1)/TrialVector(2); }
      else { TrialScalar = TrialVector(0); }
      break;
    case TotalEnergy:
      TrialScalar = CommitScalar + 0.5*( TrialVector(1) + CommitVector(1) )*( TrialVector(0) - CommitVector(0) );
      break;
    case PlasticEnergy:
      if ( TrialVector(2) > 0.0 )
	TrialScalar = CommitScalar + 0.5*( TrialVector(1) + CommitVector(1) )*( TrialVector(0) - CommitVector(0) )
	  - 0.5* TrialVector(1) * TrialVector(1) / TrialVector(2);
      break;
    }
  
  // Now calculate damage
  if ( TrialScalar >= 0.0 ) {
    TrialDmg = TrialScalar / MaxValue;
  } else {
    TrialDmg = fabs( TrialScalar / MinValue );
  }
  
  if ( fabs(TrialDmg) < CommitDmg ) TrialDmg = CommitDmg;
  
  return 0;
}


double
NormalizedPeak::getDamage (void)
{
  return TrialDmg;
}



double NormalizedPeak::getPosDamage (void)
{
  return TrialDmg;
}


double NormalizedPeak::getNegDamage (void)
{
  return TrialDmg;
}



int
NormalizedPeak::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}


int
NormalizedPeak::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  return 0;
}


void
NormalizedPeak::Print(OPS_Stream &s, int flag )
{
    s << "NormalizedPeak tag: " << this->getTag() << endln;
    s << "  MaximumValue: " << MaxValue << " MinimumValue: " << MinValue << endln;
    s << " Response type: " <<damagename << endln;
    
}

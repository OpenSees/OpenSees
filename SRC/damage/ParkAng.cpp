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
// $Source: /usr/local/cvs/OpenSees/SRC/damage/ParkAng.cpp,v $

// Written: AA,GGD 
// Created: 10/02
// Revision: AA
//
// Description: This file contains the class implementation for ParkAng 
// damage model. ParkAng damage model calculates the damage index based on
// primary half cycle and the summation of follower half cycles.
//

#include <ParkAng.h>
#include <DamageResponse.h>
#include <math.h>

ParkAng::ParkAng (int tag, double deltaU , double beta , double sigmaY )
:DamageModel(tag,DMG_TAG_ParkAng),
DeltaU( deltaU ) , Beta( beta ) , SigmaY( sigmaY )
{
  if ( DeltaU <= 0 || Beta <= 0 , SigmaY <= 0 )
    opserr << "ParkAng::ParkAng : Incorrect arguments for the damage model"<<endln;
  
  this->revertToStart();
}


ParkAng::ParkAng()
:DamageModel(0,DMG_TAG_ParkAng)
{
	// Does nothing
}

ParkAng::~ParkAng()
{
	// Does nothing
}

int
ParkAng::setTrial (double scalar, double scalarRate )
{
	opserr << "WARNING: ParkAng::setTrial Wrong Method called" << endln;
	opserr << "ParkAng Model uses vector based setTrial method" << endln;
	return -1;
}

int
ParkAng::setTrial (const Vector & trialVector )
{
	// Trial step
	double TForce, TDeformation, TUnloadingK, TEnergy, TMaxDefo, TDamage;

	// Commited state
	double CForce		= CommitInfo[0];;
	double CDeformation	= CommitInfo[1];;
	double CUnloadingK	= CommitInfo[2];;
	double CEnergy		= CommitInfo[3];;
	double CMaxDefo		= CommitInfo[4];;
	double CDamage		= CommitInfo[5];;
	
	
	// Deformation = trialVector(0);
	// Force = trialVector(1);
	//
	if ( trialVector.Size() != 3 ) {
		opserr << "WARNING: ParkAng::setTrial Wrong vector size for trial data" << endln;
		return -1;
	}

	TDeformation = trialVector(0);
	TForce = trialVector(1);
	TUnloadingK = trialVector(2);
	
	if ( TUnloadingK < 0.0 ) {
		opserr << "WARNING: ParkAng::setTrial negative unloading stiffness specified" << endln;
		return -1;
	}

	TEnergy = CEnergy + 0.5 * ( TForce + CForce ) * ( TDeformation - CDeformation );

	double PlasticEnergy;
	if (TUnloadingK != 0.0 ) {
		PlasticEnergy = TEnergy - 0.5 * TForce * TForce / TUnloadingK;
	} else {
		PlasticEnergy = TEnergy;
	}

	TMaxDefo = ( fabs( TDeformation ) > fabs( CMaxDefo ) ) ? fabs(TDeformation) : fabs(CMaxDefo);

	TDamage = ( TMaxDefo / DeltaU ) + ( Beta * PlasticEnergy / SigmaY / DeltaU );
	if ( TDamage < CDamage )  TDamage = CDamage;

	// Trial step
	TrialInfo[0] =  TForce;
	TrialInfo[1] =  TDeformation;
	TrialInfo[2] =  TUnloadingK;
	TrialInfo[3] =  TEnergy;
	TrialInfo[4] =  TMaxDefo;
	TrialInfo[5] =  TDamage;

	return 0;
}


int
ParkAng::setTrial ()
{
	opserr << "WARNING: ParkAng::setTrial Wrong Method called" << endln;
	opserr << "ParkAng Model uses vector based setTrial method" << endln;
	return -1;
}


double
ParkAng::getDamage (void)
{
	return TrialInfo[5];
}


double ParkAng::getPosDamage (void)
{
	return TrialInfo[5];
}


double ParkAng::getNegDamage (void)
{
	return TrialInfo[5];
}

    
int
ParkAng::commitState (void)
{
	for ( int i=0 ; i<6 ; i++ )
	{
		LastCommitInfo[i] = CommitInfo[i];
		CommitInfo[i] = TrialInfo[i];
	}

	return 0;
}

int
ParkAng::revertToLastCommit (void)
{
	for ( int i=0 ; i<6 ; i++ )
	{
		CommitInfo[i] = LastCommitInfo[i];
	}

	return 0;
}

int
ParkAng::revertToStart (void)
{
	for ( int i = 0 ; i< 6 ; i++ ){
		TrialInfo[i] = 0.0;
		CommitInfo[i] = 0.0;
		LastCommitInfo[i] = 0.0;
	}

	return 0;
}

DamageModel*
ParkAng::getCopy (void)
{
	ParkAng *theCopy = new ParkAng(this->getTag(), DeltaU , Beta , SigmaY );
    
 	for ( int i=0 ; i<6 ; i++ )
	{
		theCopy->TrialInfo[i] = TrialInfo[i];
		theCopy->CommitInfo[i] = CommitInfo[i];
		theCopy->LastCommitInfo[i] = LastCommitInfo[i];
	}

	return theCopy;
}

Response*
ParkAng::setResponse(const char **argv, int argc, OPS_Stream &info)
{
//
// we compare argv[0] for known response types for the Truss
//

	if ( strcmp(argv[0],"damage") == 0 || strcmp(argv[0],"damageindex") == 0 )
    return new DamageResponse( this , 1 , 0.0 );

	else if (strcmp(argv[0],"Value") == 0 || strcmp(argv[0],"Values") == 0 || strcmp(argv[0],
		"Data") == 0 )
    return new DamageResponse( this , 2 , Vector(3) );

	else if (strcmp(argv[0],"trial") == 0 || strcmp(argv[0],"trialinfo") == 0 )
    return new DamageResponse( this , 3 , Vector(6) );

	else 
		return 0;

}

int 
ParkAng::getResponse(int responseID, Information &info)
{
	switch (responseID) {
	case -1:
		return -1;
	
	case 1:
		return info.setDouble( this->getDamage() );

	case 2:
		if(info.theVector!=0)
		{
			(*(info.theVector))(0) = TrialInfo[1];
			(*(info.theVector))(1) = TrialInfo[0];
			(*(info.theVector))(2) = TrialInfo[2];
		}
		return 0;

	case 3:
		if(info.theVector!=0)
		{
			(*(info.theVector))(0) = TrialInfo[0];
			(*(info.theVector))(1) = TrialInfo[1];
			(*(info.theVector))(2) = TrialInfo[2];
			(*(info.theVector))(3) = TrialInfo[3];
			(*(info.theVector))(4) = TrialInfo[4];
			(*(info.theVector))(5) = TrialInfo[5];
		}

		return 0;

	default:
		return -1;
	}
}


int
ParkAng::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}


int
ParkAng::recvSelf(int commitTag, Channel &theChannel,
								FEM_ObjectBroker &theBroker)
{
	return 0;
}


void
ParkAng::Print(OPS_Stream &s, int flag )
{
    s << "ParkAng tag: " << this->getTag() << endln;
    s << "  DeltaU: " << DeltaU << " Beta: " << Beta << "  SigmaY: " << SigmaY << endln;
}



int ParkAng::setInputResponse ( Element *elem , const char **argv , int argc, int ndof )
{
	return -1;
}

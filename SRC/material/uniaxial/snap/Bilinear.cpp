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
// $Date: 2008-04-14 21:26:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/Bilinear.cpp,v $
//
//
// Bilinear.cpp: implementation of the Bilinear class inspried by Fortran version.
// Originally from SNAP PROGRAM by Luis Ibarra and Prof H.K. Krawinkler
//
// Written: A. Altoontash & Prof. G. Deierlein 05/03
// Revised: 05/05
//
// Purpose: This file contains the implementation for the Bilinear class.
//
//////////////////////////////////////////////////////////////////////

#include <Bilinear.h>
#include <Channel.h>
#include <Parameter.h>

#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DEBG 0

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Bilinear::Bilinear(int tag, Vector inputParam  ,DamageModel *strength,DamageModel *stiffness,DamageModel *capping)
  :UniaxialMaterial(tag,MAT_TAG_SnapBilinear),StrDamage(0),StfDamage(0),CapDamage(0)
{
	int ErrorFlag =0;
	//*********************************************************************
	//	CONSTRUCTING VARIABLES
	//
	//	elstk		=	Initial elastic stiffness
	//	fyieldPos	=	Positive yield strength
	//	fyieldNeg	=	Negative yield strength
	//	alfa		=	Strain-hardening ratio (fraction of elstk)
	//	alfaCap		=	Cap slope ratio (fraction of elstk)
	//	capDispPos	=	Cap displacement on positive side
	//  capDispNeg	=	Cap displacement on negative side
	//	flagCapenv  =	A flag to establish a cut-off limit for the force once the cap is hit
	//  Resfac		=	Residual stress ratio (fraction of yield strength)
	//	strength	=	A pointer to strength damage model
	//	stiffness	=	A pointer to stiffness damage model
	//	capping		=	A pointer to capping damage model
	//
	//  INPUT VARIABLES
	//  d			= Current Displcamenet
	//
	//
	//  OUTPUT VARIABLE
	//	f			= Calculated force response
	//	ek			= Calculated stiffness
	//	
	//
	//*********************************************************************
	
	if( (inputParam.Size()) < 9) 
		opserr << "Error: Bilinear(): inputParam, size <15\n" << "\a";

	elstk		= inputParam[0];
	fyieldPos	= inputParam[1];
	fyieldNeg	= inputParam[2];
	alfa		= inputParam[3];
	alfaCap		= inputParam[4];
	capDispPos	= inputParam[5];
	capDispNeg	= inputParam[6];
	flagCapenv  = (int)inputParam[7];
	Resfac		= inputParam[8];


	// Error check
	
	if ( fyieldPos <= 0.0 || fyieldNeg >= 0.0 )
    {
		opserr << "Error: Bilinear::Bilinear  : Incorrect yield stresse \n" << "\a";
		ErrorFlag =1;
    }

	if ( elstk <= 0.0 )
    {
		opserr << "Error: Bilinear::Bilinear  : Elastic modulus must be positive\n" << "\a";
		ErrorFlag =1;
    }
	
	if ( alfa < 0.0 || alfa > 0.8 )
    {
		opserr << "Error: Bilinear::Bilinear  : alpha is recommended to be in the range of [0.0 , 0.8]\n" << "\a";	
    }
	
	if ( alfaCap >= 0.0 || alfaCap == alfa )
    {
		opserr << "Error: Bilinear::Bilinear  : CapSlope must be negative and not equal to alfa\n" << "\a";
		ErrorFlag =1;
    }	
	
	if ( capDispPos < fyieldPos/elstk || capDispNeg > fyieldNeg/elstk )
    {
		opserr << "Error: Bilinear::Bilinear  : Capping branch must be located outside the yield criteria\n" << "\a";
		ErrorFlag =1;
    }
	
	if ( Resfac <  0.0  || Resfac > 1.0)
    {
		opserr << "Error: Bilinear::Bilinear  : Residual must be positive and less than 1.0\n" << "\a";
		ErrorFlag =1;
    }
	
	if ( DEBG ==1 )
    {
		// Open an output file for debugging
		char FileName[20];						// debugging
		sprintf(FileName, "Bilinear%d.out", tag);
		OutputFile = fopen ( FileName , "w" );	// debugging
		fprintf( OutputFile , "Constructor called\n" );	// debugging
		fprintf( OutputFile , "elstk = %f\n",inputParam[0]);
		fprintf( OutputFile , "fyieldPos = %f\n",inputParam[1]);
		fprintf( OutputFile , "fyieldNeg = %f\n",inputParam[2]);
		fprintf( OutputFile , "alfa = %f\n",inputParam[3]);
		fprintf( OutputFile , "alfaCap = %f\n",inputParam[4]);
		fprintf( OutputFile , "capDispPos = %f\n",inputParam[5]);
		fprintf( OutputFile , "capDispNeg = %f\n",inputParam[6]);
		fprintf( OutputFile , "flagCapenv = %d\n",(int) inputParam[7]);
		fprintf( OutputFile , "Resfac = %f\n",inputParam[8]);
    }
	
	
	if ( DEBG ==2 )
    {
		// Open an output file for debugging
		char FileName[20];						// debugging
		sprintf(FileName, "Bilinear%d.out", tag);
		OutputFile = fopen ( FileName , "w" );	// debugging
    }

	if ( ErrorFlag == 1 )    {
		opserr << "Error: Bilinear::Bilinear  : Error: check the input values\n" << "\a";	
		exit(-1);
	}

	
	if ( strength != NULL )
	{
		StrDamage = strength->getCopy();
		if ( StrDamage == NULL ) {
			opserr << "Error: Bilinear::Bilinear  : Can not make a copy of strength damage model\n" << "\a";	
			exit(-1);
		}
	}
	
	if ( stiffness != NULL )
	{
		StfDamage = stiffness->getCopy();
		if ( StfDamage == NULL ) {
			opserr << "Error: Bilinear::Bilinear  : Can not make a copy of stiffness damage model\n" << "\a";	
			exit(-1);
		}
	}	

	if ( capping != NULL )
	{
		CapDamage = capping->getCopy();
		if ( CapDamage == NULL ) {
			opserr << "Error: Bilinear::Bilinear  : Can not make a copy of capping damage model\n" << "\a";	
			exit(-1);
		}
	}
	
	// Initialice history data
	this->revertToStart();
}


Bilinear::Bilinear()
  :UniaxialMaterial(0,MAT_TAG_SnapBilinear)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Empty constructor called\n" );	// debugging
}


Bilinear::~Bilinear()
{
	if ( DEBG == 1 || DEBG == 2 )
    {
		fprintf( OutputFile , "Distructor called\n" );		// debugging
		fclose( OutputFile );						// debugging
    }

	if ( StrDamage != 0 ) delete StrDamage;
	if ( StfDamage != 0 ) delete StfDamage;
	if ( CapDamage != 0 ) delete CapDamage;
}


int Bilinear::revertToStart()
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Revert to start\n" );		// debugging

	hsLastCommit[0] =  0.0;													// dP
	hsLastCommit[1] =  0.0;													// fP
	hsLastCommit[2] =  elstk;												// ekP
	hsLastCommit[3] =  elstk;												// ekexcurs
	hsLastCommit[4] =  fyieldPos;											// fyPos
	hsLastCommit[5] =  fyieldNeg;											// fyNeg
	hsLastCommit[6] =  alfa*elstk;											// ekhard
	hsLastCommit[7] =  capDispPos;											// cpPos
	hsLastCommit[8] =  capDispNeg;											// cpNeg
	hsLastCommit[9] =  alfaCap*elstk;										// ekcap
	hsLastCommit[10]=  0.0;													// dmax
	hsLastCommit[11]=  0.0;													// dmin
	hsLastCommit[12]=  fyieldPos+alfa*elstk*(capDispPos-fyieldPos/elstk);	// fuPos
	hsLastCommit[13]=  fyieldNeg+alfa*elstk*(capDispNeg-fyieldNeg/elstk);	// fuNeg
	hsLastCommit[14]=  0.0;													// Enrgtot
	hsLastCommit[15]=  0.0;													// Enrgc
	hsLastCommit[16]=  0.0;													// reserved

	for( int i=0 ; i<17 ; i++) {
		hsCommit[i] = hsLastCommit[i];
		hsTrial[i] = hsLastCommit[i];
	}
	if ( StrDamage != NULL ) StrDamage->revertToStart();
	if ( StfDamage != NULL ) StfDamage->revertToStart();
	if ( CapDamage != NULL ) CapDamage->revertToStart();
	
	if ( DEBG == 1 || DEBG == 2 )
    {
		if ( StrDamage != NULL ) fprintf( OutputFile , "%d" ,StrDamage->getDamage() );		// debugging
		if ( StfDamage != NULL ) fprintf( OutputFile , "\t%d" , StfDamage->getDamage() );		// debugging
		if ( CapDamage != NULL ) fprintf( OutputFile , "\t%d\n" , CapDamage->getDamage() );		// debugging
    }

	return 0;
}


void Bilinear::Print(OPS_Stream &s, int flag)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Print\n" );	// debugging
	s << "Bilinear Tag: " << this->getTag() << endln;
	s << "d : " << hsTrial[0] << endln;
	s << "f : " << hsTrial[1] << endln;
	s << "ek: " << hsTrial[2] << endln;
	s << endln;
}


int Bilinear::revertToLastCommit()
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Revert to last commit\n" );	// debugging
	
	for(int i=0; i<17; i++) {
		hsTrial[i] = hsCommit[i];
		hsCommit[i] = hsLastCommit[i];
	}
	if ( StrDamage != NULL ) StrDamage->revertToLastCommit();
	if ( StfDamage != NULL ) StfDamage->revertToLastCommit();
	if ( CapDamage != NULL ) CapDamage->revertToLastCommit();
	
	return 0;
}


int Bilinear::commitState()
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Commit state\n" );	// debugging
	
	for(int i=0; i<17; i++) {
		hsLastCommit[i] = hsCommit[i];
		hsCommit[i] = hsTrial[i];
	}

		// Calling the damage object
	Vector InforForDamage(3);
	InforForDamage(0) = hsCommit[0];
	InforForDamage(1) = hsCommit[1];
	InforForDamage(2) = hsCommit[3];

	if ( StrDamage != NULL ) {
		StrDamage->setTrial(InforForDamage);
		StrDamage->commitState();
	}

	if ( StfDamage != NULL ) {
		StfDamage->setTrial(InforForDamage);
		StfDamage->commitState();
	}

	if ( CapDamage != NULL ) {
		CapDamage->setTrial(InforForDamage);
		CapDamage->commitState();
	}

	return 0;
}


double Bilinear::getTangent()
{
	if ( DEBG ==1 ) {
		fprintf( OutputFile , "Get tangent\n" );
		fprintf( OutputFile , "tangent = %f\n",hsTrial[2]);
	} // debugging
	return hsTrial[2];
}

double Bilinear::getInitialTangent (void)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Get initial tangent\n" );	// debugging
	return elstk;
}

double Bilinear::getStress()
{
	if ( DEBG ==1 ) {
		fprintf( OutputFile , "Get stress\n" );
		fprintf( OutputFile , "Stress = %f\n",hsTrial[1]);
	}// debugging
	return hsTrial[1];
}


double Bilinear::getStrain (void)
{
	if ( DEBG ==1 ) {
		fprintf( OutputFile , "Get strain\n" );
		fprintf( OutputFile , "Strain = %f\n",hsTrial[0]);
	}	// debugging
	return hsTrial[0];
}


int Bilinear::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Receive self\n" );	// debugging
	return 0;
}


int Bilinear::sendSelf(int cTag, Channel &theChannel)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Send self\n" );	// debugging
	return 0;
}


UniaxialMaterial *Bilinear::getCopy(void)
{
	if ( DEBG ==1 ) fprintf( OutputFile , "Get copy\n" );	// debugging
	Vector inp(9);
	
	inp[0]  = elstk;
	inp[1]  = fyieldPos;
	inp[2]  = fyieldNeg;
	inp[3]  = alfa;
	inp[4]  = alfaCap;
	inp[5]  = capDispPos;
	inp[6]  = capDispNeg;
	inp[7]  = flagCapenv;
	inp[8]  = Resfac;

	
	Bilinear *theCopy = new Bilinear(this->getTag(), inp ,StrDamage,StfDamage,CapDamage);
	
	for (int i=0; i<17; i++) {
		theCopy->hsTrial[i] = hsTrial[i];
		theCopy->hsCommit[i] = hsCommit[i];
		theCopy->hsLastCommit[i] = hsLastCommit[i];
	}

	return theCopy;
}


int Bilinear::setTrialStrain( double d, double strainRate)
{
	if ( DEBG ==1 ){
		fprintf( OutputFile , "Set trial displacement to %f\n" ,d);

	}
	// debugging
	
	// HYSTERETIC VARIABLES
	double dP,fP,ekP,ekexcurs,fyPos,fyNeg,ekhard,cpPos,cpNeg,ekcap,dmax,dmin;
	double fuPos,fuNeg,Enrgtot,Enrgc;

	// LOCAL VARIABLES
	double deltaD, f, ek, fenvPos, fenvNeg, ekenvPos, ekenvNeg;

	// Relationship between basic variables and hsTrial array --------------
	dP					= hsCommit[0];
	fP					= hsCommit[1];
	ekP					= hsCommit[2];
	ekexcurs			= hsCommit[3];
	fyPos				= hsCommit[4];
	fyNeg				= hsCommit[5];
	ekhard				= hsCommit[6];
	cpPos				= hsCommit[7];
	cpNeg				= hsCommit[8];
	ekcap				= hsCommit[9];
	dmax				= hsCommit[10];
	dmin				= hsCommit[11];
	fuPos				= hsCommit[12];
	fuNeg				= hsCommit[13];
	Enrgtot				= hsCommit[14];
	Enrgc				= hsCommit[15];

	// Initialization of variables (each time the subroutine is called) 
	deltaD = d-dP;

	// Calculate the displacement envelope for maximum and minimum strain
	if ( d > dmax ) dmax = d;
	if ( d < dmin ) dmin = d;

	// Check for a new excursion to degrade the parameters
	if ( (fP + ekexcurs * deltaD) * fP <= 0.0 )
	{		
		// degrade the model paraeters based on the damage
		// degrade the strength ( yield stress )
		if ( StrDamage != NULL )
		{
			double StrengthResidual = ( 1.0 - StrDamage->getDamage() );
			if ( StrengthResidual < 0.0 ) StrengthResidual = 0.0;
			
			fyPos = StrengthResidual * ( fyieldPos - Resfac*fyieldPos) + Resfac*fyieldPos;
			fyNeg = StrengthResidual * ( fyieldNeg - Resfac*fyieldNeg) + Resfac*fyieldNeg;
		}

		// degrade the stiffness
		if ( StfDamage != NULL )
		{
			double StiffnessResidual = ( 1.0 - StfDamage->getDamage() );
			if ( StiffnessResidual < 0.0 ) StiffnessResidual = 0.0;

			ekexcurs = StiffnessResidual * ( elstk - alfa*elstk ) + alfa*elstk;
		}
		
		// degrade the cap
		if ( CapDamage != NULL )
		{
			double CapRefPos = ( fyieldPos + (capDispPos - fyieldPos/elstk) * alfa*elstk - Resfac * fyieldPos ) / (alfaCap*elstk);
			double CapRefNeg = ( fyieldNeg + (capDispNeg - fyieldNeg/elstk) * alfa*elstk - Resfac * fyieldNeg ) / (alfaCap*elstk);

				
			double CapPosResidual = ( 1.0 - CapDamage->getPosDamage() );
			if ( CapPosResidual < 0.0 ) CapPosResidual = 0.0;
			double CapNegResidual = ( 1.0 - CapDamage->getNegDamage() );
			if ( CapNegResidual < 0.0 ) CapNegResidual = 0.0;

			cpPos = CapPosResidual * ( capDispPos - CapRefPos ) + CapRefPos;
			cpPos = ( cpPos > CapRefPos ) ? cpPos : CapRefPos;

			cpNeg = CapNegResidual * ( capDispNeg - CapRefNeg ) + CapRefNeg;
			cpNeg = ( cpNeg < CapRefNeg ) ? cpNeg : CapRefNeg;
		}
	}

	// predict the f based on excurs stiffness
	f = fP + ekexcurs * deltaD;
	ek = ekenvPos = ekenvNeg = ekexcurs;

	// Calculate the envelopes
	if ( f >= 0 )
	{
		this->envelPosCap( ekexcurs, fyPos, ekhard, cpPos, ekcap, Resfac*fyieldPos, &fuPos, d, &fenvPos, &ekenvPos );
		fenvNeg = 0.0; 
	} else
	{
		this->envelNegCap( ekexcurs, fyNeg, ekhard, cpNeg, ekcap, Resfac*fyieldNeg, &fuNeg, d, &fenvNeg, &ekenvNeg );
		fenvPos = 0.0; 
	}

	if ( DEBG == 2 ) fprintf( OutputFile , "%f  %f  %f\n" , d , cpPos, fuPos );	// debugging

	// Compare the predictor force with the envelope and correct the 
	if ( f > fenvPos )
	{
		f = fenvPos;
	} 
	else if ( f < fenvNeg )
	{
		f = fenvNeg;
	}


	if ( flagCapenv == 1 )
	{
		if ( f > fuPos )
		{
			f = fuPos;
		} 
		else if ( f < fuNeg )
		{
			f = fuNeg;
		}
	}

	if ( deltaD != 0.0 ) ek = ( f - fP ) / deltaD;

	// Relationship between basic variables and hsTrial array	for next cycle
	hsTrial[0] = d;
	hsTrial[1] = f;
	hsTrial[2] = ek;
	hsTrial[3] = ekexcurs;
	hsTrial[4] = fyPos;
	hsTrial[5] = fyNeg;
	hsTrial[6] = ekhard;
	hsTrial[7] = cpPos;
	hsTrial[8] = cpNeg;
	hsTrial[9] = ekcap;
	hsTrial[10] = dmax;
	hsTrial[11] = dmin;
	hsTrial[12] = fuPos;
	hsTrial[13] = fuNeg;
	hsTrial[14] = Enrgtot;
	hsTrial[15] = Enrgc;
	hsTrial[16] = 0.0;

	return 0;
}


Response* 
Bilinear::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{
  Response *theResponse = 0;

  if ( argv == NULL || argc == 0 ) {
    opserr << "Error: Bilinear::setResponse  : No argument specified\n" << "\a";
    return 0;
  };
  
  theOutput.tag("UniaxialMaterialOutput");
  theOutput.attr("matType", this->getClassType());
  theOutput.attr("matTag", this->getTag());
  
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"stress") == 0 ) {
    theOutput.tag("ResponseType", "sigma11");
    theResponse = new MaterialResponse(this, 1, 0.0);
  }

  else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformation") == 0 ||
	   strcmp(argv[0],"strain") == 0) {
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 2, 0.0);
  }

  else if (strcmp(argv[0],"plastic") == 0 || strcmp(argv[0],"plasticdefo") == 0 ||
	   strcmp(argv[0],"plasticdeformation") == 0 || strcmp(argv[0],"plasticstrain") == 0) {
    theOutput.tag("ResponseType", "eps1P");
    theResponse =  new MaterialResponse(this, 3, 0.0);
  }

  else if ( (strcmp(argv[0],"stiff") == 0) || (strcmp(argv[0],"stiffness") == 0) ) {
    theOutput.tag("ResponseType", "C11");
    theResponse =  new MaterialResponse(this, 4, 0.0);
  }
	
  else if ( (strcmp(argv[0],"unloading") == 0) || (strcmp(argv[0],"unloadingstiffness") == 0)
	    || (strcmp(argv[0],"unloadingstiff") == 0 ) ) {
    theOutput.tag("ResponseType", "C11_unloading");
    theResponse =  new MaterialResponse(this, 5, 0.0);
  }
  
  else if ( (strcmp(argv[0],"damage") == 0) || (strcmp(argv[0],"damages") == 0)
	    || (strcmp(argv[0],"Damage") == 0 ) || (strcmp(argv[0],"Damages") == 0 ) ) {
    theOutput.tag("ResponseType", "str_damaga");
    theOutput.tag("ResponseType", "stf_damaga");
    theOutput.tag("ResponseType", "cap_damaga");
    theResponse =  new MaterialResponse(this, 6, Vector(3));
  }

  theOutput.endTag();
  return theResponse;    
}

int Bilinear::getResponse(int responseID, Information &matInfo)
{
  switch (responseID) {
    case 1:
		return matInfo.setDouble( hsTrial[1] );
			
    case 2:
		return matInfo.setDouble( hsTrial[0] );
		
    case 3:
		return matInfo.setDouble( hsTrial[0]  - hsTrial[1] / hsTrial[3]);
		
	case 4:
		return matInfo.setDouble( hsTrial[2] );
		
	case 5:
		return matInfo.setDouble( hsTrial[3] );

	case 6:
		(*(matInfo.theVector))(0) = 0.0;
		(*(matInfo.theVector))(1) = 0.0;
		(*(matInfo.theVector))(2) = 0.0;
		if ( StrDamage != NULL ) (*(matInfo.theVector))(0) = StrDamage->getDamage();
		if ( StfDamage != NULL ) (*(matInfo.theVector))(1) = StfDamage->getDamage();
		if ( CapDamage != NULL ) (*(matInfo.theVector))(2) = CapDamage->getDamage();
	return 0;
			
    default:
		return 0;
  }
}


void Bilinear::recordInfo(int cond )
{

}


void Bilinear::envelPosCap( double ekelstk, double fy, double ekhard, double dcap,
						   double ekcap, double fRes, double *fuPos, double d, double *f, double *ek )
{
	double dy, fucap, dRes, dmin;
	
	dy = fy / ekelstk;
	dmin = dy - ( (fy-fRes) / ekhard);
//	dcap = ( dcap > dy ) ? dcap : dy;
	fucap = fRes + (dcap - dmin) * ekhard;
	dRes = dcap + ( fRes - fucap ) / ekcap;

	
	if ( DEBG ==1 )
    {		
		fprintf( OutputFile , "Positive envelope called\n" );	// debugging
		fprintf( OutputFile , "dmin = %f\n",dmin);
		fprintf( OutputFile , "dy = %f\n",dy);
		fprintf( OutputFile , "fy = %f\n",fy);
		fprintf( OutputFile , "dcap = %f\n",dcap);
		fprintf( OutputFile , "fucap = %f\n",fucap);
		fprintf( OutputFile , "dRes = %f\n",dRes);
		fprintf( OutputFile , "fRes = %f\n",fRes);
    }


	if ( d < dmin ){
		*f = fRes;
		*ek = 0.0;
	}
	else if ( d < dcap ){
		*f = ekhard * ( d - dmin ) + fRes;
		*ek = ekhard;
	}
	else if ( d < dRes ){
		*f = fucap + ekcap * (d - dcap);
		*ek = ekcap;
		if ( *f < *fuPos ) *fuPos = *f;
	}
	else{
		*f = fRes;
		*ek = 0.0;
		*fuPos = fRes;
	}
	return;
}


void Bilinear::envelNegCap( double ekelstk, double fy, double ekhard, double dcap,
						   double ekcap, double fRes, double *fuNeg,double d, double *f, double *ek )
{
	if ( fy > 0.0 || fRes > 0.0 ) {
		opserr <<" Error : Bilinear::envelNegCap wrong parameters in function call";
		exit(-1);
	} 

	double dy, fucap, dRes, dmax;
	
	dy = fy / ekelstk;
	dmax = dy - ( (fy-fRes) / ekhard);
//	dcap = ( dcap < dy ) ? dcap : dy;
	fucap = fRes + (dcap - dmax)*ekhard;
	dRes = dcap + ( fRes - fucap ) / ekcap;
	
	if ( DEBG ==1 )
    {		
		fprintf( OutputFile , "Negative envelope called\n" );	// debugging
		fprintf( OutputFile , "dRes = %f\n",dRes);
		fprintf( OutputFile , "fRes = %f\n",fRes);
		fprintf( OutputFile , "dcap = %f\n",dcap);
		fprintf( OutputFile , "fucap = %f\n",fucap);
		fprintf( OutputFile , "dy = %f\n",dy);
		fprintf( OutputFile , "fy = %f\n",fy);
		fprintf( OutputFile , "dmax = %f\n",dmax);
    }
		
	if ( d > dmax ){
		*f = fRes;
		*ek = 0.0;
	}
	else if ( d > dcap ){
		*f = ekhard * ( d - dmax ) + fRes;
		*ek = ekhard;
	}
	else if ( d > dRes ){
		*f = fucap + ekcap * (d - dcap);
		*ek = ekcap;
		if ( *f > *fuNeg ) *fuNeg = *f;
	}
	else{
		*f = fRes;
		*ek = 0.0;
		*fuNeg = fRes;
	}
	return;
}


int
Bilinear::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return 0;
  
  if (strcmp(argv[0],"elstk") == 0) 
    return param.addObject(1, this);
  
  if (strcmp(argv[0],"fyieldPos") == 0) 
    return param.addObject(2, this);
  
  if (strcmp(argv[0],"fyieldNeg") == 0) 
    return param.addObject(3, this);
  
  if (strcmp(argv[0],"alfa") == 0) 
    return param.addObject(4, this);
  
  if (strcmp(argv[0],"alfaCap") == 0) 
    return param.addObject(5, this);
  
  if (strcmp(argv[0],"capDispPos") == 0) 
    return param.addObject(6, this);
  
  if (strcmp(argv[0],"capDispNeg") == 0) 
    return param.addObject(7, this);
  
  if (strcmp(argv[0],"Resfac") == 0) 
    return param.addObject(8, this);
  
  if (strcmp(argv[0],"flagCapenv") == 0) 
    return param.addObject(9, this);
  
  else
    opserr << "WARNING: Could not set parameter in BoucWenMaterial. " << endln;
  
  return 0;
}


int
Bilinear::updateParameter(int parameterID, Information &info)
{

	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->elstk = info.theDouble;
		break;
	case 2:
		this->fyieldPos = info.theDouble;
		break;
	case 3:
		this->fyieldNeg = info.theDouble;
		break;
	case 4:
		this->alfa = info.theDouble;
		break;
	case 5:
		this->alfaCap = info.theDouble;
		break;
	case 6:
		this->capDispPos = info.theDouble;
		break;
	case 7:
		this->capDispNeg = info.theDouble;
		break;
	case 8:
		this->Resfac = info.theDouble;
		break;
	case 9:
		this->flagCapenv = info.theInt;
		break;
	default:
		return -1;
	}

	return 0;
}



int
Bilinear::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}



/*
double
Bilinear::getStressSensitivity(int gradNumber, bool conditional)
{

	// Declare output variable
	double sensitivity = 0.0;


	// Issue warning if response is zero (then an error will occur)
	if (hsTrial[0] == 0.0)  {
		opserr << "ERROR: Bilinear::getStressSensitivity() is called " << endln
			<< " with zero hysteretic deformation d." << endln;
	}

	// First set values depending on what is random
	double Delstk;
	double DfyieldPos;
	double DfyieldNeg;
	double Dalfa;
	double DalfaCap;
	double DcapDispPos;
	double DcapDispNeg;
	double DResfac;


	if (parameterID == 0) { }
	else if (parameterID == 1) {Delstk=1.0;}
	else if (parameterID == 2) {DfyieldPos=1.0;}
	else if (parameterID == 3) {DfyieldNeg=1.0;}
	else if (parameterID == 4) {Dalfa=1.0;}
	else if (parameterID == 5) {DalfaCap=1.0;}
	else if (parameterID == 6) {DcapDispPos=1.0;}
	else if (parameterID == 7) {DcapDispNeg=1.0;}
	else if (parameterID == 8) {DResfac=1.0;}



	// Pick up sensitivity history variables for this gradient number
	double DCz = 0.0;
	double DCe = 0.0;
	double DCstrain = 0.0;
	if (SHVs != 0) {
		DCz		 = (*SHVs)(0,(gradNumber-1));
		DCe		 = (*SHVs)(1,(gradNumber-1));
		DCstrain = (*SHVs)(2,(gradNumber-1));
	}

	
	// Compute sensitivity of z_{i+1} 
	// (use same equations as for the unconditional 
	// sensitivities, just set DTstrain=0.0)
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11;
	double DTstrain = 0.0; 
	double dStrain = Tstrain-Cstrain;
	double TA, Tnu, Teta, DTz, Psi, Phi, DPsi;

	c1  = DCe 
		- Dalpha*ko*dStrain*Tz
		+ (1-alpha)*Dko*dStrain*Tz
		+ (1-alpha)*ko*(DTstrain-DCstrain)*Tz;
	c2  = (1-alpha)*ko*dStrain;
	c3  = DAo - DdeltaA*Te - deltaA*c1;
	c4  = -deltaA*c2;
	c5  = DdeltaNu*Te + deltaNu*c1;
	c6  = deltaNu*c2;
	c7  = DdeltaEta*Te + deltaEta*c1;
	c8  = deltaEta*c2;
	TA = Ao - deltaA*Te;
	Tnu = 1.0 + deltaNu*Te;
	Teta = 1.0 + deltaEta*Te;
	Psi = gamma + beta*signum(dStrain*Tz);
	DPsi= Dgamma + Dbeta*signum(dStrain*Tz);
	Phi = TA - pow(fabs(Tz),n)*Psi*Tnu;
	c9  = dStrain/Teta;
	c10 = DCz + c9*c3 - c9*pow(fabs(Tz),n)*Dn*log(fabs(Tz))*Psi*Tnu
		- c9*pow(fabs(Tz),n)*DPsi*Tnu - c9*pow(fabs(Tz),n)*Psi*c5
		- Phi/(Teta*Teta)*c7*dStrain + Phi/Teta*(DTstrain-DCstrain);
	c11 = 1.0 - c9*c4 + c9*pow(fabs(Tz),n)*Psi*c6
		+ c9*pow(fabs(Tz),n)*n/fabs(Tz)*signum(Tz)*Psi*Tnu
		+ Phi/(Teta*Teta)*c8*dStrain;

	DTz = c10/c11;

	sensitivity = Dalpha*ko*Tstrain
		        + alpha*Dko*Tstrain
				- Dalpha*ko*Tz
				+ (1-alpha)*Dko*Tz
				+ (1-alpha)*ko*DTz;

	return sensitivity;
}



double
BoucWenMaterial::getTangentSensitivity(int gradNumber)
{
	return 0.0;
}

double
BoucWenMaterial::getDampTangentSensitivity(int gradNumber)
{
	return 0.0;
}

double
BoucWenMaterial::getStrainSensitivity(int gradNumber)
{
	return 0.0;
}

double
BoucWenMaterial::getRhoSensitivity(int gradNumber)
{
	return 0.0;
}


int
BoucWenMaterial::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(3,numGrads);
	}

	// First set values depending on what is random
	double Dalpha = 0.0;
	double Dko = 0.0;
	double Dn = 0.0;
	double Dgamma = 0.0;
	double Dbeta = 0.0;
	double DAo = 0.0;
	double DdeltaA = 0.0;
	double DdeltaNu = 0.0;
	double DdeltaEta = 0.0;

	if (parameterID == 1) {Dalpha=1.0;}
	else if (parameterID == 2) {Dko=1.0;}
	else if (parameterID == 3) {Dn=1.0;}
	else if (parameterID == 4) {Dgamma=1.0;}
	else if (parameterID == 5) {Dbeta=1.0;}
	else if (parameterID == 6) {DAo=1.0;}
	else if (parameterID == 7) {DdeltaA=1.0;}
	else if (parameterID == 8) {DdeltaNu=1.0;}
	else if (parameterID == 9) {DdeltaEta=1.0;}


	// Pick up sensitivity history variables for this gradient number
	double DCz = 0.0;
	double DCe = 0.0;
	double DCstrain = 0.0;
	if (SHVs != 0) {
		DCz		 = (*SHVs)(0,(gradNumber-1));
		DCe		 = (*SHVs)(1,(gradNumber-1));
		DCstrain = (*SHVs)(2,(gradNumber-1));
	}

	
	// Compute sensitivity of z_{i+1} 
	// (use same equations as for the unconditional 
	// sensitivities, just set DTstrain=0.0)
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11;
	double DTstrain = TstrainSensitivity; 
	double dStrain = Tstrain-Cstrain;
	double TA, Tnu, Teta, DTz, Psi, Phi, DPsi, DTe;

	c1  = DCe 
		- Dalpha*ko*dStrain*Tz
		+ (1-alpha)*Dko*dStrain*Tz
		+ (1-alpha)*ko*(DTstrain-DCstrain)*Tz;
	c2  = (1-alpha)*ko*dStrain;
	c3  = DAo - DdeltaA*Te - deltaA*c1;
	c4  = -deltaA*c2;
	c5  = DdeltaNu*Te + deltaNu*c1;
	c6  = deltaNu*c2;
	c7  = DdeltaEta*Te + deltaEta*c1;
	c8  = deltaEta*c2;
	TA = Ao - deltaA*Te;
	Tnu = 1.0 + deltaNu*Te;
	Teta = 1.0 + deltaEta*Te;
	Psi = gamma + beta*signum(dStrain*Tz);
	DPsi= Dgamma + Dbeta*signum(dStrain*Tz);
	Phi = TA - pow(fabs(Tz),n)*Psi*Tnu;
	c9  = dStrain/Teta;
	c10 = DCz + c9*c3 - c9*pow(fabs(Tz),n)*Dn*log(fabs(Tz))*Psi*Tnu
		- c9*pow(fabs(Tz),n)*DPsi*Tnu - c9*pow(fabs(Tz),n)*Psi*c5
		- Phi/(Teta*Teta)*c7*dStrain + Phi/Teta*(DTstrain-DCstrain);
	c11 = 1.0 - c9*c4 + c9*pow(fabs(Tz),n)*Psi*c6
		+ c9*pow(fabs(Tz),n)*n/fabs(Tz)*signum(Tz)*Psi*Tnu
		+ Phi/(Teta*Teta)*c8*dStrain;

	DTz = c10/c11;

	DTe = DCe - Dalpha*ko*dStrain*Tz
		+ (1-alpha)*Dko*dStrain*Tz
		+ (1-alpha)*ko*(DTstrain-DCstrain)*Tz
		+ (1-alpha)*ko*dStrain*DTz;


	// Save sensitivity history variables
	(*SHVs)(0,(gradNumber-1)) = DTz;
	(*SHVs)(1,(gradNumber-1)) = DTe;
	(*SHVs)(2,(gradNumber-1)) = DTstrain;

	return 0;
}

*/

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
// $Date: 2000-12-13 05:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticPPMaterial.cpp,v $
                                                                        
                                                                        
// File: ~/material/ElasticPPMaterial.C
//
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) ElasticPPMaterial.C, revA"


#include <ElasticPPMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>


ElasticPPMaterial::ElasticPPMaterial(int tag, double e, double eyp)
:UniaxialMaterial(tag,MAT_TAG_ElasticPPMaterial),
 ezero(0.0), E(e), ep(0.0)
{
	fyp = E*eyp;
	fyn = -fyp;
}

ElasticPPMaterial::ElasticPPMaterial(int tag, double e, double eyp,
	double eyn, double ez )
:UniaxialMaterial(tag,MAT_TAG_ElasticPPMaterial),
 ezero(ez), E(e), ep(0.0), eptrial(0.)
{
    if (eyp < 0) {
	cerr << "ElasticPPMaterial::ElasticPPMaterial() - eyp < 0, setting > 0\n";
	eyp *= -1.;
    }
    if (eyn > 0) {
	cerr << "ElasticPPMaterial::ElasticPPMaterial() - eyn > 0, setting < 0\n";
	eyn *= -1.;
    }    
    
    fyp = E*eyp;
    fyn = E*eyn;
}

ElasticPPMaterial::ElasticPPMaterial()
:UniaxialMaterial(0,MAT_TAG_ElasticPPMaterial),
 fyp(0.0), fyn(0.0), ezero(0.0), E(0.0), ep(0.0), eptrial(0.)
{

}

ElasticPPMaterial::~ElasticPPMaterial()
{
  // does nothing
}

int 
ElasticPPMaterial::setTrialStrain(double strain, double strainRate)
{
    trialStrain = strain;
    return 0;
}

double 
ElasticPPMaterial::getStrain(void)
{
    return trialStrain;
}

double 
ElasticPPMaterial::getStress(void)
{
    double sigtrial;	// trial stress
    double f;		// yield function
    double sig;		// stress

    // compute trial stress
    sigtrial = E * ( trialStrain - ezero - ep );

    // evaluate yield function
    if ( sigtrial >= 0.0 )
	f =  sigtrial - fyp;
    else
	f = -sigtrial + fyn;

    if ( f <= 0.0 )
	// elastic
        sig = sigtrial;

    else
    {
	// plastic
	if ( sigtrial > 0.0 )
	{
	   sig = fyp;
	   eptrial = ep + f / E;
	}
	else
	{
	   sig = fyn;
	   eptrial = ep - f / E;
	}
    }
    
    return sig;
}


double 
ElasticPPMaterial::getTangent(void)
{
    double sigtrial;	// trial stress
    double f;		// yield function

    // compute trial stress
    sigtrial = E * ( trialStrain - ezero - ep );

    // evaluate yield function
    if ( sigtrial >= 0 )
	f =  sigtrial - fyp;
    else
	f = -sigtrial + fyn;

    return ( f <= 0.0 ) ? E : 0.0;

}

double ElasticPPMaterial::getSecant ()
{
    if (trialStrain != 0.0)
	return this->getStress()/trialStrain;
    else
	return E; 
}

int 
ElasticPPMaterial::commitState(void)
{
    ep = eptrial;
    return 0;
}	


int 
ElasticPPMaterial::revertToLastCommit(void)
{
    return 0;
}


int 
ElasticPPMaterial::revertToStart(void)
{
    ep = 0.0;
    return 0;
}


UniaxialMaterial *
ElasticPPMaterial::getCopy(void)
{
    ElasticPPMaterial *theCopy =
	 new ElasticPPMaterial(this->getTag(),E,fyp/E,fyn/E,ezero);
    theCopy->ep = 0.0;
    return theCopy;
}


int 
ElasticPPMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(6);
  data(0) = this->getTag();
  data(1) = ep;
  data(2) = E;
  data(3) = ezero;
  data(4) = fyp;
  data(5) = fyn;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    cerr << "ElasticPPMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
ElasticPPMaterial::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(6);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    cerr << "ElasticPPMaterial::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    ep    = data(1);
    E     = data(2);
    ezero = data(3);
    fyp   = data(4);
    fyn   = data(5);    
  }

  return res;
}

void 
ElasticPPMaterial::Print(ostream &s, int flag)
{
    s << "ElasticPP tag: " << this->getTag() << endl;
    s << "  E: " << E << endl;
    s << "  ep: " << ep << endl;
}



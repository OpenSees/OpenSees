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

// $Revision: 1.7 $
// $Date: 2009/03/23 23:17:04 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/ElasticPPcpp.cpp,v $

// Written: fmk
//
// Description: This file contains the class implementation for
// ElasticMaterial.
//
// What: "@(#) ElasticPPcpp.C, revA"

#include <elementAPI.h>
#include "ConcreteZBH_fitted.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>


static int numConcreteZBH_fitted = 0;

void *
OPS_ConcreteZBH_fitted()
{
  // print out some KUDO's
  if (numConcreteZBH_fitted == 0) {
    opserr << "ConcreteZBH_fitted uniaxial material - Use at your Own Peril\n";
    numConcreteZBH_fitted =1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //
  int    iData[1];
  double dData[20];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConcreteZBH_fitted tag" << endln;
    return 0;
  }

  numData = 20;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid ...\n";
    return 0;
  }

  //
  // create a new material
  //

  theMaterial = new ConcreteZBH_fitted(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7],
                                dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ConcreteZBH_fitted\n";
    return 0;
  }

  // return the material
  return theMaterial;
}

ConcreteZBH_fitted::ConcreteZBH_fitted(int tag, double _fc0, double _ec0, double _Ec, double _fccs, double _eccs, double _rs,
		       double _e1, double _e2, double _e3, double _e4, double _e5, double _e6, double _e7, double _e8, double _e9,
	           double _eps_cy, double _eps_ccuf, double _sig_ccuf, double _eps_ccus, double _sig_ccus)
:UniaxialMaterial(tag, 0), 
  fc0(_fc0), ec0(_ec0), Ec(_Ec), fccs(_fccs), eccs(_eccs), rs(_rs),
  e1(_e1), e2(_e2), e3(_e3), e4(_e4), e5(_e5), e6(_e6), e7(_e7), e8(_e8), e9(_e9),
	eps_cy(_eps_cy), eps_ccuf(_eps_ccuf), sig_ccuf(_sig_ccuf), eps_ccus(_eps_ccus), sig_ccus(_sig_ccus)
  
{
  sigp   = 0.0;
  Ep     = Ec;
  elp    = 0.0;
  epsp   = 0.0;
  eminp  = 0.0;
  eunl1p = 0.0;
  eunl2p = 0.0;
  eunl3p = 0.0;
  Eunlp  = fc0/ec0;
  Eunl2p = fc0/ec0;
  Et3p   = fc0/ec0;
  sunlp  = 0.0;
  elunlp = 0.0;
  muunlp = 0.0;
  flaggp = 4;

  sig    = 0.0;
  Et     = Ec;
  eps    = 0.0;
  el     = 0.0;
  emin   = 0.0;
  eunl1  = 0.0;
  eunl2  = 0.0;
  eunl3  = 0.0;
  Eunl   = fc0/ec0;
  Eunl2  = fc0/ec0;
  Et3    = fc0/ec0;
  sunl   = 0.0;
  elunl  = 0.0;
  muunl  = 0.0;
  flagg  = 4; 

  beta = (Ec / fabs(fc0) - 1 / fabs(ec0));
  r0 = Ec / (Ec - fc0/ ec0);

}

ConcreteZBH_fitted::ConcreteZBH_fitted()
:UniaxialMaterial(0, 0),
fc0(0.0), ec0(0.0), Ec(0.0), fccs(0.0), eccs(0.0), rs(0.0),
e1(0.0), e2(0.0), e3(0.0), e4(0.0), e5(0.0), e6(0.0), e7(0.0), e8(0.0), e9(0.0),
eps_cy(0.0), eps_ccuf(0.0), sig_ccuf(0.0), eps_ccus(0.0), sig_ccus(0.0)
{
sigp   = 0.0;
  Ep     = Ec;
  elp    = 0.0;
  epsp   = 0.0;
  eminp  = 0.0;
  eunl1p = 0.0;
  eunl2p = 0.0;
  eunl3p = 0.0;
  Eunlp  = 0.0;
  Eunl2p = 0.0;
  Et3p   = 0.0;
  sunlp  = 0.0;
  elunlp = 0.0;
  muunlp = 0.0;
  flaggp = 4;

  sig    = 0.0;
  Et     = Ec;
  eps    = 0.0;
  el     = 0.0;
  emin   = 0.0;
  eunl1  = 0.0;
  eunl2  = 0.0;
  eunl3  = 0.0;
  Eunl   = 0.0;
  Eunl2  = 0.0;
  Et3    = 0.0;
  sunl   = 0.0;
  elunl  = 0.0;
  muunl  = 0.0;
  flagg  = 4;


}

ConcreteZBH_fitted::~ConcreteZBH_fitted()
{
  // does nothing
}

int
ConcreteZBH_fitted::setTrialStrain(double strain, double strainRate)
{
	// retrieve concrete hitory variables

	//emin = eminp;
	flagg = flaggp;

  //sig   = sigp;
  //Et    = Ep;
  //el    = elp;
  //eps   = epsp;
  emin  = eminp;
  eunl1 = eunl1p;
  eunl2 = eunl2p;
  eunl3 = eunl3p;
  Eunl  = Eunlp;
  Eunl2 = Eunl2p;
  Et3   = Et3p;
  sunl  = sunlp;
  elunl = elunlp;
  muunl = muunlp;

	// calculate current strain

	eps = strain;
	double deps = eps - epsp;
	emin = fmin(eps, eminp);
	//emin = eminp;

    if (fabs(deps) < 10*DBL_EPSILON){
            sig = sigp;
            Et  = Ep;
			el  = elp;
    }
    else if (eps <= 0.) {                    // compression
            if (flagg == 4) {                  // monotonic envelope
                if (deps < 0.)  {              // negative strain increment: loading
                    this->envelope(eps, deps, sig, Et, el);
                    //emin = eps;
                }
                else {                        // deps > 0 - unloading
                    elunl = elp;
					if (elunl < 0) {
						elunl = 0;
					}
					sunl = sigp;
					Eunl = Ec;
					Eunl2 = Ec / (1 + 40 * elunl);
                    eunl1 = emin-sunl/Eunl;
					//eunl1 = fmin(eunl1, 0.0);
                    eunl2 = emin-sunl/Eunl2;
					if (eunl2 > 0.0) {
						Eunl2 = sunl / emin;
						Eunl = Eunl2 * (1 + 40 * elunl);
						eunl2 = 0.0;
						eunl1 = emin - sunl / Eunl;
					}
					//eunl2 = fmin(eunl2, 0.0);
                    muunl = 0.4*elunl*Eunl2/sunl;
                    if (eps > eunl2) {        // crack has opened
                        flagg = 2;
                        sig   = 0.;
                        //Et    = 0.;
						Et    = 1e-10;
                        el    = elunl+(Ec/Eunl2)*muunl*(eunl1-emin);
                    }
                    else if (eps > eunl1) {
                        flagg = 1;
                        sig   = 0.;
                        //Et    = 0.;
						Et    = 1e-10;
                        el    = elunl+(Ec/Eunl2)*muunl*(eunl1-emin);
                    }
                    else {
                        sig   = Eunl*(eps-eunl1);
                        Et    = Eunl;
                        flagg = 0;
                        el    = elunl+(Eunl/Eunl2)*muunl*(eps-emin);
                    }
                }
            }
            else if (flagg == 0) {
                if (eps <= emin) {             // the strain increment brings back on the envelope
                    flagg = 4;
                    this->envelope(eps, deps, sig, Et, el);
                    //emin = eps;
                }
                else if (eps > eunl2) {        // crack has opened
                    flagg = 2;
                    sig   = 0.;
                    //Et    = 0.;
					Et    = 1e-10;
                    el    = elunl+(Ec/Eunl2)*muunl*(eunl1-emin);
                }
                else if (eps > eunl1) {
                    flagg = 1;
                    sig   = 0.;
                    //Et    = 0.;
					Et    = 1e-10;
                    el    = elunl+(Ec/Eunl2)*muunl*(eunl1-emin);
                }
                else {                         // eunl1 > eps > emin
                    sig   = Eunl*(eps-eunl1);
                    Et    = Eunl;
                    el    = elunl+(Eunl/Eunl2)*muunl*(eps-emin);
					flagg = 0;
                }
            }
            else if (flagg == 1) {
                if (eps <= emin) {             // the strain increment brings back on the envelope
                    flagg = 4;
                    this->envelope(eps, deps, sig, Et, el);
                    //emin = eps;
                }
                else if (eps > eunl2) {        // crack has opened
                    flagg = 2;
                    sig   = 0.;
                    //Et    = 0.;
					Et    = 1e-10;
					el    = elunl+(Ec/Eunl2)*muunl*(eunl1-emin);
                }
                else if (deps < 0.) {            // loading
                    flagg = 3;
                    eunl3 = epsp;
                    Et3   = sunl/(emin-eunl3);
                    Et    = Et3;
                    sig   = Et*(eps-eunl3);
                    el    = elunl+(Et3/Eunl2)*muunl*(eps-emin);
                }
                else {                         // keeps flag = 1
                    //Et    = 0.;
                    flagg = 1;
					Et    = 1e-10;
                    sig  = 0.;
					el   = elunl+(Ec/Eunl2)*muunl*(eunl1-emin);
                }
            }
            else if (flagg == 2) {
                if (eps <= emin) {             // the strain increment brings back on the envelope
                    flagg = 4;
                    this->envelope(eps, deps, sig, Et, el);
                    //emin = eps;
                }
                else if (eps > eunl2) {        // crack has opened
                    flagg = 2;
                    sig   = 0.;
                    //Et    = 0.;
					Et    = 1e-10;
					el    = elunl+(Ec/Eunl2)*muunl*(eunl1-emin);
                }
                else {                         // eunl2 > eps > emin
                    flagg = 2;
                    sig   = Eunl2*(eps-eunl2);
                    Et    = Eunl2;
                    el    = elunl+muunl*(eps-emin);
                }
            }
            else if (flagg == 3) {
                if (eps <= emin) {             // the strain increment brings back on the envelope
                    flagg = 4;
                    this->envelope(eps, deps, sig, Et, el);
                   // emin = eps;
                }
                else if (eps > eunl2) {         // crack has opened
                    flagg = 2;
                    sig   = 0.;
                    //Et    = 0.;
					Et    = 1e-10;
					el    = elunl+(Ec/Eunl2)*muunl*(eunl1-emin);
                }
                else if (eps > eunl3) {
                    flagg = 1;
                    sig   = 0.;
                    //Et    = 0.;
					Et    = 1e-10;
					el    = elunl+(Ec/Eunl2)*muunl*(eunl1-emin);
                }
                else {                          // eunl3 > eps > emin
                    flagg = 3;
                    Et    = Et3;
                    sig   = (eps - eunl3)*Et;
                    el    = elunl+(Et3/Eunl2)*muunl*(eps-emin);
                }
            }
    }
    else {                                  //eps > 0
            flagg = 2;
            sig   = 0.;
            //Et    = 0.;
			Et    = 1e-10;
			el    = elunl+(Ec/Eunl2)*muunl*(eunl1-emin);
    }
 return 0;

}

double
ConcreteZBH_fitted::getStrain(void)
{
  return eps;
}

double
ConcreteZBH_fitted::getStress(void)
{
  return sig;
}


double
ConcreteZBH_fitted::getTangent(void)
{
  return Et;
}

int
ConcreteZBH_fitted::commitState(void)
{
  sigp   = sig;
  Ep     = Et;
  elp    = el;
  epsp   = eps;
  eminp  = emin;
  eunl1p = eunl1;
  eunl2p = eunl2;
  eunl3p = eunl3;
  Eunlp  = Eunl;
  Eunl2p = Eunl2;
  Et3p   = Et3;
  sunlp  = sunl;
  elunlp = elunl;
  muunlp = muunl;
  flaggp = flagg;

  return 0;
}


int
ConcreteZBH_fitted::revertToLastCommit(void)
{
  sig   = sigp;
  Et    = Ep;
  el    = elp;
  eps   = epsp;
  emin  = eminp;
  eunl1  = eunl1p;
  eunl2 = eunl2p;
  eunl3 = eunl3p;
  Eunl  = Eunlp;
  Eunl2 = Eunl2p;
  Et3   = Et3p;
  sunl  = sunlp;
  elunl = elunlp;
  muunl = muunlp;
  flagg = flaggp;

  return 0;
}


int
ConcreteZBH_fitted::revertToStart(void)
{
  sigp   = 0.0;
  Ep     = Ec;
  elp    = 0.0;
  epsp   = 0.0;
  eminp  = 0.0;
  eunl1p = 0.0;
  eunl2p = 0.0;
  eunl3p = 0.0;
  Eunlp  = fc0/ec0;
  Eunl2p = fc0/ec0;
  Et3p   = fc0/ec0;
  sunlp  = 0.0;
  elunlp = 0.0;
  muunlp = 0.0;
  flaggp = 4;

  sig    = 0.0;
  Et     = Ec;
  eps    = 0.0;
  el     = 0.0;
  emin   = 0.0;
  eunl1  = 0.0;
  eunl2  = 0.0;
  eunl3  = 0.0;
  Eunl   = fc0/ec0;
  Eunl2  = fc0/ec0;
  Et3    = fc0/ec0;
  sunl   = 0.0;
  elunl  = 0.0;
  muunl  = 0.0;
  flagg  = 4;


  return 0;
}


UniaxialMaterial *
ConcreteZBH_fitted::getCopy(void)
{
  ConcreteZBH_fitted *theCopy =
    new ConcreteZBH_fitted(this->getTag(), fc0, ec0, Ec, fccs, eccs, rs,
                    e1, e2, e3, e4, e5, e6, e7, e8, e9,
		           eps_cy, eps_ccuf, sig_ccuf, eps_ccus, sig_ccus);

  return theCopy;
}


int
ConcreteZBH_fitted::sendSelf(int cTag, Channel &theChannel)
{
	int	res = 0;
	static	Vector	data(36);
	data(0) = fc0;
	data(1) = ec0;
	data(2) = Ec;
	data(3) = fccs;
	data(4) = eccs;
	data(5) = rs;
	data(6) = e1;
	data(7) = e2;
	data(8) = e3;
	data(9) = e4;
	data(10) = e5;
	data(11) = e6;
	data(12) = e7;
	data(13) = e8;
	data(14) = e9;
	data(15) = eps_cy;
	data(16) = eps_ccuf;
	data(17) = sig_ccuf;
	data(18) = eps_ccus;
	data(19) = sig_ccus;
	data(20) = sigp;
	data(21) = Ep;
	data(22) = elp;
	data(23) = epsp;
	data(24) = eminp;
	data(25) = eunl1p;
	data(26) = eunl2p;
	data(27) = eunl3p;
	data(28) = Eunlp;
	data(29) = Eunl2p;
	data(30) = Et3p;
	data(31) = sunlp;
	data(32) = elunlp;
	data(33) = muunlp;
	data(34) = flaggp;
	data(35) = this->getTag();



  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "ConcreteZBH_fitted::sendSelf() - failed to send data\n";

  return res;
}

int
ConcreteZBH_fitted::recvSelf(int cTag, Channel &theChannel,
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(36);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "ConcreteZBH_fitted::recvSelf() - failed to recv data\n";
  else {
	  fc0 = data(0);
	  ec0 = data(1);
	  Ec = data(2);
	  fccs = data(3);
	  eccs = data(4);
	  rs = data(5);
	  e1 = data(6);
	  e2 = data(7);
	  e3 = data(8);
	  e4 = data(9);
	  e5 = data(10);
	  e6 = data(11);
	  e7 = data(12);
	  e8 = data(13);
	  e9 = data(14);
	  eps_cy = data(15);
	  eps_ccuf = data(16);
	  sig_ccuf = data(17);
	  eps_ccus = data(18);
	  sig_ccus = data(19);
	  sigp = data(20);
	  Ep = data(21);
	  elp = data(22);
	  epsp = data(23);
	  eminp = data(24);
	  eunl1p = data(25);
	  eunl2p = data(26);
	  eunl3p = data(27);
	  Eunlp = data(28);
	  Eunl2p = data(29);
	  Et3p = data(30);
	  sunlp = data(31);
	  elunlp = data(32);
	  muunlp = data(33);
	  flaggp = int(data(34));
	  this->setTag(int(data(35)));

  }

  eps = epsp;
  sig = sigp;
  Et  = Ep;

  return res;
}

void
ConcreteZBH_fitted::Print(OPS_Stream &s, int flag)
{
  s << "ConcreteZBH_fitted tag: " << this->getTag() << endln;
  /*s << "  E: " << E << endln;
  s << "  ep: " << ep << endln;
  s << "  stress: " << trialStress << " tangent: " << trialTangent << endln;*/
}

void
ConcreteZBH_fitted::envelope (double eps, double deps, double &sig, double &Et, double &el)
{ 
	if ((eps / ec0) < eps_ccuf) {
		sig = ((eps / ec0) < eps_cy) ? (e1 * pow((eps / ec0), 3) + e2 * pow((eps / ec0), 2) + (Ec * ec0 / fc0) * eps / ec0) / (e3 * pow((eps / ec0), 2) + e4 * (eps / ec0) + 1.) : ((e1 * pow(eps_cy, 3) + e2 * pow(eps_cy, 2) + (Ec * ec0 / fc0) * eps_cy) / (e3 * pow(eps_cy, 2) + e4 * eps_cy + 1) + (e5 * pow(((eps / ec0) - eps_cy), 3) + e6 * pow(((eps / ec0) - eps_cy), 2) + e7 * ((eps / ec0) - eps_cy)) / (e8 * pow(((eps / ec0) - eps_cy), 2) + e9 * ((eps / ec0) - eps_cy) + 1.));
		sig = sig * fc0;
		/*Et = (eps / ec0 < eps_cy) ? ((3 * e1 * pow(eps / ec0, 2) + 2 * e2 * eps / ec0 + (Ec * ec0 / fc0)) * (e3 * pow(eps / ec0, 2) + e4 * eps / ec0 + 1) - (e1 * pow(eps / ec0, 3) + e2 * pow(eps / ec0, 2) + (Ec * ec0 / fc0) * eps / ec0) * (2 * e3 * eps / ec0 + e3)) / pow((e3 * pow(eps / ec0, 2) + e4 * eps / ec0 + 1), 2)
			: ((3 * e5 * pow((eps / ec0 - eps_cy), 2) + 2 * e6 * (eps / ec0 - eps_cy) + e7) * (e8 * pow((eps / ec0 - eps_cy), 2) + e9 * (eps / ec0 - eps_cy) + 1) - (e5 * pow((eps / ec0 - eps_cy), 3) + e6 * pow((eps / ec0 - eps_cy), 2) + e7 * (eps / ec0 - eps_cy)) * (2 * e8 * (eps / ec0 - eps_cy) + e9)) / pow((e8 * pow((eps / ec0 - eps_cy), 2) + e9 * (eps / ec0 - eps_cy) + 1), 2);
		Et = Et * fc0 / ec0;*/
		Et = (sig - sigp) / deps;
		el = (Ec * eps - sig) / (2. * beta * sig);
	}
	else if ((eps / ec0) <= eps_ccus) {
		double sig1 = fccs * (eps / eccs) *rs / (rs - 1. + pow((eps / eccs), rs));
		double sig2 = sig_ccuf*fc0 - (eps - eps_ccuf*ec0) * Ec;
		if (sig1 < sig2) {
			sig = sig1;
			Et = (fccs * rs / eccs) / (rs - 1 + pow((eps / eccs), rs)) -
				(fccs * pow(rs, 2) * (pow((eps / eccs), rs)) / eccs) / pow((rs - 1 + pow((eps / eccs), rs)), 2);
			el = (Ec * eps - sig) / (2. * beta * sig);
		}
		else {
			sig = sig2;
			Et = -Ec;
			el = (Ec * eps - sig) / (2 * beta * sig);
		}
	}
	else {
		double sig_p2 = sig_ccuf * fc0 - (eps_ccus - eps_ccuf)*ec0 * Ec;
		double sig_p = fmin(sig_ccus*fc0, sig_p2);

		double sig2 = sig_p  - (eps - eps_ccus*ec0) * Ec;
		double sig1 = fc0 * r0 * (eps / ec0) / (r0 - 1. + pow((eps / ec0), r0));
		if (sig1 < sig2) {
			sig = sig1;
			Et = (fc0 * r0 / ec0) / (r0 - 1. + pow((eps / ec0), r0)) -
				(fc0 * pow(r0, 2) * (pow((eps / ec0), r0)) / ec0) / pow((r0 - 1. + pow((eps / ec0), r0)), 2);
			el = (Ec * eps - sig) / (2 * beta * sig);
		}
		else {
			sig = sig2;
			Et = -Ec;
			el = (Ec * eps - sig) / (2. * beta * sig);
		}
	}
	return;
}

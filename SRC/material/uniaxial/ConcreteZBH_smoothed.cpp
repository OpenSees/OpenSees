#include <elementAPI.h>
#include "ConcreteZBH_smoothed.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

static int numConcreteZBH_smoothed = 0;

void *
OPS_ConcreteZBH_smoothed()
{
  // print out some KUDO's
  if (numConcreteZBH_smoothed == 0) {
    opserr << "ConcreteZBH_smoothed uniaxial material - Use at your Own Peril\n";
    numConcreteZBH_smoothed =1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //
  int    iData[1];
  double dData[18];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConcreteZBH_smoothed tag" << endln;
    return 0;
  }

  numData = 18;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid ...\n";
    return 0;
  }

  //
  // create a new material
  //

  theMaterial = new ConcreteZBH_smoothed(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7],
                                dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], dData[15],
                                dData[16], dData[17]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ConcreteZBH_original\n";
    return 0;
  }

  // return the material
  return theMaterial;
}

ConcreteZBH_smoothed::ConcreteZBH_smoothed(int tag, double _fc0, double _ec0, double _Ec,
		       double _Es, double _fy, double _eults, double _s, double _As_t,
		       double _Ef, double _eultf, double _tf, double _D, double _Ds,
		       double _As_l, double _kg_f, double _ks_s, double _ks_f, double _type_reinf)
:UniaxialMaterial(tag, 0),
 fc0(_fc0), ec0(_ec0), Ec(_Ec), Es(_Es), fy(_fy), eults(_eults),
  s(_s), As_t(_As_t), Ef(_Ef), eultf(_eultf), tf(_tf),
  D(_D), Ds(_Ds), As_l(_As_l), kg_f(_kg_f), ks_s(_ks_s),
  ks_f(_ks_f), type_reinf(_type_reinf)
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
  flp    = 0.0;
  flunlp = 0.0;
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
  fl     = 0.0;
  flunl  = 0.0;
  elunl  = 0.0;
  muunl  = 0.0;
  flagg  = 4;

  roj_f  = 4 * tf / D;
  roj_s  = 4 * As_t / (s*Ds);
  roj_sl = As_l / (0.25*3.1416*pow(Ds, 2));
  kg_s = (1 - 0.5*(s - 2 * pow((As_t / 3.1416), 0.5)) / Ds);
  kg_s = (kg_s > 0) ? (pow((1 - 0.5*(s - 2 * pow((As_t / 3.1416), 0.5)) / Ds), type_reinf) / (1 - roj_sl)) : 0;
  kg_s = fmin(kg_s, 1.0);
  //kg_s = pow((1 - 0.5*(s - 2 * pow((As_t / 3.1416), 0.5)) / Ds), 2) / (1 - roj_sl);
  
  beta = Ec/fabs(fc0)-1/fabs(ec0);

  fls = 0.5*ks_s*kg_s*roj_s*fy;
  fccs = (2.254 * pow((1 + 7.94 * fls / fabs(fc0)), 0.5) - 2 * fls / fabs(fc0) - 1.254) * fc0;
  eccs = ec0 * (1 + 5 * (fccs / fc0 - 1));
  eps_ccus = -0.004 - 1.4 * roj_s * fy * eults / fabs(fccs);
  double x = eps_ccus / eccs;
  double Esecs = fccs / eccs;
  rs = Ec / (Ec - Esecs);
  sig_ccus = fccs * x * rs / (rs - 1 + pow(x, rs));

  double flf = 0.5 * ks_f * kg_f * roj_f * Ef * eultf + 0.5 * ks_s * kg_s * roj_s * fy;
  double fccf = (2.254 * pow((1 + 7.94 * flf / fabs(fc0)), 0.5) - 2 * flf / fabs(fc0) - 1.254) * fc0;
  double eccf = ec0 * (1 + 5 * (fccf / fc0 - 1));
  double Esecf = fccf / eccf;
  double Esecu = Ec / (1 + 2 * beta * eultf);
  eps_ccuf = eccf * pow(((Esecf * (Ec - Esecu)) / (Esecu * (Ec - Esecf))), 1 - (Esecf / Ec));
  sig_ccuf = Esecu * eps_ccuf;

  r0 = Ec / (Ec - fc0 / ec0);

  if (fabs(eps_ccuf) > fabs(eps_ccus)) {
	  eps_ccus = eps_ccuf;
      sig_ccus = sig_ccuf;
  }

}

ConcreteZBH_smoothed::ConcreteZBH_smoothed()
:UniaxialMaterial(0, 0),
 fc0(0.0), ec0(0.0), Ec(0.0), Es(0.0), fy(0.0), eults(0.0),
  s(0.0), As_t(0.0), Ef(0.0), eultf(0.0), tf(0.0),
  D(0.0), Ds(0.0), As_l(0.0), kg_f(0.0), ks_s(0.0),
  ks_f(0.0), type_reinf(0.0)
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
  flp    = 0.0;
  flunlp = 0.0;
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
  fl     = 0.0;
  flunl  = 0.0;
  elunl  = 0.0;
  muunl  = 0.0;
  flagg  = 4;
}

ConcreteZBH_smoothed::~ConcreteZBH_smoothed()
{
  // does nothing
}

int
ConcreteZBH_smoothed::setTrialStrain(double strain, double strainRate)
{
	// retrieve concrete hitory variables

	//emin = eminp;
	flagg = flaggp;

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
  fl    = flp;
  flunl = flunlp;
  elunl = elunlp;
  muunl = muunlp;

	// calculate current strain

	eps = strain;
	double deps = eps - epsp;
	emin = fmin(eps, eminp);

    if (fabs(deps) < 10*DBL_EPSILON){
            sig = sigp;
            Et  = Ep;
			el  = elp;
    }
    else if (eps <= 0.) {                    // compression
            if (flagg == 4) {                  // monotonic envelope
                if (deps < 0.)  {              // negative strain increment: loading
                    this->envelope(eps, deps, fl, sig, Et, el);
                   // emin = eps;
                }
                else {                        // deps > 0 - unloading
                    elunl = elp;
                    Eunl2 = Ec/(1+2*20*elunl);
                    sunl  = sigp;
                    eunl1 = emin-sigp/Ec;
                    eunl2 = emin-sigp/Eunl2;
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
                        sig   = Ec*(eps-eunl1);
                        Et    = Ec;
                        flagg = 0;
                        el    = elunl+(Ec/Eunl2)*muunl*(eps-emin);
                    }
                }
            }
            else if (flagg == 0) {
                if (eps <= emin) {             // the strain increment brings back on the envelope
                    flagg = 4;
                    this->envelope(eps, deps, fl, sig, Et, el);
                   // emin = eps;
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
                    sig   = Ec*(eps-eunl1);
                    Et    = Ec;
                    el    = elunl+(Ec/Eunl2)*muunl*(eps-emin);
					flagg = 0;
                }
            }
            else if (flagg == 1) {
                if (eps <= emin) {             // the strain increment brings back on the envelope
                    flagg = 4;
                    this->envelope(eps, deps, fl, sig, Et, el);
                   // emin = eps;
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
                    this->envelope(eps, deps, fl, sig, Et, el);
                  //  emin = eps;
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
                    this->envelope(eps, deps, fl, sig, Et, el);
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
ConcreteZBH_smoothed::getStrain(void)
{
  return eps;
}

double
ConcreteZBH_smoothed::getStress(void)
{
  return sig;
}


double
ConcreteZBH_smoothed::getTangent(void)
{
  return Et;
}

int
ConcreteZBH_smoothed::commitState(void)
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
  flp    = fl;
  flunlp = flunl;
  elunlp = elunl;
  muunlp = muunl;
  flaggp = flagg;

  return 0;
}


int
ConcreteZBH_smoothed::revertToLastCommit(void)
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
  fl    = flp;
  flunl = flunlp;
  elunl = elunlp;
  muunl = muunlp;
  flagg = flaggp;

  return 0;
}


int
ConcreteZBH_smoothed::revertToStart(void)
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
  flp    = 0.0;
  flunlp = 0.0;
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
  fl     = 0.0;
  flunl  = 0.0;
  elunl  = 0.0;
  muunl  = 0.0;
  flagg  = 4;
  return 0;
}


UniaxialMaterial *
ConcreteZBH_smoothed::getCopy(void)
{
  ConcreteZBH_smoothed *theCopy =
    new ConcreteZBH_smoothed(this->getTag(), fc0, ec0, Ec, Es, fy, eults,
									s, As_t, Ef, eultf, tf, D, Ds, As_l, kg_f, ks_s, ks_f, type_reinf);

  return theCopy;
}


int
ConcreteZBH_smoothed::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(46);

  data(0) = fc0;
  data(1) = ec0;
  data(2) = Ec;
  data(3) = Es;
  data(4) = fy;
  data(5) = eults;
  data(6) = s;
  data(7) = As_t;
  data(8) = Ef;
  data(9) = eultf;
  data(10) = tf;
  data(11) = D;
  data(12) = Ds;
  data(13) = As_l;
  data(14) = kg_f;
  data(15) = ks_s;
  data(16) = ks_f;
  data(17) = type_reinf;
  data(18) = beta;
  data(19) = fls;
  data(20) = eps_ccuf;
  data(21) = eps_ccus;
  data(22) = sig_ccuf;
  data(23) = sig_ccus;
  data(24) = fccs;
  data(25) = eccs;
  data(26) = rs;
  data(27) = r0;
  data(28) = sigp;
  data(29) = Ep;
  data(30) = elp;
  data(31) = epsp;
  data(32) = eminp;
  data(33) = eunl1p;
  data(34) = eunl2p;
  data(35) = eunl3p;
  data(36) = Eunlp;
  data(37) = Eunl2p;
  data(38) = Et3p;
  data(39) = sunlp;
  data(40) = flp;
  data(41) = flunlp;
  data(42) = elunlp;
  data(43) = muunlp;
  data(44) = flaggp;

  data(45) = this->getTag();

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "ConcreteZBH_original::sendSelf() - failed to send data\n";

  return res;
}

int
ConcreteZBH_smoothed::recvSelf(int cTag, Channel &theChannel,
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(46);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "ConcreteZBH_smoothed::recvSelf() - failed to recv data\n";
  else {
	  fc0 = data(0);
	  ec0 = data(1);
	  Ec = data(2);
	  Es = data(3);
	  fy = data(4);
	  eults = data(5);
	  s = data(6);
	  As_t = data(7);
	  Ef = data(8);
	  eultf = data(9);
	  tf = data(10);
	  D = data(11);
	  Ds = data(12);
	  As_l = data(13);
	  kg_f = data(14);
	  ks_s = data(15);
	  ks_f = data(16);
	  type_reinf = data(17);
	  beta = data(18);
	  fls = data(19);
	  eps_ccuf = data(20);
	  eps_ccus = data(21);
	  sig_ccuf = data(22);
	  sig_ccus = data(23);
	  fccs = data(24);
	  eccs = data(25);
	  rs = data(26);
	  r0 = data(27);
	  sigp = data(28);
	  Ep = data(29);
	  elp = data(30);
	  epsp = data(31);
	  eminp = data(32);
	  eunl1p = data(33);
	  eunl2p = data(34);
	  eunl3p = data(35);
	  Eunlp = data(36);
	  Eunl2p = data(37);
	  Et3p = data(38);
	  sunlp = data(39);
	  flp = data(40);
	  flunlp = data(41);
	  elunlp = data(42);
	  muunlp = data(43);
	  flaggp = int(data(44));
this->setTag(int(data(45)));
  }

  eps = epsp;
  sig = sigp;
  Et  = Ep;
  el  = elp;
  fl  = flp;

  return res;
}

void
ConcreteZBH_smoothed::Print(OPS_Stream &s, int flag)
{
  s << "ConcreteZBH_smoothed tag: " << this->getTag() << endln;
  /*s << "  E: " << E << endln;
  s << "  ep: " << ep << endln;
  s << "  stress: " << trialStress << " tangent: " << trialTangent << endln;*/
}

void
ConcreteZBH_smoothed::Conf_Pressure (double eps, double flp, double &fc, double &fl, double &el)
{

double fcc   = (2.254*pow((1+7.94*flp/fabs(fc0)),0.5)-2*flp/fabs(fc0)-1.254)*fc0;
double ecc   = ec0*(1+5*(fcc/fc0-1));
double x     = eps/ecc;
double Esecc = fcc/ecc;
double r     = Ec/(Ec-Esecc);
fc    = fcc*x*r/(r-1+pow(x,r));
el    = (Ec*eps-fc)/(2*beta*fc);


fl = 0.5*ks_f*kg_f*roj_f*Ef*el + 0.5*ks_s*kg_s*roj_s*fmin(Es*fabs(el), fy);
  return;
}

void                   
ConcreteZBH_smoothed::envelope(double eps, double deps, double &fl, double &sig, double &Et, double &el)
{
	if (fls == 0.0) {
		if (eps >= eps_ccuf) {
			this->Conf_Pressure(eps, flp, sig, fl, el);
			int count = 0;
			while (fabs(fl - flp) > fmax(fl / 100000000, 0.00000000001))
			{
				count = count + 1;
				flp = fl;
				this->Conf_Pressure(eps, flp, sig, fl, el);
				if (count > 20) { break; }
			}
			Et = (sig - sigp) / deps;
		}
		else {
			double sig2 = sig_ccuf - (eps - eps_ccuf) * Ec;
			double sig1 = fc0 * r0 * (eps / ec0) / (r0 - 1 + pow((eps / ec0), r0));
			if (sig1 <= sig2) {
				sig = sig1;
				Et = (fc0 * r0 / ec0) / (r0 - 1 + pow((eps / ec0), r0)) -
					(fc0 * pow(r0, 2) * (pow((eps / ec0), r0)) / ec0) / pow((r0 - 1 + pow((eps / ec0), r0)), 2);
				el = (Ec * eps - sig) / (2 * beta * sig);
			}
			else {
				sig = sig2;
				Et = -Ec;
				el = (Ec * eps - sig) / (2 * beta * sig);
			}
		}
	}
	else {
		if (eps >= eps_ccuf) {
			this->Conf_Pressure(eps, flp, sig, fl, el);
			int count = 0;
			while (fabs(fl - flp) > fmax(fl / 100000000, 0.00000000001))
			{
				count = count + 1;
				flp = fl;
				this->Conf_Pressure(eps, flp, sig, fl, el);
				if (count > 20) { break; }
			}
			Et = (sig - sigp) / deps;
		}
		else if (eps >= eps_ccus)
		{
			double sig1 = fccs * (eps / eccs) *rs / (rs - 1 + pow((eps / eccs), rs));
			double sig2 = sig_ccuf - (eps - eps_ccuf) * Ec;
			if (sig1 <= sig2) {
				sig = sig1;
				Et = (fccs * rs / eccs) / (rs - 1 + pow((eps / eccs), rs)) -
					(fccs * pow(rs, 2) * (pow((eps / eccs), rs)) / eccs) / pow((rs - 1 + pow((eps / eccs), rs)), 2);
				el = (Ec * eps - sig) / (2 * beta * sig);
			}
			else {
				sig = sig2;
				Et = -Ec;
				el = (Ec * eps - sig) / (2 * beta * sig);
			}
		}
		else {
			double sig_p2 = sig_ccuf - (eps_ccus - eps_ccuf) * Ec;
			double sig_p = fmin(sig_ccus, sig_p2);

			double sig2 = sig_p - (eps - eps_ccus) * Ec;
			double sig1 = fc0 * r0 * (eps / ec0) / (r0 - 1 + pow((eps / ec0), r0));
			if (sig1 <= sig2) {
				sig = sig1;
				Et = (fc0 * r0 / ec0) / (r0 - 1 + pow((eps / ec0), r0)) -
					(fc0 * pow(r0, 2) * (pow((eps / ec0), r0)) / ec0) / pow((r0 - 1 + pow((eps / ec0), r0)), 2);
				el = (Ec * eps - sig) / (2 * beta * sig);
			}
			else {
				sig = sig2;
				Et = -Ec;
				el = (Ec * eps - sig) / (2 * beta * sig);
			}
		}
	}
	
	return;
}

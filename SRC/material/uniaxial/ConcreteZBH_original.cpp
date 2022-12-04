#include <elementAPI.h>
#include "ConcreteZBH_original.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

static int numConcreteZBH_original = 0;

void *
OPS_ConcreteZBH_original()
{
  // print out some KUDO's
  if (numConcreteZBH_original == 0) {
    opserr << "ConcreteZBH_original uniaxial material - Use at your Own Peril\n";
    numConcreteZBH_original =1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  // parse the input line for the material parameters
  int    iData[1];
  double dData[18];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConcreteZBH_original tag" << endln;
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

  theMaterial = new ConcreteZBH_original(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7],
                                dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], dData[15],
                                dData[16], dData[17]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ConcreteZBH_original\n";
    return 0;
  }

  // return the material
  return theMaterial;
}

ConcreteZBH_original::ConcreteZBH_original(int tag, double _fc0, double _ec0, double _Ec,
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

  roj_f  = 4 * tf / D;   //reinforcement ratio
  roj_s  = 4 * As_t / (s*Ds);
  roj_sl = As_l / (0.25*3.1416*pow(Ds, 2));
  kg_s = (1 - 0.5*(s - 2 * pow((As_t / 3.1416), 0.5)) / Ds);
  kg_s = (kg_s > 0) ? (pow((1 - 0.5*(s - 2 * pow((As_t / 3.1416), 0.5)) / Ds), type_reinf) / (1 - roj_sl)) : 0;
  kg_s = fmin(kg_s, 1.0);   //effectiveness coefficient by Mander et al. (1988)
  
  beta = Ec/fabs(fc0)-1/fabs(ec0);

  double fls = 0.5 * ks_s * kg_s * roj_s * fy;
  double fccs = (2.254 * pow((1 + 7.94 * fls / fabs(fc0)), 0.5) - 2 * fls / fabs(fc0) - 1.254) * fc0;
  eccus = -0.004 - 1.4 * roj_s * fy * eults / fabs(fccs);

}

ConcreteZBH_original::ConcreteZBH_original()
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

ConcreteZBH_original::~ConcreteZBH_original()
{
  // does nothing
}

int
ConcreteZBH_original::setTrialStrain(double strain, double strainRate)
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
    else if (eps <= 0.) {                      // compression
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
ConcreteZBH_original::getStrain(void)
{
  return eps;
}

double
ConcreteZBH_original::getStress(void)
{
  return sig;
}


double
ConcreteZBH_original::getTangent(void)
{
  return Et;
}

int
ConcreteZBH_original::commitState(void)
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
ConcreteZBH_original::revertToLastCommit(void)
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
ConcreteZBH_original::revertToStart(void)
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
ConcreteZBH_original::getCopy(void)
{
  ConcreteZBH_original *theCopy =
    new ConcreteZBH_original(this->getTag(), fc0, ec0, Ec, Es, fy, eults,
									s, As_t, Ef, eultf, tf, D, Ds, As_l, kg_f, ks_s, ks_f, type_reinf);

  return theCopy;
}


int
ConcreteZBH_original::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(37);

data(0) =fc0;
data(1) =ec0;
data(2) =Ec;
data(3) =Es;
data(4) =fy;
data(5) =eults;
data(6) =s;
data(7) =As_t;
data(8) =Ef;
data(9) =eultf;
data(10) =tf;
data(11) =D;
data(12) =Ds;
data(13) =As_l;
data(14) =kg_f;
data(15) =ks_s;
data(16) =ks_f;
data(17) = type_reinf;
data(18) =sigp;
data(19) =Ep;
data(20) =elp;
data(21) =epsp;
data(22) =eminp;
data(23) =eunl1p;
data(24) =eunl2p;
data(25) =eunl3p;
data(26) =Eunlp;
data(27) =Eunl2p;
data(28) =Et3p;
data(29) =sunlp;
data(30) =flp;
data(31) =flunlp;
data(32) =elunlp;
data(33) =muunlp;
data(34) =flaggp;
data(35) = eccus;
  data(36) = this->getTag();

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "ConcreteZBH_original::sendSelf() - failed to send data\n";

  return res;
}

int
ConcreteZBH_original::recvSelf(int cTag, Channel &theChannel,
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(37);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "ConcreteZBH_original::recvSelf() - failed to recv data\n";
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
	  sigp = data(18);
	  Ep = data(19);
	  elp = data(20);
	  epsp = data(21);
	  eminp = data(22);
	  eunl1p = data(23);
	  eunl2p = data(24);
	  eunl3p = data(25);
	  Eunlp = data(26);
	  Eunl2p = data(27);
	  Et3p = data(28);
	  sunlp = data(29);
	  flp = data(30);
	  flunlp = data(31);
	  elunlp = data(32);
	  muunlp = data(33);
	  flaggp=	  int(data(34));
	  eccus = data(35);
	  this->setTag(int(data(36)));
  }

  eps = epsp;
  sig = sigp;
  Et  = Ep;
  el  = elp;
  fl  = flp;

  return res;
}

void
ConcreteZBH_original::Print(OPS_Stream &s, int flag)
{
  s << "ConcreteZBH_original tag: " << this->getTag() << endln;

}

void
ConcreteZBH_original::Conf_Pressure (double eps, double flp, double &fc, double &fl, double &el)
{

double fcc   = (2.254*pow((1+7.94*flp/fabs(fc0)),0.5)-2*flp/fabs(fc0)-1.254)*fc0;
double ecc   = ec0*(1+5*(fcc/fc0-1));
double x     = eps/ecc;
double Esecc = fcc/ecc;
double r     = Ec/(Ec-Esecc);
fc    = fcc*x*r/(r-1+pow(x,r));
el    = (Ec*eps-fc)/(2*beta*fc);
double el0   = 0;
double el1   = 0;

if (el >= eultf) {
    el0 = 0;
}
else {
    el0 = fabs(el);
}
if (eps <= eccus) {
	el1 = 0;
}
else {
	el1 = fabs(el);
}

fl    = 0.5*ks_f*kg_f*roj_f*Ef*el0+0.5*ks_s*kg_s*roj_s*fmin(Es*fabs(el1),fy);

  return;
}

void                   
ConcreteZBH_original::envelope (double eps, double deps, double &fl, double &sig, double &Et, double &el)
{
	this->Conf_Pressure(eps, flp, sig, fl, el);
	int count = 0;
	while(fabs(fl-flp)>fmax(fl/10000,0.0000001))
	{
		count = count+1;
		flp = fl;
		this->Conf_Pressure(eps, flp, sig, fl, el);
		if (count>20){break;}
	}
	Et  = (sig-sigp)/deps;

	return;
}

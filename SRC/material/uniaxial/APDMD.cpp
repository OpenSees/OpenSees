#include <elementAPI.h>
#include "APDMD.h"
#include <iostream>
#include <Vector.h>
#include <float.h>
#include <Channel.h>
#include <math.h>

static int numAPDMD = 0;

void *
OPS_APDMD(void)
{
    if (numAPDMD == 0) {
        opserr << "APDMD unaxial material - Written by BUCEA 2024; \n";
        numAPDMD++;
    }

    opserr << "Due to known issues and unreliable results, this material has been" << endln;
    opserr << "temporarily removed from the compiled versions of OpenSees (Tcl and Py)" << endln;
    opserr << "The material source code remains available. Compile at your own risk." << endln;
    return 0;
  
    UniaxialMaterial* theMaterial = 0;

    int    iData[1];
    double dData[11];
    double mData[4];
    int numData;

    numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial APDMD tag" << endln;
        return 0;
    }

    numData = 11;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial APDMD\n";
        return 0;
    }

    theMaterial = new APDMD(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9],dData[10]);
 
    if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type APDMD\n";
    return 0;
  }
  return theMaterial;
}


APDMD::APDMD(int tag,double _Fy1, double _E1,double _ad2,double _Fy2, double _E2, double _b1, double _b2, double _rezaA1, double _rezaN1, double _rezaA2, double _rezaN2):
  UniaxialMaterial(tag, MAT_TAG_APDMD),
  Fy1(_Fy1), E1(_E1), ad2(_ad2),Fy2(_Fy2), E2(_E2), b1(_b1),b2(_b2),rezaAA1(_rezaA1), rezaNN1(_rezaN1), rezaAA2(_rezaA2), rezaNN2(_rezaN2)
{
  konP = 0;
  kon = 0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;

  epsmaxP1 = Fy1/E1;
  epsminP1 = -epsmaxP1;
  epsmaxP2 = Fy2 / E2+ ad2;
  epsminP2= -epsmaxP2;

  epsplP = 0.0;
  epssAP = 0.0;
  sigsAP = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;
  CrotMax= 0.0 ;
} 


APDMD::APDMD(void):
  UniaxialMaterial(0, MAT_TAG_APDMD)
{
  konP = 0;
}

APDMD::~APDMD(void)
{
  // Does nothing
}

UniaxialMaterial*
APDMD::getCopy(void)
{
  APDMD *theCopy = new APDMD(this->getTag(), Fy1, E1,ad2, Fy2, E2,b1,b2, rezaAA1, rezaNN1, rezaAA2, rezaNN2);
  return theCopy;
}

double
APDMD::getInitialTangent(void)
{
  return E1;
}

int
APDMD::setTrialStrain(double strain, double strainRate)
{
    double epsy1 = Fy1 / E1;
    eps = strain;
    double deps = eps - epsP;
   
    TrotMax = CrotMax;
    if (strain >= CrotMax) {
        TrotMax = strain;
    }
    if (TrotMax < ad2) {
        epsmax1 = epsmaxP1;
        epsmin1 = epsminP1;
        epspl = epsplP;
        epssA = epssAP;
        sigsA = sigsAP;
      
        sigr = sigsrP;
        kon = konP;
        if (kon == 0 || kon == 3) {
            if (fabs(deps) < 10.0 * DBL_EPSILON) {
                sig = 0;
                kon = 3;
                return 0;
            }
            else {
                epsmax1 = epsy1;
                epsmin1 = -epsy1;
                if (deps < 0.0) {
                    kon = 2;
                    epssA = epsmin1;
                    sigsA = Fy1;
                    epspl = epsmin1;
                }
                else {
                    kon = 1;
                    epssA = epsmax1;
                    sigsA = Fy1;
                    epspl = epsmax1;
                }
            }

        }
        if (kon == 2 && deps > 0.0) {
            kon = 20;
        }
        else if (kon == 1 && deps < 0.0) {
            kon = 10;
            epsr = epsP;
            sigr = sigP;
        }
        else if (kon == 10 && deps > 0.0) {
            kon = 20;
        }
        else if (kon == 20 && eps > CrotMax) {
            kon = 1;
        }   
        if (kon == 2 || kon == 1) {
            double stress1 = findstress1(eps, E1, Fy1, 0, rezaAA1, rezaNN1);
            sig = stress1;
        }
        else if (kon == 10) {
            double stress1n = findstressn(eps, sigr, E1, Fy1, rezaAA1, rezaNN1);
            sig = stress1n;
        }
        else if (kon == 20) {
            double stress1p = findstressp(eps, sigr, E1, Fy1, rezaAA1, rezaNN1);
            sig = stress1p;
        }
       

    }

    else if (TrotMax >= ad2) {
        double epsy2 = Fy2/ E2 + ad2;
        epsmax2 = epsmaxP2;
        epsmin2 = epsminP2;
        epssC = epssCP;
        sigsC = sigsCP;
        epsr = epssrP;
        sigr = sigsrP;
        kon = konP;

        if (kon == 0 || kon == 3) {
            if (fabs(deps) < 10.0 * DBL_EPSILON) {
                sig = 0;
                kon = 3;
                return 0;
            }
            else {
                epsmax2 = epsy2;
                epsmin2 = -epsy2;
                if (deps < 0.0) {
                    kon = 2;
                    epssC = epsmin2;
                    sigsC = Fy2 + (Fy1 + E1 * (epssC - ad2));
                    epspl = epsmin2;
                }
                else {
                    kon = 1;
                    epssC = epsmax2;
                    sigsC = Fy2 +(Fy1+E1*(epssC-ad2));
                    epspl = epsmax2;
                }
            }
         
        }

        if (epsr == 0) {
            stap = ad2;
            stan = -ad2;
        }

        double Fyb = findstress1(ad2, E1, Fy1, 0, rezaAA1, rezaNN1);
        double momz = sigr - findstress1(epsr, E1, Fy1, 0, rezaAA1, rezaNN1);
        double rotch = fabs(momz) / E2 * (1 + rezaAA2 * pow(fabs(momz) /(Fy2 * 2), rezaNN2)) - TrotMax;
        double rotc = -fabs(momz) / E2 * (1 + rezaAA2 * pow(fabs(momz) / (Fy2 * 2), rezaNN2)) + TrotMax;
        stan = -ad2 + (rotc - ad2);
        stap = rotc - (stan - rotch);

        if (stap > ad2) {
            stap = ad2;
        }
        if (stan < -ad2) {
            stan = -ad2;
        }

        if (kon == 2 && deps > 0.0) {
            kon = 20;
        }
        else if (kon == 1 && deps < 0.0) {
            kon = 10;
        }
        else if (kon == 10 && eps < stan) {
            kon = 2;
        }
        else if (kon == 20 && eps > stap) {
            kon = 1;
        }

      

        if (kon == 1) {
            if (eps < stap) {
                double stress1 = findstress1(eps, E1, Fy1, 0, rezaAA1, rezaNN1);
                sig = stress1;
            }
            else {
                if (stap < 0) {
                    double stress1 = findstress1(eps, E1, Fy1,2*stap, rezaAA1, rezaNN1);
                    double stress2 = findstress1(eps, E2, Fy2,stap, rezaAA2, rezaNN2);
                    sig = stress1 + stress2;
              
                }
                else if (stap > 0) {
                    double stress1 = findstress1(eps, E1, Fy1, 0, rezaAA1, rezaNN1);
                    double stress2 = findstress1(eps, E2, Fy2, stap, rezaAA2, rezaNN2);
                    sig = stress1 + stress2;
              
                }
              
            }
            double stress1p = findstressp(eps, Fyb, E1, Fy1, rezaAA1, rezaNN1);
            if (sig < stress1p) {
                sig = stress1p;
            }
        }
        else if (kon == 2) {
            if (eps > stan) {
                double stress1 = findstress1(eps, E1, Fy1, 0, rezaAA1, rezaNN1);
                sig = stress1;
            }
            else {

                if (stan < 0) {
                    double stress1 = findstress1(eps, E1, Fy1, 0, rezaAA1, rezaNN1);
                    double stress3 = findstress1(eps, E2, Fy2, stan, rezaAA2, rezaNN2);
                    sig= stress1 + stress3;
                }
                else if (stan > 0) {

                    double stress1 = findstress1(eps, E1, Fy1, 2*stan, rezaAA1, rezaNN1);
                    double stress3 = findstress1(eps, E2, Fy2, stan, rezaAA2, rezaNN2);
                    sig = stress1 +stress3;
                   
                }
            }
           double stress1n = findstressn(eps, Fyb, E1, Fy1, rezaAA1, rezaNN1);
            if (sig > stress1n) {
                sig = stress1n;
            }
        }
        else if (kon == 20) {

            double rotch = fabs(momz) / E2 * (1 + rezaAA2 * pow(fabs(momz) / (Fy2 * 2), rezaNN2)) - TrotMax;
           

            if (eps < rotch)
            {
                double stress1p = findstressp(eps, Fyb, E1, Fy1, rezaAA1, rezaNN1);
                double stress2p = findstressp(eps, momz, E2, Fy2, rezaAA2, rezaNN2);
                sig = stress1p + stress2p;
              
            }
            else
            {
                double stress1p = findstressp(eps, Fyb, E1, Fy1, rezaAA1, rezaNN1);
                sig = stress1p;
               
            }
        }
        else if (kon == 10) {
            double rotc = -fabs(momz) / E2 * (1 + rezaAA2 * pow(fabs(momz) / (Fy2 * 2), rezaNN2)) + TrotMax;
            if (eps > rotc)
            {
                double stress1n = findstressn(eps, Fyb, E1, Fy1, rezaAA1, rezaNN1);  
                double stress2n = findstressn(eps, momz, E2, Fy2, rezaAA2, rezaNN2);
                sig = stress1n + stress2n;
            }
            else
            {
                double stress1n = findstressn(eps, Fyb, E1, Fy1, rezaAA1, rezaNN1); 
                sig = stress1n;
            }
        }

    }

    if (strain == TrotMax) {
        epsr = eps;
        sigr = sig;
    }

}
double APDMD::findstress1(double targetStrain, double K, double mq, double yi, double rezaAA, double rezaNN)
{
    double low = -1000.0;
    double high = 1000.0;
    double tolerance = 1e-8;
    while (high - low > tolerance) {
        double mid = (low + high) / 2;
        double currentStrain = (mid / K) * (1 + rezaAA * pow(fabs(mid / mq), rezaNN-1)) + yi;

        if (currentStrain < targetStrain) {
            low = mid;
        }
        else {
            high = mid;
        }
    }
    return (low + high) / 2;
}
double APDMD::findstressp(double targetStrain, double mz, double K, double  mq, double rezaAA, double rezaNN)
{
    double low = -1000.0;
    double high = 1000.0;
    double tolerance = 1e-8;
    while (high - low > tolerance) {
        double mid = (low + high) / 2;
        double currentStrain = (mid + mz) / K * (1 + rezaAA * pow(fabs(mid + mz) / mq / 2, rezaNN-1)) - TrotMax;
        if (currentStrain < targetStrain) {
            low = mid;
        }
        else {
            high = mid;
        }
    }
    return (low + high) / 2;
}
double APDMD::findstressn(double targetStrain, double mz, double K, double  mq, double rezaAA, double rezaNN)
{
    double low = -1000.0;
    double high = 1000.0;
    double tolerance = 1e-6;
    while (high - low > tolerance) {
        double mid = (low + high) / 2;
        double currentStrain = ((mid - mz) / K) * (1 + rezaAA * pow(fabs(mid - mz) / mq / 2, rezaNN-1)) + TrotMax;

        if (currentStrain < targetStrain) {
            low = mid;
        }
        else {
            high = mid;
        }
    }
    return (low + high) / 2;
}

double
APDMD::getStrain(void)
{
  return eps;
}

double
APDMD::getStress(void)
{
  return sig;
}

double
APDMD::getTangent(void)
{
    return E1;
}

int
APDMD::commitState(void)
{
  epsminP2 = epsmin2;
  epsmaxP2 = epsmax2;
  epsminP2 = epsmin2;
  epsmaxP2 = epsmax2;
  epsplP = epspl;
  epssAP = epssA;
  sigsAP = sigsA;
  epssCP = epssC;
  sigsCP = sigsC;
  epssrP = epsr;
  sigsrP = sigr;
  konP = kon;
  CrotMax = TrotMax;

  sigP = sig;
  epsP = eps;

  return 0;
}

int
APDMD::revertToLastCommit(void)
{
  epsmin2 = epsminP2;
  epsmax2 = epsmaxP2;
  epsmin2 = epsminP2;
  epsmax2 = epsmaxP2;
  epspl = epsplP;
  epssA = epssAP;
  sigsA = sigsAP;
  epssC = epssCP;
  sigsC = sigsCP;
  epsr = epssrP;
  sigr = sigsrP;
  kon = konP;
  TrotMax = CrotMax;

  sig = sigP;
  eps = epsP;

  return 0;
}

int
APDMD::revertToStart(void)
{
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;

  konP = 0;
  epsmaxP1 = Fy1/E1;
  epsminP1 = -epsmaxP1;
  epsmaxP2 = Fy2 / E2 + ad2;
  epsminP2 = -epsmaxP2;
  CrotMax = 0.0;
  epsplP = 0.0;
  epssAP = 0.0;
  sigsAP = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

  return 0;
}

int
APDMD::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(20);
  data(0) = Fy1;
  data(1) = E1;
  data(2) = Fy2;
  data(3) = E2;
  data(4) = epsminP1;
  data(5) = epsmaxP1;
  data(6) = epsminP2;
  data(7) = epsmaxP2;
  data(8) = epsplP;
  data(9) = epssAP;
  data(10) = sigsAP;
  data(11) = epssAP;
  data(12) = sigsCP;
  data(13) = epssCP;
  data(14) = sigsrP;
  data(15) = konP;  
  data(16) = epsP;  
  data(17) = sigP;  
  data(18) = this->getTag();
  data(19) = sigini;

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "APDMD::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int
APDMD::recvSelf(int commitTag, Channel &theChannel,
             FEM_ObjectBroker &theBroker)
{
  static Vector data(20);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "APDMD::recvSelf() - failed to recvSelf\n";
    return -1;
  }

  Fy1 = data(0);
  E1 = data(1);
  Fy2 = data(2);
  E2 =data(3)  ;
  epsminP1 = data(4);
  epsmaxP1 = data(5);
  epsminP2 = data(6);
  epsmaxP2 = data(7);
  epsplP = data(8);
  epssAP = data(9);
  sigsAP = data(10);
  epssAP = data(11);
  sigsCP = data(12);
  epssCP = data(13);
  sigsrP = data(14);
  konP = data(15);  
  epsP = data(16);  
  sigP = data(17);  
  this->setTag(data(18));
  sigini = data(19);


  sig = sigP;
  eps = epsP;
 
  return 0;
}

void
APDMD::Print(OPS_Stream &s, int flag)
{
  s << "APDMD:(strain, stress, tangent) " << eps << " " << sig  << endln;
}

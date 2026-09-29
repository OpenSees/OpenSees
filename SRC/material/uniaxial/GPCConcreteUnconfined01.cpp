#include <GPCConcreteUnconfined01.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Channel.h>
#include <Vector.h>
#include <math.h>

const double GPCConcreteUnconfined01::TENS_STIFF_RATIO = 1.0e-6;

// uniaxialMaterial GPCUnconfined tag? fck? ecu?
void* OPS_GPCUnconfined()
{
    int iData[1];
    double dData[2];
    int numData = 1;

    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial GPCUnconfined tag" << endln;
        return 0;
    }
    numData = OPS_GetNumRemainingInputArgs();
    if (numData != 2) {
        opserr << "Invalid #args, want: uniaxialMaterial GPCUnconfined tag? fck? ecu?" << endln;
        return 0;
    }
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid args, want: uniaxialMaterial GPCUnconfined " << iData[0]
               << " fck? ecu?" << endln;
        return 0;
    }
    if (dData[0] == 0.0 || dData[1] == 0.0) {
        opserr << "GPCUnconfined material: fck and ecu must be nonzero" << endln;
        return 0;
    }

    UniaxialMaterial* theMaterial =
        new GPCConcreteUnconfined01(iData[0], dData[0], dData[1]);
    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type GPCUnconfined" << endln;
        return 0;
    }
    return theMaterial;
}

GPCConcreteUnconfined01::GPCConcreteUnconfined01(int tag, double FCK, double ECU)
 : UniaxialMaterial(tag, MAT_TAG_GPCConcreteUnconfined01),
   fck(fabs(FCK)), ecu(fabs(ECU))
{
    if (fck <= 0.0 || ecu <= 0.0) {
        opserr << "GPCConcreteUnconfined01: fck and ecu must be nonzero" << endln;
    }
    E0 = 6.75 * fck / ecu;   // d(sigma)/d(eps) at eps = 0

    ceps = 0.0; csig = 0.0; ctang = E0; ceps_max = 0.0;
    teps = 0.0; tsig = 0.0; ttang = E0; teps_max = 0.0;
}

GPCConcreteUnconfined01::GPCConcreteUnconfined01()
 : UniaxialMaterial(0, MAT_TAG_GPCConcreteUnconfined01),
   fck(0.0), ecu(0.0), E0(0.0)
{
    ceps = 0; csig = 0; ctang = 0; ceps_max = 0;
    teps = 0; tsig = 0; ttang = 0; teps_max = 0;
}

GPCConcreteUnconfined01::~GPCConcreteUnconfined01() {}

// Internal envelope, COMPRESSION POSITIVE (paper convention)
void GPCConcreteUnconfined01::envelope(double eps, double &sig, double &tang)
{
    if (eps < 0.0) {              // tension: no resistance
        sig = 0.0;
        tang = TENS_STIFF_RATIO * E0;
        return;
    }
    if (eps >= ecu) {             // beyond ultimate strain: crushed
        sig = 0.0;
        tang = TENS_STIFF_RATIO * E0;
        return;
    }
    double a = 6.75 * fck / (ecu * ecu * ecu);
    double d = eps - ecu;
    sig  = a * eps * d * d;
    tang = a * (d * d + 2.0 * eps * d);
}

int GPCConcreteUnconfined01::setTrialStrain(double strain, double strainRate)
{
    teps = strain;               // OpenSees convention: compression negative
    double e = -strain;          // internal: compression positive

    double sComp, tComp;
    if (e >= ceps_max) {
        // virgin loading: on the envelope
        envelope(e, sComp, tComp);
        teps_max = e;
    } else {
        // unload/reload: single straight line (slope E0) through the
        // reversal point; same line for both directions
        double sMax, tMax;
        envelope(ceps_max, sMax, tMax);
        double sLin = sMax - E0 * (ceps_max - e);
        if (sLin < 0.0) {
            sComp = 0.0;
            tComp = TENS_STIFF_RATIO * E0;
        } else {
            sComp = sLin;
            tComp = E0;
        }
        teps_max = ceps_max;
    }

    tsig  = -sComp;   // back to OpenSees sign convention
    ttang = tComp;    // d(-sigma)/d(-eps) = d(sigma)/d(eps)
    return 0;
}

double GPCConcreteUnconfined01::getStrain(void)  { return teps; }
double GPCConcreteUnconfined01::getStress(void)  { return tsig; }
double GPCConcreteUnconfined01::getTangent(void) { return ttang; }
double GPCConcreteUnconfined01::getInitialTangent(void) { return E0; }

int GPCConcreteUnconfined01::commitState(void)
{
    ceps = teps; csig = tsig; ctang = ttang; ceps_max = teps_max;
    return 0;
}

int GPCConcreteUnconfined01::revertToLastCommit(void)
{
    teps = ceps; tsig = csig; ttang = ctang; teps_max = ceps_max;
    return 0;
}

int GPCConcreteUnconfined01::revertToStart(void)
{
    ceps = 0; csig = 0; ctang = E0; ceps_max = 0;
    teps = 0; tsig = 0; ttang = E0; teps_max = 0;
    return 0;
}

UniaxialMaterial *GPCConcreteUnconfined01::getCopy(void)
{
    GPCConcreteUnconfined01 *theCopy =
        new GPCConcreteUnconfined01(this->getTag(), fck, ecu);
    theCopy->ceps = ceps; theCopy->csig = csig;
    theCopy->ctang = ctang; theCopy->ceps_max = ceps_max;
    theCopy->teps = teps; theCopy->tsig = tsig;
    theCopy->ttang = ttang; theCopy->teps_max = teps_max;
    return theCopy;
}

int GPCConcreteUnconfined01::sendSelf(int commitTag, Channel &theChannel)
{
    static Vector data(6);
    data(0) = this->getTag();
    data(1) = fck; data(2) = ecu;
    data(3) = ceps; data(4) = csig; data(5) = ceps_max;
    if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "GPCConcreteUnconfined01::sendSelf() - failed to send data" << endln;
        return -1;
    }
    return 0;
}

int GPCConcreteUnconfined01::recvSelf(int commitTag, Channel &theChannel,
                                      FEM_ObjectBroker &theBroker)
{
    static Vector data(6);
    if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "GPCConcreteUnconfined01::recvSelf() - failed to receive data" << endln;
        return -1;
    }
    this->setTag(int(data(0)));
    fck = data(1); ecu = data(2);
    ceps = data(3); csig = data(4); ceps_max = data(5);
    E0 = 6.75 * fck / ecu;
    ctang = E0;
    teps = ceps; tsig = csig; ttang = ctang; teps_max = ceps_max;
    return 0;
}

void GPCConcreteUnconfined01::Print(OPS_Stream &os, int flag)
{
    os << "GPCConcreteUnconfined01, tag: " << this->getTag()
       << ", fck: " << fck << ", ecu: " << ecu << endln;
}
#include <GPCConcrete01.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Channel.h>
#include <Vector.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double GPCConcrete01::TENS_STIFF_RATIO = 1.0e-6;

// Unconfined:
//   uniaxialMaterial GPCConcrete01 tag? fck? ecu?
// Confined (rectangular section):
//   uniaxialMaterial GPCConcrete01 tag? fck? ecu? fyh? diaLat? bk? hk? s? <Mfck?>
void* OPS_GPCConcrete01()
{
    int iData[1];
    double dData[8];
    int numData = 1;

    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial GPCConcrete01 tag" << endln;
        return 0;
    }

    numData = OPS_GetNumRemainingInputArgs();
    if (numData != 2 && numData != 7 && numData != 8) {
        opserr << "Invalid #args, want: uniaxialMaterial GPCConcrete01 tag? fck? ecu? "
               << "<fyh? diaLat? bk? hk? s? <Mfck?>>" << endln;
        return 0;
    }

    dData[7] = 0.7;   // default Mfck
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid args for uniaxialMaterial GPCConcrete01 " << iData[0] << endln;
        return 0;
    }
    // dData: 0 fck, 1 ecu, 2 fyh, 3 diaLat, 4 bk, 5 hk, 6 s, 7 Mfck

    if (dData[0] == 0.0 || dData[1] == 0.0) {
        opserr << "GPCConcrete01 material: fck and ecu must be nonzero" << endln;
        return 0;
    }

    UniaxialMaterial* theMaterial = 0;
    if (numData == 2) {
        theMaterial = new GPCConcrete01(iData[0], dData[0], dData[1]);
    } else {
        if (dData[4] == 0.0 || dData[5] == 0.0 || dData[6] == 0.0) {
            opserr << "GPCConcrete01 material: bk, hk and s must be nonzero" << endln;
            return 0;
        }
        theMaterial = new GPCConcrete01(iData[0], dData[0], dData[1],
                                        dData[2], dData[3], dData[4],
                                        dData[5], dData[6], dData[7]);
    }

    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type GPCConcrete01" << endln;
        return 0;
    }
    return theMaterial;
}

// Unconfined constructor
GPCConcrete01::GPCConcrete01(int tag, double FCK, double ECU)
 : UniaxialMaterial(tag, MAT_TAG_GPCConcrete01),
   fck(fabs(FCK)), ecu(fabs(ECU)), confined(false),
   fyh(0.0), diaLat(0.0), bk(0.0), hk(0.0), s(0.0), Mfck(0.7)
{
    computeEnvelopeParams();
    ceps = 0.0; csig = 0.0; ctang = E0; ceps_max = 0.0;
    teps = 0.0; tsig = 0.0; ttang = E0; teps_max = 0.0;
}

// Confined constructor
GPCConcrete01::GPCConcrete01(int tag, double FCK, double ECU,
                             double FYH, double DIALAT,
                             double BK, double HK, double S, double MFCK)
 : UniaxialMaterial(tag, MAT_TAG_GPCConcrete01),
   fck(fabs(FCK)), ecu(fabs(ECU)), confined(true),
   fyh(fabs(FYH)), diaLat(fabs(DIALAT)), bk(fabs(BK)), hk(fabs(HK)),
   s(fabs(S)), Mfck(MFCK)
{
    computeEnvelopeParams();
    ceps = 0.0; csig = 0.0; ctang = E0; ceps_max = 0.0;
    teps = 0.0; tsig = 0.0; ttang = E0; teps_max = 0.0;
}

GPCConcrete01::GPCConcrete01()
 : UniaxialMaterial(0, MAT_TAG_GPCConcrete01),
   fck(0.0), ecu(0.0), confined(false),
   fyh(0.0), diaLat(0.0), bk(0.0), hk(0.0), s(0.0), Mfck(0.7),
   fcc(0.0), eccu(0.0), E0(0.0)
{
    ceps = 0.0; csig = 0.0; ctang = 0.0; ceps_max = 0.0;
    teps = 0.0; tsig = 0.0; ttang = 0.0; teps_max = 0.0;
}

GPCConcrete01::~GPCConcrete01() {}

// Sets fcc, eccu and E0 from the input parameters.
void GPCConcrete01::computeEnvelopeParams(void)
{
    if (fck <= 0.0 || ecu <= 0.0) {
        opserr << "GPCConcrete01: fck and ecu must be nonzero" << endln;
    }

    fcc  = fck;
    eccu = ecu;

    if (confined) {
        if (s <= 0.0 || bk <= 0.0 || hk <= 0.0 || fck <= 0.0) {
            opserr << "GPCConcrete01: confined parameters require nonzero "
                   << "fck, bk, hk and s; using unconfined envelope" << endln;
        } else {
            double ls    = 2.0 * (bk + hk);
            double rho_s = diaLat * diaLat * (M_PI / 4.0) * ls / (s * bk * hk);
            double K     = 1.0 + rho_s * fyh / fck;

            double e50u = 2.0 * ecu / 3.0;
            double e50h = 2.0 * rho_s * sqrt(bk / s) / 3.0;

            eccu = 4.8 * (e50u + e50h) / 3.6;
            fcc  = Mfck * K * fck;
        }
    }

    if (eccu <= 0.0 || fcc <= 0.0) {
        opserr << "GPCConcrete01: computed fcc/eccu are non-positive; "
               << "check input values" << endln;
    }

    E0 = 6.75 * fcc / eccu;   // d(sigma)/d(eps) at eps = 0
}

// Internal envelope, COMPRESSION POSITIVE
void GPCConcrete01::envelope(double eps, double &sig, double &tang)
{
    if (eps < 0.0 || eps >= eccu) {   // tension, or beyond ultimate strain
        sig  = 0.0;
        tang = TENS_STIFF_RATIO * E0;
        return;
    }
    double a = 6.75 * fcc / (eccu * eccu * eccu);
    double d = eps - eccu;
    sig  = a * eps * d * d;
    tang = a * (d * d + 2.0 * eps * d);
}

int GPCConcrete01::setTrialStrain(double strain, double strainRate)
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

double GPCConcrete01::getStrain(void)  { return teps; }
double GPCConcrete01::getStress(void)  { return tsig; }
double GPCConcrete01::getTangent(void) { return ttang; }
double GPCConcrete01::getInitialTangent(void) { return E0; }

int GPCConcrete01::commitState(void)
{
    ceps = teps; csig = tsig; ctang = ttang; ceps_max = teps_max;
    return 0;
}

int GPCConcrete01::revertToLastCommit(void)
{
    teps = ceps; tsig = csig; ttang = ctang; teps_max = ceps_max;
    return 0;
}

int GPCConcrete01::revertToStart(void)
{
    ceps = 0.0; csig = 0.0; ctang = E0; ceps_max = 0.0;
    teps = 0.0; tsig = 0.0; ttang = E0; teps_max = 0.0;
    return 0;
}

UniaxialMaterial *GPCConcrete01::getCopy(void)
{
    GPCConcrete01 *theCopy;
    if (confined)
        theCopy = new GPCConcrete01(this->getTag(), fck, ecu,
                                    fyh, diaLat, bk, hk, s, Mfck);
    else
        theCopy = new GPCConcrete01(this->getTag(), fck, ecu);

    theCopy->ceps = ceps; theCopy->csig = csig;
    theCopy->ctang = ctang; theCopy->ceps_max = ceps_max;
    theCopy->teps = teps; theCopy->tsig = tsig;
    theCopy->ttang = ttang; theCopy->teps_max = teps_max;
    return theCopy;
}

int GPCConcrete01::sendSelf(int commitTag, Channel &theChannel)
{
    static Vector data(13);
    data(0)  = this->getTag();
    data(1)  = fck;
    data(2)  = ecu;
    data(3)  = confined ? 1.0 : 0.0;
    data(4)  = fyh;
    data(5)  = diaLat;
    data(6)  = bk;
    data(7)  = hk;
    data(8)  = s;
    data(9)  = Mfck;
    data(10) = ceps;
    data(11) = csig;
    data(12) = ceps_max;
    if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "GPCConcrete01::sendSelf() - failed to send data" << endln;
        return -1;
    }
    return 0;
}

int GPCConcrete01::recvSelf(int commitTag, Channel &theChannel,
                            FEM_ObjectBroker &theBroker)
{
    static Vector data(13);
    if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "GPCConcrete01::recvSelf() - failed to receive data" << endln;
        return -1;
    }
    this->setTag(int(data(0)));
    fck      = data(1);
    ecu      = data(2);
    confined = (data(3) != 0.0);
    fyh      = data(4);
    diaLat   = data(5);
    bk       = data(6);
    hk       = data(7);
    s        = data(8);
    Mfck     = data(9);
    ceps     = data(10);
    csig     = data(11);
    ceps_max = data(12);

    computeEnvelopeParams();
    ctang = E0;
    teps = ceps; tsig = csig; ttang = ctang; teps_max = ceps_max;
    return 0;
}

void GPCConcrete01::Print(OPS_Stream &os, int flag)
{
    os << "GPCConcrete01, tag: " << this->getTag()
       << ", fck: " << fck << ", ecu: " << ecu;
    if (confined) {
        os << ", confined: yes, fcc(derived): " << fcc
           << ", eccu(derived): " << eccu;
    } else {
        os << ", confined: no";
    }
    os << endln;
}

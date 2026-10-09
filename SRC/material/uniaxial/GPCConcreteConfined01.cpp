#include <GPCConcreteConfined01.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Channel.h>
#include <Vector.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double GPCConcreteConfined01::TENS_STIFF_RATIO = 1.0e-6;

// uniaxialMaterial GPCConfined tag? fck? ecu? fyh? diaLat? bk? hk? s? <Mfck?>
void* OPS_GPCConfined()
{
    int iData[1];
    double dData[8];
    int numData = 1;

    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial GPCConfined tag" << endln;
        return 0;
    }
    numData = OPS_GetNumRemainingInputArgs();
    if (numData != 7 && numData != 8) {
        opserr << "Invalid #args, want: uniaxialMaterial GPCConfined tag? fck? ecu? "
               << "fyh? diaLat? bk? hk? s? <Mfck?>" << endln;
        return 0;
    }
    dData[7] = 0.7;   // default Mfck
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid args for uniaxialMaterial GPCConfined " << iData[0] << endln;
        return 0;
    }
    // dData: 0 fck, 1 ecu, 2 fyh, 3 diaLat, 4 bk, 5 hk, 6 s, 7 Mfck
    if (dData[0] == 0.0 || dData[1] == 0.0 || dData[4] == 0.0 ||
        dData[5] == 0.0 || dData[6] == 0.0) {
        opserr << "GPCConfined material: fck, ecu, bk, hk and s must be nonzero" << endln;
        return 0;
    }

    UniaxialMaterial* theMaterial = new GPCConcreteConfined01(
        iData[0], dData[0], dData[1], dData[2], dData[3],
        dData[4], dData[5], dData[6], dData[7]);
    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type GPCConfined" << endln;
        return 0;
    }
    return theMaterial;
}

GPCConcreteConfined01::GPCConcreteConfined01(int tag,
    double FCK, double ECU, double FYH, double DIALAT,
    double BK, double HK, double S, double MFCK)
 : UniaxialMaterial(tag, MAT_TAG_GPCConcreteConfined01),
   fck(fabs(FCK)), ecu(fabs(ECU)), fyh(fabs(FYH)), diaLat(fabs(DIALAT)),
   bk(fabs(BK)), hk(fabs(HK)), s(fabs(S)), Mfck(MFCK)
{
    computeConfinedParams();
    E0 = 6.75 * fcc / eccu;
    ceps = 0; csig = 0; ctang = E0; ceps_max = 0;
    teps = 0; tsig = 0; ttang = E0; teps_max = 0;
}

GPCConcreteConfined01::GPCConcreteConfined01()
 : UniaxialMaterial(0, MAT_TAG_GPCConcreteConfined01),
   fck(0), ecu(0), fyh(0), diaLat(0), bk(0), hk(0), s(0), Mfck(0.7),
   fcc(0), eccu(0), E0(0)
{
    ceps = 0; csig = 0; ctang = 0; ceps_max = 0;
    teps = 0; tsig = 0; ttang = 0; teps_max = 0;
}

GPCConcreteConfined01::~GPCConcreteConfined01() {}

void GPCConcreteConfined01::computeConfinedParams(void)
{
    double ls    = 2.0 * (bk + hk);
    double rho_s = diaLat * diaLat * (M_PI / 4.0) * ls / (s * bk * hk);
    double K     = 1.0 + rho_s * fyh / fck;

    double e50u = 2.0 * ecu / 3.0;
    double e50h = 2.0 * rho_s * sqrt(bk / s) / 3.0;

    eccu = 4.8 * (e50u + e50h) / 3.6;
    fcc  = Mfck * K * fck;

    if (eccu <= 0.0 || fcc <= 0.0) {
        opserr << "GPCConcreteConfined01: computed fcc/eccu are non-positive; "
               << "check input values" << endln;
    }
}

// Internal envelope, COMPRESSION POSITIVE
void GPCConcreteConfined01::envelope(double eps, double &sig, double &tang)
{
    if (eps < 0.0 || eps >= eccu) {
        sig = 0.0;
        tang = TENS_STIFF_RATIO * E0;
        return;
    }
    double a = 6.75 * fcc / (eccu * eccu * eccu);
    double d = eps - eccu;
    sig  = a * eps * d * d;
    tang = a * (d * d + 2.0 * eps * d);
}

int GPCConcreteConfined01::setTrialStrain(double strain, double strainRate)
{
    teps = strain;               // OpenSees convention: compression negative
    double e = -strain;          // internal: compression positive

    double sComp, tComp;
    if (e >= ceps_max) {
        envelope(e, sComp, tComp);
        teps_max = e;
    } else {
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

    tsig  = -sComp;
    ttang = tComp;
    return 0;
}

double GPCConcreteConfined01::getStrain(void)  { return teps; }
double GPCConcreteConfined01::getStress(void)  { return tsig; }
double GPCConcreteConfined01::getTangent(void) { return ttang; }
double GPCConcreteConfined01::getInitialTangent(void) { return E0; }

int GPCConcreteConfined01::commitState(void)
{
    ceps = teps; csig = tsig; ctang = ttang; ceps_max = teps_max;
    return 0;
}

int GPCConcreteConfined01::revertToLastCommit(void)
{
    teps = ceps; tsig = csig; ttang = ctang; teps_max = ceps_max;
    return 0;
}

int GPCConcreteConfined01::revertToStart(void)
{
    ceps = 0; csig = 0; ctang = E0; ceps_max = 0;
    teps = 0; tsig = 0; ttang = E0; teps_max = 0;
    return 0;
}

UniaxialMaterial *GPCConcreteConfined01::getCopy(void)
{
    GPCConcreteConfined01 *theCopy = new GPCConcreteConfined01(
        this->getTag(), fck, ecu, fyh, diaLat, bk, hk, s, Mfck);
    theCopy->ceps = ceps; theCopy->csig = csig;
    theCopy->ctang = ctang; theCopy->ceps_max = ceps_max;
    theCopy->teps = teps; theCopy->tsig = tsig;
    theCopy->ttang = ttang; theCopy->teps_max = teps_max;
    return theCopy;
}

int GPCConcreteConfined01::sendSelf(int commitTag, Channel &theChannel)
{
    static Vector data(12);
    data(0) = this->getTag();
    data(1) = fck; data(2) = ecu; data(3) = fyh; data(4) = diaLat;
    data(5) = bk; data(6) = hk; data(7) = s; data(8) = Mfck;
    data(9) = ceps; data(10) = csig; data(11) = ceps_max;
    if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "GPCConcreteConfined01::sendSelf() - failed to send data" << endln;
        return -1;
    }
    return 0;
}

int GPCConcreteConfined01::recvSelf(int commitTag, Channel &theChannel,
                                    FEM_ObjectBroker &theBroker)
{
    static Vector data(12);
    if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "GPCConcreteConfined01::recvSelf() - failed to receive data" << endln;
        return -1;
    }
    this->setTag(int(data(0)));
    fck = data(1); ecu = data(2); fyh = data(3); diaLat = data(4);
    bk = data(5); hk = data(6); s = data(7); Mfck = data(8);
    ceps = data(9); csig = data(10); ceps_max = data(11);
    computeConfinedParams();
    E0 = 6.75 * fcc / eccu;
    ctang = E0;
    teps = ceps; tsig = csig; ttang = ctang; teps_max = ceps_max;
    return 0;
}

void GPCConcreteConfined01::Print(OPS_Stream &os, int flag)
{
    os << "GPCConcreteConfined01, tag: " << this->getTag()
       << ", fck: " << fck << ", fcc(derived): " << fcc
       << ", eccu(derived): " << eccu << endln;
}
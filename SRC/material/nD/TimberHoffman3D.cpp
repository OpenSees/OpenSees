
#include "TimberHoffman3D.h"

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <cmath>
#include <cstring>
#include <algorithm>

namespace {

double vecNormInf(const Vector &v)
{
    double m = 0.0;
    for (int i = 0; i < v.Size(); ++i)
        m = std::max(m, std::fabs(v(i)));
    return m;
}

double quadForm(const Vector &v, const Matrix &A)
{
    Vector Av(v.Size());
    Av.addMatrixVector(0.0, A, v, 1.0);
    return v ^ Av;
}

}

void *OPS_TimberHoffman3D()
{
    if (OPS_GetNumRemainingInputArgs() < 28) {
        opserr
            << "WARNING insufficient arguments for TimberHoffman3D\n"
            << "Want: nDMaterial TimberHoffman3D tag "
            << "E1 E2 E3 nu12 nu13 nu23 G12 G13 G23 "
            << "fc1 fc2 fc3 ft1 ft2 ft3 f12 f13 f23 "
            << "h sigmaE0 Acomp Bcomp Gf1t Gf2t Gf3t eta Lc <dt>\n";
        return nullptr;
    }

    int tag;
    int nInt = 1;
    if (OPS_GetIntInput(&nInt, &tag) < 0) {
        opserr << "WARNING invalid TimberHoffman3D tag\n";
        return nullptr;
    }

    double data[27];
    int nDouble = 27;

    if (OPS_GetDoubleInput(&nDouble, data) < 0) {
        opserr << "WARNING invalid TimberHoffman3D material parameters\n";
        return nullptr;
    }

    double dt = 1.0;
    if (OPS_GetNumRemainingInputArgs() > 0) {
        int one = 1;
        if (OPS_GetDoubleInput(&one, &dt) < 0) {
            opserr << "WARNING invalid TimberHoffman3D dt\n";
            return nullptr;
        }
    }

    if (data[0] <= 0.0 || data[1] <= 0.0 || data[2] <= 0.0 ||
        data[6] <= 0.0 || data[7] <= 0.0 || data[8] <= 0.0) {
        opserr << "WARNING TimberHoffman3D elastic moduli must be positive\n";
        return nullptr;
    }

    if (data[9] <= 0.0 || data[10] <= 0.0 || data[11] <= 0.0 ||
        data[12] <= 0.0 || data[13] <= 0.0 || data[14] <= 0.0 ||
        data[15] <= 0.0 || data[16] <= 0.0 || data[17] <= 0.0) {
        opserr << "WARNING TimberHoffman3D strength parameters must be positive\n";
        return nullptr;
    }

    if (data[19] <= 0.0) {
        opserr << "WARNING TimberHoffman3D sigmaE0 must be positive\n";
        return nullptr;
    }

    if (data[22] <= 0.0 || data[23] <= 0.0 || data[24] <= 0.0) {
        opserr << "WARNING TimberHoffman3D fracture energies must be positive\n";
        return nullptr;
    }

    if (data[25] < 0.0 || data[26] <= 0.0 || dt <= 0.0) {
        opserr << "WARNING TimberHoffman3D requires eta >= 0, Lc > 0, and dt > 0\n";
        return nullptr;
    }
    return new TimberHoffman3D(
        tag,
        data[0], data[1], data[2],
        data[3], data[4], data[5],
        data[6], data[7], data[8],
        data[9], data[10], data[11],
        data[12], data[13], data[14],
        data[15], data[16], data[17],
        data[18],
        data[19],
        data[20], data[21],
        data[22], data[23], data[24],
        data[25],
        data[26],
        dt
    );
}

TimberHoffman3D::TimberHoffman3D()
: NDMaterial(0, ND_TAG_TimberHoffman3D),
  E1(0), E2(0), E3(0),
  nu12(0), nu13(0), nu23(0),
  G12(0), G13(0), G23(0),
  fc1(0), fc2(0), fc3(0),
  ft1(0), ft2(0), ft3(0),
  f12(0), f13(0), f23(0),
  hardeningModulus(0), sigmaE0(1.0),
  Acomp(0), Bcomp(0),
  Gf1t(0), Gf2t(0), Gf3t(0),
  eta(1e-4), Lc(1.0), dt(1.0),
  C0(6,6), Cdam(6,6), P(6,6), Q0(6), Z(6,6),
  strainTrial(6), stressTrial(6), effectiveStressTrial(6),
  plasticStrainTrial(6),
  eqpTrial(0.0), damageTTrial(3), damageC1Trial(0.0),
  damageVTrial(3), FmaxTTrial(3), FmaxC1Trial(1.0),
  strainCommit(6), stressCommit(6), effectiveStressCommit(6),
  plasticStrainCommit(6),
  eqpCommit(0.0), damageTCommit(3), damageC1Commit(0.0),
  damageVCommit(3), FmaxTCommit(3), FmaxC1Commit(1.0),
  tol(1e-8), maxIter(60), innerTol(1e-10), innerMaxIter(30)
{
    revertToStart();
}

TimberHoffman3D::TimberHoffman3D(
    int tag,
    double E1_, double E2_, double E3_,
    double nu12_, double nu13_, double nu23_,
    double G12_, double G13_, double G23_,
    double fc1_, double fc2_, double fc3_,
    double ft1_, double ft2_, double ft3_,
    double f12_, double f13_, double f23_,
    double h_,
    double sigmaE0_,
    double Acomp_, double Bcomp_,
    double Gf1t_, double Gf2t_, double Gf3t_,
    double eta_,
    double Lc_,
    double dt_
)
: NDMaterial(tag, ND_TAG_TimberHoffman3D),
  E1(E1_), E2(E2_), E3(E3_),
  nu12(nu12_), nu13(nu13_), nu23(nu23_),
  G12(G12_), G13(G13_), G23(G23_),
  fc1(fc1_), fc2(fc2_), fc3(fc3_),
  ft1(ft1_), ft2(ft2_), ft3(ft3_),
  f12(f12_), f13(f13_), f23(f23_),
  hardeningModulus(h_), sigmaE0(sigmaE0_),
  Acomp(Acomp_), Bcomp(Bcomp_),
  Gf1t(Gf1t_), Gf2t(Gf2t_), Gf3t(Gf3t_),
  eta(eta_), Lc(Lc_), dt(dt_),
  C0(6,6), Cdam(6,6), P(6,6), Q0(6), Z(6,6),
  strainTrial(6), stressTrial(6), effectiveStressTrial(6),
  plasticStrainTrial(6),
  eqpTrial(0.0), damageTTrial(3), damageC1Trial(0.0),
  damageVTrial(3), FmaxTTrial(3), FmaxC1Trial(1.0),
  strainCommit(6), stressCommit(6), effectiveStressCommit(6),
  plasticStrainCommit(6),
  eqpCommit(0.0), damageTCommit(3), damageC1Commit(0.0),
  damageVCommit(3), FmaxTCommit(3), FmaxC1Commit(1.0),
  tol(1e-8), maxIter(60), innerTol(1e-10), innerMaxIter(30)
{
    buildElasticStiffness();
    buildHoffmanPQ();

    Z.Zero();
    Z(0,0)=1.0; Z(1,1)=1.0; Z(2,2)=1.0;
    Z(3,3)=0.5; Z(4,4)=0.5; Z(5,5)=0.5;

    revertToStart();
}

TimberHoffman3D::~TimberHoffman3D() {}

void TimberHoffman3D::buildElasticStiffness()
{
    Matrix S(6,6);
    S.Zero();

    const double nu21 = nu12 * E2 / E1;
    const double nu31 = nu13 * E3 / E1;
    const double nu32 = nu23 * E3 / E2;

    S(0,0)=1.0/E1;
    S(1,1)=1.0/E2;
    S(2,2)=1.0/E3;

    S(0,1)=-nu21/E2;
    S(1,0)=-nu12/E1;

    S(0,2)=-nu31/E3;
    S(2,0)=-nu13/E1;

    S(1,2)=-nu32/E3;
    S(2,1)=-nu23/E2;

    S(3,3)=1.0/G12;
    S(4,4)=1.0/G13;
    S(5,5)=1.0/G23;

    if (S.Invert(C0) < 0) {
        opserr << "WARNING TimberHoffman3D failed to invert elastic compliance matrix\n";
        C0.Zero();
    }
    Cdam = C0;
}

void TimberHoffman3D::buildHoffmanPQ()
{
    const double se2 = sigmaE0*sigmaE0;

    const double a12 =
        0.5*se2*(1.0/(fc1*ft1) + 1.0/(fc2*ft2) - 1.0/(fc3*ft3));

    const double a23 =
        0.5*se2*(-1.0/(fc1*ft1) + 1.0/(fc2*ft2) + 1.0/(fc3*ft3));

    const double a13 =
        0.5*se2*(1.0/(fc1*ft1) - 1.0/(fc2*ft2) + 1.0/(fc3*ft3));

    const double a11 = se2*(fc1-ft1)/(fc1*ft1);
    const double a22 = se2*(fc2-ft2)/(fc2*ft2);
    const double a33 = se2*(fc3-ft3)/(fc3*ft3);

    const double a44 = se2/(3.0*f12*f12);
    const double a55 = se2/(3.0*f13*f13);
    const double a66 = se2/(3.0*f23*f23);

    P.Zero();
    P(0,0)=2.0*(a13+a12);
    P(1,1)=2.0*(a23+a12);
    P(2,2)=2.0*(a13+a23);

    P(0,1)=P(1,0)=-2.0*a12;
    P(0,2)=P(2,0)=-2.0*a13;
    P(1,2)=P(2,1)=-2.0*a23;

    P(3,3)=6.0*a44;
    P(4,4)=6.0*a55;
    P(5,5)=6.0*a66;

    Q0.Zero();
    Q0(0)=a11;
    Q0(1)=a22;
    Q0(2)=a33;
}

double TimberHoffman3D::sigmaEk(double eqp) const
{
    return sigmaE0 + hardeningModulus*eqp;
}

double TimberHoffman3D::yieldFunction(const Vector &s, double eqp) const
{
    const double sek = sigmaEk(eqp);
    return 0.5*quadForm(s,P) + (s ^ Q0) - sek*sek;
}

void TimberHoffman3D::flowDirection(
    const Vector &s, double eqp, Vector &g) const
{
    g.Zero();
    g.addMatrixVector(0.0, P, s, 1.0);
    g += Q0; // fixed-Q interpretation
}

double TimberHoffman3D::equivalentPlasticIncrement(
    const Vector &depsP) const
{
    return std::sqrt(std::max(0.0, (2.0/3.0)*quadForm(depsP,Z)));
}

int TimberHoffman3D::stateForLambda(
    const Vector &sigmaTrial,
    double dLambda,
    Vector &sigma,
    double &eqp,
    Vector &depsP
)
{
    if (dLambda <= 0.0) {
        sigma = sigmaTrial;
        eqp = eqpCommit;
        depsP.Zero();
        return 0;
    }

    Matrix A(6,6);
    A.Zero();
    for (int i=0;i<6;++i) A(i,i)=1.0;

    Matrix CP(6,6);
    CP.addMatrixProduct(0.0, C0, P, 1.0);
    A.addMatrix(1.0, CP, dLambda);

    Vector CQ(6);
    CQ.addMatrixVector(0.0, C0, Q0, 1.0);

    Vector rhs(sigmaTrial);
    rhs.addVector(1.0, CQ, -dLambda);

    Matrix Ainv(6,6);
    if (A.Invert(Ainv) < 0)
        return -1;

    sigma.addMatrixVector(0.0, Ainv, rhs, 1.0);

    Vector g(6);
    flowDirection(sigma, eqpCommit, g);

    depsP = g;
    depsP *= dLambda;

    eqp = eqpCommit + equivalentPlasticIncrement(depsP);
    return 0;
}

int TimberHoffman3D::returnMapHardening(
    const Vector &sigmaTrial,
    Vector &sigmaNew,
    double &eqpNew,
    Vector &depsP,
    double &dLambda
)
{
    bool anyCompression = false;
    for (int i=0;i<3;++i)
        if (sigmaTrial(i) < 0.0) anyCompression = true;

    if (!anyCompression || yieldFunction(sigmaTrial,eqpCommit) <= tol) {
        sigmaNew = sigmaTrial;
        eqpNew = eqpCommit;
        depsP.Zero();
        dLambda = 0.0;
        return 0;
    }

    double lo = 0.0;
    double hi = 1.0e-12;

    Vector sHi(6), dpHi(6);
    double kHi = eqpCommit;

    bool bracketed = false;
    for (int i=0;i<140;++i) {
        if (stateForLambda(sigmaTrial,hi,sHi,kHi,dpHi) < 0)
            return -1;

        if (yieldFunction(sHi,kHi) <= 0.0) {
            bracketed = true;
            break;
        }
        hi *= 10.0;
    }

    if (!bracketed) return -2;

    Vector sMid(6), dpMid(6);
    double kMid = eqpCommit;

    for (int i=0;i<maxIter;++i) {
        const double mid = 0.5*(lo+hi);

        if (stateForLambda(sigmaTrial,mid,sMid,kMid,dpMid) < 0)
            return -3;

        const double f = yieldFunction(sMid,kMid);

        if (std::fabs(f) <= tol) {
            sigmaNew = sMid;
            eqpNew = kMid;
            depsP = dpMid;
            dLambda = mid;
            return 0;
        }

        if (f > 0.0) lo = mid;
        else hi = mid;
    }

    dLambda = 0.5*(lo+hi);
    if (stateForLambda(sigmaTrial,dLambda,sigmaNew,eqpNew,depsP) < 0)
        return -4;

    return 0;
}

void TimberHoffman3D::failureIndices(
    const Vector &s,
    Vector &Ft,
    double &F1c
) const
{
    Ft.Zero();

    F1c = (s(0) < 0.0) ? (-s(0)/fc1) : 0.0;

    Ft(0) = (s(0) > 0.0) ? (s(0)/ft1) : 0.0;

    if (s(1) > 0.0) {
        Ft(1) =
            std::pow(s(1)/ft2,2.0) +
            std::pow(s(3)/f12,2.0) +
            std::pow(s(5)/f23,2.0);
    }

    if (s(2) > 0.0) {
        Ft(2) =
            std::pow(s(2)/ft3,2.0) +
            std::pow(s(4)/f13,2.0) +
            std::pow(s(5)/f23,2.0);
    }
}

double TimberHoffman3D::tensileDamage(double F, int dir) const
{
    if (F <= 1.0) return 0.0;

    double Ei=E1, ft=ft1, Gf=Gf1t;

    if (dir==2) { Ei=E2; ft=ft2; Gf=Gf2t; }
    if (dir==3) { Ei=E3; ft=ft3; Gf=Gf3t; }

    const double exponent =
        (1.0-F)*Lc*ft*ft/(Ei*Gf);

    double d = 1.0 - (1.0/F)*std::exp(exponent);
    return std::max(0.0,std::min(0.999999,d));
}

double TimberHoffman3D::compressionDamage(double F1c) const
{
    if (F1c <= 1.0) return 0.0;

    const double d =
        1.0
        - (1.0/F1c)*(1.0-Acomp)
        - Acomp*std::exp(Bcomp*(1.0-F1c));

    return std::max(0.0,std::min(0.999999,d));
}

void TimberHoffman3D::updateDamage(const Vector &stressEff)
{
    Vector Ft(3);
    double F1c = 0.0;
    failureIndices(stressEff,Ft,F1c);

    for (int i=0;i<3;++i)
        FmaxTTrial(i)=std::max(FmaxTCommit(i),Ft(i));

    FmaxC1Trial=std::max(FmaxC1Commit,F1c);

    for (int i=0;i<3;++i)
        damageTTrial(i)=std::max(
            damageTCommit(i),
            tensileDamage(FmaxTTrial(i),i+1)
        );

    damageC1Trial=std::max(
        damageC1Commit,
        compressionDamage(FmaxC1Trial)
    );

    Vector dRaw(3);
    dRaw(0)=std::max(damageTTrial(0),damageC1Trial);
    dRaw(1)=damageTTrial(1);
    dRaw(2)=damageTTrial(2);

    const double dtEff = (ops_Dt > 0.0) ? ops_Dt : dt;
    const double aOld = eta/(eta + dtEff);
    const double aNew = dtEff/(eta + dtEff);

    for (int i=0;i<3;++i) {
        damageVTrial(i)=
            aOld*damageVCommit(i)+aNew*dRaw(i);

        damageVTrial(i)=std::max(
            damageVCommit(i),
            damageVTrial(i)
        );

        damageVTrial(i)=std::max(
            0.0,
            std::min(0.999999,damageVTrial(i))
        );
    }
}

void TimberHoffman3D::buildDamagedStiffness(const Vector &d)
{
    const double d1=d(0);
    const double d2=d(1);
    const double d3=d(2);

    Cdam.Zero();

    Cdam(0,0)=(1.0-d1)*C0(0,0);
    Cdam(1,1)=(1.0-d2)*C0(1,1);
    Cdam(2,2)=(1.0-d3)*C0(2,2);

    Cdam(0,1)=Cdam(1,0)=(1.0-d1)*(1.0-d2)*C0(0,1);
    Cdam(0,2)=Cdam(2,0)=(1.0-d1)*(1.0-d3)*C0(0,2);
    Cdam(1,2)=Cdam(2,1)=(1.0-d2)*(1.0-d3)*C0(1,2);

    Cdam(3,3)=(1.0-d1)*(1.0-d2)*C0(3,3);
    Cdam(4,4)=(1.0-d1)*(1.0-d3)*C0(4,4);
    Cdam(5,5)=(1.0-d2)*(1.0-d3)*C0(5,5);
}

int TimberHoffman3D::setTrialStrain(const Vector &strain)
{
    if (strain.Size()!=6) return -1;

    strainTrial=strain;

    Vector deps(6);
    deps = strainTrial;
    deps.addVector(1.0,strainCommit,-1.0);

    Vector sigmaTrialEff(effectiveStressCommit);
    sigmaTrialEff.addMatrixVector(
        1.0,C0,deps,1.0
    );

    Vector depsP(6);
    double dLambda=0.0;

    if (returnMapHardening(
        sigmaTrialEff,
        effectiveStressTrial,
        eqpTrial,
        depsP,
        dLambda
    ) < 0)
        return -2;

    plasticStrainTrial=plasticStrainCommit;
    plasticStrainTrial += depsP;

    Vector elasticStrain(strainTrial);
    elasticStrain.addVector(
        1.0,
        plasticStrainTrial,
        -1.0
    );

    updateDamage(effectiveStressTrial);
    buildDamagedStiffness(damageVTrial);

    stressTrial.addMatrixVector(
        0.0,
        Cdam,
        elasticStrain,
        1.0
    );

    return 0;
}

int TimberHoffman3D::setTrialStrain(
    const Vector &strain,
    const Vector &rate)
{
    return setTrialStrain(strain);
}

const Vector &TimberHoffman3D::getStress(void)
{
    return stressTrial;
}

const Vector &TimberHoffman3D::getStrain(void)
{
    return strainTrial;
}

const Matrix &TimberHoffman3D::getTangent(void)
{
    // IMPORTANT:
    // This is the damaged elastic tangent, NOT yet the
    // fully consistent elastoplastic-damage algorithmic tangent.
    return Cdam;
}

const Matrix &TimberHoffman3D::getInitialTangent(void)
{
    return C0;
}

int TimberHoffman3D::commitState(void)
{
    strainCommit=strainTrial;
    stressCommit=stressTrial;
    effectiveStressCommit=effectiveStressTrial;
    plasticStrainCommit=plasticStrainTrial;

    eqpCommit=eqpTrial;

    damageTCommit=damageTTrial;
    damageC1Commit=damageC1Trial;
    damageVCommit=damageVTrial;

    FmaxTCommit=FmaxTTrial;
    FmaxC1Commit=FmaxC1Trial;

    return 0;
}

int TimberHoffman3D::revertToLastCommit(void)
{
    strainTrial=strainCommit;
    stressTrial=stressCommit;
    effectiveStressTrial=effectiveStressCommit;
    plasticStrainTrial=plasticStrainCommit;

    eqpTrial=eqpCommit;

    damageTTrial=damageTCommit;
    damageC1Trial=damageC1Commit;
    damageVTrial=damageVCommit;

    FmaxTTrial=FmaxTCommit;
    FmaxC1Trial=FmaxC1Commit;

    buildDamagedStiffness(damageVTrial);

    return 0;
}

int TimberHoffman3D::revertToStart(void)
{
    strainTrial.Zero();
    stressTrial.Zero();
    effectiveStressTrial.Zero();
    plasticStrainTrial.Zero();

    strainCommit.Zero();
    stressCommit.Zero();
    effectiveStressCommit.Zero();
    plasticStrainCommit.Zero();

    eqpTrial=0.0;
    eqpCommit=0.0;

    damageTTrial.Zero();
    damageTCommit.Zero();

    damageVTrial.Zero();
    damageVCommit.Zero();

    damageC1Trial=0.0;
    damageC1Commit=0.0;

    FmaxTTrial.Zero();
    FmaxTCommit.Zero();

    for (int i=0;i<3;++i) {
        FmaxTTrial(i)=1.0;
        FmaxTCommit(i)=1.0;
    }

    FmaxC1Trial=1.0;
    FmaxC1Commit=1.0;

    Cdam=C0;
    return 0;
}

NDMaterial *TimberHoffman3D::getCopy(void)
{
    return new TimberHoffman3D(*this);
}

NDMaterial *TimberHoffman3D::getCopy(const char *type)
{
    if (strcmp(type,"ThreeDimensional")==0 ||
        strcmp(type,"3D")==0)
        return getCopy();

    return nullptr;
}

const char *TimberHoffman3D::getType(void) const
{
    return "ThreeDimensional";
}

int TimberHoffman3D::getOrder(void) const
{
    return 6;
}

int TimberHoffman3D::sendSelf(int commitTag, Channel &theChannel)
{
    if (this->getDbTag() == 0)
        this->setDbTag(theChannel.getDbTag());

    static Vector data(69);
    int c = 0;
    data(c++) = this->getTag();
    data(c++) = E1; data(c++) = E2; data(c++) = E3;
    data(c++) = nu12; data(c++) = nu13; data(c++) = nu23;
    data(c++) = G12; data(c++) = G13; data(c++) = G23;
    data(c++) = fc1; data(c++) = fc2; data(c++) = fc3;
    data(c++) = ft1; data(c++) = ft2; data(c++) = ft3;
    data(c++) = f12; data(c++) = f13; data(c++) = f23;
    data(c++) = hardeningModulus; data(c++) = sigmaE0;
    data(c++) = Acomp; data(c++) = Bcomp;
    data(c++) = Gf1t; data(c++) = Gf2t; data(c++) = Gf3t;
    data(c++) = eta; data(c++) = Lc; data(c++) = dt;
    data(c++) = tol; data(c++) = maxIter; data(c++) = innerTol; data(c++) = innerMaxIter;

    for (int i=0;i<6;++i) data(c++) = strainCommit(i);
    for (int i=0;i<6;++i) data(c++) = stressCommit(i);
    for (int i=0;i<6;++i) data(c++) = effectiveStressCommit(i);
    for (int i=0;i<6;++i) data(c++) = plasticStrainCommit(i);
    data(c++) = eqpCommit;
    for (int i=0;i<3;++i) data(c++) = damageTCommit(i);
    data(c++) = damageC1Commit;
    for (int i=0;i<3;++i) data(c++) = damageVCommit(i);
    for (int i=0;i<3;++i) data(c++) = FmaxTCommit(i);
    data(c++) = FmaxC1Commit;

    if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "TimberHoffman3D::sendSelf - failed to send Vector\n";
        return -1;
    }
    return 0;
}

int TimberHoffman3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    if (this->getDbTag() == 0)
        this->setDbTag(theChannel.getDbTag());

    static Vector data(69);
    if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "TimberHoffman3D::recvSelf - failed to recv Vector\n";
        return -1;
    }

    int c = 0;
    this->setTag((int)data(c++));
    E1 = data(c++); E2 = data(c++); E3 = data(c++);
    nu12 = data(c++); nu13 = data(c++); nu23 = data(c++);
    G12 = data(c++); G13 = data(c++); G23 = data(c++);
    fc1 = data(c++); fc2 = data(c++); fc3 = data(c++);
    ft1 = data(c++); ft2 = data(c++); ft3 = data(c++);
    f12 = data(c++); f13 = data(c++); f23 = data(c++);
    hardeningModulus = data(c++); sigmaE0 = data(c++);
    Acomp = data(c++); Bcomp = data(c++);
    Gf1t = data(c++); Gf2t = data(c++); Gf3t = data(c++);
    eta = data(c++); Lc = data(c++); dt = data(c++);
    tol = data(c++); maxIter = (int)data(c++); innerTol = data(c++); innerMaxIter = (int)data(c++);

    buildElasticStiffness();
    buildHoffmanPQ();
    Z.Zero();
    Z(0,0)=1.0; Z(1,1)=1.0; Z(2,2)=1.0;
    Z(3,3)=0.5; Z(4,4)=0.5; Z(5,5)=0.5;

    for (int i=0;i<6;++i) strainCommit(i) = data(c++);
    for (int i=0;i<6;++i) stressCommit(i) = data(c++);
    for (int i=0;i<6;++i) effectiveStressCommit(i) = data(c++);
    for (int i=0;i<6;++i) plasticStrainCommit(i) = data(c++);
    eqpCommit = data(c++);
    for (int i=0;i<3;++i) damageTCommit(i) = data(c++);
    damageC1Commit = data(c++);
    for (int i=0;i<3;++i) damageVCommit(i) = data(c++);
    for (int i=0;i<3;++i) FmaxTCommit(i) = data(c++);
    FmaxC1Commit = data(c++);

    strainTrial = strainCommit;
    stressTrial = stressCommit;
    effectiveStressTrial = effectiveStressCommit;
    plasticStrainTrial = plasticStrainCommit;
    eqpTrial = eqpCommit;
    damageTTrial = damageTCommit;
    damageC1Trial = damageC1Commit;
    damageVTrial = damageVCommit;
    FmaxTTrial = FmaxTCommit;
    FmaxC1Trial = FmaxC1Commit;

    buildDamagedStiffness(damageVCommit);
    return 0;
}

void TimberHoffman3D::Print(
    OPS_Stream &s,
    int flag)
{
    s << "TimberHoffman3D tag: " << this->getTag() << endln;
    s << "  E1 E2 E3: "
      << E1 << " " << E2 << " " << E3 << endln;
    s << "  sigmaE0 h: "
      << sigmaE0 << " " << hardeningModulus << endln;
    s << "  Acomp Bcomp: "
      << Acomp << " " << Bcomp << endln;
    s << "  Lc eta dt: "
      << Lc << " " << eta << " " << dt << endln;
}

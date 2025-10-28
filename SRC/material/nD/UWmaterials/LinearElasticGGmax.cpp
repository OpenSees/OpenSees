/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** LinearElasticGGmax: linear elastic nD with G/Gmax degradation      **
** ****************************************************************** */

#include "LinearElasticGGmax.h"
#include <classTags.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <cmath>
#include <algorithm>

Vector LinearElasticGGmax::sigma(6);
Matrix LinearElasticGGmax::D(6,6);

/* ----------------------- OPS Factory Function ---------------------- */
void* OPS_LinearElasticGGmaxMaterial()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 5) {
        opserr << "nDMaterial LinearElasticGGmax tag G K|nu rho curveType <params|userCurve>\n";
        return 0;
    }

    int tag; int num = 1;
    if (OPS_GetIntInput(&num, &tag) < 0) {
        opserr << "LinearElasticGGmax: invalid tag\n";
        return 0;
    }

    double props[3]; num = 3;
    if (OPS_GetDoubleInput(&num, props) < 0) {
        opserr << "LinearElasticGGmax: need G, K|nu, rho\n";
        return 0;
    }
    double G = props[0];
    double K_or_nu = props[1];
    double rho = props[2];

    int curveType; num = 1;
    if (OPS_GetIntInput(&num, &curveType) < 0) {
        opserr << "LinearElasticGGmax: invalid curveType\n";
        return 0;
    }

    if (curveType == 0) {
        int npts; num = 1;
        if (OPS_GetIntInput(&num, &npts) < 0 || npts < 2) {
            opserr << "LinearElasticGGmax: user curve needs numPoints >= 2\n";
            return 0;
        }
        std::vector<double> strains(npts), gg(npts);
        num = npts;
        if (OPS_GetDoubleInput(&num, strains.data()) < 0) {
            opserr << "LinearElasticGGmax: bad strain list\n";
            return 0;
        }
        num = npts;
        if (OPS_GetDoubleInput(&num, gg.data()) < 0) {
            opserr << "LinearElasticGGmax: bad G/Gmax list\n";
            return 0;
        }
        return new LinearElasticGGmax(tag, G, K_or_nu, rho, strains, gg);
    } else {
        double params[3] = {0.,0.,0.};
        int rem = OPS_GetNumRemainingInputArgs();
        if (rem > 0) {
            if (rem > 3) rem = 3;
            num = rem;
            if (OPS_GetDoubleInput(&num, params) < 0) {
                opserr << "LinearElasticGGmax: invalid curve params\n";
                return 0;
            }
        }
        return new LinearElasticGGmax(tag, G, K_or_nu, rho, curveType, params[0], params[1], params[2]);
    }
}

/* --------------------------- Constructors -------------------------- */
LinearElasticGGmax::LinearElasticGGmax(int tag, double G_in, double K_or_nu_in, double r,
                                       int cType, double p1, double p2, double p3)
: NDMaterial(tag, ND_TAG_LinearElasticGGmax),
  G0(G_in), K0(0.0), nu(0.0), hasK(false), rho0(r),
  mu_c(0.0), lambda_c(0.0),
  curveType(cType), param1(p1), param2(p2), param3(p3),
  epsilon(6), Cepsilon(6), nDim(3)
{
    // detect whether second is nu or K
    if (K_or_nu_in > -0.999 && K_or_nu_in < 0.5) { nu = K_or_nu_in; hasK = false; }
    else { K0 = K_or_nu_in; hasK = true; }

    epsilon.Zero(); Cepsilon.Zero();
    sigma.Zero(); D.Zero();
}

LinearElasticGGmax::LinearElasticGGmax(int tag, double G_in, double K_or_nu_in, double r,
                                       const std::vector<double>& strains,
                                       const std::vector<double>& ggmax)
: NDMaterial(tag, ND_TAG_LinearElasticGGmax),
  G0(G_in), K0(0.0), nu(0.0), hasK(false), rho0(r),
  mu_c(0.0), lambda_c(0.0),
  curveType(0), param1(0.0), param2(0.0), param3(0.0),
  userStrains(strains), userGGmax(ggmax),
  epsilon(6), Cepsilon(6), nDim(3)
{
    if (K_or_nu_in > -0.999 && K_or_nu_in < 0.5) { nu = K_or_nu_in; hasK = false; }
    else { K0 = K_or_nu_in; hasK = true; }

    epsilon.Zero(); Cepsilon.Zero();
    sigma.Zero(); D.Zero();
}

LinearElasticGGmax::LinearElasticGGmax()
: NDMaterial(0, ND_TAG_LinearElasticGGmax),
  G0(0.0), K0(0.0), nu(0.0), hasK(true), rho0(0.0),
  mu_c(0.0), lambda_c(0.0),
  curveType(1), param1(0.0), param2(0.0), param3(0.0),
  epsilon(6), Cepsilon(6), nDim(3)
{
    epsilon.Zero(); Cepsilon.Zero();
    sigma.Zero(); D.Zero();
}

LinearElasticGGmax::~LinearElasticGGmax() {}

/* ------------------------ NDMaterial methods ----------------------- */

int LinearElasticGGmax::setTrialStrain(const Vector &strn)
{
    if (epsilon.Size() != strn.Size()) {
        // detect dimensionality (3 for 2D plane strain/plane stress vectors)
        if (strn.Size() == 3) { nDim = 2; epsilon.resize(3); Cepsilon.resize(3); }
        else { nDim = 3; epsilon.resize(6); Cepsilon.resize(6); }
    }
    epsilon = strn;

    // compute G/Gmax from equivalent shear strain (uses normals too)
    const double gamma_eq = computeShearStrain(epsilon);
    const double gg = computeGGmax(gamma_eq);
    computeTangent(gg);   // sets mu_c, lambda_c and updates D

    // FAST path: sigma = 2*mu*eps + lambda*tr(eps)*I in Voigt form
    if (nDim == 3) {
        const double exx = epsilon(0), eyy = epsilon(1), ezz = epsilon(2);
        const double gxy = epsilon(3), gyz = epsilon(4), gxz = epsilon(5);
        const double tr  = exx + eyy + ezz;

        sigma(0) = lambda_c*tr + 2.0*mu_c*exx;
        sigma(1) = lambda_c*tr + 2.0*mu_c*eyy;
        sigma(2) = lambda_c*tr + 2.0*mu_c*ezz;
        sigma(3) = mu_c * gxy;   // engineering shear
        sigma(4) = mu_c * gyz;
        sigma(5) = mu_c * gxz;
    } else {
        // plane strain: epsilon = [exx, eyy, gxy], with ezz = 0.0
        const double exx = epsilon(0), eyy = epsilon(1), gxy = epsilon(2);
        const double tr  = exx + eyy;

        sigma(0) = lambda_c*tr + 2.0*mu_c*exx;
        sigma(1) = lambda_c*tr + 2.0*mu_c*eyy;
        sigma(2) = mu_c * gxy;   // tau_xy
    }

    return 0;
}

int LinearElasticGGmax::setTrialStrain(const Vector &strain, const Vector &rate)
{ return setTrialStrain(strain); }

int LinearElasticGGmax::setTrialStrainIncr(const Vector &strain)
{
    Vector trial = Cepsilon;
    trial.addVector(1.0, strain, 1.0);
    return setTrialStrain(trial);
}

int LinearElasticGGmax::setTrialStrainIncr(const Vector &strain, const Vector &rate)
{ return setTrialStrainIncr(strain); }

const Vector &LinearElasticGGmax::getStrain(void) { return epsilon; }
const Vector &LinearElasticGGmax::getStress(void) { return sigma; }
const Matrix &LinearElasticGGmax::getTangent(void) { return D; }
const Matrix &LinearElasticGGmax::getInitialTangent(void)
{
    // initial tangent: use ggmax=1
    computeTangent(1.0);
    return D;
}

int LinearElasticGGmax::commitState(void)
{
    Cepsilon = epsilon;
    return 0;
}

int LinearElasticGGmax::revertToLastCommit(void)
{
    // Restore trial strain to last committed state
    epsilon = Cepsilon;

    // Recompute tangent and stress consistent with restored strain
    const double gamma_eq = computeShearStrain(epsilon);
    const double gg = computeGGmax(gamma_eq);
    computeTangent(gg);

    if (nDim == 3) {
        const double exx = epsilon(0), eyy = epsilon(1), ezz = epsilon(2);
        const double gxy = epsilon(3), gyz = epsilon(4), gxz = epsilon(5);
        const double tr  = exx + eyy + ezz;

        sigma(0) = lambda_c*tr + 2.0*mu_c*exx;
        sigma(1) = lambda_c*tr + 2.0*mu_c*eyy;
        sigma(2) = lambda_c*tr + 2.0*mu_c*ezz;
        sigma(3) = mu_c * gxy;
        sigma(4) = mu_c * gyz;
        sigma(5) = mu_c * gxz;
    } else {
        const double exx = epsilon(0), eyy = epsilon(1), gxy = epsilon(2);
        const double tr  = exx + eyy;

        sigma(0) = lambda_c*tr + 2.0*mu_c*exx;
        sigma(1) = lambda_c*tr + 2.0*mu_c*eyy;
        sigma(2) = mu_c * gxy;
    }

    return 0;
}

int LinearElasticGGmax::revertToStart(void)
{
    epsilon.Zero(); Cepsilon.Zero(); sigma.Zero();
    computeTangent(1.0);
    return 0;
}

NDMaterial* LinearElasticGGmax::getCopy(void)
{
    double second = hasK ? K0 : nu;
    LinearElasticGGmax *theCopy;
    if (curveType == 0) {
        theCopy = new LinearElasticGGmax(this->getTag(), G0, second, rho0, userStrains, userGGmax);
    } else {
        theCopy = new LinearElasticGGmax(this->getTag(), G0, second, rho0, curveType, param1, param2, param3);
    }
    theCopy->hasK = hasK;
    theCopy->epsilon = epsilon;
    theCopy->Cepsilon = Cepsilon;
    theCopy->nDim = nDim;
    return theCopy;
}

NDMaterial* LinearElasticGGmax::getCopy(const char *type)
{
    return getCopy();
}

const char *LinearElasticGGmax::getType(void) const
{
    return (nDim==2) ? "PlaneStrain" : "ThreeDimensional";
}

int LinearElasticGGmax::getOrder(void) const
{
    return (nDim==2) ? 3 : 6;
}

/* ------------------------- Channel I/O ------------------------------ */

int LinearElasticGGmax::sendSelf(int commitTag, Channel &theChannel)
{
    static Vector data(12);
    data(0) = this->getTag();
    data(1) = G0;
    data(2) = hasK ? K0 : nu;
    data(3) = rho0;
    data(4) = curveType;
    data(5) = param1;
    data(6) = param2;
    data(7) = param3;
    data(8) = nDim;
    data(9) = (double)userStrains.size();
    data(10)= hasK ? 1.0 : 0.0;
    data(11)= 0.0;

    int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0) return res;

    if (curveType == 0 && !userStrains.empty()) {
        Vector s(userStrains.size()), g(userGGmax.size());
        for (size_t i=0;i<userStrains.size();++i){ s(i)=userStrains[i]; g(i)=userGGmax[i]; }
        if ((res = theChannel.sendVector(this->getDbTag(), commitTag, s)) < 0) return res;
        if ((res = theChannel.sendVector(this->getDbTag(), commitTag, g)) < 0) return res;
    }
    return 0;
}

int LinearElasticGGmax::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    static Vector data(12);
    int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) return res;

    this->setTag((int)data(0));
    G0 = data(1);
    double second = data(2);
    rho0 = data(3);
    curveType = (int)data(4);
    param1 = data(5); param2 = data(6); param3 = data(7);
    nDim = (int)data(8);
    int nPts = (int)data(9);
    hasK = (data(10) > 0.5);

    if (hasK) { K0 = second; nu = 0.0; }
    else      { nu = second; K0 = 0.0; }

    if (curveType == 0 && nPts > 0) {
        Vector s(nPts), g(nPts);
        if ((res = theChannel.recvVector(this->getDbTag(), commitTag, s)) < 0) return res;
        if ((res = theChannel.recvVector(this->getDbTag(), commitTag, g)) < 0) return res;
        userStrains.resize(nPts); userGGmax.resize(nPts);
        for (int i=0;i<nPts;i++){ userStrains[i]=s(i); userGGmax[i]=g(i); }
    }
    return 0;
}

/* ----------------------------- Print ------------------------------- */
void LinearElasticGGmax::Print(OPS_Stream &s, int flag)
{
    s << "LinearElasticGGmax, tag: " << this->getTag() << endln;
    s << "  G0: " << G0 << endln;
    if (hasK) s << "  K0: " << K0 << endln;
    else      s << "  nu: " << nu << endln;
    s << "  rho: " << rho0 << endln;
    s << "  curveType: " << curveType << endln;
}

/* ------------------------ Private helpers -------------------------- */

double LinearElasticGGmax::computeShearStrain(const Vector& strain) const
{
    double gamma = 0.0;

    if (nDim == 3) {
        // Voigt: [eps11, eps22, eps33, g12, g23, g31] with engineering shear
        double eps11 = strain(0);
        double eps22 = strain(1);
        double eps33 = strain(2);
        double eps12 = strain(3) * 0.5;  // eng -> tensorial
        double eps23 = strain(4) * 0.5;
        double eps31 = strain(5) * 0.5;

        double epsm = (eps11 + eps22 + eps33) / 3.0;
        double e11 = eps11 - epsm;
        double e22 = eps22 - epsm;
        double e33 = eps33 - epsm;

        double J2 = 0.5 * (e11*e11 + e22*e22 + e33*e33)
                    + (eps12*eps12 + eps23*eps23 + eps31*eps31);

        gamma = std::sqrt(2.0 * J2);
    } else {
        // 2D (plane strain/stress vector: [eps11, eps22, g12])
        double eps11 = strain(0);
        double eps22 = strain(1);
        double eps12 = strain(2) * 0.5;  // eng -> tensorial

        // For pure 2D we use mean over in-plane normals
        double epsm = (eps11 + eps22) / 2.0;
        double e11 = eps11 - epsm;
        double e22 = eps22 - epsm;

        double J2 = 0.5 * (e11*e11 + e22*e22) + eps12*eps12;
        gamma = std::sqrt(2.0 * J2);
    }

    return std::fabs(gamma);
}

double LinearElasticGGmax::computeGGmax(double gamma)
{
    gamma = std::abs(gamma);
    if (curveType == 0) return interpolateUserCurve(gamma);
    if (curveType == 1) return hardinDrnevich(gamma);
    if (curveType == 2) return vuceticDobry(gamma);
    if (curveType == 3) return darendeli(gamma);
    return 1.0;
}

void LinearElasticGGmax::computeTangent(double gg)
{
    // Degraded shear modulus
    double mu = std::max(0.0, G0 * gg);
    double lambda;
    if (hasK) {
        double K = K0; // keep bulk modulus constant
        lambda = K - (2.0/3.0)*mu;
    } else {
        lambda = (2.0*mu*nu) / (1.0 - 2.0*nu);
    }

    // Cache values for fast stress computation
    mu_c = mu;
    lambda_c = lambda;

    // Fill D for getTangent/getInitialTangent
    if (nDim == 3) {
        D.Zero();
        double mup2 = lambda + 2.0*mu;
        D(0,0)=D(1,1)=D(2,2)=mup2;
        D(0,1)=D(1,0)=lambda;
        D(0,2)=D(2,0)=lambda;
        D(1,2)=D(2,1)=lambda;
        D(3,3)=mu; D(4,4)=mu; D(5,5)=mu;
    } else {
        // plane strain: 3x3 [exx, eyy, gxy]
        D.Zero();
        double mup2 = lambda + 2.0*mu;
        D(0,0)=D(1,1)=mup2;
        D(0,1)=D(1,0)=lambda;
        D(2,2)=mu;
    }
}

/* -------------------------- Curve models --------------------------- */

double LinearElasticGGmax::hardinDrnevich(double gamma) const
{
    // gamma_ref in param1
    double gref = (param1>0.0? param1 : 1e-4);
    return 1.0 / (1.0 + gamma / gref);
}

double LinearElasticGGmax::vuceticDobry(double gamma) const
{
    // Simplified: PI in param1 -> adjust reference strain
    double PI = std::max(0.0, param1);
    double gref = 1e-4 * std::pow(10.0, -0.014*PI); // crude trend
    return 1.0 / (1.0 + gamma / gref);
}

double LinearElasticGGmax::darendeli(double gamma) const
{
    // Very simplified surrogate based on param1=PI, param2=p', param3=OCR
    double PI  = std::max(0.0, param1);
    double p   = (param2>0.0? param2:100.0);
    double OCR = (param3>0.0? param3:1.0);
    double gref = 1e-4 * std::pow(p/100.0, 0.3) * std::pow(OCR, 0.1) * std::pow(10.0, -0.01*PI);
    return 1.0 / (1.0 + gamma / gref);
}

double LinearElasticGGmax::interpolateUserCurve(double gamma) const
{
    if (userStrains.empty() || userGGmax.empty()) return 1.0;
    // assume userStrains strictly increasing (positive gammas)
    size_t n = userStrains.size();
    if (gamma <= userStrains.front()) return userGGmax.front();
    if (gamma >= userStrains.back())  return userGGmax.back();

    // log-log linear interpolation
    double lg = std::log10(gamma);
    size_t i=1;
    while (i<n && gamma > userStrains[i]) ++i;
    // interpolate between i-1 and i
    double x1=std::log10(userStrains[i-1]), x2=std::log10(userStrains[i]);
    double y1=std::log10(std::max(1e-12, userGGmax[i-1])), y2=std::log10(std::max(1e-12, userGGmax[i]));
    double y = y1 + (y2-y1)*( (lg-x1)/(x2-x1) );
    return std::pow(10.0, y);
}

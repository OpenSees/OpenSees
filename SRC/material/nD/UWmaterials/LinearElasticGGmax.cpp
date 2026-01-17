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
#include <Parameter.h>
#include <Information.h>
#include <cmath>
#include <algorithm>
#include <cstring>   // for strcmp

// Instance members are initialized in constructors; remove static shared state.

/* OPS Factory Function */
void* OPS_LinearElasticGGmaxMaterial()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 5) {
        opserr << "nDMaterial LinearElasticGGmax tag G K|nu rho curveType [chi] <params|userCurve>\n";
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
    // chi is optional; if present, read it, else default to 0.0
    double chi = 0.0;
    int remAfterCurveType = OPS_GetNumRemainingInputArgs();
    if (remAfterCurveType > 0) {
        // Peek next arg by trying to read one double; if it fails and curveType==0, we assume user provided pairs next
        // We can't truly peek; instead, we will branch based on curveType cases below.
    }

    if (curveType == 0) {
        // Optional chi before pairs: detect whether first remaining token is chi or part of pairs
        int rem = OPS_GetNumRemainingInputArgs();
        if (rem < 4) {
            opserr << "LinearElasticGGmax: curveType=0 expects interleaved (gamma, G/Gmax) pairs (>= 2 pairs)" << endln;
            return nullptr;
        }

        // Heuristic: if rem is odd, treat first as chi then pairs; if even, no chi and all are pairs
        bool hasChiInline = (rem % 2) == 1;
        if (hasChiInline) {
            int nread = 1;
            if (OPS_GetDoubleInput(&nread, &chi) < 0) {
                opserr << "LinearElasticGGmax: invalid chi for user curve" << endln;
                return nullptr;
            }
            rem = OPS_GetNumRemainingInputArgs();
        }

        if (rem < 4 || (rem % 2) != 0) {
            opserr << "LinearElasticGGmax: curveType=0 expects interleaved (gamma, G/Gmax) pairs (>= 2 pairs)" << endln;
            return nullptr;
        }

        Vector flat(rem);
        if (OPS_GetDoubleInput(&rem, &flat(0)) < 0) {
            opserr << "LinearElasticGGmax: unable to read (gamma, G/Gmax) pairs" << endln;
            return nullptr;
        }

        const int nPts = rem / 2;
        std::vector<double> strains(nPts), gg(nPts);
        for (int i = 0; i < nPts; ++i) {
            strains[i] = flat(2*i);
            gg[i]      = flat(2*i + 1);
        }

        // Basic validation/clamping
        for (int i = 0; i < nPts; ++i) {
            if (strains[i] <= 0.0) strains[i] = (i == 0 ? 1e-12 : strains[i-1] * 1.000000000001);
            if (i > 0 && strains[i] <= strains[i-1]) {
                strains[i] = strains[i-1] * 1.000000000001; // enforce strictly increasing
            }
            if (gg[i] <= 0.0) gg[i] = 1e-6;
            if (gg[i] > 1.0)  gg[i] = 1.0;
        }
        return new LinearElasticGGmax(tag, G, K_or_nu, rho, strains, gg, chi);
    } else {
        // Predefined curve types: read optional chi then up to three optional parameters (p1, p2, p3)
        int rem = OPS_GetNumRemainingInputArgs();
        // If there are 1..4 remaining args, first may be chi, followed by up to 3 params
        double params[3] = {0.0, 0.0, 0.0};
        if (rem > 0) {
            // Try reading chi first
            int nread = 1;
            double chiCandidate = 0.0;
            if (OPS_GetDoubleInput(&nread, &chiCandidate) == 0) {
                chi = chiCandidate;
                rem = OPS_GetNumRemainingInputArgs();
            }
        }
        int toRead = std::min(rem, 3);
        if (toRead > 0) {
            if (OPS_GetDoubleInput(&toRead, params) < 0) {
                opserr << "LinearElasticGGmax: unable to read curve parameters (p1 p2 p3)" << endln;
                return nullptr;
            }
        }
        return new LinearElasticGGmax(tag, G, K_or_nu, rho, curveType, params[0], params[1], params[2], chi);
    }

    return nullptr;
}

/* --------------------------- Constructors -------------------------- */
LinearElasticGGmax::LinearElasticGGmax(int tag, double G_in, double K_or_nu_in, double r,
                                                                             int cType, double p1, double p2, double p3,
                                                                             double chi_in)
: NDMaterial(tag, ND_TAG_LinearElasticGGmax),
    G0(G_in), K0(0.0), nu(0.0), hasK(false), rho0(r), chi(chi_in),
        mu_c(0.0), lambda_c(0.0),
    curveType(cType), param1(p1), param2(p2), param3(p3),
    epsilon(6), Cepsilon(6), nDim(3)
{
    // detect whether second is nu or K
    if (K_or_nu_in > -0.999 && K_or_nu_in < 0.5) { nu = K_or_nu_in; hasK = false; }
    else { K0 = K_or_nu_in; hasK = true; }

    epsilon.Zero(); Cepsilon.Zero();
    sigma.resize(6); sigma_vis.resize(6); D.resize(6,6); D0.resize(6,6); Ddamp.resize(6,6); Cep.resize(6,6);
    sigma.Zero(); sigma_vis.Zero(); D.Zero(); D0.Zero(); Ddamp.Zero();
    Cep.Zero();

    gammaMaxTrial  = 0.0;
    gammaMaxCommit = 0.0;
    // G update policy defaults (preserve legacy behavior)
    updateStride = 1;
    stepCounter = 0;
    updateOnDemand = false;
    pendingOneShotUpdate = false;
    lastGG = 1.0;
    strideByStep = false;
    commitCounter = 0;
    lastUpdateCommit = -1;

    // Initialize initial and damping tangents so Ddamp is available immediately
    computeInitialElasticTangent();
    computeDampingTangent();
}

LinearElasticGGmax::LinearElasticGGmax(int tag, double G_in, double K_or_nu_in, double r,
                                                                             const std::vector<double>& strains,
                                                                             const std::vector<double>& ggmax,
                                                                             double chi_in)
: NDMaterial(tag, ND_TAG_LinearElasticGGmax),
    G0(G_in), K0(0.0), nu(0.0), hasK(false), rho0(r), chi(chi_in),
        mu_c(0.0), lambda_c(0.0),
    curveType(0), param1(0.0), param2(0.0), param3(0.0),
    userStrains(strains), userGGmax(ggmax),
    epsilon(6), Cepsilon(6), nDim(3)
{
    if (K_or_nu_in > -0.999 && K_or_nu_in < 0.5) { nu = K_or_nu_in; hasK = false; }
    else { K0 = K_or_nu_in; hasK = true; }

    epsilon.Zero(); Cepsilon.Zero();
    sigma.resize(6); sigma_vis.resize(6); D.resize(6,6); D0.resize(6,6); Ddamp.resize(6,6); Cep.resize(6,6);
    sigma.Zero(); sigma_vis.Zero(); D.Zero(); D0.Zero(); Ddamp.Zero();
    Cep.Zero();

    gammaMaxTrial  = 0.0;
    gammaMaxCommit = 0.0;
    updateStride = 1;
    stepCounter = 0;
    updateOnDemand = false;
    pendingOneShotUpdate = false;
    lastGG = 1.0;
    strideByStep = false;
    commitCounter = 0;
    lastUpdateCommit = -1;   

    // Initialize initial and damping tangents so Ddamp is available immediately
    computeInitialElasticTangent();
    computeDampingTangent();
}

LinearElasticGGmax::LinearElasticGGmax()
: NDMaterial(0, ND_TAG_LinearElasticGGmax),
    G0(0.0), K0(0.0), nu(0.0), hasK(true), rho0(0.0), chi(0.0),
    mu_c(0.0), lambda_c(0.0),
  curveType(1), param1(0.0), param2(0.0), param3(0.0),
  epsilon(6), Cepsilon(6), nDim(3)
{
    epsilon.Zero(); Cepsilon.Zero();
    sigma.resize(6); sigma_vis.resize(6); D.resize(6,6); D0.resize(6,6); Ddamp.resize(6,6); Cep.resize(6,6);
    sigma.Zero(); sigma_vis.Zero(); D.Zero(); D0.Zero(); Ddamp.Zero();
    Cep.Zero();

    gammaMaxTrial = 0.0;
    gammaMaxCommit = 0.0;
    updateStride = 1;
    stepCounter = 0;
    updateOnDemand = false;
    pendingOneShotUpdate = false;
    lastGG = 1.0;
    strideByStep = false;
    commitCounter = 0;
    lastUpdateCommit = -1;    

    // Initialize initial and damping tangents so Ddamp is available immediately
    computeInitialElasticTangent();
    computeDampingTangent();
}

LinearElasticGGmax::~LinearElasticGGmax() {}
// double LinearElasticGGmax::getRho(void) { return rho0; }

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
    const double gamma_used = computeShearStrain(epsilon);
    const double gg = computeGGmax(gamma_used);

    // Decide whether to update tangent (and thus mu_c)
    bool doUpdate = false;
    if (updateOnDemand) {
        if (pendingOneShotUpdate) {
            doUpdate = true;
            pendingOneShotUpdate = false;
            lastUpdateCommit = commitCounter;
        }
    } else if (strideByStep) {
        // At most one update per converged time step, every N steps
        const bool atNewStep = (commitCounter != lastUpdateCommit);
        if (atNewStep) {
            if (updateStride <= 1 || (commitCounter > 0 && (commitCounter % updateStride) == 0)) {
                doUpdate = true;
                lastUpdateCommit = commitCounter;
            }
        }
    } else {
        // Count raw setTrialStrain calls
        ++stepCounter;
        int effectiveCounter = (stepCounter > 0) ? (stepCounter - 1) : 0;
        if (updateStride <= 1 || (effectiveCounter % updateStride) == 0) {
            doUpdate = true;
        }
    }

    // Display step index for debug (time step or trial count)
    const int   displayStep = strideByStep ? (commitCounter + 1)
                                           : ((stepCounter > 0) ? (stepCounter - 1) : 0);
    const char* label       = strideByStep ? "timeStep " : "trial ";
    const int   stepIndexEst = commitCounter + 1;

    if (doUpdate) {
        if (debugUpdate) {
            opserr << "LinearElasticGGmax(tag=" << this->getTag()
                   << "): UPDATE at " << label << displayStep
                   << " (gg=" << gg << ", gamma_used=" << gamma_used << ")"
                   << " [step=" << stepIndexEst << "]" << endln;
        }
        computeTangent(gg);
        lastGG = gg;
    } else if (debugUpdate) {
        opserr << "LinearElasticGGmax(tag=" << this->getTag()
           << "): SKIP   at " << label << displayStep
           << " (using lastGG=" << lastGG << ")"
           << " [step=" << stepIndexEst << "]" << endln;
    }

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

    // Add viscous stress contribution by approximating strain rate as dStrain/ops_Dt (consistent with J2CyclicBoundingSurface)
    sigma_vis.Zero();
    extern double ops_Dt;

    if (ops_Dt > 0.0) {
        if (nDim == 3) {
            Vector dE(6);
            for (int i=0;i<6;++i) dE(i) = epsilon(i) - Cepsilon(i);
            for (int i=0;i<6;++i) {
                double acc = 0.0;
                for (int j=0;j<6;++j) acc += Ddamp(i,j) * dE(j);
                sigma_vis(i) = acc / ops_Dt;
            }
            // opserr << "Ddamp " << Ddamp << endln;

        } else {
            Vector dE(3);
            for (int i=0;i<3;++i) dE(i) = epsilon(i) - Cepsilon(i);
            for (int i=0;i<3;++i) {
                double acc = 0.0;
                for (int j=0;j<3;++j) acc += Ddamp(i,j) * dE(j);
                sigma_vis(i) = acc / ops_Dt;
            }
        }
        sigma.addVector(1.0, sigma_vis, 1.0);
    }

    return 0;
}

int LinearElasticGGmax::setTrialStrain(const Vector &strain, const Vector &rate)
{
    // Elements typically do not provide rate; ignore and rely on dStrain/ops_Dt approximation
    return setTrialStrain(strain);
}

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
const Matrix &LinearElasticGGmax::getTangent(void) {
    // Return elastic tangent plus viscous contribution if time step is positive
    extern double ops_Dt; // from OPS_Globals
    Cep.Zero();
    if (ops_Dt > 0.0 && chi != 0.0) {
        if (nDim == 3) {
            for (int i=0;i<6;++i)
                for (int j=0;j<6;++j)
                    Cep(i,j) = D(i,j) + (1.0/ops_Dt) * Ddamp(i,j);
        } else {
            // plane strain: use top-left 3x3 block
            for (int i=0;i<3;++i)
                for (int j=0;j<3;++j)
                    Cep(i,j) = D(i,j) + (1.0/ops_Dt) * Ddamp(i,j);
        }
        return Cep;
    }
    // No damping term; return elastic tangent in Cep for consistency
    if (nDim == 3) {
        for (int i=0;i<6;++i)
            for (int j=0;j<6;++j)
                Cep(i,j) = D(i,j);
    } else {
        for (int i=0;i<3;++i)
            for (int j=0;j<3;++j)
                Cep(i,j) = D(i,j);
    }
    return Cep;
}
const Matrix &LinearElasticGGmax::getInitialTangent(void)
{
    // Return initial elastic tangent D0 (no degradation)
    computeInitialElasticTangent();
    return D0;
}

const Matrix &LinearElasticGGmax::getDampTangent(void)
{
    // Ensure damping tangent is up to date with current elastic tangent
    computeDampingTangent();
    return Ddamp;
}

int LinearElasticGGmax::commitState(void)
{
    Cepsilon = epsilon;
    gammaMaxCommit = gammaMaxTrial;
    // Count converged time steps
    ++commitCounter;    
    return 0;
}

int LinearElasticGGmax::revertToLastCommit(void)
{
    // Restore trial strain to last committed state
    epsilon = Cepsilon;
    gammaMaxTrial = gammaMaxCommit;

    // Recompute tangent and stress consistent with restored strain
    const double gamma_used = computeShearStrain(epsilon);
    const double gg = computeGGmax(gamma_used);
    computeTangent(gg);
    lastGG = gg;

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
    epsilon.Zero(); Cepsilon.Zero(); sigma.Zero(); sigma_vis.Zero();
    gammaMaxTrial = 0.0;
    gammaMaxCommit = 0.0;
    computeTangent(1.0);
    computeDampingTangent();
    lastGG = 1.0;
    return 0;
}

NDMaterial* LinearElasticGGmax::getCopy(void)
{
    double second = hasK ? K0 : nu;
    LinearElasticGGmax *theCopy;
    if (curveType == 0) {
        theCopy = new LinearElasticGGmax(this->getTag(), G0, second, rho0, userStrains, userGGmax, chi);
    } else {
        theCopy = new LinearElasticGGmax(this->getTag(), G0, second, rho0, curveType, param1, param2, param3, chi);
    }
    theCopy->hasK = hasK;
    theCopy->epsilon = epsilon;
    theCopy->Cepsilon = Cepsilon;
    theCopy->nDim = nDim;

    theCopy->gammaMaxTrial  = gammaMaxTrial;
    theCopy->gammaMaxCommit = gammaMaxCommit;

    // copy update policy
    theCopy->updateStride = updateStride;
    theCopy->stepCounter = 0; // fresh counter for the copy
    theCopy->updateOnDemand = updateOnDemand;
    theCopy->pendingOneShotUpdate = false;
    theCopy->lastGG = lastGG;   
    theCopy->debugUpdate = debugUpdate;    
    theCopy->dampUsesInitial = dampUsesInitial;

    // ensure damping tangent is consistent on the copy
    theCopy->computeInitialElasticTangent();
    theCopy->computeDampingTangent();

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
    Vector data(15);
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
    data(11)= gammaMaxCommit;
    data(12)= chi; // persist damping coefficient
    data(13)= gammaScale; // optional: persist eq-linear scale
    data(14)= dampUsesInitial ? 1.0 : 0.0; // persist damping mode

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
    Vector data(15);
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
    gammaMaxCommit = data(11);
    chi = data(12);
    gammaScale = data(13);
    dampUsesInitial = (data(14) > 0.5);
    gammaMaxTrial  = gammaMaxCommit;

    if (hasK) { K0 = second; nu = 0.0; }
    else      { nu = second; K0 = 0.0; }

    if (curveType == 0 && nPts > 0) {
        Vector s(nPts), g(nPts);
        if ((res = theChannel.recvVector(this->getDbTag(), commitTag, s)) < 0) return res;
        if ((res = theChannel.recvVector(this->getDbTag(), commitTag, g)) < 0) return res;
        userStrains.resize(nPts); userGGmax.resize(nPts);
        for (int i=0;i<nPts;i++){ userStrains[i]=s(i); userGGmax[i]=g(i); }
    }

    // Recompute tangents on restored state
    computeInitialElasticTangent();
    computeTangent(1.0);
    computeDampingTangent();
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

double LinearElasticGGmax::computeShearStrain(const Vector& strain)
{
    double gamma = 0.0;

    // Octahedral-equivalent engineering shear gamma_oct = 2*sqrt(2/3 * e_dev:e_dev)
    if (nDim == 3) {
        double exx = strain(0), eyy = strain(1), ezz = strain(2);
        double gxy = strain(3), gyz = strain(4), gxz = strain(5); // engineering shear
        // Convert to tensorial shear halves for octahedral measure
        double exy = 0.5*gxy, eyz = 0.5*gyz, exz = 0.5*gxz;
        double em  = (exx+eyy+ezz)/3.0;
        double e11=exx-em, e22=eyy-em, e33=ezz-em;
        // gamma_oct (engineering) = 2*sqrt(2/3*(e11^2+e22^2+e33^2+2(exy^2+eyz^2+exz^2)))
        double val = 2.0*std::sqrt( (2.0/3.0)*( e11*e11 + e22*e22 + e33*e33 +
                                               2.0*(exy*exy + eyz*eyz + exz*exz) ) );

	gamma = val;
    } else {
        // 2D (plane) vector is [exx, eyy, gxy]
        double exx = strain(0), eyy = strain(1), gxy = strain(2);
        double exy = 0.5*gxy;
        double em  = (exx+eyy)/2.0; // 2D mean
        double e11 = exx-em, e22=eyy-em;
        double val = 2.0*std::sqrt( (2.0/3.0)*( e11*e11 + e22*e22 + 2.0*(exy*exy) ) );
        gamma = val;
    }

    gamma = std::fabs(gamma);
    // Apply equivalent-linear scaling factor
    gamma *= gammaScale;

    // --- keep a running maximum over the analysis history (trial state) ---
    if (gamma > gammaMaxTrial) {
        gammaMaxTrial = gamma;
    }

    return gammaMaxTrial;
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

    // Update damping tangent consistently (proportional to elastic tangent)
    computeDampingTangent();
}

void LinearElasticGGmax::computeDampingTangent()
{
    // Damping tangent: either proportional to CURRENT elastic tangent (D)
    // or to the INITIAL elastic tangent (D0) depending on dampUsesInitial flag.
    Ddamp.Zero();
    if (dampUsesInitial) {
        if (nDim == 3) {
            for (int i=0;i<6;++i)
                for (int j=0;j<6;++j)
                    Ddamp(i,j) = chi * D0(i,j);
        } else {
            for (int i=0;i<3;++i)
                for (int j=0;j<3;++j)
                    Ddamp(i,j) = chi * D0(i,j);
        }
    } else {
        // proportional to current tangent (default)
        if (nDim == 3) {
            for (int i=0;i<6;++i)
                for (int j=0;j<6;++j)
                    Ddamp(i,j) = chi * D(i,j);
        } else {
            for (int i=0;i<3;++i)
                for (int j=0;j<3;++j)
                    Ddamp(i,j) = chi * D(i,j);
        }
    }
}

void LinearElasticGGmax::computeInitialElasticTangent()
{
    // Build D0 using initial moduli (gg=1)
    double mu = std::max(0.0, G0);
    double lambda;
    if (hasK) {
        double K = K0;
        lambda = K - (2.0/3.0)*mu;
    } else {
        lambda = (2.0*mu*nu) / (1.0 - 2.0*nu);
    }

    if (nDim == 3) {
        D0.Zero();
        double mup2 = lambda + 2.0*mu;
        D0(0,0)=D0(1,1)=D0(2,2)=mup2;
        D0(0,1)=D0(1,0)=lambda;
        D0(0,2)=D0(2,0)=lambda;
        D0(1,2)=D0(2,1)=lambda;
        D0(3,3)=mu; D0(4,4)=mu; D0(5,5)=mu;
    } else {
        D0.Zero();
        double mup2 = lambda + 2.0*mu;
        D0(0,0)=D0(1,1)=mup2;
        D0(0,1)=D0(1,0)=lambda;
        D0(2,2)=mu;
    }
}

/* --------------------- Parameter interface  ------------------- */
// Supported names: "updateStride", "updateOnDemand", "requestUpdate"
int LinearElasticGGmax::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1 || argv == nullptr) return -1;

    // --- Passthrough router for element->material forwarding ---
    // Accept forms: -ndMaterial <tag> <name>, ndMaterial <tag> <name>,
    //               -material <tag> <name>, material <tag> <name>
    if ((strcmp(argv[0], "-ndMaterial") == 0 || strcmp(argv[0], "ndMaterial") == 0 ||
         strcmp(argv[0], "-material")    == 0 || strcmp(argv[0], "material")    == 0)) {
        if (argc >= 3) {
            int matTag = atoi(argv[1]);
            if (matTag == this->getTag()) {
                return this->setParameter(&argv[2], argc - 2, param);
            } else {
                return -1; // not our tag
            }
        } else {
            return -1; // not enough tokens
        }
    }

    // Supported names
    if (strcmp(argv[0], "updateStride") == 0) {
        return param.addObject(1, this);
    } else if (strcmp(argv[0], "updateOnDemand") == 0) {
        return param.addObject(2, this);
    } else if (strcmp(argv[0], "requestUpdate") == 0) {
        return param.addObject(3, this);
    } else if (strcmp(argv[0], "debugUpdate") == 0) {
        return param.addObject(4, this);
    } else if (strcmp(argv[0], "strideByStep") == 0 || strcmp(argv[0], "updateStrideByStep") == 0) {
        return param.addObject(5, this);
    } else if (strcmp(argv[0], "gammaScale") == 0 || strcmp(argv[0], "equivLinearFactor") == 0) {
        return param.addObject(6, this);
    } else if (strcmp(argv[0], "chi") == 0 || strcmp(argv[0], "dampingCoeff") == 0) {
        return param.addObject(7, this);
    } else if (strcmp(argv[0], "dampUseInitial") == 0 || strcmp(argv[0], "dampingUseInitial") == 0) {
        return param.addObject(8, this);
    }
    return -1;
}

int LinearElasticGGmax::updateParameter(int paramID, Information &info)
{
    switch (paramID) {
    case 1: {
     // cast to int safely; treat <1 as 1
        int n = static_cast<int>(std::lround(info.theDouble));
        updateStride = (n < 1) ? 1 : n;
        opserr << "LinearElasticGGmax(tag=" << this->getTag()
               << "): updateStride=" << updateStride << endln;
        return 0;
    }
    case 2: {
        int flag = static_cast<int>(std::lround(info.theDouble));
        updateOnDemand = (flag != 0);
        opserr << "LinearElasticGGmax(tag=" << this->getTag()
               << "): updateOnDemand=" << updateOnDemand << endln;        
        return 0;
    }
    case 3: {
     // arm a one-shot update on next trial step
        pendingOneShotUpdate = true;
        opserr << "LinearElasticGGmax(tag=" << this->getTag()
               << "): requestUpdate armed" << endln;        
        return 0;
    }
    case 4: {
        int flag = static_cast<int>(std::lround(info.theDouble));
        debugUpdate = (flag != 0);
        opserr << "LinearElasticGGmax(tag=" << this->getTag()
               << "): debugUpdate=" << debugUpdate << endln;
        return 0;
    }
    case 5: {
        int flag = static_cast<int>(std::lround(info.theDouble));
        strideByStep = (flag != 0);
        opserr << "LinearElasticGGmax(tag=" << this->getTag()
               << "): strideByStep=" << strideByStep << endln;
        return 0;    
    }
    case 6: {
        gammaScale = info.theDouble;
        if (gammaScale <= 0.0) gammaScale = 1.0; // guard
        opserr << "LinearElasticGGmax(tag=" << this->getTag()
               << "): gammaScale=" << gammaScale << endln;
        return 0;
    }
    case 7: {
        // viscous damping coefficient chi
        chi = info.theDouble;
        computeDampingTangent();
        opserr << "LinearElasticGGmax(tag=" << this->getTag() << "): chi=" << chi << endln;
        return 0;
    }
    case 8: {
        int flag = static_cast<int>(std::lround(info.theDouble));
        dampUsesInitial = (flag != 0);
        computeDampingTangent();
        opserr << "LinearElasticGGmax(tag=" << this->getTag() << "): dampUseInitial=" << dampUsesInitial << endln;
        return 0;
    }
    default:
        return -1;
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
    // Inputs (same semantics as your code):
    //   param1 = PI (Plasticity Index, unitless)
    //   param2 = p' (mean effective stress in kPa)
    //   param3 = OCR (Overconsolidation Ratio, unitless)
    double PI  = std::max(0.0, param1);
    double p   = (param2 > 0.0 ? param2 : 100.0);  // kPa
    double OCR = (param3 > 0.0 ? param3 : 1.0);

    // Reference strain gref:
    // Slightly stronger dependence than the simple surrogate to match Darendeli/Menq tables.
    // Tune exponents minimally if you observe consistent offsets in your GGmax_org.
    // Typical range: gref ~ 1e-6 to 1e-3 depending on PI, p', OCR.
    double gref = 1.0e-4
                  * std::pow(p / 100.0, 0.35)   // confining effect
                  * std::pow(OCR,        0.12)  // overconsolidation effect
                  * std::pow(10.0,      -0.012 * PI); // plasticity effect

    // Shape exponent beta:
    // Allows matching curvature across mid/high strains (GG ≈ 0.2–0.8).
    // Keep beta in a reasonable band to avoid pathological slopes.
    double beta = 1.0
                  + 0.006 * PI                     // more plastic → steeper drop
                  + 0.04  * (OCR - 1.0)            // OCR increases steepness
                  + 0.02  * std::log10(std::max(1e-6, p / 100.0)); // weak confining effect
    if (beta < 0.8) beta = 0.8;
    if (beta > 1.6) beta = 1.6;

    // Compute G/Gmax using the two-parameter law.
    double g = std::fabs(gamma);
    double denom = 1.0 + std::pow(g / std::max(1e-12, gref), beta);
    double GG = 1.0 / denom;

    // Clamp to (0, 1] for safety
    if (GG > 1.0) GG = 1.0;
    if (GG < 1.0e-12) GG = 1.0e-12;

    return GG;
}


double LinearElasticGGmax::interpolateUserCurve(double gamma) const
{
    if (userStrains.empty() || userGGmax.empty()) return 1.0;

    // Assume userStrains strictly increasing (positive gammas)
    const size_t n = userStrains.size();
    const double g = std::abs(gamma);

    if (g <= userStrains.front())
        return userGGmax.front();

    // --- Beyond last point: hold shear stress constant at tau_max ---
    if (g >= userStrains.back()) {
        const double gammaLast = userStrains.back();
        const double ggLast    = std::max(1e-12, userGGmax.back());
        const double G0pos     = std::max(1e-12, G0);

        // tau_max = G0 * gg_last * gamma_last
        const double tauMax = G0pos * ggLast * gammaLast;

        // gg(gamma) = tau_max / (G0 * gamma)
        double gg = tauMax / (G0pos * g);

        // clamp to (0, 1]
        if (gg > 1.0) gg = 1.0;
        if (gg < 1e-12) gg = 1e-12;

        return gg;
    }

    // --- Log-log linear interpolation inside the user-defined range ---
    const double lg = std::log10(g);
    size_t i = 1;
    while (i < n && g > userStrains[i]) ++i;

    // interpolate between i-1 and i
    const double x1 = std::log10(userStrains[i-1]);
    const double x2 = std::log10(userStrains[i]);
    const double y1 = std::log10(std::max(1e-12, userGGmax[i-1]));
    const double y2 = std::log10(std::max(1e-12, userGGmax[i]));
    const double y  = y1 + (y2 - y1) * ((lg - x1) / (x2 - x1));

    double gg = std::pow(10.0, y);
    if (gg > 1.0) gg = 1.0;
    if (gg < 1e-12) gg = 1e-12;
    return gg;
}




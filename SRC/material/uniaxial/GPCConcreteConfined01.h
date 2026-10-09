// GPCConcreteConfined01.h
//
// Confined geopolymer concrete (GPC) uniaxial material based on:
//   Kocaer O, Aldemir A. "Confined compressive stress-strain model
//   for rectangular geopolymer reinforced concrete members."
//   Structural Concrete. 2024. https://doi.org/10.1002/suco.202300973
//
// Envelope (internal convention: compression positive), same cubic
// shape as the unconfined model with fcc/eccu in place of fck/ecu:
//   sigma(eps) = 6.75*fcc*eps*(eps-eccu)^2/eccu^3,  0 <= eps <= eccu
//
// Confinement (modified Kent & Park 1971):
//   K     = 1 + rho_s*fyh/fck
//   e50u  = 2*ecu/3
//   ls    = 2*(bk+hk)
//   rho_s = diaLat^2*(pi/4)*ls/(s*bk*hk)
//   e50h  = 2*rho_s*sqrt(bk/s)/3
//   eccu  = 4.8*(e50u+e50h)/3.6
//   fcc   = Mfck*K*fck     (Mfck = 0.7 calibrated in the source paper)
//
// Units: any consistent set (fck, fyh same stress unit; diaLat, bk, hk, s
// same length unit; ecu dimensionless).
//
// LIMITATIONS: rectangular sections only;
// pure compression model (tensile stress always zero)

#ifndef GPCConcreteConfined01_h
#define GPCConcreteConfined01_h

#include <UniaxialMaterial.h>

class GPCConcreteConfined01 : public UniaxialMaterial
{
public:
    GPCConcreteConfined01(int tag,
                          double fck, double ecu,
                          double fyh, double diaLat,
                          double bk, double hk, double s,
                          double Mfck = 0.7);
    GPCConcreteConfined01();
    ~GPCConcreteConfined01();

    const char *getClassType(void) const { return "GPCConcreteConfined01"; };

    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    UniaxialMaterial *getCopy(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &os, int flag = 0);

protected:
    // inputs
    double fck, ecu, fyh, diaLat, bk, hk, s, Mfck;

    // derived confined envelope parameters
    double fcc, eccu, E0;

    // committed state (ceps_max: largest compressive strain reached)
    double ceps, csig, ctang;
    double ceps_max;

    // trial state
    double teps, tsig, ttang;
    double teps_max;

    static const double TENS_STIFF_RATIO;

    void computeConfinedParams(void);
    void envelope(double eps, double &sig, double &tang);
};

#endif
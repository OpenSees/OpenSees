// GPCConcrete01.h
//
// Uniaxial compressive stress-strain material for geopolymer concrete (GPC),
// with optional confinement. Based on:
//
//   Unconfined model:
//     Kocaer O, Aldemir A. "Compressive stress-strain model for the
//     estimation of the flexural capacity of reinforced geopolymer
//     concrete members." Structural Concrete. 2023;24(4):5102-5121.
//     https://doi.org/10.1002/suco.202200914
//
//   Confined model (rectangular sections):
//     Kocaer O, Aldemir A. "Confined compressive stress-strain model
//     for rectangular geopolymer reinforced concrete members."
//     Structural Concrete. 2024. https://doi.org/10.1002/suco.202300973
//
// Envelope (internal convention: compression positive):
//   sigma(eps) = 6.75*fcc*eps*(eps-eccu)^2/eccu^3,   0 <= eps <= eccu
//   sigma(eps) = 0 otherwise (no tension, no strength beyond eccu)
//
// Unconfined case: fcc = fck, eccu = ecu.
// Confined case (modified Kent & Park 1971), computed once at construction:
//   K     = 1 + rho_s*fyh/fck
//   e50u  = 2*ecu/3
//   ls    = 2*(bk+hk)
//   rho_s = diaLat^2*(pi/4)*ls/(s*bk*hk)
//   e50h  = 2*rho_s*sqrt(bk/s)/3
//   eccu  = 4.8*(e50u+e50h)/3.6
//   fcc   = Mfck*K*fck       (Mfck = 0.7 calibrated in the source paper)
//
// After fcc/eccu are known, the stress-strain rules are identical for the
// confined and unconfined cases.
//
// Units: any consistent set (fck, fyh in the same stress unit; diaLat, bk,
// hk, s in the same length unit; ecu dimensionless).
//
// bk, hk: width and height of the rectangular confined core, measured
// center-to-center of the transverse reinforcement.
//
// Interface convention: OpenSees sign convention (compression negative).
// fck and ecu may be given as positive or negative values (magnitudes used).
//
// Cyclic rule (NOT given in the source papers; added assumption):
//   Loading follows the envelope. On reversal, unloading follows a straight
//   line of slope E0 (initial tangent) anchored at the reversal point,
//   clipped at zero stress. Reloading retraces the same line until the
//   reversal strain is reached, then rejoins the envelope.
//
// LIMITATIONS: pure compression model (tensile stress is always zero);
// confinement formulation is for rectangular sections only; Mfck was
// calibrated by the source paper on a limited number of column tests.

#ifndef GPCConcrete01_h
#define GPCConcrete01_h

#include <UniaxialMaterial.h>

class GPCConcrete01 : public UniaxialMaterial
{
public:
    // Unconfined
    GPCConcrete01(int tag, double fck, double ecu);

    // Confined (rectangular section)
    GPCConcrete01(int tag, double fck, double ecu,
                  double fyh, double diaLat,
                  double bk, double hk, double s,
                  double Mfck = 0.7);

    GPCConcrete01();
    ~GPCConcrete01();

    const char *getClassType(void) const { return "GPCConcrete01"; };

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

private:
    // inputs
    double fck, ecu;
    bool   confined;
    double fyh, diaLat, bk, hk, s, Mfck;   // used only if confined

    // envelope parameters actually used (fcc = fck, eccu = ecu if unconfined)
    double fcc, eccu, E0;

    // committed state (ceps_max: largest compressive strain reached)
    double ceps, csig, ctang;
    double ceps_max;

    // trial state
    double teps, tsig, ttang;
    double teps_max;

    static const double TENS_STIFF_RATIO;

    void computeEnvelopeParams(void);
    void envelope(double eps, double &sig, double &tang);
};

#endif

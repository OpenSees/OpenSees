// GPCConcreteUnconfined01.h
//
// Unconfined geopolymer concrete (GPC) uniaxial material based on:
//   Kocaer O, Aldemir A. "Compressive stress-strain model for the
//   estimation of the flexural capacity of reinforced geopolymer
//   concrete members." Structural Concrete. 2023;24(4):5102-5121.
//   https://doi.org/10.1002/suco.202200914
//
// Monotonic envelope (internal convention: compression positive):
//   sigma(eps) = 6.75*fck*eps*(eps-ecu)^2/ecu^3,  0 <= eps <= ecu
//   sigma(eps) = 0 otherwise (no tension, no strength beyond ecu)
//
// Interface convention: OpenSees sign convention (compression negative).
// fck and ecu may be given as positive or negative values (magnitudes used).
//
// Cyclic rule:
//   Loading follows the envelope. On reversal, unloading follows a
//   straight line of slope E0 (initial tangent) anchored at the reversal
//   point, clipped at zero stress. Reloading retraces the same line
//   until the reversal strain is reached, then rejoins the envelope.
//

#ifndef GPCConcreteUnconfined01_h
#define GPCConcreteUnconfined01_h

#include <UniaxialMaterial.h>

class GPCConcreteUnconfined01 : public UniaxialMaterial
{
public:
    GPCConcreteUnconfined01(int tag, double fck, double ecu);
    GPCConcreteUnconfined01();
    ~GPCConcreteUnconfined01();

    const char *getClassType(void) const { return "GPCConcreteUnconfined01"; };

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
    double fck;   // peak compressive strength (magnitude)
    double ecu;   // ultimate strain (magnitude)
    double E0;    // initial tangent = 6.75*fck/ecu

    // committed state (ceps_max: largest compressive strain reached)
    double ceps, csig, ctang;
    double ceps_max;

    // trial state
    double teps, tsig, ttang;
    double teps_max;

    static const double TENS_STIFF_RATIO;

    // envelope in compression-positive convention
    void envelope(double eps, double &sig, double &tang);
};

#endif
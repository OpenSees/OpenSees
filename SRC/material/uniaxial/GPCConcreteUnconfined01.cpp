// GPCConcreteUnconfined01.h
//
// Unconfined GPC uniaxial material based on:
//   Kocaer O, Aldemir A. "Compressive stress-strain model for the
//   estimation of the flexural capacity of reinforced geopolymer
//   concrete members." Structural Concrete. 2023;24(4):5102-5121.
//   https://doi.org/10.1002/suco.202200914
//
// Monotonic envelope (compression positive):
//   sigma(eps) = 6.75 * fck * eps * (eps - ecu)^2 / ecu^3,  0 <= eps <= ecu
//   sigma(eps) = 0                                          otherwise (no tension)
//
// Cyclic rule:
//   - Loading follows the envelope curve above.
//   - On any strain reversal, unloading follows a straight line with
//     slope E0 (initial tangent) anchored at the reversal point
//     (eps_max, envelope(eps_max)), clipped to sigma >= 0.
//   - Reloading retraces the SAME straight line (it is a pure function
//     of trial strain and the stored reversal point, not of direction),
//     until strain reaches eps_max again, at which point the material
//     rejoins the envelope.
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

    void Print(OPS_Stream &s, int flag = 0);

protected:
    // input parameters
    double fck;   // unconfined compressive strength (positive)
    double ecu;   // ultimate strain (positive)
    double E0;    // initial tangent = 6.75*fck/ecu

    // committed history
    double ceps, csig, ctang;
    double ceps_max;   // largest strain ever reached (reversal point)

    // trial state
    double teps, tsig, ttang;
    double teps_max;

    static const double TENS_STIFF_RATIO; // small stiffness used when sigma clipped to 0

    void envelope(double eps, double &sig, double &tang);
};

#endif
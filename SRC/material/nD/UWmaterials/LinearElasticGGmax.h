/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** LinearElasticGGmax: linear elastic nD material with G/Gmax curves  **
** Accepts (G,K) or (G,nu) and degrades only G by G/Gmax(gamma)       **
**                                                                    **
** ****************************************************************** */
#ifndef LinearElasticGGmax_h
#define LinearElasticGGmax_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>

class LinearElasticGGmax : public NDMaterial
{
public:
    // Predefined curves: curveType in {1=Hardin-Drnevich, 2=Vucetic-Dobry, 3=Darendeli}
    LinearElasticGGmax(int tag, double G_in, double K_or_nu, double rho,
                       int curveType, double p1 = 0.0, double p2 = 0.0, double p3 = 0.0,
                       double chi_in = 0.0);

    // User curve (curveType==0)
    LinearElasticGGmax(int tag, double G_in, double K_or_nu, double rho,
                       const std::vector<double>& strains, const std::vector<double>& ggmax,
                       double chi_in = 0.0);

    LinearElasticGGmax();
    ~LinearElasticGGmax();

    const char *getClassType(void) const override { return "LinearElasticGGmax"; };

    // Mandatory NDMaterial interface
    int setTrialStrain(const Vector &strain) override;
    int setTrialStrain(const Vector &strain, const Vector &rate) override;
    int setTrialStrainIncr(const Vector &strain) override;
    int setTrialStrainIncr(const Vector &strain, const Vector &rate) override;

    const Vector &getStrain(void) override;
    const Vector &getStress(void) override;
    const Matrix &getTangent(void) override;
    const Matrix &getInitialTangent(void) override;
    // Damping tangent (viscous dashpot contribution)
    const Matrix &getDampTangent(void);

    int commitState(void) override;
    int revertToLastCommit(void) override;
    int revertToStart(void) override;
    NDMaterial *getCopy(void) override;
    NDMaterial *getCopy(const char *type) override;

    const char *getType(void) const override;
    int getOrder(void) const override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;

    void Print(OPS_Stream &s, int flag = 0) override;

    // NDMaterial::getRho() defaults to 0.0; override to return material density
    double getRho(void) {return rho0;}

    // Parameter control
    int setParameter(const char **argv, int argc, Parameter &param) override;
    int updateParameter(int paramID, Information &info) override;

private:
    // Material properties
    double G0;       // initial shear modulus
    double K0;       // bulk modulus (valid if hasK==true)
    double nu;       // Poisson's ratio (valid if hasK==false)
    bool   hasK;     // true if user gave K; false if user gave nu
    double rho0;     // density (not used in constitutive law here)
    double chi;      // viscous damping coefficient (dashpot scaling)
    double mu_c;     // cached shear modulus for current trial state
    double lambda_c; // cached first Lamé parameter for current trial state

    // Running maximum of equivalent shear strain
    double gammaMaxTrial{0.0};
    double gammaMaxCommit{0.0};
    // Scale factor for equivalent-linear adjustments to gamma (e.g., 0.85)
    double gammaScale{1.0};

    // G/Gmax curve definition
    int    curveType;               // 0=user, 1=Hardin-Drnevich, 2=Vucetic-Dobry, 3=Darendeli
    double param1, param2, param3;  // model params (gamma_ref, PI, p', OCR, etc.)
    std::vector<double> userStrains; // engineering gamma values (>0)
    std::vector<double> userGGmax;   // ratio values in [0,1]

    // State
    Vector epsilon;    // trial strain (6 for 3D Voigt, 3 for 2D)
    Vector Cepsilon;   // committed strain
    Vector sigma;      // trial stress (resized by ctor)
    Vector sigma_vis;  // viscous stress contribution
    Matrix D;          // current elastic tangent
    Matrix D0;         // initial elastic tangent (gg=1)
    Matrix Ddamp;      // damping tangent (chi * D)
    Matrix Cep;        // consistent tangent including viscous term
    int nDim;           // 2 or 3

    // --- G update policy ---
    int  updateStride{1};          // recompute G every N trial steps (default=1 -> current behavior)
    int  stepCounter{0};           // counts setTrialStrain calls
    bool updateOnDemand{false};    // only update when explicitly requested
    bool pendingOneShotUpdate{false}; // armed by "requestUpdate"
    bool debugUpdate{false};      // <-- debug switch (off by default)

    // Last used gg
    double lastGG{1.0}; 

    // --- Optional: stride by *committed time steps* instead of trial calls ---
    bool   strideByStep{false};   // default false -> legacy behavior
    int    commitCounter{0};      // increments once per commitState()
    int    lastUpdateCommit{-1};  // the commit index when we last updated

    // Helpers
    double computeGGmax(double gamma);
    double computeShearStrain(const Vector& strain); // octahedral-equivalent
    void   computeTangent(double gg_ratio);
    void   computeInitialElasticTangent();
    void   computeDampingTangent();

    // Curve models
    double hardinDrnevich(double gamma) const;
    double vuceticDobry(double gamma) const;
    double darendeli(double gamma) const;
    double interpolateUserCurve(double gamma) const;
};

#endif // LINEAR_ELASTIC_GGMAX_H

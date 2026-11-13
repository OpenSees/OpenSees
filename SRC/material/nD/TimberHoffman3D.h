
#ifndef TimberHoffman3D_h
#define TimberHoffman3D_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <classTags.h>



class TimberHoffman3D : public NDMaterial
{
public:
    TimberHoffman3D();

    TimberHoffman3D(
        int tag,
        double E1, double E2, double E3,
        double nu12, double nu13, double nu23,
        double G12, double G13, double G23,
        double fc1, double fc2, double fc3,
        double ft1, double ft2, double ft3,
        double f12, double f13, double f23,
        double h,
        double sigmaE0,
        double Acomp, double Bcomp,
        double Gf1t, double Gf2t, double Gf3t,
        double eta,
        double Lc,
        double dt = 1.0
    );

    ~TimberHoffman3D();

    int setTrialStrain(const Vector &strain) override;
    int setTrialStrain(const Vector &strain, const Vector &rate) override;

    const Vector &getStress(void) override;
    const Vector &getStrain(void) override;
    const Matrix &getTangent(void) override;
    const Matrix &getInitialTangent(void) override;

    int commitState(void) override;
    int revertToLastCommit(void) override;
    int revertToStart(void) override;

    NDMaterial *getCopy(void) override;
    NDMaterial *getCopy(const char *type) override;

    const char *getType(void) const override;
    int getOrder(void) const override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker) override;

    void Print(OPS_Stream &s, int flag = 0) override;

private:
    void buildElasticStiffness();
    void buildHoffmanPQ();

    double sigmaEk(double eqp) const;
    double yieldFunction(const Vector &stressEff, double eqp) const;
    void flowDirection(const Vector &stressEff, double eqp, Vector &g) const;

    double equivalentPlasticIncrement(const Vector &depsP) const;

    int returnMapHardening(
        const Vector &sigmaTrial,
        Vector &sigmaNew,
        double &eqpNew,
        Vector &depsP,
        double &dLambda
    );

    int stateForLambda(
        const Vector &sigmaTrial,
        double dLambda,
        Vector &sigma,
        double &eqp,
        Vector &depsP
    );

    void failureIndices(
        const Vector &stressEff,
        Vector &Ft,
        double &F1c
    ) const;

    double tensileDamage(double F, int direction) const;
    double compressionDamage(double F1c) const;

    void updateDamage(const Vector &stressEff);
    void buildDamagedStiffness(const Vector &damage);

private:
    // Material parameters
    double E1, E2, E3;
    double nu12, nu13, nu23;
    double G12, G13, G23;

    double fc1, fc2, fc3;
    double ft1, ft2, ft3;
    double f12, f13, f23;

    double hardeningModulus;
    double sigmaE0;

    double Acomp;
    double Bcomp;

    double Gf1t;
    double Gf2t;
    double Gf3t;

    double eta;
    double Lc;
    double dt;

    // Elastic and Hoffman tensors
    Matrix C0;
    Matrix Cdam;
    Matrix P;
    Vector Q0;
    Matrix Z;

    // Trial state
    Vector strainTrial;
    Vector stressTrial;
    Vector effectiveStressTrial;
    Vector plasticStrainTrial;

    double eqpTrial;
    Vector damageTTrial;
    double damageC1Trial;
    Vector damageVTrial;
    Vector FmaxTTrial;
    double FmaxC1Trial;

    // Committed state
    Vector strainCommit;
    Vector stressCommit;
    Vector effectiveStressCommit;
    Vector plasticStrainCommit;

    double eqpCommit;
    Vector damageTCommit;
    double damageC1Commit;
    Vector damageVCommit;
    Vector FmaxTCommit;
    double FmaxC1Commit;

    // Solver controls
    double tol;
    int maxIter;
    double innerTol;
    int innerMaxIter;
};

void *OPS_TimberHoffman3D();

#endif

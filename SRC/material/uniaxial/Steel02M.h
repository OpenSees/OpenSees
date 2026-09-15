/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

#ifndef Steel02M_h
#define Steel02M_h

#include <UniaxialMaterial.h>
#include <tuple> 

class Steel02M : public UniaxialMaterial
{
  public:
    Steel02M(int tag, double E0,
             double Fy01, double Fy02, 
             double alphaP, double alphaN,
             double R0, double cR1, double cR2,
             double a1, double a2,
             double Fysfy01, double Fysfy02,
             double b1, double b2);
    
    Steel02M();
    virtual ~Steel02M();

    const char *getClassType(void) const {return "Steel02M";};

    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);      
    double getStress(void);
    double getTangent(void);
    
    int commitState(void) override;
    int revertToLastCommit(void);    
    int revertToStart(void);         
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag = 0);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    static int globalMaxIter;

  private:
  
    static std::tuple<double, double> Steel02M_Intersection(double E0, double epsr, double Esh, double fy, double sigr);


    static std::tuple<double, double, double> Steel02M_IntersectionBasic(double alpha, double epsb0, double epsbr, 
                                                                  double epsr, double R, double sigb0, 
                                                                  double sigbr, double sigr, double eps_m, double sig_m);
  
    static std::tuple<double, double, double> Steel02M_IntersectionMCL1(double alphaL[3][2], double E0, double eps0L[3][2], 
                                                                 double epsrL[3][2], double epsr, 
                                                                 double RL[3][2], double sig0L[3][2], double sigrL[3][2], 
                                                                 double sigr, int idx, double eps_n, double sig_n, double sig_nb, double Et);

    static std::tuple<double, double> Steel02M_GMP_Stress(double alpha, double eps, double eps0, double epsr, 
                                                   double R, double sig0, double sigr);
												   

    static std::tuple<double, double> Steel02M_MCL1_Stress(double aL[3][2], double e, double e0L[3][2], 
                                                    double erL[3][2], double rL[3][2], double s0L[3][2], 
                                                    double srL[3][2], int idx);
												

    static std::tuple<double, double> Steel02M_MCL2_Stress(double aL[3][2], double e, double e0L[3][2], 
                                                    double erL[3][2], double rL[3][2], double s0L[3][2], 
                                                    double srL[3][2], int idx);
																			

    // --- Material Parameters (Input) ---
    double E0;
    double R0;
    double alpha[2];    // [0]: Positive, [1]: Negative
    double Fy0[2];      // [0]: Positive, [1]: Negative
    double cR[2];
    double a[2];
    double b[2];
    double Fysfy0[2];

    // --- State Variables (Trial) ---
    double eps;
    double sig;
    double Et;
    double epsr;
    double sigr;
    double sig0;
    int    Flag[3];     // [0]: Direction, [1]: Curve Level, [2]: Yield trigger
    
    // 3x2 History Matrices (Rows: 0=Basic, 1=MCL1, 2=MCL2 | Columns: 0=Pos, 1=Neg)
    double alphaL[3][2];
    double epsPL[3][2];

    double eps0L[3][2];
    double epsrL[3][2];
    double sig0L[3][2];
    double sigrL[3][2];
    double RL[3][2];

    double epsy[2];
    double Fy[2];
    double Fys[2];

    // --- Committed Variables (History) ---
    double epsP;
    double sigP;
    double eP;
    double epsr_P;
    double sigr_P;
    double sig0P;
    int    FlagP[3];

    double alphaLP[3][2];
    double epsPL_P[3][2];

    double eps0LP[3][2];
    double epsrLP[3][2];
    double sig0LP[3][2];
    double sigrLP[3][2];
    double RLP[3][2];

    double FyP[2];
    double FysP[2];
};

#endif
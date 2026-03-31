/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        

// $LGap: gap length to simulate the gap length due to the pin tolerance
// $NM:	Employed adaptive numerical algorithm (default value NM = 1; 1 = Dormand-Prince54, 2=6th order Adams-Bashforth-Moulton, 3=modified Rosenbrock Triple)
// $RelTol:	Tolerance for absolute relative error control of the adaptive iterative algorithm (default value 10^-6)
// $AbsTol:	Tolerance for absolute error control of adaptive iterative algorithm (default value 10^-6)
// $MaxHalf: Maximum number of sub-step iterations within an integration step, h=dt*(0.5)^MaxHalf (default value 15)

#ifndef APDVFD_h
#define APDVFD_h

#include <UniaxialMaterial.h>


class APDVFD : public UniaxialMaterial
{
  public: 
      APDVFD(int tag, double K, double G1,double G2,double Alpha, double L ,double LC, double DP, double DG, double N1,double N2, double DO1, double DO2, double DC, double S, double HP, double HC, double LGap, double NM, double RelTol, double AbsTol, double MaxHalf);
      APDVFD();
    ~APDVFD();

    const char *getClassType(void) const {return "APDVFD";};

    int setTrialStrain(double strain, double strainRate); 
    double getStrain(void); 
    double getStrainRate(void);
    double getStress(void);

    double getTangent(void);
    double getInitialTangent(void);
    double getDampTangent(void);


    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    double sgn(double dVariable);
    int DormandPrince(double vel0, double vel1, double y0, double h, double& yt, double& eps, double& error);
    int ABM6(double vel0, double vel1, double y0, double h, double& yt, double& eps, double& error);
    int ROS(double vel0, double vel1, double y0, double h, double& y2, double& eps, double& error);
    double f(double v, double fd);


        
    UniaxialMaterial *getCopy(void);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
                 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
  protected:
    
  private:
        // Fixed Input Material Variables
    double K;
    double G1;
    double G2;
    double Alpha;

    double L;
    double LC;
    double DP;
    double DG;
    double N1;
    double N2;
    double DO1;
    double DO2;
    double DC;
    double S;
    double HP;
    double HC;
    double LGap;
    double NM;
    double RelTol;
    double AbsTol;
    double MaxHalf;
    
    // Trial State Variables
    double C;
    
    double Tstrain; // Trial Strain
    double Tstress; // Trial Stress
    double Ttangent; // Trial Tangent
    double TVel;   // Trial Velocity
    double Tpugr;   // Trial gap initiation displacement
    double Tnugr;   // Trial gap initiation displacement

    
    // Committeed State Variables
    double Cstrain;  // Committed Strain
    double Cstress;  // Committed Stress
    double Ctangent; // Committed Tangent
    double CVel;    // Committed velocity
    double Cpugr;	// Committed gap initiation displacement
    double Cnugr;   // Trial gap initiation displacement
};

#endif

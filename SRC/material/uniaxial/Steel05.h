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
                                                                        
// $Revision: 1.0 $                                                                     
// Written: S. A. Jalalli 10/2019
// Adding Cyclic and in-cycle deterioration modes to steel02 UniaxialMaterial


#ifndef Steel05_h
#define Steel05_h

#include <UniaxialMaterial.h>

class Steel05 : public UniaxialMaterial
{
public:
    Steel05(int tag,
        double fy, double E0, double b,
        double ductC, double pcEFac, double gama, double _c, double _resFac,
        double R0, double cR1, double cR2,
        double a1, double a2, double a3, double a4, double sigInit = 0.0);

    Steel05(void);
    virtual ~Steel05();


    const char* getClassType(void) const { return "Steel05"; };

    double getInitialTangent(void);
    UniaxialMaterial* getCopy(void);

    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStress(void);
    double getTangent(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel,
        FEM_ObjectBroker& theBroker);

    void Print(OPS_Stream& s, int flag = 0);

    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

    virtual double getEnergy() { return EnergyP; };

    double getInitYieldStrain() { return Fy / E0; }
    virtual void resetEnergy(void) { EnergyP = 0; }

protected:

private:
    double epsPCFac;		//ratio between post cap strain and yield strain
    double pstcpEFac;		//ratio between post cap stiffness and initial stiffness
    double FailEnerg, c, gama;			//damage parameters
    double resFac;
    double FydP, FydN;		//Pos and Neg Fy's affected by damage
    double ExcurEnergy;
    void updateDamage();
    double EnergyP; //by SAJalali
    // matpar : STEEL FIXED PROPERTIES
    double Fy;  //  = matpar(1)  : yield stress
    double E0;  //  = matpar(2)  : initial stiffness
    double b;   //  = matpar(3)  : hardening ratio (Esh/E0)
    double R0;  //  = matpar(4)  : exp transition elastic-plastic
    double cR1; //  = matpar(5)  : coefficient for changing R0 to R
    double cR2; //  = matpar(6)  : coefficient for changing R0 to R
    double a1;  //  = matpar(7)  : coefficient for isotropic hardening in compression
    double a2;  //  = matpar(8)  : coefficient for isotropic hardening in compression
    double a3;  //  = matpar(9)  : coefficient for isotropic hardening in tension
    double a4;  //  = matpar(10) : coefficient for isotropic hardening in tension
    double sigini; // initial 
    // hstvP : STEEL HISTORY VARIABLES
    double epsminP; //  = hstvP(1) : max eps in compression
    double epsmaxP; //  = hstvP(2) : max eps in tension
    double epsplP;  //  = hstvP(3) : plastic excursion
    double epss0P;  //  = hstvP(4) : eps at asymptotes intersection
    double sigs0P;  //  = hstvP(5) : sig at asymptotes intersection
    double epssrP;  //  = hstvP(6) : eps at last inversion point
    double sigsrP;  //  = hstvP(7) : sig at last inversion point
    int    konP;    //  = hstvP(8) : index for loading/unloading
    // hstv : STEEL HISTORY VARIABLES   
    double epsP;  //  = strain at previous converged step
    double sigP;  //  = stress at previous converged step
    double eP;    //   stiffness modulus at last converged step;

    double epsmin;
    double epsmax;
    double epspl;
    double epss0;
    double sigs0;
    double epsr;
    double sigr;
    int    kon;
    double sig;
    double e;
    double eps;   //  = strain at current step
};


#endif


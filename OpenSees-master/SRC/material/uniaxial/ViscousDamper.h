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
                                                                        
// $Revision: C $
// $Date: May 2015 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ViscousDamper.cpp,v $
                                                                        
// Written: Sarven Akcelyan and Dimitrios G. Lignos, PhD, McGill University
// Created: January 2013
// Updated: May 2015
// Revision: C
//
// Description: This file contains the class interface for 
// Viscous Damper Model Relationship of the form F = K*u_s = C*pow(V_d,alpha)
// Reference: 
// Akcelyan, S., and Lignos, D.G. (2015), “Adaptive Numerical Method Algorithms for Nonlinear Viscous and Bilinear Oil Damper Models Under Random Vibrations”, ASCE Journal of Engineering Mechanics, (under review)
// Kasai K, Oohara K. (2001). “Algorithm and Computer Code To Simulate Response of Nonlinear Viscous Damper”. Proceedings Passively Controlled Structure Symposium 2001, Yokohama, Japan.

// Variables:
// $K: Elastic stiffness of linear spring (to model the axial flexibility of a viscous damper (brace and damper portion)
// $C: Viscous damping coefficient of the damper
// $Alpha: Viscous damper exponent
// $LGap: gap length to simulate the gap length due to the pin tolerance
// $NM:	Employed adaptive numerical algorithm (default value NM = 1; 1 = Dormand-Prince54, 2=6th order Adams-Bashforth-Moulton, 3=modified Rosenbrock Triple)
// $RelTol:	Tolerance for absolute relative error control of the adaptive iterative algorithm (default value 10^-6)
// $AbsTol:	Tolerance for absolute error control of adaptive iterative algorithm (default value 10^-6)
// $MaxHalf: Maximum number of sub-step iterations within an integration step, h=dt*(0.5)^MaxHalf (default value 15)

#ifndef ViscousDamper_h
#define ViscousDamper_h

#include <UniaxialMaterial.h>

class ViscousDamper : public UniaxialMaterial
{
  public: 
    ViscousDamper(int tag, double K, double C, double Alpha, double LGap, double NM, double RelTol, double AbsTol, double MaxHalf);   
    ViscousDamper(); 
    ~ViscousDamper();

    const char *getClassType(void) const {return "ViscousDamper";};

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
    double C;
    double Alpha; 
    double LGap;
    double NM;
    double RelTol;
    double AbsTol;
    double MaxHalf;
    
    // Trial State Variables
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

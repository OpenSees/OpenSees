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
                                                                        
// $Revision: 1 $
// $Date: 2011/02/01 12:35:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ViscousDamper.h,v $
                                                                        
// Written: Sarven Akcelyan and Dimitrios G. Lignos, PhD, McGill University
// Created: January 2013
// Revision: A
//
// Description: This file contains the class interface for 
// Viscous Damper Model Relationship of the form F = K*u + C*pow(V,a)
// Reference: Kasai K, Oohara K. (2001). “Algorithm and Computer Code To Simulate Response of Nonlinear Viscous Damper”. 
// Proceedings Passively Controlled Structure Symposium 2001, Yokohama, Japan.
// Variables:
// K: axial stiffness of a damper
// C: Velocity constant of a damper
// Alpha: Exponent of velocity of a damper

#ifndef ViscousDamper_h
#define ViscousDamper_h

#include <UniaxialMaterial.h>

class ViscousDamper : public UniaxialMaterial
{
  public:
	 ViscousDamper(int tag, double K, double C, double Alpha, double NM, double Tol, double MaxHalf);   
    ViscousDamper(); 
    ~ViscousDamper();

   const char *getClassType(void) const {return "ViscousDamper";};

    int setTrialStrain(double velocity, double strainRate = 0.0); 
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
	int DormandPrince(double vel0, double vel1, double y0, double h, double& yt, double& eps);
	int ABM6(double vel0, double vel1, double y0, double h, double& yt, double& eps);
	int ROS(double vel0, double vel1, double y0, double h, double& y2, double& eps);
	double f(double v, double fd);

        
    UniaxialMaterial *getCopy(void);
    
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
	double NM;
    double Tol;
	double MaxHalf;

        // Trial State Variables
        double Tstrain; // Trial Strain
        double Tstress; // Trial Stress
        double Ttangent; // Trial Tangent
        double TdVel;   // Trial Incremental Velocity

        
        // Committeed State Variables
        double Cstrain;  // Committed Strain
        double Cstress;  // Committed Stress
        double Ctangent; // Committed Tangent
        double CdVel;    // Committed incremental velocity

};


#endif
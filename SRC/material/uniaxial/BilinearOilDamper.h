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
                                                                        
// $Revision: 0 $
// $Date: May 2015 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BilinearOilDamper.h,v $
                                                                        
// Written: Sarven Akcelyan and Dimitrios G. Lignos, PhD, McGill University
// Created: May 2015
// Updated: May 2015
// Revision: A

// Description: This file contains the class interface for 
// Oil Damper Model Relationship of the form  before relief valve ==> F = K*u_s = C*V_d  after relief valve ==> F= K*u_s = Fr + p*C*(V_d-Fr/C)
//
// References: 
// Akcelyan, S., and Lignos, D.G. (2015), “Adaptive Numerical Method Algorithms for Nonlinear Viscous and Bilinear Oil Damper Models Under Random Vibrations”, ASCE Journal of Engineering Mechanics, (under review)
// Kasai, K., Takahashi, O., and Sekiguchi, Y. (2004). "JSSI manual for building passive control technology part-10 time-history analysis model for nonlinear oil dampers." Proc., The 13th World Conference on Earthquake Engineering, Vancouver, B.C., Canada.

// Variables:
// $K: Elastic stiffness of linear spring (to model the axial flexibility of a viscous damper (brace and damper portion)
// $C: Viscous damping coefficient of the damper
// $Fr: Relief load
// $p: post-relief viscous damping coefficient ratio, (p=C2/C)
// $LGap: gap length to simulate the gap length due to the pin tolerance
// $NM:	Employed adaptive numerical algorithm (default value NM = 1; 1 = Dormand-Prince54, 2 = Finite differences)
// $RelTol: Tolerance for absolute relative error control of the adaptive iterative algorithm (default value 10^-6)
// $AbsTol: Tolerance for absolute error control of adaptive iterative algorithm (default value 10^-10)
// $MaxHalf: Maximum number of sub-step iterations within an integration step (default value 15)

#ifndef BilinearOilDamper_h
#define BilinearOilDamper_h

#include <UniaxialMaterial.h>

class BilinearOilDamper : public UniaxialMaterial
{
  public:

    BilinearOilDamper(int tag, double K, double C, double Fr, double p, double LGap, double NM, double RelTol, double AbsTol, double MaxHalf);    
    BilinearOilDamper(); 
    ~BilinearOilDamper();

    const char *getClassType(void) const {return "BilinearOilDamper";}

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void); 
    double getStrainRate(void);
    double getStress(void);

    double getTangent(void);
    double getInitialTangent(void);
    double getDampTangent(void);


    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
        
    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
                 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    
  protected:
    
  private:
    // Fixed Input Material Variables
    double K;
    double C;
    double Fr;   
    double p;   
    double LGap;
    double NM;
    double RelTol;
    double AbsTol;
    double MaxHalf;
    // Trial State Variables
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

    double sgn(double dVariable);
    int DormandPrince(double vel0, double vel1, double y0, double h, double& yt, double& eps, double& error);
    double f(double v, double fd);    
};


#endif

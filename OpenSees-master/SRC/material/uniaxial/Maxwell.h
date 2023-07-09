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
                                                                        
// $Revision$
// $Date$
// $Source$
                                                                        
// Written: Dimitrios G. Lignos, PhD, Assistant Professor, McGill University
// Created: February 2011
// Revision: A
//
// Description: This file contains the class interface for 
// Maxwell Model.  Maxwell defines a force(F)-velocity(v)
// relationship of the form F = K*u + C*pow(V,a)
// Variables:
// K: axial stiffness of a damper
// C: Velocity constant of a damper
// Alpha: Exponent of velocity of a damper
// L: Length of the Damper

// Olsson, A.K., and Austrell, P-E., (2001), "A fitting procedure for viscoelastic-elastoplastic material models," 
// Proceedings of the Second European Conference on Constitutive Models for Rubber, Germany, 2001.
// Ottosen, N.S., and Ristinmaa, M., (1999). "The mechanics of constitutive modelling, (Numerical and thermodynamical topics)," 
// Lund University,Division of Solid Mechanics, Sweden, 1999. 

#ifndef Maxwell_h
#define Maxwell_h

#include <UniaxialMaterial.h>

class Maxwell : public UniaxialMaterial
{
  public:
  Maxwell(int tag, double K, double C, double Alpha, double L, int returnD = 0);    
    Maxwell(); 
    ~Maxwell();

    const char *getClassType(void) const {return "Maxwell";};

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
    double L;
    int returnD;
    
    // Trial State Variables
    double Tstrain; // Trial Strain
    double Tstress; // Trial Stress
    double Ttangent; // Not a state variable but declared for convenience
    
    // Committeed State Variables
    double Cstrain; // Committed Strain
    double Cstress; // Committed Stress
    double Ctangent; 
};


#endif


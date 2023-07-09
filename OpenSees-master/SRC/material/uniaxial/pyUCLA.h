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

// $Date: 2004/09/22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/pyUCLA.h,v $

#ifndef pyUCLA_h
#define pyUCLA_h

// Written: HyungSuk Shin
// Created: Sept. 2004
//
// Description: This file contains the class definition for 
// pyUCLA.  pyUCLA provides the abstraction
// for a one-dimensional p-y contact material. 
// This material is based on Tarciroglu's model

#include <UniaxialMaterial.h>
#include <Matrix.h>

class pyUCLA : public UniaxialMaterial
{
  public:
    pyUCLA(int tag, int soilType, double pult, double y50, double Cd);
    pyUCLA();
    ~pyUCLA();

	const char *getClassType(void) const {return "pyUCLA";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return E;};

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
    // Material parameters
	int    soilType;	// soilType (1 for soft clay)
	double pult;		// pult
	double y50;			// y50 = 2.5 * epsilon_c * D
	double Cd;			// drag resistance ratio (Cd = dragStress/pult)

	double E;			// Elastic modulus
	double Ed;			// Elastic modulus for drag element
	double pult50;		// pult50 = 0.5 * Pu * Ltr = 0.5 * pult
	double n;			// n = 1/3
	double theta;		// theta
	double epsilonY;	// yield strain
	double limitStress; // limit stress
	double dragStress;	// drag stress
	
    // Committed history variables
    double CplasticStrain1;	// Committed plastic strain 1
    double CplasticStrain2;	// Committed plastic strain 2
    double CplasticStrain3;	// Committed plastic strain 3
    double Chardening1;		// Committed internal hardening variable 1
    double Chardening2;		// Committed internal hardening variable 2

	// Trial history variables
    double TplasticStrain1;	// Trial plastic strain 1
    double TplasticStrain2;	// Trial plastic strain 2
    double TplasticStrain3;	// Trial plastic strain 3
    double Thardening1;		// Trial internal hardening variable 1 
    double Thardening2;		// Trial internal hardening variable 2

    // Trial state variables
    double Tstrain;			// Trial strain
	double Tstress;			// Trial stress
	double Ttangent;		// Trial tangent

    double Tstress1;		// Trial stress 1
	double Tstress2;		// Trial stress 2
	double Tstress3;		// Trial stress 3
	double Tstrain1;		// Trial strain 1
	double Tstrain2;		// Trial strain 2
	double Tstrain3;		// Trial strain 3
	double Ttangent1;		// Trial tangent 1
    double Ttangent2;		// Trial tangent 2
    double Ttangent3;		// Trial tangent 3

	// functions
	void projectStressTangent();
 
};


#endif


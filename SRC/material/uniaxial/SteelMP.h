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

// $Revision: 1.2 $
// $Date: 2009-03-27 19:19:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SteelMP.h,v $
                                                                        
// $Function contributed by  Quan Gu & Michele Barbato
#ifndef SteelMP_h
#define SteelMP_h

#define MAT_TAG_SteelMP 1989

#include <UniaxialMaterial.h>

// Default values for isotropic hardening parameters a1, a2, a3, and a4
#define STEEL_03_DEFAULT_A1        0.0
#define STEEL_03_DEFAULT_A2       55.0
#define STEEL_03_DEFAULT_A3        0.0
#define STEEL_03_DEFAULT_A4       55.0

class SteelMP : public UniaxialMaterial
{
  public:
    SteelMP(int tag, double FY, double E, double B, double R,  double CR1, double CR2, double A1, double A2);
    
    ~SteelMP();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {
		return E0;
	};

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);


// AddingSensitivity:BEGIN //////////////////////////////////////////
    int    setParameter             (const char **argv, int argc, Parameter &param);
    int    updateParameter          (int parameterID, Information &info);
	int    activateParameter        (int parameterID);
	double getStressSensitivity     (int gradNumber, bool conditional);
	double getStrainSensitivity		(int gradNumber);
	double getInitialTangentSensitivity(int gradNumber);
	int    commitSensitivity        (double strainGradient, int gradNumber, int numGrads);
// AddingSensitivity:END ///////////////////////////////////////////
    
  protected:
    
  private:
    /*** Material Properties ***/
    double fy;  // Yield stress
    double E0;  // stiffness
    double b;   // Hardening ratio (b = Esh/E0)
    double r0;	// radius of rounded corners  20
    
	double coeffR1;  //18.5
    double coeffR2;  //0.15
    double a1;            //0
    double a2;  // 0

    /*** CONVERGED History Variables ***/
    double CminStrain;  // Minimum strain in compression
    double CmaxStrain;  // Maximum strain in tension

    int Cloading;       // Flag for loading/unloading
                        // 1 = loading (positive strain increment)
                        // -1 = unloading (negative strain increment)
                        // 0 initially
    double CYieldStrain;
    double CYieldStress;
    double CReverStrain;
    double CReverStress;
    double CPlasticExcursion;


    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cstress;
    double Ctangent;


    /*** TRIAL History Variables ***/
    double TminStrain;   // abs of minimum strain
    double TmaxStrain;   // abs of minimum strain

    int Tloading;
    double TYieldStrain;
    double TYieldStress;
    double TReverStrain;
    double TReverStress;
    double TPlasticExcursion;
    
    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience


    // Calculates the trial state variables based on the trial strain
    void determineTrialState (double dStrain);
    double getR (double x_in);



// AddingSensitivity:BEGIN //////////////////////////////////////////


    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////

};

#endif

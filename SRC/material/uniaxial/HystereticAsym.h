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
                                                                        
// $Revision: 1.00 $
// $Date: 2022-Apr-22 12:15:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HystereticAsym.cpp,v $


                                                                        
                                                                        
#ifndef HystereticAsym_h
#define HystereticAsym_h

// Written: Salvatore Sessa Mail: salvatore.sessa2@unina.it
// Created: 04/2022
// Revision: A
//
// Description: This file contains the class implementation for 
// HystereticAsym providing asymmetric hysteresis loops
//
// What: "@(#) HystereticAsym.C, revA"


#include <UniaxialMaterial.h>


class HystereticAsym : public UniaxialMaterial
{
  public:
    HystereticAsym(int tag, double ka, double kb, double fbar, double b1, double b2, double g);
    HystereticAsym();
    ~HystereticAsym();

    const char *getClassType(void) const {return "HystereticAsym";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    // int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);  


//	Response *setResponse(const char **argv, int argc,
    void Print(OPS_Stream &s, int flag =0);


// AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    int    updateParameter          (int parameterID, Information &info);
    int    activateParameter        (int parameterID);
	double getStrainSensitivity(int gradNumber);
    double getStressSensitivity     (int gradIndex, bool conditional);
    double getInitialTangentSensitivity(int gradIndex);
    int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);
// AddingSensitivity:END ///////////////////////////////////////////

    
 protected:
    
 private:
    /*** Material Properties ***/
    double ka;  // Initial stiffness
    double kb;  // Infinite stiffness
    double fo;   // Transition parameter 
    double b1;	// Hardening parameter
    double b2;   // Hardening parameter
	double g;	// modulating parameter
  //double tol; // Tolerance on u00 (default 1.0e-20)

  //double uo; // Auxiliary parameter depending on (k1, k2 and a).
	double a; // Auxiliary parameter
	double InitTangent;		//initial tangent;

	// constant aux vars:
	//double kab;
	//double delta1;
	//double expauo;

	// aux vars:
	double delta2;
	double delta3; 
  //double delta4;
  //double exputujuo;
	double fe;

    
    /*** CONVERGED History Variables ***/
  //double Cuj;


	double st;			// Sign of the trial velocity
                        // 1 = loading (positive strain increment)
                        // -1 = unloading (negative strain increment)
	double uj;

    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cstress;
    double Ctangent;
	

    /*** TRIAL History Variables ***/
	double dStrain;
    
    
    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience

	double signum(double value);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;

  //double Dka;
  //double Dkb;
  //	double Dfo;
  //	double Db1;
  //	double Db2;
  //	double Dg;
  //	double Dtol;
	

	
  //	double Dut;
  //	double Dfc = 0.0;
  //	double Duc;

// AddingSensitivity:END ///////////////////////////////////////////

};

#endif

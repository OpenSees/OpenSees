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
// Added by Liming Jiang, Mian Zhou (07/2016)
// Created: 07/16
//
// Description: This file contains the class definition for 
// StainlessECThermal. StainlessECThermal is modified on the basis of SteelECThermal
// and SteelECthermal is developed for modelling steel material 
// which strictly satisfies Eurocode regarding the temperature dependent properties.
//StainlessECThermal is developed for modelling strainless steel material based on BS EN 1993 1-2:2005 and Annex C


                                                                                                                                                
#ifndef StainlessECThermal_h
#define StainlessECThermal_h



#include <UniaxialMaterial.h>
//
class StainlessECThermal : public UniaxialMaterial
{
  public:
	 StainlessECThermal(int tag, int gradeTag, double fy, double E0, double fu, double sigInit=0.0);
    StainlessECThermal();
    ~StainlessECThermal();

    const char *getClassType(void) const {return "StainlessECThermal";};


    double getThermalElongation(void); //return ThermalElongation
    double getElongTangent(double, double&, double&, double);//Added for temperature dependent tangent and thermalElongation


    int setTrialStrain(double strain, double strainRate =0); 
    int setTrialStrain(double strain, double FiberTemperature, double strainRate); //Added for Temperature-strain-stress

    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return E0;};

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

    int getVariable(const char *variable, Information &);
    
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    int    updateParameter          (int parameterID, Information &info);
    int    activateParameter        (int parameterID);
    double getStressSensitivity     (int gradIndex, bool conditional);
    double getInitialTangentSensitivity(int gradIndex);
    int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////
    
 protected:
    
 private:
	 //
    /////For Input   properties///////////////////////////////start
    int gradeTag;
    //double Temp;  // material temp  
    double fyT;
    double E0T;
	double fuT;  // Ultimate tensile strength, added by MianZhou
    double TemperautreC;
	double EctT;
	double EpsiUT;
	double sigini; // initial  stresses
	double epsini; // initial strain
//
/////For Input properties///////////////////////////////////end    

    /*** Temperature-dependent Material Properties ***/
    double fy;  // proof strength at 0.2% plastic strain
    double E0;  // Initial stiffness
	double fu;
    double b;   // Hardening ratio (b = Esh/E0)
	double Ect;
	double EpsiU;
	double ThermalElongation; // 
	//

//
//
    /*** CONVERGED History Variables ***/
    double CminStrain;  // Minimum strain in compression
    double CmaxStrain;  // Maximum strain in tension
    double CshiftP;     // Shift in hysteresis loop for positive loading
    double CshiftN;     // Shift in hysteresis loop for negative loading
    int Cloading;       // Flag for loading/unloading
                        // 1 = loading (positive strain increment)
                        // -1 = unloading (negative strain increment)
                        // 0 initially
    
    /*** CONVERGED State Variables ***/    
	double Ctemp;
    double Cstrain;
    double Cstress;
    double Ctangent;    
	

    /*** TRIAL History Variables ***/
    double TminStrain;
    double TmaxStrain;
    double TshiftP;
    double TshiftN;
    int Tloading;
    
    /*** TRIAL State Variables ***/
	double Ttemp;
    double Tstrain;
    double Tstress;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience

    // Calculates the trial state variables based on the trial strain
    void determineTrialState (double dStrain);
	//
	// determine the yield surface for the residual/initial stress, convert initial stress to initial strain
	double determineYieldSurface(double sigini);
	//
    // Determines if a load reversal has occurred based on the trial strain
    void detectLoadReversal (double dStrain);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif

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
                                                                        
// $Revision: 1.7 $
// $Date: 2008-08-26 16:30:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticMaterialThermal.h,v $
                                                                        
                                                                        
#ifndef ElasticMaterial_h
#define ElasticMaterial_h

// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for 
// ElasticMaterialThermal. ElasticMaterialThermal provides the abstraction
// of an viscoelastic uniaxial material,
// i.e. stress = E*strain + eta*strainrate
// Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io] 


#include <UniaxialMaterial.h>

class ElasticMaterialThermal : public UniaxialMaterial
{
  public:
    ElasticMaterialThermal(int tag, double Epos, double alpha =0.0, double et=0.0, double Eneg=0.0, int softindex = 0);
    ElasticMaterialThermal();
    ~ElasticMaterialThermal();

    const char *getClassType(void) const {return "ElasticMaterialThermal";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0); 
    double getStrain(void) {return trialStrain;};
    double getStrainRate(void) {return trialStrainRate;};
    double getStress(void);
    double getTangent(void);
    double getDampTangent(void) {return eta;};
    double getInitialTangent(void);

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

	double getThermalElongation(void);//returning the Thermal Elongation
	double getElongTangent(double, double&, double&, double);
	int getVariable(const char *variable, Information &);  //For providing the universal function


    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int activateParameter(int parameterID);
    double getStressSensitivity(int gradIndex, bool conditional);
    double getTangentSensitivity(int gradIndex);
    double getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(double strainGradient, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    double trialStrain;
    double trialStrainRate;
    double committedStrain;
    double committedStrainRate;
    double Epos;
    double Eneg;
    double eta;
	double Alpha;
	double E0;
	double Eneg0;
	int softIndex;
	double ThermalElongation;
	double Temp;
	
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////

    static double ConcRedFactors[12];
    static double SteelRedFactors[12];
};


#endif


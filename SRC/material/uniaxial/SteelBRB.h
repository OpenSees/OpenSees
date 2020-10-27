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
                                                                        

// Written: Q. Gu & Y.Peng
// Created: May 2012
//
// Description: This file contains the class definition for steel buckling restrained braces
// refer to: Alessandro Zona, Andrea Dall'Asta, 2012, "Elastoplastic model for 
// steel buckling restrained braces", Journal of Constructional steel research 68(2012), 118-125
 

#ifndef SteelBRB_h
#define SteelBRB_h

#include <UniaxialMaterial.h>
#include <Information.h>

class SteelBRB : public UniaxialMaterial
{
 public:
  SteelBRB(int tag, 
	double E,
	double sigmaY0,  
	double sigmaY_T,
	double alpha_T,
	double alpha_C,
	double sigmaY_C,
	double beta_T,
	double beta_C,
	double delta_T,
	double delta_C,
	double Tol);

   
  ~SteelBRB();
  
  int setTrialStrain(double strain, double strainRate = 0.0); 
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
  
  void Print(OPS_Stream &s, int flag =0);


 // AddingSensitivity:BEGIN //////////////////////////////////////////
    int    setParameter             (const char **argv, int argc, Information &info);
    int    updateParameter          (int parameterID, Information &info);
	int    activateParameter        (int parameterID);
	double getStressSensitivity     (int gradNumber, bool conditional);
	double getStrainSensitivity		(int gradNumber);
	double getInitialTangentSensitivity(int gradNumber);
	int    commitSensitivity        (double strainGradient, int gradNumber, int numGrads);
	Response* setResponse(const char **argv, int argc, OPS_Stream & theOutput);
	int getResponse(int responseID, Information &matInfo);
// AddingSensitivity:END ///////////////////////////////////////////
	

 protected:
  
 private:


  double  Newton_BRB(double CStress, double beta, double CPlastStrain, double sigmaY, double cumPlastStrain,
	         double delta, double alpha, double strainInc, double x0, double Tol, int N0);

  double PlastStrainIncRes(double CStress, double beta, double CPlastStrain, double sigmaY,
	       double cumPlastStrain, double delta, double alpha, double strainInc, double plastStrainInc);

  double PlastStrainIncResDev(double CStress, double beta, double CPlastStrain, double sigmaY, 
	       double cumPlastStrain, double delta, double alpha, double strainInc, double plastStrainInc );

 
	double tangent;  // trial tangent


// ---- committed values in the last converged step---- 
	double CStress;
	double CPlastStrain;
	double CCumPlastStrain;
	double CYieldStress;
	double CStrain;


	double stress;
	double plastStrain;
	double cumPlastStrain;
	double yieldStress;
	double strain;

// -------- model data ----------
	double E;

	double sigmaY0;  
	double sigmaY_T;
	double alpha_T;
	double alpha_C;
	double sigmaY_C;
	double beta_T;
	double beta_C;
	double delta_T;
	double delta_C;

// -- 
	double CDissipatedEnergy;
	double dissipatedEnergy;

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////


	double Tol;
 
#if !_DLL
	ofstream* debug1;
#endif


};

#endif


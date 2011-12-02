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

// $Revision: 1.1 $
// $Date: 2009-05-29 19:10:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/UniaxialJ2Plasticity.h,v $

#ifndef UniaxialJ2Plasticity_h
#define UniaxialJ2Plasticity_h

// Written: QGU. UCSD Group for response sensitivity analysis
//  Contact: Quan Gu(guquan2005@gmail.com)  
//           Joel P. Conte (jpconte@ucsd.edu)
//           Michele Barbato (mbarbato@lsu.edu)
//
// Created: May 2009
// 
// refer to paper:
// Conte, J. P., Vijalapura, P., and Meghella, M., "Consistent Finite Element Response Sensitivities Analysis," 
// Journal of Engineering Mechanics, ASCE, Vol. 129, No. 12, pp. 1380-1393, 2003.
//

#include <UniaxialMaterial.h>
#include <Matrix.h>


class UniaxialJ2Plasticity : public UniaxialMaterial
{
  public:
    UniaxialJ2Plasticity(int tag, double E, double sigmaY,
		      double Hkin, double Hiso);
    ~UniaxialJ2Plasticity();

    const char *getClassType(void) const {return "UniaxialJ2Plasticity";};

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
    
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter (const char **argv, int argc, Parameter &param);
    int    updateParameter          (int parameterID, Information &info);
    int    activateParameter        (int parameterID);
    double getStressSensitivity     (int gradIndex, bool conditional);
    double getInitialTangentSensitivity    (int gradIndex);
    int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);
	double getStrainSensitivity(int gradIndex);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    // Material parameters
    double E;	// Elastic modulus
    double sigmaY;	// Yield stress
    double Hiso;	// Isotropic hardening parameter
    double Hkin;	// Kinematic hardening parameter
    //double eta;
	
    // Committed history variables
    double CPlasticStrain;	// Committed plastic strain
    double CBackStress;		// Committed back stress;
    double CAccumulatedPlasticStrain;		// Committed accumulated plastic strain

	// Trial history variables
    double TPlasticStrain;	// Trial plastic strain
    double TBackStress;		// Trial back stress
    double TAccumulatedPlasticStrain;		// Trial accumulated plastic strain

    // Trial state variables
    double TStrain;		// Trial strain
    double TStress;		// Trial stress
    double TTangent;	// Trial tangent

    double CStrain;		// Committed strain
    double CStress;		// Committed stress
    double CTangent;	// Committed tangent


// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////
};


#endif


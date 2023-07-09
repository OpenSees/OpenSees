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
// $Date: 2008-08-26 16:30:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticBDMaterial.h,v $
                                                                        
#ifndef ElasticBDMaterial_h
#define ElasticBDMaterial_h

#include <UniaxialMaterial.h>

#define MAT_TAG_ElasticBDMaterial 1976

class ElasticBDMaterial : public UniaxialMaterial
{
  public:
    ElasticBDMaterial(int tag, double E, double b, double d, double a, double eta = 0.0);    
    ElasticBDMaterial();    
    ~ElasticBDMaterial();

    const char *getClassType(void) const {return "ElasticBDMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0); 
    double getStrain(void) {return trialStrain;};
    double getStrainRate(void) {return trialStrainRate;};
    double getStress(void);
    double getTangent(void) {return a*b*d*E;};
    double getDampTangent(void) {return eta;};
    double getInitialTangent(void) {return a*b*d*E;};

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

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int activateParameter(int parameterID);
    double getStressSensitivity(int gradIndex, bool conditional);
    double getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(double strainGradient, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    double trialStrain;
    double trialStrainRate;
    double E;
    double b;
    double d;
    double a;
    double eta;

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////
};


#endif


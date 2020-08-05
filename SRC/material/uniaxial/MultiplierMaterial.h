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
                                                      
// Description: This file contains the class definition for 
// MultiplierMaterial.  MultiplierMaterial wraps a UniaxialMaterial
// and multiplies the stress and tangent of the wrapped UniaxialMaterial
// object by a factor. This wrapper can be used to apply overstrength
// factors to materials and p-y multipliers for shadowing effects in
// pile groups. An example is shown below with a multiplier of 0.8.

#ifndef MultiplierMaterial_h
#define MultiplierMaterial_h

#include <UniaxialMaterial.h>

class MultiplierMaterial : public UniaxialMaterial
{
  public:
    MultiplierMaterial(int tag, UniaxialMaterial &material, double mult); 
    MultiplierMaterial();
    ~MultiplierMaterial();
    
    const char *getClassType(void) const {return "MultiplierMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrialStrain(double strain, double FiberTemperature, double strainRate); 
    double getStrain(void);          
    double getStrainRate(void);
    double getStress(void);
    double getTangent(void);
    double getDampTangent(void);
    double getInitialTangent(void) {return theMaterial->getInitialTangent();}

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
	UniaxialMaterial *theMaterial;

	double multiplier;
};


#endif


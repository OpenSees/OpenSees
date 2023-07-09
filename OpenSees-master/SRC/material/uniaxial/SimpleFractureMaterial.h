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
                                                                        
// Written: fmk
// Created: Feb 2015
//
// Description: This file contains the class definition for 
// SimpleFractureMaterial.  SimpleFractureMaterial wraps a UniaxialMaterial
// and implements a fracture of the material. if material reaches a certain
// tensile strain, the material will no longer provide any tensile resistance
// but will provide compressive resistance once the material reaches a strain
// at which unloading from the max strain would cause the original material to
// again provide compressive stress (an estimate of the elastic strain)

#ifndef SimpleFractureMaterial_h
#define SimpleFractureMaterial_h

#include <UniaxialMaterial.h>

class SimpleFractureMaterial : public UniaxialMaterial
{
  public:
    SimpleFractureMaterial(int tag, UniaxialMaterial &material, double max); 
    SimpleFractureMaterial();
    ~SimpleFractureMaterial();
    
    const char *getClassType(void) const {return "SimpleFractureMaterial";}

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrialStrain(double strain, double FiberTemperature, double strainRate); 
    double getStrain(void);          
    double getStrainRate(void);
    double getStress(void);
    double getTangent(void);
    double getDampTangent(void);
  double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    bool hasFailed(void) {return Cfailed;}

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
  
  protected:
    
  private:
	UniaxialMaterial *theMaterial;

	double maxStrain;

	double TstartCompStrain;
	double CstartCompStrain;

	bool Tfailed;
	bool Cfailed;

	double Tstress;
	double Tstrain;
	double Ttangent;

	double Cstress;
	double Cstrain;
	double Ctangent;
};


#endif


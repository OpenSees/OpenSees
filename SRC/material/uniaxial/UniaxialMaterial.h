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
                                                                        
// $Revision: 1.12 $
// $Date: 2008-08-26 16:36:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/UniaxialMaterial.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 05/98
// Revision: A
//
// Description: This file contains the class definition for 
// UniaxialMaterial. UniaxialMaterial is a base class and 
// thus no objects of it's type can be instantiated. It has pure virtual 
// functions which must be implemented in it's derived classes. 
//
// What: "@(#) UniaxialMaterial.h, revA"

#ifndef UniaxialMaterial_h
#define UniaxialMaterial_h

#define POS_INF_STRAIN        1.0e16
#define NEG_INF_STRAIN       -1.0e16

#include <Material.h>
class ID;
class Vector;
class Matrix;
class Information;
class Response;

class SectionForceDeformation;

class UniaxialMaterial : public Material
{
  public:
    UniaxialMaterial(int tag, int classTag);    
    UniaxialMaterial();
    virtual ~UniaxialMaterial();

    virtual int setTrialStrain (double strain, double strainRate =0) =0;
    virtual int setTrialStrain (double strain, double temperature, double strainRate);
    virtual int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    virtual int setTrial (double strain, double temperature, double &stress, double &tangent, double &thermalElongation, double strainRate = 0.0);

    virtual double getStrain (void) = 0;
    virtual double getStrainRate (void);
    virtual double getStress (void) = 0;
    virtual double getTangent (void) = 0;
    virtual double getInitialTangent (void) = 0;
    virtual double getDampTangent (void);
    virtual double getRho(void);
    
    virtual int commitState (void) = 0;
    virtual int revertToLastCommit (void) = 0;    
    virtual int revertToStart (void) = 0;        
    
    virtual UniaxialMaterial *getCopy (void) = 0;
    virtual UniaxialMaterial *getCopy(SectionForceDeformation *s);
    
    virtual Response *setResponse (const char **argv, int argc, 
				   OPS_Stream &theOutputStream);
    virtual int getResponse (int responseID, Information &matInformation);    
    virtual bool hasFailed(void) {return false;}

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    virtual double getStressSensitivity     (int gradIndex, bool conditional);
    virtual double getStrainSensitivity     (int gradIndex);
    virtual double getTangentSensitivity(int gradIndex);
    virtual double getInitialTangentSensitivity(int gradIndex);
    virtual double getDampTangentSensitivity(int gradIndex);
    virtual double getRhoSensitivity        (int gradIndex);
    virtual int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////
    
 protected:
    
 private:
};

extern bool OPS_addUniaxialMaterial(UniaxialMaterial *newComponent);
extern UniaxialMaterial *OPS_getUniaxialMaterial(int tag);
extern void OPS_clearAllUniaxialMaterial(void);
extern void OPS_printUniaxialMaterial(OPS_Stream &s, int flag = 0);

#endif

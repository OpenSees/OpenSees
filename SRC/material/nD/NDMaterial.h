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
**                                                                    **
** ****************************************************************** */

// $Revision: 1.23 $
// $Date: 2010-09-13 21:29:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/NDMaterial.h,v $


#ifndef NDMaterial_h
#define NDMaterial_h

// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class definition for NDMaterial.
// NDMaterial is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes.
//
// What: "@(#) NDMaterial.h, revA"

#include <Material.h>

class Matrix;
class ID;
class Vector;
class Information;
class Response;

class NDMaterial : public Material
{
  public:
    NDMaterial(int tag, int classTag);
    NDMaterial();
    virtual ~NDMaterial();

    // methods to set state and retrieve state using Matrix and Vector classes
    virtual double getRho(void);

    virtual int setTrialStrain(const Vector &v);
    virtual int setTrialStrain(const Vector &v, const Vector &r);
    virtual int setTrialStrainIncr(const Vector &v);
    virtual int setTrialStrainIncr(const Vector &v, const Vector &r);
    virtual const Matrix &getTangent(void);
    virtual const Matrix &getInitialTangent(void) {return this->getTangent();};

	//Added by L.Jiang, [SIF]
	virtual double getThermalTangentAndElongation(double &TempT, double &, double &);
	virtual double setThermalTangentAndElongation(double &TempT, double &, double &);
	virtual const Vector& getTempAndElong(void);
	//Added by L.Jiang, [SIF]

    virtual const Vector &getStress(void);
    virtual const Vector &getStrain(void);

    virtual int commitState(void) = 0;
    virtual int revertToLastCommit(void) = 0;
    virtual int revertToStart(void) = 0;

    virtual NDMaterial *getCopy(void) = 0;
    virtual NDMaterial *getCopy(const char *code);

    virtual const char *getType(void) const = 0;
    virtual int getOrder(void) const {return 0;};  //??

    virtual Response *setResponse (const char **argv, int argc, OPS_Stream &s);
    virtual int getResponse (int responseID, Information &matInformation);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    virtual const Vector & getStressSensitivity         (int gradIndex, bool conditional);
    virtual const Vector & getStrainSensitivity         (int gradIndex);
    virtual const Matrix & getTangentSensitivity        (int gradIndex);
    virtual const Matrix & getInitialTangentSensitivity (int gradIndex);
    virtual const Matrix & getDampTangentSensitivity    (int gradIndex);
    virtual double         getRhoSensitivity            (int gradIndex);
    virtual int            commitSensitivity            (const Vector & strainGradient, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:

  private:
    static Matrix errMatrix;
    static Vector errVector;
};

extern bool OPS_addNDMaterial(NDMaterial *newComponent);
extern NDMaterial *OPS_getNDMaterial(int tag);
extern bool OPS_removeNDMaterial(int tag);
extern void OPS_clearAllNDMaterial(void);
extern void OPS_printNDMaterial(OPS_Stream &s, int flag = 0);

#endif

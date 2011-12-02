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
** Additions and changes by:                                          **
**   Boris Jeremic (@ucdavis.edu)                                     **
**                                                                    **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.22 $
// $Date: 2009-10-13 21:11:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/NDMaterial.h,v $


#ifndef NDMaterial_h
#define NDMaterial_h

// File: ~/material/NDMaterial.h
//
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
#include <Tensor.h>

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

    virtual const Vector &getStress(void);
    virtual const Vector &getStrain(void);

    // methods to set and retrieve state using the Tensor class
    virtual int setTrialStrain(const Tensor &v);
    virtual int setTrialStrain(const Tensor &v, const Tensor &r);
    virtual int setTrialStrainIncr(const Tensor &v);
    virtual int setTrialStrainIncr(const Tensor &v, const Tensor &r);
    virtual const Tensor& getTangentTensor(void);
    virtual const stresstensor& getStressTensor(void);
    virtual const straintensor& getStrainTensor(void);
    //Added Joey Aug. 13, 2001
    virtual const straintensor& getPlasticStrainTensor(void);

    // added Sept 22 2003 for Large Deformation, F is the Deformation Gradient
    virtual int setTrialF(const straintensor &f);
    virtual int setTrialFIncr(const straintensor &df);
    virtual int setTrialC(const straintensor &c);
    virtual int setTrialCIncr(const straintensor &dc);
    virtual const stresstensor& getPK1StressTensor(void);
    virtual const stresstensor& getCauchyStressTensor(void);
    virtual const straintensor& getF(void);
    virtual const straintensor& getC(void);
    virtual const straintensor& getFp(void);
    // Only For Large Deformation, END////////////////////////////////////////

    virtual int commitState(void) = 0;
    virtual int revertToLastCommit(void) = 0;
    virtual int revertToStart(void) = 0;

    virtual NDMaterial *getCopy(void) = 0;
    virtual NDMaterial *getCopy(const char *code);

    virtual const char *getType(void) const = 0;
    virtual int getOrder(void) const {return 0;};  //??

    virtual Response *setResponse (const char **argv, int argc, 
				   OPS_Stream &s);
    virtual int getResponse (int responseID, Information &matInformation);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    virtual const Vector & getStressSensitivity     (int gradIndex, bool conditional);
    virtual const Vector & getStrainSensitivity     (int gradIndex);
    virtual const Matrix & getTangentSensitivity    (int gradIndex);
    virtual const Matrix & getDampTangentSensitivity(int gradIndex);
    virtual double         getRhoSensitivity        (int gradIndex);
    virtual int            commitSensitivity        (Vector & strainGradient, int gradIndex, int numGrads);
// AddingSensitivity:END ///////////////////////////////////////////

  protected:

  private:
    static Matrix errMatrix;
    static Vector errVector;
    static Tensor errTensor;
    static stresstensor errstresstensor;
    static straintensor errstraintensor;
};


#endif

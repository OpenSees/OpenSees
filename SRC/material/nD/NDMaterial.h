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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
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

class NDMaterial : public Material
{
  public:
    NDMaterial(int tag, int classTag);
    NDMaterial();
    virtual ~NDMaterial();

    // methods to set state and retrieve state using Matrix and Vector classes
    virtual int setTrialStrain(const Vector &v) = 0;
    virtual int setTrialStrain(const Vector &v, const Vector &r) = 0;
    virtual int setTrialStrainIncr(const Vector &v) = 0;
    virtual int setTrialStrainIncr(const Vector &v, const Vector &r) = 0;
    virtual const Matrix &getTangent(void) = 0;
    virtual const Vector &getStress(void) = 0;
	virtual const Vector &getStrain(void) = 0;

    // methods to set and retrieve state using the Tensor class    
    virtual int setTrialStrain(const Tensor &v);
    virtual int setTrialStrain(const Tensor &v, const Tensor &r);
    virtual int setTrialStrainIncr(const Tensor &v);
    virtual int setTrialStrainIncr(const Tensor &v, const Tensor &r);
    virtual const Tensor &getTangentTensor(void);
    virtual const Tensor &getStressTensor(void);
    virtual const Tensor &getStrainTensor(void);

    virtual int commitState(void) = 0;
    virtual int revertToLastCommit(void) = 0;
    virtual int revertToStart(void) = 0;
    
    virtual NDMaterial *getCopy(void) = 0;
    virtual NDMaterial *getCopy(const char *code) = 0;

    virtual const char *getType(void) const = 0;
    virtual int getOrder(void) const = 0;

    virtual int setResponse (char **argv, int argc, Information &matInformation);
    virtual int getResponse (int responseID, Information &matInformation);

  protected:
    
  private:
};


#endif

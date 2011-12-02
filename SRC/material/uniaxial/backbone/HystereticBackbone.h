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
// $Date: 2008-11-09 06:05:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/HystereticBackbone.h,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the interface for HystereticBackbone,
// which represents a backbone curve for hysteretic models.

#ifndef HystereticBackbone_h
#define HystereticBackbone_h

#include <TaggedObject.h>
#include <MovableObject.h>

class Information;

class HystereticBackbone : public TaggedObject, public MovableObject
{
 public:
  HystereticBackbone(int tag, int classTag);
  virtual ~HystereticBackbone();
  
  virtual double getStress(double strain) = 0;
  virtual double getTangent(double strain) = 0;
  virtual double getEnergy(double strain) = 0;
  
  virtual double getYieldStrain(void) = 0;
  
  virtual HystereticBackbone *getCopy(void) = 0;
  
  virtual int setVariable(char *argv);
  virtual int getVariable(int varID, double &theValue);
  
  virtual int setParameter(char **argv, int argc, Information &eleInformation);
  virtual int updateParameter(int responseID, Information &eleInformation);	
  
 protected:
  
 private:
  
};

extern bool OPS_addHystereticBackbone(HystereticBackbone *newComponent);
extern HystereticBackbone *OPS_getHystereticBackbone(int tag);
extern void OPS_clearAllHystereticBackbone(void);

#endif

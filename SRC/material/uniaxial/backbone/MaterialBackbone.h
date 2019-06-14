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

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the implementation of 
// MaterialBackbone, which treats a uniaxial material as
// a hysteretic backbone by removing path dependency, i.e.
// commitState is never called on the uniaxial material.

#ifndef MaterialBackbone_h
#define MaterialBackbone_h

#include <HystereticBackbone.h>
#include <UniaxialMaterial.h>

class MaterialBackbone : public HystereticBackbone
{
 public:
  MaterialBackbone(int tag, UniaxialMaterial &material);
  MaterialBackbone();
  ~MaterialBackbone();
  
  double getStress(double strain);
  double getTangent(double strain);
  double getEnergy(double strain);
  
  double getYieldStrain(void);
  
  HystereticBackbone *getCopy(void);
  
  void Print(OPS_Stream &s, int flag = 0);
  
  int setVariable(char *argv);
  int getVariable(int varID, double &theValue);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
 protected:
  
 private:
  UniaxialMaterial *theMaterial;
};

#endif

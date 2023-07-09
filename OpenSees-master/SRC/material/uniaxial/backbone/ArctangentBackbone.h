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
// $Date: 2008-11-09 06:03:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/ArctangentBackbone.h,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the implementation of 
// ArctangentBackbone, which is a continuous function given
// by K1*atan(K2*strain); as developed by Ranzo and Petrangeli (1998)

#ifndef ArctangentBackbone_h
#define ArctangentBackbone_h

#include <HystereticBackbone.h>
#include <Vector.h>

class ArctangentBackbone : public HystereticBackbone
{
 public:
  ArctangentBackbone(int tag, double K1, double gammaY, double alpha);
  ArctangentBackbone();
  ~ArctangentBackbone();
  
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
  double K1;
  double K2;
  double gammaY;
  double alpha;
};

#endif

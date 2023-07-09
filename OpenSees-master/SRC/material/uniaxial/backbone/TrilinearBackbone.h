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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/TrilinearBackbone.h,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the implementation of 
// TrilinearBackbone, which is a trilinear backbone

#ifndef TrilinearBackbone_h
#define TrilinearBackbone_h

#include <HystereticBackbone.h>
#include <Vector.h>

class TrilinearBackbone : public HystereticBackbone
{
 public:
  TrilinearBackbone(int tag, double e1, double s1, 
		    double e2, double s2, double e3, double s3);
  TrilinearBackbone(int tag, double e1, double s1, 
		    double e2, double s2);
  TrilinearBackbone();
  ~TrilinearBackbone();
  
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
  double E1, E2, E3;
  double e1, e2, e3;
  double s1, s2, s3;
};

#endif

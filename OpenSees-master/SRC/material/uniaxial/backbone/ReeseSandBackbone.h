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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-11-24 17:12:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/ReeseSandBackbone.h,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the implementation of 
// ReeseSandBackbone.

#ifndef ReeseSandBackbone_h
#define ReeseSandBackbone_h

#include <HystereticBackbone.h>
#include <Vector.h>

class ReeseSandBackbone : public HystereticBackbone
{
 public:
  ReeseSandBackbone(int tag, double kx, double ym, double pm,
		    double yu, double pu);
  ReeseSandBackbone();
  ~ReeseSandBackbone();
  
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
  double kx;
  double ym;
  double pm;
  double yu;
  double pu;
};

#endif

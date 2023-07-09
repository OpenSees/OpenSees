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
// $Date: 2003-02-14 23:00:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/eigenAlgo/StandardEigenAlgo.h,v $
                                                                        
// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition of StandardEigenAlgo.
// StandardEigenAlgo is a class which performs a eigen solution algorithm
// to solve standard eigenvalue equations. It is not expected that 
// this class will have subclasses.

#ifndef StandardEigenAlgo_h
#define StandardEigenAlgo_h

#include <EigenAlgorithm.h>

class StandardEigenAlgo : public EigenAlgorithm
{
 public:
  StandardEigenAlgo();
  virtual ~StandardEigenAlgo();
  
  virtual int solveCurrentStep(int numModes);
  
  virtual int sendSelf(int commitTag, Channel &theChannel);
  virtual int recvSelf(int commitTag, Channel &theChannel,
		       FEM_ObjectBroker &theBroker);
  
  virtual void Print(OPS_Stream &s, int flag = 0);
  
 protected:
  
 private:
  
};

#endif

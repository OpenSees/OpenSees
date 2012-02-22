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
// $Date: 2008-08-26 17:07:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/AcceleratedNewton.h,v $
                                                                        
#ifndef AcceleratedNewton_h
#define AcceleratedNewton_h

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for 
// AcceleratedNewton.  AcceleratedNewton is a class which uses a Krylov
// subspace accelerator on the modified Newton method.
// The accelerator is described by Carlson and Miller in
// "Design and Application of a 1D GWMFE Code"
// from SIAM Journal of Scientific Computing (Vol. 19, No. 3,
// pp. 728-765, May 1998)

#include <EquiSolnAlgo.h>
#include <Vector.h>
//#include <Timer.h>

class Accelerator;

class AcceleratedNewton: public EquiSolnAlgo
{
 public:
  AcceleratedNewton(int tangent = CURRENT_TANGENT);
  AcceleratedNewton(ConvergenceTest &theTest, Accelerator *theAccel,
		    int tangent = CURRENT_TANGENT);
  ~AcceleratedNewton();
  
  int solveCurrentStep(void);    
  int setConvergenceTest(ConvergenceTest *theNewTest);
  ConvergenceTest *getTest(void);     
  
  int getNumFactorizations(void) {return numFactorizations;}
  int getNumIterations(void) {return numIterations;}
  //double getTotalTimeCPU(void)   {return totalTimeCPU;}
  //double getTotalTimeReal(void)  {return totalTimeReal;}
  //double getSolveTimeCPU(void)   {return solveTimeCPU;}
  //double getSolveTimeReal(void)  {return solveTimeReal;}
  //double getAccelTimeCPU(void)   {return accelTimeCPU;}
  //double getAccelTimeReal(void)  {return accelTimeReal;}
  
  virtual int sendSelf(int commitTag, Channel &theChannel);
  virtual int recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);    
  
 protected:
  
 private:
  ConvergenceTest *theTest;
  int tangent;
  
  Accelerator *theAccelerator;
  
  // Storate for accelerated mod-Newton prediction
  Vector *vAccel;
  
  int numFactorizations;
  int numIterations;

  //Timer totalTimer;
  //double totalTimeReal;
  //double totalTimeCPU;

  //Timer solveTimer;
  //double solveTimeReal;
  //double solveTimeCPU;

  //Timer accelTimer;
  //double accelTimeReal;
  //double accelTimeCPU;

  bool firstTangent;
};

#endif

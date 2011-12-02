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
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/MinUnbalDispNorm.h,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/MinUnbalDispNorm.h
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the class definition for MinUnbalDispNorm.
// MinUnbalDispNorm is an algorithmic class for perfroming a static analysis
// using the minimum unbalanced displacement norm (Chan IJNME 26(2657:2669)1988
//
// What: "@(#) MinUnbalDispNorm.h, revA"

#ifndef MinUnbalDispNorm_h
#define MinUnbalDispNorm_h

#include <StaticIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;

class MinUnbalDispNorm : public StaticIntegrator
{
  public:
    MinUnbalDispNorm(double lambda1, int specNumIterStep, 
		     double dlambda1min, double dlambda1max);

    ~MinUnbalDispNorm();

    int newStep(void);    
    int update(const Vector &deltaU);
    int domainChanged(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(ostream &s, int flag =0);    
    
  protected:
    
  private:
    double dLambda1LastStep;                  // dLambda1 at step (i-1)
    double specNumIncrStep, numIncrLastStep;    // Jd & J(i-1) 

    Vector *deltaUhat, *deltaUbar, *deltaU, *deltaUstep; // vectors for disp measures
    Vector *phat; 	                                 // the reference load vector

    double deltaLambdaStep, currentLambda; // dLambda(i) & current value of lambda  
    int signLastDeltaLambdaStep;           // sign of dLambda(i-1)
    double dLambda1min, dLambda1max;       // min & max values for dlambda1 at step (i) 
};

#endif


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
                                                                        
                                                                        
                                                                        
// File: ~/analysis/integrator/DisplacementPath.h
// 
// Written: JZhong 
// Created: 01/04

//
// Description: This file contains the class definition for DisplacementPath.
// DisplacementPath is an algorithmic class for performing a static analysis
// using the user-defined displacement path.
// Revised based on DisplacementControl
// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.

#ifndef DisplacementPath_h
#define DisplacementPath_h

#include <StaticIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class Domain;

class DisplacementPath : public StaticIntegrator
{
  public:
    DisplacementPath(int node, int dof, Vector &incrementVector, Domain *Domain); // change from *theDomain
	
    ~DisplacementPath();

    int newStep(void);    
    int update(const Vector &deltaU);
    int domainChanged(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
    
  private:
    int theNode;          // the node that is being followed
    int theDof;           // the dof at the node being followed
    Vector *theIncrementVector; // the vector containing the increment for each step
    Domain *theDomain;    // the domain containing the node being followed
    int theDofID;         // the system level id of the dof being followed
    
    Vector *deltaUhat, *deltaUbar, *deltaU, *deltaUstep;
    Vector *phat;                           // the reference load vector
    double deltaLambdaStep, currentLambda;  // dLambda(i) & current value of lambda  

	double theCurrentIncrement;  // deltaU at step (i)
    int currentStep;
};

#endif


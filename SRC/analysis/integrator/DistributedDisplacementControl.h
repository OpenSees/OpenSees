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
// $Date: 2008-04-11 23:37:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/DistributedDisplacementControl.h,v $

// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for DistributedDisplacementControl.
// DistributedDisplacementControl is an algorithmic class for performing a static analysis
// using the arc length scheme, that is within a load step the following
// constraint is enforced: 
//  i=1        delta U^T delta U + alpha^2 delta lambda^2 = delta s^2
//  i>1        dU^T delta U + alpha^2 dLambda delta lambda = 0
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and DistributedDisplacementControl is a control parameter.
//
// What: "@(#) DistributedDisplacementControl.h, revA"

#ifndef DistributedDisplacementControl_h
#define DistributedDisplacementControl_h

#include <StaticIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class Domain;

class DistributedDisplacementControl : public StaticIntegrator
{
  public:
    DistributedDisplacementControl(int node, int dof, double increment,
				   int numIncrStep, double minIncrement, double maxIncrement);
    DistributedDisplacementControl();

    ~DistributedDisplacementControl();

    int newStep(void);    
    int update(const Vector &deltaU);
    int domainChanged(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    

    int setProcessID(int processTag);
    int setChannels(int numChannels, Channel **theChannels);
    
  protected:
    
  private:
    int processID;         // processID
    Channel **theChannels; // Channel array
    int numChannels;       // numChannels in theChannel array

    int theNode;          // the node that is being followed
    int theDof;           // the dof at the node being followed
    double theIncrement;  // deltaU at step (i)
    int theDofID;         // the system level id of the dof being followed
    bool allHaveDofID;
    
    Vector *deltaUhat, *deltaUbar, *deltaU, *deltaUstep;
    Vector *phat;                           // the reference load vector
    double deltaLambdaStep, currentLambda;  // dLambda(i) & current value of lambda  

    double specNumIncrStep, numIncrLastStep; // Jd & J(i-1) 
    double minIncrement, maxIncrement; // min/max values of deltaU at (i)
};

#endif


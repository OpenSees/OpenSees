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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/DisplacementControl.h,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/DisplacementControl.h
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for DisplacementControl.
// DisplacementControl is an algorithmic class for perfroming a static analysis
// using the arc length scheme, that is within a load step the follwing
// constraint is enforced: 
//  i=1        delta U^T delta U + alpha^2 delta lambda^2 = delta s^2
//  i>1        dU^T delta U + alpha^2 dLambda delta lambda = 0
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and DisplacementControl is a control parameter.
//
// What: "@(#) DisplacementControl.h, revA"

#ifndef DisplacementControl_h
#define DisplacementControl_h

#include <StaticIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class Domain;

class DisplacementControl : public StaticIntegrator
{
  public:
    DisplacementControl(int node, int dof, double increment, Domain *theDomain,
			int numIncrStep, double minIncrement, double maxIncrement);

    ~DisplacementControl();

    int newStep(void);    
    int update(const Vector &deltaU);
    int domainChanged(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(ostream &s, int flag =0);    
    
  protected:
    
  private:
    int theNode;          // the node that is being followed
    int theDof;           // the dof at the node being followed
    double theIncrement;  // deltaU at step (i)
    Domain *theDomain;    // the domain containg the noe being followed
    int theDofID;         // the syste level id of the dof being followed
    
    Vector *deltaUhat, *deltaUbar, *deltaU, *deltaUstep;
    Vector *phat;                           // the reference load vector
    double deltaLambdaStep, currentLambda;  // dLambda(i) & current value of lambda  

    double specNumIncrStep, numIncrLastStep; // Jd & J(i-1) 
    double minIncrement, maxIncrement; // min/max values of deltaU at (i)
};

#endif


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

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/PerformanceFunctionCoeff.h,v $
                                                                        
#ifndef PerformanceFunctionCoeff_h
#define PerformanceFunctionCoeff_h

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class interface for Node.
// A Node provides the abstraction for a node in the FEM.
// Nodes have original position, trial displacement, velocity and 
// acceleration and committed displacement, velocity and acceleration 
// (the last committed() trial quantities).
//
// What: "@(#) Node.h, revA"

#include <DomainComponent.h>

class Channel;
class FEM_ObjectBroker;
class OPS_Stream;

class PerformanceFunctionCoeff : public DomainComponent
{
  public:
    // constructors
    PerformanceFunctionCoeff
		(int Tag, int NodeID, int Direction, double Coefficient);
    ~PerformanceFunctionCoeff();

    // public methods dealing with the DOF at the node
    int  getTag(void) const;    
    int  getNodeID(void) const;    
    int  getDirection(void) const;    
    double  getCoefficient(void) const;    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

    
  private:
    int NodeID;
    int Direction;
    int Tag; 
	double Coefficient;
};

#endif

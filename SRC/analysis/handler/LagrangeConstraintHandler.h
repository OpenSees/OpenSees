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
// $Date: 2005-11-28 21:37:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/handler/LagrangeConstraintHandler.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: May 1998.
// Revision: A
//
// Description: This file contains the class definition for 
// LagrangeConstraintHandler. LagrangeConstraintHandler is a 
// constraint handler for handling constraints using the Lagrange method.
// for each element and degree-of-freedom at a node it constructs regular
// FE_Element and DOF_Groups; for each SP_Constraint and MP_Constraint
// LagrangeSP_FE and LagrangeMP_FE elements are created.
//
// What: "@(#) LagrangeConstraintHandler.h, revA"

#ifndef LagrangeConstraintHandler_h
#define LagrangeConstraintHandler_h

#include <ConstraintHandler.h>

class FE_Element;
class DOF_Group;

class LagrangeConstraintHandler : public ConstraintHandler
{
  public:
    LagrangeConstraintHandler(double alphaSP, double alphaMP);
    ~LagrangeConstraintHandler();

    int handle(const ID *nodesNumberedLast =0);
    void clearAll(void);    

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    
  protected:
    
  private:
    double alphaSP;
    double alphaMP;
};

#endif





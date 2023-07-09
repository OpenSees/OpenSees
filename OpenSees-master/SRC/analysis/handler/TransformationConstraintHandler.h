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
                                                                        
// $Revision: 1.4 $
// $Date: 2005-11-28 21:35:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/handler/TransformationConstraintHandler.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: May 1998.
// Revision: A
//
// Description: This file contains the class definition for 
// TransformationConstraintHandler. TransformationConstraintHandler is a 
// constraint handler for handling constraints using the Transformation method.
// for each element and degree-of-freedom at a node it constructs regular
// FE_Element and DOF_Groups if there is no SP_Constraint or MP_Constraint
// constraining an elements node or the node; otherwise a TransformationFE
// element and a TransformationDOF_Group are created. 
//
// What: "@(#) TransformationConstraintHandler.h, revA"

#ifndef TransformationConstraintHandler_h
#define TransformationConstraintHandler_h

#include <ConstraintHandler.h>

class FE_Element;
class DOF_Group;

class TransformationConstraintHandler : public ConstraintHandler
{
  public:
    TransformationConstraintHandler();
    ~TransformationConstraintHandler();

    int handle(const ID *nodesNumberedLast =0);
    int applyLoad();
    void clearAll(void);    
    int enforceSPs(void);    
    int doneNumberingDOF(void);        

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

  protected:
    
  private:
    FE_Element 	**theFEs;
    DOF_Group 	**theDOFs;

    int 	numFE;
    int 	numDOF;
    Domain  *theDomain;
    int numConstrainedNodes;

    int numTransformationFEs;
};

#endif





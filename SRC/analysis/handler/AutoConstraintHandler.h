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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/handler/AutoConstraintHandler.h,v $
                                                                        
                                                                        
// Written: Massimo Petracca
// Created: June 2024.
// Revision: A
//
// Description: This file contains the class definition for 
// AutoConstraintHandler. AutoConstraintHandler is a 
// constraint handler for handling constraints using the different methods.
//
// For each element and degree-of-freedom at a node it constructs:
//
// 1) regular FE_Element and DOF_Groups if there is no SP_Constraint or MP_Constraint;
// 2) TransformationDOF_Group for SP constraints (as in the Transformation method)
// 3) PenaltyMP_FE for MP constraints
//
// Notes:
// 1) For each PenaltyMP_FE, by default it automatically selects a proper penalty
//    value based on the stiffness values found on the DOFs involved in the
//    MP constraint
//
// What: "@(#) AutoConstraintHandler.h, revA"

#ifndef AutoConstraintHandler_h
#define AutoConstraintHandler_h

#include <ConstraintHandler.h>
#include <vector>

class FE_Element;
class TransformationDOF_Group;

class AutoConstraintHandler : public ConstraintHandler
{
public:
    AutoConstraintHandler();
    AutoConstraintHandler(
        bool _verbose,
        bool _auto_penalty,
        double _auto_penalty_oom,
        double _user_penalty);
    ~AutoConstraintHandler();

    int handle(const ID* nodesNumberedLast = 0);
    int applyLoad();
    void clearAll(void);
    int doneNumberingDOF(void);

    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel,
        FEM_ObjectBroker& theBroker);

private:
    bool verbose = false;
    bool auto_penalty = true;
    double auto_penalty_oom = 3.0;
    double user_penalty = 0.0;
    std::vector<TransformationDOF_Group*> theDOFs;
};

#endif





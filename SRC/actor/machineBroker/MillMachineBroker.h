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
// $Date: 2003-08-29 07:19:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/machineBroker/MillMachineBroker.h,v $
                                                                        
                                                                        
// File: ~/actor/machineBroker/MillMachineBroker.h
//
// Written: fmk
// Created: Mar 1998
// Revision: A
//
// Purpose: This file contains the class definition for MillMachineBroker.
// MillMachineBroker is a machine broker for use with the cluster of
// millenium machines in Soda Hall.
//
// What: "@(#) MillMachineBroker.h, revA"

#ifndef MillMachineBroker_h
#define MillMachineBroker_h

#include <MachineBroker.h>

class MillMachineBroker : public MachineBroker
{
  public:
    MillMachineBroker(FEM_ObjectBroker *theBroker);
    virtual ~MillMachineBroker();

    virtual int startActor(char *actorProgram, 
			   Channel &theChannel,
			   int compDemand =0);

  protected:
    
  private:
    int currentMachine;
    int maxNumMachines;
    char **machines;
};

#endif

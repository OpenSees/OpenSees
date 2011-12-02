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
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/machineBroker/AlphaMachineBroker.h,v $
                                                                        
                                                                        
// File: ~/actor/broker/AlphaMachineBroker.h
//
// Written: fmk
// Created: Fri Sept 20 12:27:47: 1996
// Revision: A
//
// Purpose: This file contains the class definition for AlphaMachineBroker.
// AlphaMachineBroker is a machine broker for use with the cluster of
// alpha machines in 509 Davis Hall.
//
// What: "@(#) AlphaMachineBroker.h, revA"

#ifndef AlphaMachineBroker_h
#define AlphaMachineBroker_h

#include <MachineBroker.h>

class AlphaMachineBroker : public MachineBroker
{
  public:
    AlphaMachineBroker();
    virtual ~AlphaMachineBroker();

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

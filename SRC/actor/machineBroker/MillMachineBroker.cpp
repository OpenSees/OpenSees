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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-08-29 07:19:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/machineBroker/MillMachineBroker.cpp,v $
                                                                        
                                                                        
// File: ~/actor/broker/MillMachineBroker.h
//
// Written: fmk
// Created: March 1998
// Revision: A
//
// Purpose: This file contains the class definition for MillMachineBroker.
// MillMachineBroker is a machine broker for use with the cluster of
// Millennium machines in Sida Hall.
//
// What: "@(#) MillMachineBroker.h, revA"

#include <MillMachineBroker.h>
#include <stdlib.h>

#include <string.h>
#include <remote.h>
#include <Channel.h>

MillMachineBroker::MillMachineBroker(FEM_ObjectBroker *theBroker)
  :MachineBroker(theBroker), currentMachine(0),maxNumMachines(5)
{
    char *mill0 = "mill0";
    char *mill1 = "mill1";
    char *mill2 = "mill2";
    char *mill3 = "mill3";
    char *mill4 = "mill4";
    
    char **theMachines = (char **)malloc(5*sizeof(char *));
    theMachines[0] = mill0;
    theMachines[1] = mill1;
    theMachines[2] = mill2;
    theMachines[3] = mill3;
    theMachines[4] = mill4;
    
    machines = theMachines;
}
MillMachineBroker::~MillMachineBroker()
{
}


int 
MillMachineBroker::startActor(char *actorProgram, 
			       Channel &theChannel,
			       int compDemand)
{ 
    char  remotecmd[400];

    // get the next machine, a round-robin approach
    char *machine;
    if (currentMachine < maxNumMachines)
	machine = machines[currentMachine];
    else {
	currentMachine = 0;
	machine = machines[currentMachine];
    }

    currentMachine++;
    
    strcpy(remotecmd,REMOTE);
    strcat(remotecmd," ");          
    strcat(remotecmd,machine);
    strcat(remotecmd," ");
    strcat(remotecmd,actorProgram);
    strcat(remotecmd," ");
    strcat(remotecmd,theChannel.addToProgram());    
    strcat(remotecmd,"\n");

    // opserr << "MillMachineBroker::Constructor - command\n"<< remotecmd;
    system(remotecmd);

    return 0;
}


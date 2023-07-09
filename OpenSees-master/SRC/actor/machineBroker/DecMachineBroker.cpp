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
// $Source: /usr/local/cvs/OpenSees/SRC/actor/machineBroker/DecMachineBroker.cpp,v $
                                                                        
                                                                        
// File: ~/actor/broker/DecMachineBroker.h
//
// Written: fmk
// Created: Fri Sept 20 12:27:47: 1996
// Revision: A
//
// Purpose: This file contains the class definition for DecMachineBroker.
// DecMachineBroker is a machine broker for use with the cluster of
// dec machines in 509 Davis Hall.
//
// What: "@(#) DecMachineBroker.h, revA"

#include <DecMachineBroker.h>
#include <stdlib.h>

#include <string.h>
#include <remote.h>
#include <Channel.h>


DecMachineBroker::DecMachineBroker(FEM_ObjectBroker *theBroker)
:MachineBroker(theBroker), currentMachine(0), maxNumMachines(9)
{
    char *dec1 = "dec-1";
    char *dec2 = "dec-2";
    char *dec3 = "dec-3";
    char *dec4 = "dec-4";
    char *dec5 = "dec-5";
    char *dec6 = "dec-6";    
    char *dec7 = "dec-7";    
    char *dec8 = "dec-8";        
    char *dec9 = "dec-9";            
    
    char **theMachines = (char **)malloc(9*sizeof(char *));
    theMachines[0] = dec9;
    theMachines[1] = dec8;
    theMachines[2] = dec7;
    theMachines[3] = dec6;
    theMachines[4] = dec5;
    theMachines[5] = dec4;
    theMachines[6] = dec3;    
    theMachines[7] = dec2;    
    theMachines[8] = dec1;        
    
    machines = theMachines;
}
DecMachineBroker::~DecMachineBroker()
{
}


int 
DecMachineBroker::startActor(char *actorProgram, 
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

    system(remotecmd);

    return 0;
}


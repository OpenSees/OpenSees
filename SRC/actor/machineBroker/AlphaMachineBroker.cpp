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
// $Source: /usr/local/cvs/OpenSees/SRC/actor/machineBroker/AlphaMachineBroker.cpp,v $
                                                                        
                                                                        
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

#include <AlphaMachineBroker.h>
#include <stdlib.h>

#include <string.h>
#include <remote.h>
#include <Channel.h>

AlphaMachineBroker::AlphaMachineBroker(FEM_ObjectBroker *theBroker)
  :MachineBroker(theBroker), currentMachine(0),maxNumMachines(8)
{
    char *alpha1 = "alpha-1";
    char *alpha2 = "alpha-2";
    char *alpha3 = "alpha-3";
    char *alpha4 = "alpha-4";
    char *alpha5 = "alpha-5";
    char *alpha6 = "alpha-6";    
    char *alpha7 = "alpha-7";    
    char *alpha8 = "alpha-8";        
    
    char **theMachines = (char **)malloc(8*sizeof(char *));
    theMachines[0] = alpha5;
    theMachines[1] = alpha2;
    theMachines[2] = alpha8;
    theMachines[3] = alpha1;
    theMachines[4] = alpha4;
    theMachines[5] = alpha6;
    theMachines[6] = alpha3;    
    theMachines[7] = alpha7;    
    
    machines = theMachines;
}
AlphaMachineBroker::~AlphaMachineBroker()
{
}


int 
AlphaMachineBroker::startActor(char *actorProgram, 
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
    //opserr << "CurrentMachine : " << machine <<  endln;
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


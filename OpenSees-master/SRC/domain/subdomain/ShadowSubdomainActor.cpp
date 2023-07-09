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
// $Date: 2003-02-14 23:01:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/subdomain/ShadowSubdomainActor.cpp,v $
                                                                        
                                                                        
// Written: fmk
// Created: 12/96
// Revision: A
//
// Purpose: This file contains the code for a node process. 
// What: "@(#) ShadowSubdomainActor.C, revA"

#include <ActorSubdomain.h>
#include <FEM_ObjectBroker.h>
#include <TCP_Socket.h>
// #include <TCP_SocketNoDelay.h>
#include <UDP_Socket.h>
#include <SocketAddress.h>
#include <Message.h>

int main(int argv, char **argc)
{
    int channelType = atoi(argc[1]);
    

    //
    // create a channel
    //
    
    Channel *theChannel = 0;
    if (channelType == 1) {
	char *machine = argc[2];    	
	int port = atoi(argc[3]);	
	theChannel = new TCP_Socket(port,machine);
    }
    //    else if (channelType == 2) {
    //	char *machine = argc[2];    	
    //	int port = atoi(argc[3]);	    
    //	theChannel = new UDP_Socket(port,machine);
    //    {    
    // else if (channelType == 3) {
    // char *machine = argc[2];    	
    // int port = atoi(argc[3]);	    
    // theChannel = new TCP_SocketNoDelay(port,machine);
    // }   
    else {
	opserr << "ACTOR PROCESS: ShadowSubdomainActor Exiting ";	
	opserr << "- invalid channel type " << channelType << endln;
	exit(-1);
    }

    if (theChannel == 0) {
	opserr << "ACTOR PROCESS: ShadowSubdomainActor Exiting ";	
	opserr << "- could not create the channel " << endln;
	exit(-1);
    }	

    //
    // create an object broker
    //
    
    FEM_ObjectBroker theObjectBroker;

    //
    // create the actor object and start it running
    //
    
    ActorSubdomain theActor(*theChannel,theObjectBroker);
    theActor.run();

    // exit normally
    opserr << "ShadowSubdomainActor:: ACTOR PROCESS EXITING\n\n";
    // opserr << theActor;
    exit(0);
}



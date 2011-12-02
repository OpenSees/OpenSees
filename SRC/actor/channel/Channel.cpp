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
// $Date: 2002-06-18 01:05:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/Channel.cpp,v $
                                                                        
                                                                        
// File: ~/actor/channel/Channel.C
//
// Written: fmk
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementation of Channel.
//
// What: "@(#) Channel.C, revA"

#include <Channel.h>
#include <Message.h>
#include <MovableObject.h>
#include <FEM_ObjectBroker.h>

Channel::Channel ()
{

}

Channel::~Channel()
{
    
}    

int
Channel::isDatastore(void)
{
  return 0;
}

int
Channel::getDbTag(void)
{
  return 0;
}


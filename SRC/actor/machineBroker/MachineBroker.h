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
// $Source: /usr/local/cvs/OpenSees/SRC/actor/machineBroker/MachineBroker.h,v $
                                                                        
                                                                        
// File: ~/actor/broker/MachineBroker.h
//
// Written: fmk
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class definition for MachineBroker.
// MachineBroker is an abstract base class, a subclass of which must
// be written for each parallel machine. A MachineBroker is responsible
// for getting an actor process running on the parallel machine.
//
// What: "@(#) MachineBroker.h, revA"

#ifndef MachineBroker_h
#define MachineBroker_h

class Channel;

class MachineBroker
{
  public:
    MachineBroker() {};
    virtual ~MachineBroker() {};

    virtual int startActor(char *actorProgram, 
			   Channel &theChannel,
			   int compDemand =0) =0;

  protected:
    
  private:

};

#endif

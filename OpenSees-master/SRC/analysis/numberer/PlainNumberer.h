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
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/numberer/PlainNumberer.h,v $
                                                                        
                                                                        
// File: ~/analysis/numberer/PlainNumberer.h
// 
// Written: fmk 
// Created: 9/96
// Revision: A
//
// Description: This file contains the class definition for PlainNumberer.
// PlainNumberer is a subclass of DOF_Numberer. The PlainNumberer assigns
// equation numbers to the DOFs on a first come first serve basis; that is 
// it gets the DOF_GrpIter and assigns the equation numbers to the DOFs
// as it iterates through the iter.
//
// What: "@(#) PlainNumberer.h, revA"

#ifndef PlainNumberer_h
#define PlainNumberer_h

#define START_EQN_NUMBER 0

#include <DOF_Numberer.h>

class PlainNumberer: public DOF_Numberer
{
  public:
    PlainNumberer();
    ~PlainNumberer();

    int numberDOF(int lastDOF = -1);
    int numberDOF(ID &lastDOFs);    

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

  protected:

  private:
    
};

#endif


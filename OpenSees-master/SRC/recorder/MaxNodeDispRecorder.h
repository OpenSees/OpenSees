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
// $Date: 2004-11-24 22:45:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/MaxNodeDispRecorder.h,v $
                                                                        
                                                                        
#ifndef MaxNodeDispRecorder_h
#define MaxNodeDispRecorder_h

// File: ~/recorder/MaxNodeDispRecorder.h
//
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for MaxNodeDispRecorder.
// A MaxNodeDispRecorder is used to determine the max nodal displacement
// at a collection of nodes over an analysis. (between commitTag of 0 and
// last commitTag).
//
// What: "@(#) MaxNodeDispRecorder.h, revA"

#include <Recorder.h>
#include <ID.h>
#include <Vector.h>

class Domain;

class MaxNodeDispRecorder: public Recorder
{
  public:
    MaxNodeDispRecorder(int dof, const ID &theNodes, Domain &theDomain);
    ~MaxNodeDispRecorder();
    int record(int commitTag, double timeStamp);
    int playback(int commitTag);

    int restart(void);    
    
  protected:
    
  private:	
    ID theNodes;
    Vector maxDisp;
    int dof;
    Domain *theDomain;
};


#endif

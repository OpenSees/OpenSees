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
// $Date: 2003-02-14 23:01:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/MaxNodeDispRecorder.cpp,v $
                                                                        
                                                                        
// File: ~/recorder/MaxNodeDispRecorder.C
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
// What: "@(#) MaxNodeDispRecorder.C, revA"

#include <MaxNodeDispRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <Vector.h>
#include <ID.h>

MaxNodeDispRecorder::MaxNodeDispRecorder(int theDof, 
					 const ID &nodes, 
					 Domain &theDom)
:theNodes(nodes), maxDisp(nodes.Size()), 
 dof(theDof), theDomain(&theDom)
{
    if (dof < 0) dof = 0;
}

MaxNodeDispRecorder::~MaxNodeDispRecorder()
{
    
}

int 
MaxNodeDispRecorder::record(int commitTag, double timeStamp)
{
    for (int i=0; i<theNodes.Size(); i++) {
	Node *theNode = theDomain->getNode(theNodes(i));
	if (theNode != 0) {
	    const Vector &theDisp = theNode->getTrialDisp();
	    if (theDisp.Size() > dof) {
		double disp = theDisp(dof);
		if (disp > 0 && disp > maxDisp(i))
		    maxDisp(i) = disp;
		else if (disp < 0 && -disp > maxDisp(i))
		    maxDisp(i) = -disp;
	    }
	}
    }
    return 0;
}


int 
MaxNodeDispRecorder::playback(int commitTag)
{
    opserr << "Max Recorded Displacement: " << maxDisp << endln;
    return 0;
}


void
MaxNodeDispRecorder::restart(void)
{
    maxDisp.Zero();
}


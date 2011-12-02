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

// $Revision: 1.5 $
// $Date: 2003-04-02 22:02:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/DriftRecorder.cpp,v $

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for DriftRecorder.

#include <math.h>

#include <DriftRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>

#include <iomanip>
using std::ios;

#include <string.h>

DriftRecorder::DriftRecorder(int ni, int nj, int df, int dirn,
			     Domain &theDom, const char *fileName, int sflag)
  :ndI(ni), ndJ(nj), dof(df), perpDirn(dirn), oneOverL(0.0),
   theDomain(&theDom), flag(sflag)
 
{
  theFileName = new char[strlen(fileName)+1];
  if (theFileName == 0) {
    opserr << "DriftRecorder::DriftRecorder -- out of memory copying fileName " << endln;
    exit(-1);
  }
  
  strcpy(theFileName, fileName);    
  
  theFile.open(fileName, ios::out);
  if (theFile.bad()) {
    opserr << "WARNING - FileNodeDispRecorder::FileNodeDispRecorder()";
    opserr << " - could not open file " << fileName << endln;
  }    
  
  Node *nodeI = theDomain->getNode(ndI);
  Node *nodeJ = theDomain->getNode(ndJ);

  const Vector &crdI = nodeI->getCrds();
  const Vector &crdJ = nodeJ->getCrds();

  if (crdI(dirn) == crdJ(dirn)) {
    opserr << "DriftRecorder::DriftRecorder-- Nodal projection has zero component along chosen direction\n";
    oneOverL = 0.0;
  }
  else 
    oneOverL = 1.0/fabs(crdJ(dirn) - crdI(dirn));
}

DriftRecorder::~DriftRecorder()
{
  if (!theFile)
    theFile.close();
  if (theFileName != 0)
    delete [] theFileName;
}

int 
DriftRecorder::record(int commitTag, double timeStamp)
{
  Node *nodeI = theDomain->getNode(ndI);
  Node *nodeJ = theDomain->getNode(ndJ);

  const Vector &dispI = nodeI->getTrialDisp();
  const Vector &dispJ = nodeJ->getTrialDisp();

  double dx = dispJ(dof)-dispI(dof);

  if (flag == 1)
    theFile << theDomain->getCurrentTime() << " ";
  else if (flag == 2)
    theFile << theDomain->getCurrentTime() << " ";
  
  theFile << dx*oneOverL;
  
  theFile << endln;
  theFile.flush();
  
  return 0;
}

int 
DriftRecorder::playback(int commitTag)
{
  opserr << "WARNING -- DriftRecorder::playback() -- not implemented" << endln;

  return 0;
}

void
DriftRecorder::restart(void)
{
  theFile.close();
  theFile.open(theFileName, ios::out);
  if (theFile.bad()) {
    opserr << "WARNING - DriftRecorder::restart() - could not open file ";
    opserr << theFileName << endln;
  }
}

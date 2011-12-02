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
// $Date: 2000-09-15 08:23:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/NodeRecorder.cpp,v $
                                                                        

// File: ~/recorder/NodeRecorder.C
//
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for NodeRecorder.
// A NodeRecorder is used to record the specified dof responses 
// at a collection of nodes over an analysis. (between commitTag of 0 and
// last commitTag).
//
// What: "@(#) NodeRecorder.C, revA"

#include <NodeRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <Vector.h>
#include <ID.h>

#include <iostream.h>
#include <fstream.h>
#include <string.h>

NodeRecorder::NodeRecorder(const ID &dofs, 
			   const ID &nodes, 
			   Domain &theDom,
			   char *fileName,
			   char *dataToStore,
			   int startFlag)
:theDofs(dofs), theNodes(nodes), disp(nodes.Size()*dofs.Size()), 
 theDomain(&theDom), flag(startFlag)
{
  if (strlen(fileName) > MAX_FILENAMELENGTH) 
    g3ErrorHandler->fatal("FileNodeDispRecorder::FileNodeDispREcorder - fileName too long %d max\n",
		     MAX_FILENAMELENGTH);

  strcpy(theFileName, fileName);    
    
  theFile.open(fileName, ios::out);
  if (!theFile) {
    cerr << "WARNING - FileNodeDispRecorder::FileNodeDispRecorder()";
    cerr << " - could not open file " << fileName << endl;
  }    

  disp.Zero();

  if ((strcmp(dataToStore, "disp") == 0)) {
    dataFlag = 0;
  } else if ((strcmp(dataToStore, "vel") == 0)) {
    dataFlag = 1;
  } else if ((strcmp(dataToStore, "accel") == 0)) {
    dataFlag = 2;
  } else if ((strcmp(dataToStore, "incrDisp") == 0)) {
    dataFlag = 3;
  } else if ((strcmp(dataToStore, "incrDeltaDisp") == 0)) {
    dataFlag = 4;
  } else {
    dataFlag = 6;
    cerr << "NodeRecorder::NodeRecorder - dataToStore " << *dataToStore;
    cerr << "not recognized (disp, vel, accel, incrDisp, incrDeltaDisp)\n";
  }
}

NodeRecorder::~NodeRecorder()
{
    if (!theFile)
	theFile.close();
}

int 
NodeRecorder::record(int commitTag)
{
    // now we go get the displacements from the nodes
    int numDOF = theDofs.Size();
    int numNodes = theNodes.Size();

    for (int i=0; i<numNodes; i++) {
	int cnt = i*numDOF;
	Node *theNode = theDomain->getNode(theNodes(i));
	if (theNode != 0) {
	  if (dataFlag == 0) {
	    const Vector &theDisp = theNode->getTrialDisp();
	    for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		if (theDisp.Size() > dof) {
		    disp(cnt) = theDisp(dof);
		}
		cnt++;
	    }
	  } else if (dataFlag == 1) {
	    const Vector &theDisp = theNode->getTrialVel();
	    for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		if (theDisp.Size() > dof) {
		    disp(cnt) = theDisp(dof);
		}
		cnt++;
	    }
	  } else if (dataFlag == 2) {
	    const Vector &theDisp = theNode->getTrialAccel();
	    for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		if (theDisp.Size() > dof) {
		    disp(cnt) = theDisp(dof);
		}
		cnt++;
	    }
	  } else if (dataFlag == 3) {
	    const Vector &theDisp = theNode->getIncrDisp();
	    for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		if (theDisp.Size() > dof) {
		    disp(cnt) = theDisp(dof);
		}
		cnt++;
	    }
	  } else if (dataFlag == 4) {
	    const Vector &theDisp = theNode->getIncrDeltaDisp();
	    for (int j=0; j<numDOF; j++) {
		int dof = theDofs(j);
		if (theDisp.Size() > dof) {
		    disp(cnt) = theDisp(dof);
		}
		cnt++;
	    }
	  }
	}
    }
    
    // now we write them to the file
    if (flag == 1)
	theFile << theDomain->getCurrentTime() << " ";
    else if (flag == 2)
	theFile << theDomain->getCurrentTime() << " ";


    for (int j=0; j<numNodes*numDOF; j++)
	theFile << disp(j) << " ";
    
    theFile << endl;
    
    return 0;
}

int 
NodeRecorder::playback(int commitTag)
{
  if (!theFile)
    return 0;
    
  // close o/p file to ensure all buffered data gets written to file
  theFile.close(); 

  int numDOF = theDofs.Size();
  int numNodes = theNodes.Size();  
  
  // open a stream for reading from the file
  ifstream inputFile;
  inputFile.open(theFileName, ios::in);
  if (!inputFile) {
    cerr << "WARNING - FileNodeDispRecorder::playback() - could not open file ";
    cerr << theFileName << endl;
    return -1;
  }   

  double data;
  // read file up until line we want
  for (int i=0; i<(commitTag-1); i++)
    // now read in a line
    if (flag == 1 || flag == 2) {
	inputFile >> data;
	for (int i=0; i<numNodes*numDOF; i++)
	    inputFile >> data;
    }

  // now read in our line and print out
  if (flag == 1 || flag == 2) {
      inputFile >> data;
      cerr << data << " ";
      for (int i=0; i<numNodes*numDOF; i++) {
	  inputFile >> data;
	  cerr << data << " ";
      }	
      cerr << endl;
  }
  inputFile.close();
      
    // open file again for writing
    theFile.open(theFileName, ios::app);
    if (!theFile) {
      cerr << "WARNING - FileNodeDispRecorder::playback() - could not open file ";
      cerr << theFileName << endl;
      return -1;
    }    
  
  // does nothing
  return 0;
}

void
NodeRecorder::restart(void)
{
  theFile.close();
  theFile.open(theFileName, ios::out);
  if (!theFile) {
    cerr << "WARNING - FileNodeDispRecorder::restart() - could not open file ";
    cerr << theFileName << endl;
  }
}






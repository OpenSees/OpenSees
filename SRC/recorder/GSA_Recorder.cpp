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
// $Date: 2003-02-25 23:34:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/GSA_Recorder.cpp,v $

// Written: fmk 
// Created: 02/03
//
// What: "@(#) GSA_Recorder.C, revA"

#include <GSA_Recorder.h>
#include <Domain.h>
#include <Node.h>
#include <Element.h>
#include <SP_Constraint.h>
#include <NodeIter.h>
#include <ElementIter.h>
#include <SP_ConstraintIter.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>

GSA_Recorder::GSA_Recorder(Domain &theDom, 
			   const char *fileName, 
			   const char *title1,
			   const char *title2,
			   const char *title3,
			   const char *jobno,
			   const char *initials,
			   const char *spec,
			   const char *currency,
			   const char *length,
			   const char *force,
			   const char *temp,
			   double dT)
: theDomain(&theDom), ndm(3), ndf(6), counter(0), deltaT(dT), nextTimeStampToRecord(0.0)
{
  // open file 
  if (theFile.setFile(fileName, OVERWRITE) < 0) {
    opserr << "WARNING - GSA_Recorder::GSA_Recorder()";
    opserr << " - could not open file " << fileName << endln;
    exit(-1);
  } 

  // spit out header data
  if (title1 != 0)
    theFile << "TITLE\t" << title1;
  else
    theFile << "TITLE\t" << "No Title";
  theFile << "\t" << title2;
  theFile << "\t" << title3;
  theFile << "\t" << jobno << "\t" << initials << endln;
  theFile << "SPEC\t" << spec << endln;
  theFile << "CURRENCY\t" << currency << endln;
  theFile << "UNIT_DATA\tLENGTH\t" << length << endln; 
  theFile << "UNIT_DATA\tFORCE\t" << force << endln;
  theFile << "UNIT_DATA\tTEMP\t" << temp << endln;

  // spit out nodal data
  NodeIter &theNodes = theDomain->getNodes();
  Node *theNode;
  while ((theNode=theNodes()) != 0) {
    int nodeTag = theNode->getTag();
    theFile << "NODE\t" << nodeTag;
    const Vector &crds = theNode->getCrds();
    if (crds.Size() != ndm) {
      opserr << "WARNING - GSA_Recorder::GSA_Recorder() - node: " <<  nodeTag ;
      opserr << " has invalid number of coordinates, expecting: " << ndm << " got: " << crds.Size() << endln;
      exit(-1);
    }
    const Vector &disp = theNode->getTrialDisp();
    if (disp.Size() != ndf) {
      opserr << "WARNING - GSA_Recorder::GSA_Recorder() - node: " <<  nodeTag ;
      opserr << " has invalid number of dof, expecting: " << ndf << " got: " << disp.Size() << endln;
      exit(-1);
    }
    for (int i=0; i<ndm; i++)
      theFile << "\t" << crds(i);
    theFile << endln;
  }
  
  
  // open file and spit out the initial data
  SP_ConstraintIter &theSPs = theDomain->getSPs();
  SP_Constraint *theSP;
  ID theConstrainedNodes(0,6);
  ID theSpMatrix(0, 6*ndf);
  int numNodesWithSP = 0;
  while ((theSP=theSPs()) != 0) {
    int nodeTag =  theSP->getNodeTag();
    int location = theConstrainedNodes.getLocation(nodeTag);
    if (location < 0) {
      theConstrainedNodes[numNodesWithSP] = nodeTag;
      for (int i=0; i<ndf; i++)
	theSpMatrix[numNodesWithSP*ndf+i] = 0;	  
      location = numNodesWithSP++;
    }
    int id = theSP->getDOF_Number();
    theSpMatrix[location*ndf + id] = 1;
  }
  
  for (int j=0; j<numNodesWithSP; j++) {
    theFile << "SPC\t" <<  theConstrainedNodes[j] << "\t0";
    for (int i=0; i<ndf; i++)
      theFile << "\t" << theSpMatrix[j*ndf+i];
    theFile << endln;
  }
  
  ElementIter &theElements = theDomain->getElements();
  Element *theElement;
  while ((theElement=theElements()) != 0) {
    theElement->Print(theFile, -1);
  }

}

GSA_Recorder::~GSA_Recorder()
{
  theFile << "END\n";
  theFile.close();
}

int 
GSA_Recorder::record(int commitTag, double timeStamp)
{
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {
      
    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;

    counter++;

    theFile << "ANAL_CASE\t" << counter << "\tStep" << counter << "\tL" << counter << 
      "\tGSS\tSTATIC\tPOST\t" << counter << "\topensees\t" << "20030204165318	0" << endln;

    theFile << "!\n!RESULTS FOR ANALYSIS CASE\t" << counter << "\n!\n";

    // spit out nodal displacements
    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    while ((theNode=theNodes()) != 0) {
      int nodeTag = theNode->getTag();
      const Vector &disp = theNode->getTrialDisp();
      if (ndm == 3 && ndf == 6) {
	theFile << "DISP\t" << nodeTag << "\t" << counter;
	for (int i=0; i<ndm; i++)
	  theFile << "\t" << disp(i);
	theFile << endln;
	theFile << "ROTN\t" << nodeTag << "\t" << counter;
	for (int j=0; j<ndm; j++)
	  theFile << "\t" << disp(ndm+j);
	theFile << endln;
      }
    }

    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement=theElements()) != 0) {
      theElement->Print(theFile, (counter+1)*-1);  // starts at -2, as already using -1

    }
  }
  
  return 0;
}

int 
GSA_Recorder::playback(int commitTag)
{
  // does nothing
  return 0;
}

void
GSA_Recorder::restart(void)
{

}






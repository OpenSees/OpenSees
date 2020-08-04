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
// $Date: 2009-12-23 22:56:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/RemoveRecorder.h,v $
                                                                        
#ifndef RemoveRecorder_h
#define RemoveRecorder_h
#define RECORDER_TAGS_RemoveRecorder 15


// Written: M Talaat, 06/07
// Created: 06/07
// Revision: M Talaat
//
// Description: This file contains the class definition for RemoveRecorder.
// A RemoveRecorder is used to check the removal criteria  from a node or element
// section/material during the analysis remove the element or node if applicable.
//
// What: "@(#) RemoveRecorder.h, revA"

#include <Recorder.h>
#include <OPS_Stream.h>
#include <FileStream.h>
#include <fstream>
using std::ofstream;

#include <Information.h>
#include <ID.h>
#include <Vector.h>


class Domain;
class Vector;
class Matrix;
class Element;
class Node;
class Response;
class FE_Datastore;
class Response;

class RemoveRecorder: public Recorder
{
  public:
   RemoveRecorder(int nodeID, 
		  ID &eleIDs,
		  ID &secIDs,
		  ID &secondaryEleIDs,
		  const Vector remCriteria,
		  Domain &theDomainPtr,
		  OPS_Stream &s,
		  bool echotimeflag, 
		  double deltat , 
		  const char *filename, 
		  Vector eleMass, 
		  double gAcc, 
		  int gDir, 
		  int gPat, 
		  int nTagbotn =0,
		  int nTagmidn =0,
		  int nTagtopn =0,
		  int globgrav =0,
		  const char *thefileNameinf =0);
   
   ~RemoveRecorder();
   int record(int commitTag, double timeStamp);
   int playback(int commitTag);
   
   int restart(void);    
   // changed
   int checkEleRemoval(Element* theEle, Response *eleResponse, int &theComponent,const Vector &Criteria);
   //	int checkNodeRemoval(Element* theEle, int &theComponent,const Vector Criteria);
   int elimElem(int theDeadEleTag, double timeStamp = 0);
   int elimNode(int theDeadNodeTag, double timeStamp = 0);
   int elimSecondaries(double timeStamp = 0);
   int updateNodalMasses(int theEleTag, double theEleMass);
   
   static int numRecs;
   static ID remEleList;
   static ID remNodeList;
   static int numRemEles;
   static int numRemNodes;
   static Element** remEles;
   static Node** remNodes;
   
 protected:
   
 private:	
   
   // Removable node (optional), number of elements, sections, and rules to check  
   int nodeTag;
   int numEles;
   int numSecs;
   int numRules;  
   
   // Tags for primary elements, sections to check (if applicable), secondary element (optional)
   ID eleTags, secTags, secondaryEleTags;	
   
   Vector criteria;
   bool secondaryFlag; // flag indiacting if secondary elements should be removed if all primary elements collaspe
   Vector eleMasses, eleWeights;
   double gAcc;
   int gDir, gPat;
  
   Domain *theDomain;
   bool echoTimeFlag;   // flag indicating if pseudo time also printed
   static char *fileName;  // file name  
   
   //static FileStream theFile;	// output stream
   static ofstream theFile;
   
   double deltaT;
   double nextTimeStampToRecord;
   
   // new
   static char *fileNameinf;  // file name  
   int nTagbotn, nTagmidn, nTagtopn;
   int globgrav;

   Response **eleResponses;
};


#endif

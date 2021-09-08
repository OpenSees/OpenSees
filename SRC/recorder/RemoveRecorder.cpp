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
// $Date: 2010-02-09 21:29:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/RemoveRecorder.cpp,v $
                                                                        
// Written: M Talaat
// Created: 06/07
// Revision: A
//
// Description: This file contains the class implementatation of ElementRecorder.
//

// What: "@(#) RemoveRecorder.C, revA"

#include <stdlib.h>

#include <stdio.h>
#include <RemoveRecorder.h>
#include <Domain.h>
#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <string.h>
#include <Response.h>
#include <FE_Datastore.h>
//#include <DamageModel.h>
#include <ElementIter.h>
#include <NodeIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <ElementalLoad.h>
#include <ElementalLoadIter.h>
#include <NodalLoad.h>
#include <NodalLoadIter.h>
//#include <G3string.h>

#include <SP_Constraint.h> //Joey UC Davis
#include <SP_ConstraintIter.h> //Joey UC Davis

#include <ElementResponse.h>

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <MeshRegion.h>
#include <DummyStream.h>

#include <iomanip>
using std::ios;

void *OPS_RemoveRecorder() {
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING recorder Collapse -ele eleID <eleID2? ...>  "
                  "-node nodeID <-time> <-file fileName?> ? "
               << "\n or recorder Collapse -ele eleID1 <eleID2? ...>? "
                  "<-sec secID1? secID2? ...> -crit crit1? value1?"
               << " <-crit crit2? value2?> <-time> <-file fileName?> "
                  "<-mass mass1? mass2? ...> <-g gAcc gDir? gPat?>?\n";
        return 0;
    }

    // current maximum number of either-or criteria for an element
    int maxNumCriteria = 2;
    double dT = 0.0;
    bool echoTime = false;
    int endEleIDs = 2;
    int numEle = endEleIDs - 2;
    int loc = endEleIDs;
    int flags = 0;
    int eleData = 0;
    int nodeTag = 0;
    //  new
    int nTagbotn = 0;
    int nTagmidn = 0;
    int nTagtopn = 0;
    int globgrav = 0;
    const char *fileNameinf = 0;
    //  end of new
    int numSecondaryEle = 0;

    // create an ID to hold ele tags
    // ID eleIDs(numEle, numEle+1);
    ID eleIDs;
    eleIDs = ID(1);
    ID secondaryEleIDs = ID(1);
    secondaryEleIDs[0] = 0;
    bool secondaryFlag = false;
    ID secIDs;
    //	secIDs = 0;

    // optional mass and weight definition
    Vector eleMass(1);
    eleMass.Zero();
    double gAcc = 0;
    int gDir = 0, gPat = 0;

    Vector remCriteria(2 * maxNumCriteria);
    remCriteria.Zero();
    int numCrit = 0;
    const char *fileName = 0;

    OPS_Stream *theOutputStream = 0;

    Domain *domain = OPS_GetDomain();
    while (flags == 0 && OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (strcmp(opt, "-node") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING need nodeTag\n";
                return 0;
            }
            int num = 1;
            if (OPS_GetIntInput(&num, &nodeTag) < 0) {
                opserr << "WARNING recorder Collapse -node - invalid node "
                          "tag\n";
                return 0;
            }

            Node *theNode = domain->getNode(nodeTag);
            if (theNode == 0) {
                opserr
                    << "WARNING recorder Collapse -node - invalid node \n";
                return 0;
            }

        } else if (strcmp(opt, "-file_infill") == 0 &&
                   OPS_GetNumRemainingInputArgs() > 0) {
            fileNameinf = OPS_GetString();

        } else if (strcmp(opt, "-checknodes") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "WARNING need nTagbotn nTagmidn, nTagtopn\n";
                return 0;
            }
            int data[3], num = 3;
            if (OPS_GetIntInput(&num, &data[0]) < 0) {
                opserr << "WARNING recorder Collapse -node - invalid "
                          "nTagbotn "
                          "nTagmidn, or nTagtopn \n";
                return 0;
            }
            nTagbotn = data[0];
            nTagmidn = data[1];
            nTagtopn = data[2];

        } else if (strcmp(opt, "-global_gravaxis") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "WARNING need globgrav\n";
                return 0;
            }
            int num = 1;
            if (OPS_GetIntInput(&num, &globgrav) < 0) {
                opserr << "WARNING recorder Collapse -global_gravaxis - "
                          "invalid global axis for gravity \n";
                return 0;
            }

        } else if ((strcmp(opt, "-slave") == 0) ||
                   (strcmp(opt, "-secondary") == 0)) {
            secondaryFlag = true;

        } else if ((strcmp(opt, "-ele") == 0) ||
                   (strcmp(opt, "-eles") == 0) ||
                   (strcmp(opt, "-element") == 0)) {
            // ensure no segmentation fault if user messes up
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING recorder Collapse .. -ele tag1? .. - "
                          "no ele tags specified\n";
                return 0;
            }

            const char *all = OPS_GetString();
            if (strcmp(all, "all") == 0) {
                ElementIter &theEleIter = domain->getElements();
                Element *theEle;
                while ((theEle = theEleIter()) != 0) {
                    eleIDs[numEle++] = theEle->getTag();
                }
            } else {
                OPS_ResetCurrentInputArg(-1);
                //
                // read in a list of ele until end of command or other flag
                //
                while (OPS_GetNumRemainingInputArgs() > 0) {
                    int num = 1;
                    int eleTag;
                    if (OPS_GetIntInput(&num, &eleTag) < 0) {
                        OPS_ResetCurrentInputArg(-1);
                        break;
                    }
                    if (!secondaryFlag) {
                        eleIDs[numEle++] = eleTag;
                    } else {
                        secondaryEleIDs[numSecondaryEle++] = eleTag;
                    }
                }
                secondaryFlag = false;
            }

        } else if (strcmp(opt, "-eleRange") == 0) {
            // ensure no segmentation fault if user messes up
            if (OPS_GetNumRemainingInputArgs() < 2) {
                opserr << "WARNING recorder Element .. -eleRange start? "
                          "end?  .. - no ele tags specified\n";
                return 0;
            }

            //
            // read in start and end tags of two elements & add set
            // [start,end]
            //
            int startend[2];
            int num = 2;
            if (OPS_GetIntInput(&num, &startend[0]) < 0) {
                opserr << "WARNING recorder Element -eleRange start? end? "
                          "- invalid start or end\n";
                return 0;
            }
            if (startend[0] > startend[1]) {
                int swap = startend[1];
                startend[1] = startend[0];
                startend[0] = swap;
            }
            for (int i = startend[0]; i <= startend[1]; ++i) {
                if (!secondaryFlag) {
                    eleIDs[numEle++] = i;
                } else {
                    secondaryEleIDs[numSecondaryEle++] = i;
                }
            }
            secondaryFlag = false;

        } else if (strcmp(opt, "-region") == 0) {
            // allow user to specif elements via a region

            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING recorder Element .. -region tag?  .. - "
                          "no region specified\n";
                return 0;
            }
            int tag;
            int num = 1;
            if (OPS_GetIntInput(&num, &tag) < 0) {
                opserr << "WARNING recorder Element -region tag? - "
                          "invalid tag \n";
                return 0;
            }
            MeshRegion *theRegion = domain->getRegion(tag);
            if (theRegion == 0) {
                opserr << "WARNING recorder Element -region " << tag
                       << " - region does not exist\n";
                return 0;
            }
            const ID &eleRegion = theRegion->getElements();
            for (int i = 0; i < eleRegion.Size(); i++) {
                if (secondaryFlag == false) {
                    eleIDs[numEle++] = eleRegion(i);
                } else {
                    secondaryEleIDs[numSecondaryEle++] = eleRegion(i);
                }
            }
            secondaryFlag = false;

        } else if ((strcmp(opt, "-time") == 0) ||
                   (strcmp(opt, "-load") == 0)) {
            // allow user to specify const load
            echoTime = true;
        }

        else if (strcmp(opt, "-dT") == 0) {
            // allow user to specify time step size for recording
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: no dT is given\n";
                return 0;
            }
            int num = 1;
            if (OPS_GetDoubleInput(&num, &dT) < 0) {
                opserr << "WARNING: failed to get dT\n";
                return 0;
            }
        }

        else if (strcmp(opt, "-file") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: no fileName is given\n";
                return 0;
            }
            fileName = OPS_GetString();
            //	    simulationInfo.addWriteFile(fileName);
        }

        else if ((strcmp(opt, "-mass") == 0)) {
            if (OPS_GetNumRemainingInputArgs() < numEle) {
                opserr << "WARNING: need numEle mass values\n";
                return 0;
            }
            eleMass.resize(numEle);
            eleMass.Zero();
            if (OPS_GetDoubleInput(&numEle, &eleMass(0)) < 0) {
                opserr << "WARNING: failed to get mass\n";
                return 0;
            }
        }

        else if ((strcmp(opt, "-g") == 0)) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "WARNING: need gAcc, gDir, gPat\n";
                return 0;
            }
            int num = 3;
            double data[3];
            if (OPS_GetDoubleInput(&num, &data[0]) < 0) {
                opserr << "WARNING: failed to read gAcc, gDir, gPat\n";
                return 0;
            }

            gAcc = data[0];
            gDir = data[1];
            gPat = data[2];
        }

        else if ((strcmp(opt, "-section") == 0) ||
                 (strcmp(opt, "-sec") == 0) ||
                 (strcmp(opt, "-comp") == 0)) {
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                int secID;
                if (OPS_GetIntInput(&num, &secID) < 0) {
                    OPS_ResetCurrentInputArg(-1);
                    break;
                }
                secIDs.insert(secID);
            }
        }

        else if (strcmp(opt, "-criteria") == 0 ||
                 strcmp(opt, "-crit") == 0) {
            int critTag = 0;
            double critValue = 0;

            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: need criteria type\n";
                return 0;
            }

            const char *opt1 = OPS_GetString();
            if (strcmp(opt1, "minStrain") == 0) {
                critTag = 1;
            } else if (strcmp(opt1, "maxStrain") == 0) {
                critTag = 2;
            } else if (strcmp(opt1, "axialDI") == 0) {
                critTag = 3;

            } else if (strcmp(opt1, "flexureDI") == 0) {
                critTag = 4;

            } else if (strcmp(opt1, "axialLS") == 0) {
                critTag = 5;
            } else if (strcmp(opt1, "shearLS") == 0) {
                critTag = 6;

            } else if (strcmp(opt1, "INFILLWALL") == 0) {
                critTag = 7;

            } else {
                opserr << "Error: RemoveRecorder - Removal Criteria "
                       << opt1 << " not recognized\n";
                return 0;
            }

            if (critTag != 7) {
                if (OPS_GetNumRemainingInputArgs() < 1) {
                    opserr << "WARNING: need criteria value\n";
                    return 0;
                }
                int num = 1;
                if (OPS_GetDoubleInput(&num, &critValue) < 0) {
                    opserr << "WARNING recorder Remove -crit critTag? "
                              "critValue?... invalid critValue \n";
                    return 0;
                }
            }

            remCriteria[2 * numCrit] = critTag;
            //      new
            if (critTag != 7) {
                remCriteria[2 * numCrit + 1] = critValue;
            } else {
                remCriteria[2 * numCrit + 1] = 100.0;
            }
            numCrit++;
            if (critTag != 7) {
                loc += 3;
            } else {
                loc += 2;
                secIDs = ID(1);
                secIDs[1] = 1;  // this is not used directly, gets rid of
                                // the "-sec" for infillwall
            }
        } else {
            flags = 1;
        }

        // if user has specified no element tags lets assume he wants them
        // all
        if (numEle == 0) {
            ElementIter &theEleIter = domain->getElements();
            Element *theEle;
            while ((theEle = theEleIter()) != 0) {
                eleIDs[numEle++] = theEle->getTag();
            }
        }
    }
    theOutputStream = new DummyStream();

    return new RemoveRecorder(
        nodeTag, eleIDs, secIDs, secondaryEleIDs, remCriteria, *domain,
        *theOutputStream, echoTime, dT, fileName, eleMass, gAcc, gDir,
        gPat, nTagbotn, nTagmidn, nTagtopn, globgrav, fileNameinf);
}

//#define MMTDEBUG
//#define MMTDEBUGIO

// initiatie class-wide static members to keep track of removed components
int RemoveRecorder::numRecs = 0;
ID RemoveRecorder::remEleList(0);
ID RemoveRecorder::remNodeList(0);
int RemoveRecorder::numRemEles = 0;
int RemoveRecorder::numRemNodes = 0;

Element** RemoveRecorder::remEles = 0;
Node** RemoveRecorder::remNodes = 0;


char* RemoveRecorder::fileName = 0;
// new
char* RemoveRecorder::fileNameinf = 0;
//
//FileStream* theStream = new FileStream();
//FileStream RemoveRecorder::theFile = *theStream;
ofstream RemoveRecorder::theFile;
// changed
RemoveRecorder::RemoveRecorder(int nodeID, 
			       ID &eleIDs, 
			       ID &secIDs, 
			       ID &secondaryTags, 
			       Vector remCriteria, 
			       Domain &theDomainPtr, 
			       OPS_Stream &s,
			       bool echotimeflag, 
			       double deltat, 
			       const char *theFileName ,
			       Vector eleMass, 
			       double gAcc, 
			       int gDir, 
			       int gPat, 
			       int nTagbotn, 
			       int nTagmidn, 
			       int nTagtopn, 
			       int globgrav, 
			       const char *thefileNameinf)
  :Recorder(RECORDER_TAGS_RemoveRecorder),
   nodeTag(nodeID), 
   numEles(eleIDs.Size()), 
   eleTags(eleIDs.Size()), 
   secTags(secIDs.Size()), 
   numSecs(secIDs.Size()), 
   criteria(remCriteria),
   theDomain(&theDomainPtr), 
   secondaryEleTags(secondaryTags.Size()), 
   secondaryFlag(false),
   echoTimeFlag(echotimeflag), 
   deltaT(deltat), 
   nextTimeStampToRecord(0.0), 
   gAcc(gAcc), 
   gDir(gDir), 
   gPat(gPat),
   nTagbotn(nTagbotn), 
   nTagmidn(nTagmidn), 
   nTagtopn(nTagtopn), 
   globgrav(globgrav),
   eleResponses(0)
{
  numRecs++;
#ifdef MMTDEBUGIO
  opserr<<"RemoveRecorder, constructor called"<<endln;
#endif
  
  numRules = criteria.Size()/2;

#ifdef MMTDEBUG
  for (int h=0 ; h<2 ; h++) {
    opserr<<"remCriteria["<<h<<"] = "<<criteria[h]<<endln;
  }
#endif

  eleResponses = new Response *[numEles];  
  for (int l=0 ; l<numEles ; l++) {
    eleTags(l) = eleIDs(l);
    eleResponses[l] = 0;
  }

  if (secIDs[0] != 0 || secIDs.Size() != 1) {
    for (int l=0 ; l<numSecs ; l++) {
      secTags(l) = secIDs(l);
#ifdef MMTDEBUG
      opserr<<"storing secID = "<<secIDs[l]<<endln;
#endif
    }
  } else
    secTags[0] = 0;
  
  if (secondaryTags[0] != 0 || secondaryTags.Size() != 1) {
    secondaryFlag = true;
    for (int l=0 ; l<secondaryTags.Size() ; l++) {
      secondaryEleTags(l) = secondaryTags(l);
#ifdef MMTDEBUG
      opserr<<"storing secondaryEleID = "<<secondaryTags[l]<<endln;
#endif
    }
  } else
    secondaryEleTags[0] = 0;


  if (thefileNameinf != 0)  { 
    int fileNameLength2 = strlen(thefileNameinf) + 1;
    //   fileName = new char[fileNameLength + 8];
    fileNameinf = new char[fileNameLength2];
    strcpy(fileNameinf, thefileNameinf);
  }
  
  Element *theEle = 0;
  const char **argv = new const char *[1];
  if (fileNameinf == 0)
    argv[0] = "getRemCriteria1";
  else
    argv[0] = "getRemCriteria2";


  // Get the element	
  for (int j= 0; j<numEles; j++) {
    Element *theEle = theDomainPtr.getElement(eleTags[j]);
    if ( theEle == NULL ) {
      opserr << "WARNING RemoveRecorder::RemoveRecorder() - no element with tag: "
	     << eleTags[j] << " exists in Domain\n";
      eleResponses[j] = 0;
    } else {

      // set up the element responses
      eleResponses[j] = theEle->setResponse(argv, 1, s);
      if (eleResponses[j] == 0) {
	opserr << "WARNING :: getRemCriteria - not a response quantity of element\n";
      } else {
	if (fileNameinf != 0) {
	  Information &eleInfo = eleResponses[j]->getInformation();
	  eleInfo.setString(fileNameinf);
	}
      }  
    }
  }

  delete [] argv;


  if (secondaryEleTags[0] != 0) {
    for (int k= 0; k<secondaryTags.Size(); k++) {
      Element *theEle = theDomainPtr.getElement(secondaryTags[k]);
      if ( theEle == NULL ) {
	opserr << "WARNING RemoveRecorder::RemoveRecorder() - no element with tag: "
	       << secondaryTags[k] << " exists in Domain\n";
	exit(-1);
      }
    }
  }
  
  // check if the recorder is for node removal
  if (nodeTag != 0) {
    Node *theNode = theDomainPtr.getNode(nodeTag);
    if ( theNode == NULL ) {
      opserr << "WARNING RemoveRecorder::RemoveRecorder() - no node with tag: "
	     << nodeTag << " exists in Domain\n";
      exit(-1);
    }
  }
  
  
  // extract masses
  eleMasses = eleMass;
  
  // if file is specified, copy name and open the file if it hasn't been opened
  if (theFileName != 0 && fileName == 0) {
    
    // create char array to store file name
    int fileNameLength = strlen(theFileName) + 1;
    //   fileName = new char[fileNameLength + 8];
    fileName = new char[fileNameLength];
    
    if (fileName == 0) {
      opserr << "RemoveRecorder::RemoveRecorder - out of memory creating string " <<
	fileNameLength << " long\n";
      exit(-1);
    }
    
    // compose and copy file name string
    strcpy(fileName, theFileName);
    //	char* dumName = new char(fileNameLength) + 1;
    //	strcpy(dumName, theFileName); 
    
    //	char *temp = "Collapse";
    //	fileName = strcat(temp, dumName);	
#ifdef MMTDEBUGIO
    //	opserr<<"theFileName "<<theFileName<<"len "<<strlen(theFileName)<<" fileNameLength "<<fileNameLength<<endln;
    //	opserr<<"filename "<<fileName<<"len "<<strlen(fileName)<<" temp "<<temp<<" len "<<strlen(temp)<<endln;
    
    //	opserr<<"filename "<<fileName*<<"len "<<strlen(fileName)<<" temp "<<temp*<<" len "<<strlen(temp)<<endln;
#endif
    
#ifdef MMTDEBUGIO
    opserr<<"filename "<<fileName<<" len "<<strlen(fileName)<<endln;
#endif
    //	delete [] temp;
    
    // open the file
    theFile.open(fileName, ios::out);
    if ((&theFile) == NULL) {
      opserr << "WARNING - RemoveRecorder::RemoveRecorder()";
      opserr << " - could not open file " << fileName << endln;
    }    
  }
  
#ifdef MMTDEBUG
  opserr<<"RemoveRecorder, constructor finished"<<endln;
#endif
}

RemoveRecorder::~RemoveRecorder()
{
  numRecs--;

  // if last recorder to be deleted, check if any elements or nodes had been "removed" and destroy them
  if (numRecs == 0) {
   
    for (int i=0; i<numRemEles; i++)
		if (remEles[i] != 0)
          delete remEles[i];
 
    for (int i=0; i<numRemNodes; i++)
		if (remNodes[i] != 0)
           delete remNodes[i];
 
    if (remEles != 0)
      delete [] remEles;
 
    if (remNodes != 0)
      delete [] remNodes;
 
    numRemEles = 0;
    numRemNodes = 0;
    remEles = 0;
    remNodes = 0;

    if (fileName != 0)
      delete [] fileName;
    fileName = 0;

    // close the file
    if (!theFile.bad()) 
      theFile.close();
  }

#ifdef MMTDEBUG	
  opserr<<"RemoveRecorder, destructor finished"<<endln;
#endif
}


int
// changed
RemoveRecorder::record(int commitTag, double timeStamp)
{
  //	currentTime = timeStamp;	// store so that elimElem() and elimNode() can write to the 
#ifdef MMTDEBUG	
  opserr<<"entering record()"<<endln;
#endif
  int result = 0;
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {
    
    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;
    
    if (int(nodeTag) != 0) {
      
      // check if node had already been removed
      int remFlag = 0;
      for (int m =0; m<RemoveRecorder::numRemNodes; m++) {
	if (nodeTag == RemoveRecorder::remNodeList[m]) {
	  remFlag = 1;
#ifdef MMTDEBUG
	  opserr<<" node "<<nodeTag<<" already removed "<<endln;
#endif
	}
      }
      
      if (remFlag == 0) {
	
	// go over connected elements to check connectivity
#ifdef MMTDEBUG
	opserr<<" checking node "<<nodeTag<<endln;
#endif
	int numConEles = numEles;
	for (int j=0; j<numEles; j++) {
	  for (int m =0; m<RemoveRecorder::numRemEles; m++) {
	    if (eleTags[j] == RemoveRecorder::remEleList[m]) {
	      numConEles --;
#ifdef MMTDEBUG
	      opserr<<" element "<<eleTags[j]<<" already removed from node "<<nodeTag<<" "<<numConEles<<" elements left"<<endln;
#endif
	    }
	  }
	}
	
	
	// if all surrounding elements have been removed, eliminate node
	if (numConEles == 0) {
	  Node* theNode = theDomain->getNode(nodeTag);
	  if (theNode != 0)
	    elimNode(nodeTag, timeStamp);
	  //					remNodeList[numRemNodes] = nodeTag;
	  //					remNodes[numRemNodes] = theNode;
	  //					numRemNodes ++;
#ifdef MMTDEBUG
	  opserr<<"Node removal # "<<numRemNodes<<" of node # "<<remNodeList[numRemNodes-1]<<endln;
	  
#endif      
	}			}
      
    } else {
#ifdef MMTDEBUG	
      opserr<<"checking elements"<<endln;
#endif	
      int eleCount = 0; // a counter to see if all elements of this specific recorder were removed 
      for (int j=0; j<numEles; j++) {
	// check if the element has been already removed
	int remFlag = 0;
	for (int m =0; m<RemoveRecorder::numRemEles; m++) {
	  if (eleTags[j] == RemoveRecorder::remEleList[m]) {
	    remFlag = 1;
#ifdef MMTDEBUG			
	    opserr<<" element "<<eleTags[j]<<" already removed "<<endln;
#endif
	  }
	}
	
	if (remFlag == 0) {
	  
	  Element* theEle = theDomain->getElement(eleTags[j]);
	  if (theEle != 0) {
	    int eleFlag = 0;
	    
	    for (int i=0; i< numSecs; i++) {
#ifdef MMTDEBUG			
	      opserr<<" section to be checked "<<secTags[i]<<endln;
#endif
	      // changed
	      int secFlag = this->checkEleRemoval(theEle, eleResponses[j], secTags[i], criteria);
	      if (secFlag == -1)
		return -1;
	      eleFlag += secFlag;
	    }
#ifdef MMTDEBUG	
	    opserr<<" finished checking element, flag = "<<eleFlag<<endln;
#endif
	    if (eleFlag >0) {
#ifdef MMTDEBUG
	      opserr<<" about to eliminate element "<<eleTags[j]<<endln;
#endif
	      this->updateNodalMasses(eleTags[j], eleMasses[j]);
	      this->elimElem(eleTags[j], timeStamp);
#ifdef MMTDEBUG
	      opserr<<" eliminated element "<<eleTags[j]<<endln;
#endif		
	      //						RemoveRecorder::remEleList[RemoveRecorder::numRemEles] = eleTags[j];
	      //						RemoveRecorder::remEles[RemoveRecorder::numRemEles] = theEle;
	      //						RemoveRecorder::numRemEles ++;
#ifdef MMTDEBUG
	      opserr<<"Element removal # "<<numRemEles<<" of element # "<<remEleList[numRemEles-1]<<endln;
#endif
	      remFlag = 1;
	    }
	  }
	}
	
	eleCount += remFlag;
      }
      
      // now check if the secondary elements need to be removed
#ifdef MMTDEBUG
      opserr<<"eleCount = "<<eleCount<<" numEles = "<<numEles<<" secondaryFlag = "<<int(secondaryFlag)<<endln;
#endif
      if (eleCount == numEles && secondaryFlag == true) {
	if(this->elimSecondaries(timeStamp) != 0) {
	  opserr<<"Error: Collapse Recorder - failed to remove secondary components to element "<<eleTags[0]<<endln;
	  return -1;
	}
	secondaryFlag = false;
      }
    }
    
    
    
    
    if (fileName != 0) {
      //		theFile << " \n";
      theFile.flush();
    }
  }
  
  // succesfull completion - return 0
  return result;
}
int 
RemoveRecorder::playback(int commitTag)
{
  return 0;
}

int
RemoveRecorder::restart(void)
{
  theFile.close();
  theFile.open(fileName, ios::out);
  if (!(&theFile)) {
    opserr << "WARNING - RemoveRecorder::restart() - could not open file ";
    opserr << fileName << endln;
  }    
  return 0;
}

int
// changed
RemoveRecorder::checkEleRemoval(Element* theEle, Response *eleResponse, int &theComponent,const Vector &criteria)
{
  if (eleResponse == 0)
    return 0;

  int secFlag = 0;

  Information &theInfo = eleResponse->getInformation();
  
  for (int k =0; k<numRules; k++) {
    double currentValue = 0;
    double checkwallvalue1=0;
    // changed
    //		if (int(criteria[2*k] != 0 && criteria[2*k] != 7)) {
    //			if (theEle->getRemCriteria(theComponent, int(criteria[2*k]), currentValue) != 0) {
    //				opserr<<" RemoveRecorder::checkRemoval() - failed to retrieve information from element "<<theEle->getTag()<<endln;
    //				return -1;
    //			}
    //		}
    // new

    if (int(criteria[2*k] == 7)) {

      static ID idData(6);
      idData(0) = theComponent;
      idData(1) = int(criteria[2*k]);
      idData(2) = nTagbotn;
      idData(3) = nTagmidn;
      idData(4) = nTagtopn;
      idData(5) = globgrav;;
      theInfo.setID(idData);

      eleResponse->getResponse();

      Vector *result = theInfo.theVector;
      currentValue = (*result)(0);
      checkwallvalue1 = (*result)(1);
    } 
    //
    
    switch (int(criteria[2*k])) {
      
    case 0:
      // no removal,for comparison;
      break;
      
    case 1:
      // min strain (maximum compressive), notice the inequality direction
      if (currentValue <= criteria [2*k+1]) 
	secFlag ++;
      break;
      
    case 2:
      // max strain
      if (currentValue >= criteria [2*k+1]) 
	secFlag ++;
      break;
      
    case 3:
      // axial failure damage index
      if (currentValue >= criteria [2*k+1]) 
	secFlag ++;
      break;
      
    case 4:
      // flexural failure damage index
      if (currentValue >= criteria [2*k+1]) 
	secFlag ++;
      break;
      
    case 5:
      // axial limit state
      if (currentValue >= criteria [2*k+1]) 
	secFlag ++;
      break;
      
    case 6:
      // shear limit state
      if (currentValue >= criteria [2*k+1]) 
	secFlag ++;
      break;
      //          new					
    case 7:
      // infill wall element
      if (currentValue >= checkwallvalue1) 
	secFlag ++;
      break;
      //				
    default:
      // unknown removal criteria switch
      opserr<<"RemoveRecorder::checkRemoval() - Failed to identify removal criterion "<<criteria [2*k]<<endln;
      return -1;
      break;
      
    }
  }
#ifdef MMTDEBUG
  opserr<<"secFlag = "<<secFlag<<endln;
#endif
  return secFlag;
}


int
RemoveRecorder::elimElem(int theEleTag, double timeStamp)
{

#ifdef MMTDEBUG
  opserr << "RemoveRecorder::elimElem() remving ele: " << theEleTag << " at timeStamp: " << timeStamp << endln;;
#endif

  Element *theEle = theDomain->removeElement(theEleTag);
  if (theEle != 0) {
    // we also have to remove any elemental loads from the domain
    LoadPatternIter &theLoadPatterns = theDomain->getLoadPatterns();
    LoadPattern *thePattern;

    // go through all load patterns
    while ((thePattern = theLoadPatterns()) != 0) {

      ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
      ElementalLoad *theLoad;
      
      // go through all elemental loads in the pattern
      while ((theLoad = theEleLoads()) != 0) {
	
	// remove & destroy elemental from elemental load if there
	// Note - if last element in load, remove the load and delete it

	int loadEleTag = theLoad->getElementTag();
	if (loadEleTag == theEleTag) {
	  opserr << "RemoveRecorder::elimElem() -3 removing  eleLoad\n";

	  ElementalLoad *theElementalLoad = thePattern->removeElementalLoad(theLoad->getTag());
	  if (theElementalLoad != 0) {
	    delete theElementalLoad;
	  }
	}
      }
    }

    // finally invoke the destructor on the element 
    //	delete theEle;
    /////////////// M.Talaat : Avoid recorder trouble at element removal and just set it to zero after removing it from the domain!
    theEle->revertToStart();

    RemoveRecorder::remEleList[RemoveRecorder::numRemEles] = theEle->getTag();

    Element **newRemEles = new Element *[numRemEles+1];
    for (int ii=0; ii<numRemEles; ii++)
      newRemEles[ii] = remEles[ii];

    newRemEles[numRemEles] = theEle;
    if (remEles != 0)
      delete [] remEles;

    remEles = newRemEles;

    numRemEles ++;

    // now give us some notice of what happened
    if (fileName != 0)
      theFile<<timeStamp<<" Elem "<<theEle->getTag()<<"\n";
    if (echoTimeFlag == true)
#ifdef MMTDEBUG
    opserr<< " element "<<theEle->getTag()<<" removed automatically"<<endln; 
#endif
    ;
  }
  return 0;
}


int
RemoveRecorder::elimNode(int theNodeTag, double timeStamp)
{  
  // remove from domain but do not delete yet!
  Node *theNode = theDomain->removeNode(theNodeTag);
  
  // go through all load patterns and remove associated loads and constraints
  LoadPatternIter &theLoadPatterns = theDomain->getLoadPatterns();
  LoadPattern *thePattern;
  
  while ((thePattern = theLoadPatterns()) != 0) {
    
    // start with nodal laods
    NodalLoadIter theLoads = thePattern->getNodalLoads();
    NodalLoad *theLoad;
    
    //		ID theLoadTags(0,12); 
    //		int cnt=0;
    while ((theLoad = theLoads()) != 0) {
      
      int theNode = theLoad->getNodeTag();
      if (theNode == theNodeTag) {
	//				theLoadTags[cnt] = theLoad->getTag();
	//				cnt++;
#ifdef MMTDEBUG
	opserr<<"identified load pattern "<<theLoad->getTag()<<" acting on node "<<theNode<<endln;
#endif
	NodalLoad *theNodalLoad = thePattern->removeNodalLoad(theLoad->getTag());
	if (theNodalLoad != 0) {
#ifdef MMTDEBUG
	  opserr<<"deleting nodal load pattern "<<theLoad->getTag()<<endln;
#endif
	  delete theNodalLoad;
	}	
      }
    }
    
    //		for (int i=0; i<cnt; i++) {
    //			NodalLoad *theNodalLoad = thePattern->removeNodalLoad(theLoadTags[i]);
    //			if (theNodalLoad != 0) {
    //				delete theNodalLoad;
    //																			#ifdef MMTDEBUG
    //				opserr<<"deleting nodal load pattern "<<theLoadTags(i)<<endln;
    //																			#endif
    //			}	
    //		}
    
    // follow with sp constraints
    SP_ConstraintIter &theSPs = thePattern->getSPs();
    SP_Constraint *theSP;
    
    while ((theSP = theSPs()) != 0) {
      
      int spNode = theSP->getNodeTag();
      if (spNode == theNodeTag) {
//				theSPTags[cnt] = theSP->getTag();
//				cnt++;
#ifdef MMTDEBUG
	opserr<<"identified SP_Constraint "<<theSP->getTag()<<" acting on node "<<spNode<<endln;
#endif
	SP_Constraint *theSPConstraint = thePattern->removeSP_Constraint(theSP->getTag());
	if (theSPConstraint != 0) {
#ifdef MMTDEBUG
	  opserr<<"deleting SP_Constraint "<<theSP->getTag()<<endln;
#endif
	  delete theSPConstraint;
	}	
      }
    }
  }
  
  
  // we also have to remove any sp constraints from the domain that do not belong to load patterns (support fixity)
  SP_ConstraintIter &theSPs = theDomain->getSPs();
  SP_Constraint *theSP;
  
  //	  ID theSPTags(0,12); 
  //	  int cnt=0;
  while ((theSP = theSPs()) != 0) {
    
    int spNode = theSP->getNodeTag();
    if (spNode == theNodeTag) {
      //				theSPTags[cnt] = theSP->getTag();
      //				cnt++;
#ifdef MMTDEBUG
      opserr<<"identified SP_Constraint "<<theSP->getTag()<<" acting on node "<<spNode<<endln;
#endif
      SP_Constraint *theSPConstraint = theDomain->removeSP_Constraint(theSP->getTag());
      if (theSPConstraint != 0) {
#ifdef MMTDEBUG
	opserr<<"deleting SP_Constraint "<<theSP->getTag()<<endln;
#endif
	delete theSPConstraint;
      }	
    }
  }
  
  //		for (int i=0; i<cnt; i++) {
  //		  SP_Constraint *theSPconstraint = theDomain->removeSP_Constraint(theSPTags[i]);
  //		  if (theSPconstraint != 0) {
  //		    delete theSPconstraint;
  //		  }	
  //		}
  
  if (theNode != 0) {
    // delete theNode;
    /////////////////// M.Talaat : Again, avoid recorder trouble
    theNode->revertToStart();
  }
  
  RemoveRecorder::remNodeList[numRemNodes] = theNode->getTag();
  //  RemoveRecorder::remNodes[numRemNodes] = theNode;
  //  RemoveRecorder::numRemNodes ++;

  Node **newRemNodes = new Node *[numRemNodes+1];
  for (int ii=0; ii<numRemNodes; ii++)
    newRemNodes[ii] = remNodes[ii];
  newRemNodes[numRemNodes] = theNode;
  if (remNodes != 0)
    delete [] remNodes;
  remNodes = newRemNodes;
  
  numRemNodes++;

	       
  
  // now give us some notice of what happened
  if (fileName != 0)
    theFile<<timeStamp<<" Node "<<theNode->getTag()<<"\n";
  if (echoTimeFlag == true)
    opserr<<"Node "<<theNode->getTag()<<" removed, Time/Load Factor = " <<timeStamp<<endln;
  
  return 0;
}



int 
RemoveRecorder::elimSecondaries(double timeStamp)
{
																			#ifdef MMTDEBUG	
	opserr<<"entering elimSecondaries()"<<endln;
																			#endif
	int result = 0;


	
	// remove secondary elements and nodes
	for (int i =0; i<secondaryEleTags.Size(); i++) {

		int remFlag = 0;
		for (int m =0; m<RemoveRecorder::numRemEles; m++) {
			if (secondaryEleTags[i] == RemoveRecorder::remEleList[m]) {
				remFlag = 1;
																			#ifdef MMTDEBUG			
			opserr<<" element "<<eleTags[i]<<" already removed "<<endln;
																			#endif
			}
		}

		if (remFlag == 0) {
			Element *theEle = theDomain->getElement(secondaryEleTags[i]);
			ID secondaryNodes = theEle->getExternalNodes();
			
			for (int k = 0; k<theEle->getNumExternalNodes(); k++) {
				
				int nodeRemFlag = 0;
				for (int m =0; m<RemoveRecorder::numRemNodes; m++) {
						
					if (secondaryNodes[k] == RemoveRecorder::remNodeList[m]) {
						nodeRemFlag = 1;
																			#ifdef MMTDEBUG			
					opserr<<" node "<<secondaryNodes[k]<<" already removed "<<endln;
																			#endif
					}
				}

				if (nodeRemFlag == 0) {
				this->elimNode(secondaryNodes[k], timeStamp);
																			#ifdef MMTDEBUG
				opserr<<" eliminated node "<<secondaryNodes[k]<<endln;
				opserr<<"Node removal # "<<numRemNodes<<" of node # "<<remNodeList[numRemNodes-1]<<endln;
																			#endif
				}
			}
				
			this->elimElem(secondaryEleTags[i], timeStamp);
																			#ifdef MMTDEBUG
			opserr<<" eliminated element "<<secondaryEleTags[i]<<endln;
			opserr<<"Element removal # "<<numRemEles<<" of element # "<<remEleList[numRemEles-1]<<endln;
																			#endif
		}

    }

  return result;
}


int 
RemoveRecorder::updateNodalMasses(int theEleTag, double theEleMass)
{
 // update the nodal masses and loads if requested
	  if (theEleMass != 0) {

		  Element* theEle = theDomain->getElement(theEleTag);
		  ID endNodes = theEle->getExternalNodes();
	  	  for (int k = 0; k<theEle->getNumExternalNodes(); k++) {

			  Node* theNode=theDomain->getNode(endNodes(k));
			  Matrix theNodalMass = theNode->getMass();
			  for (int l = 0; l<theNodalMass.noRows(); l++) {
				  if (theNodalMass(l,l) != 0)
				  theNodalMass(l,l) -= 0.5*theEleMass;
			  }
			  if (theDomain->setMass(theNodalMass, theNode->getTag()) != 0) {
				  opserr << "Remove Recorder::WARNING failed to set mass at node " << theNode->getTag() << endln;
			  }


			  if (this->gAcc != 0) {
			  	  double theEleWeight = theEleMass * gAcc;
			  	  Vector theNodalForces(theNode->getNumberDOF());
				  theNodalForces.Zero();
				  theNodalForces(gDir-1) = 0.5 * theEleWeight;
				  bool isLoadConst = true;
				  opserr<<"eleWeight "<<theEleWeight<<" NodalForces "<<theNodalForces(0)<<" "<<theNodalForces.Size()<<endln;
				  int nodeLoadTag = 987654 + theEleTag * 30 + k; // big number + many element nodes + current node
				  NodalLoad* theNodalLoad = new NodalLoad(nodeLoadTag, theNode->getTag(), theNodalForces, isLoadConst);
				 if (theNodalLoad == 0) {
					opserr << "RemoveRecorder::WARNING ran out of memory while updating loads on node  " << theNode->getTag();
    			 }
				 if (theDomain->addNodalLoad(theNodalLoad, gPat) == false) {
					 opserr << "RemoveRecorder::WARNING could not add updated nodal load to domain\n";
					 delete theNodalLoad;
				 }
			  }	
		  }
	  }
	  return 0;
}

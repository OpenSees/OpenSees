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
                                                                        
// $Revision: 1.34 $
// $Date: 2005-03-30 20:12:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/TclRecorderCommands.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 04/98
// Revision: AA
//
// Description: This file contains the function that is invoked
// by the interpreter when the comand 'record' is invoked by the 
// user.
//
// What: "@(#) commands.C, revA"


#include <tcl.h>
#include <tk.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <EquiSolnAlgo.h>

// recorders
#include <NodeRecorder.h>
#include <EnvelopeNodeRecorder.h>
#include <EnvelopeElementRecorder.h>
#include <PatternRecorder.h>
#include <DriftRecorder.h>
#include <ElementRecorder.h>


#include <NodeIter.h>
#include <ElementIter.h>
#include <Node.h>
#include <Element.h>
#include <DamageModel.h>
#include <DamageRecorder.h>
#include <MeshRegion.h>
#include <GSA_Recorder.h>
#include <TclModelBuilder.h>

#include <NEESData.h>

#include <DataOutputFileHandler.h>
#include <DataOutputDatabaseHandler.h>
#include <DataOutputStreamHandler.h>

extern TclModelBuilder *theDamageTclModelBuilder;

#include <EquiSolnAlgo.h>

#ifdef _NOGRAPHICS

#else
#include <TclFeViewer.h>
#include <FilePlotter.h>
#include <YsVisual.h> //!!
#include <AlgorithmIncrements.h>
#endif



static EquiSolnAlgo *theAlgorithm =0;
extern FE_Datastore *theDatabase;
extern FEM_ObjectBroker theBroker;

int
TclCreateRecorder(ClientData clientData, Tcl_Interp *interp, int argc,
		  TCL_Char **argv, Domain &theDomain, Recorder **theRecorder)
{
    // make sure at least one other argument to contain integrator
    if (argc < 2) {
	opserr << "WARNING need to specify a Recorder type\n";
	return TCL_ERROR;
    }

    //
    // check argv[1] for type of Recorder, parse in rest of arguments
    // needed for the type of Recorder, create the object and add to Domain
    //
    (*theRecorder) = 0;
    (*theRecorder) = 0;
    FE_Datastore *theRecorderDatabase = 0;
    DataOutputHandler *theDataOutputHandler = 0;
    TCL_Char *fileName = 0;
    TCL_Char *tableName = 0;

    // an Element Recorder or ElementEnvelope Recorder
    if ((strcmp(argv[1],"Element") == 0) || (strcmp(argv[1],"EnvelopeElement") == 0)
	|| (strcmp(argv[1],"ElementEnvelope") == 0)) {

        /* KEEP - FOR LEGACY REASONS NEED TO KEEP THE FOLLOWING UGLY STUFF */
        int eleID;
        if (argc < 4) {
	    opserr << "WARNING recorder Element eleID1? eleID2? ...  <-time> "
		<< "<-file fileName?> parameters";
	    return TCL_ERROR;
	}    

	int endEleIDs = 2;
	int allFlag = 0;
	while (Tcl_GetInt(interp, argv[endEleIDs], &eleID) == TCL_OK) {
	  endEleIDs++;
	} 
	Tcl_ResetResult(interp);

	// determine the number of elements
	int numEle = 0;
	if (strcmp(argv[endEleIDs],"all") == 0) {
	  endEleIDs += 1;
	  allFlag = 1;
	  numEle = theDomain.getNumElements();
	} else if (Tcl_GetInt(interp, argv[endEleIDs], &eleID) != TCL_OK) 
	  ;
	else
	  endEleIDs++;
	  
	numEle = endEleIDs-2;
	
	// create an ID to hold ele tags
        ID eleIDs(numEle, numEle+1); 

	// read in the ele tags to the ID
	if (allFlag == 1) {
	  int loc = 0;
	  ElementIter &theEleIter = theDomain.getElements();
	  Element *theEle;
	  while ((theEle = theEleIter()) != 0)
	    eleIDs[loc++] = theEle->getTag();
	} else {
	  for (int i=2; i<endEleIDs; i++) {
	    if (Tcl_GetInt(interp, argv[i], &eleID) != TCL_OK)	
	      return TCL_ERROR;	
	    eleIDs[i-2] = eleID;	  
	  }
	}
        /* ********************* END OF KEEP ****************************** */

	double dT = 0.0;
	bool echoTime = false;
	echoMode eMode = NONE;       // enum found in DataOutputFileHandler.h
	int loc = endEleIDs;
	int flags = 0;
	int eleData = 0;

	while (flags == 0 && loc < argc) {

	  if ((strcmp(argv[loc],"-ele") == 0) ||
	      (strcmp(argv[loc],"-eles") == 0) ||
	      (strcmp(argv[loc],"-element") == 0)) {
      
	    // ensure no segmentation fault if user messes up
	    if (argc < loc+2) {
	      opserr << "WARNING recorder Element .. -ele tag1? .. - no ele tags specified\n";
	      return TCL_ERROR;
	    }
	    
	    //
	    // read in a list of ele until end of command or other flag
	    //
	    loc++;
	    int eleTag;
	    while (loc < argc && Tcl_GetInt(interp, argv[loc], &eleTag) == TCL_OK) {
	      eleIDs[numEle++] = eleTag;
	      loc++;
	    }
	    Tcl_ResetResult(interp);
	    
	    if (loc == argc) {
	      opserr << "ERROR: No response type specified for element recorder. " << endln;
	      return TCL_ERROR;
	    }

	    if (strcmp(argv[loc],"all") == 0) {
	      ElementIter &theEleIter = theDomain.getElements();
	      Element *theEle;
	      while ((theEle = theEleIter()) != 0)
		eleIDs[numEle++] = theEle->getTag();
	      loc++;
	    }

	  } else if (strcmp(argv[loc],"-eleRange") == 0) {
	    
	    // ensure no segmentation fault if user messes up
	    if (argc < loc+3) {
	      opserr << "WARNING recorder Element .. -eleRange start? end?  .. - no ele tags specified\n";
	      return TCL_ERROR;
	    }
	    
	    //
	    // read in start and end tags of two elements & add set [start,end]
	    //
	    
	    int start, end;
	    if (Tcl_GetInt(interp, argv[loc+1], &start) != TCL_OK) {
	      opserr << "WARNING recorder Element -eleRange start? end? - invalid start " << argv[loc+1] << endln;
	      return TCL_ERROR;
	    }      
	    if (Tcl_GetInt(interp, argv[loc+2], &end) != TCL_OK) {
	      opserr << "WARNING recorder Element -eleRange start? end? - invalid end " << argv[loc+2] << endln;
	      return TCL_ERROR;
	    }      
	    if (start > end) {
	      int swap = end;
	      end = start;
	      start = swap;
	    }

	    for (int i=start; i<=end; i++)
	      eleIDs[numEle++] = i;	    

	    loc += 3;
	  } 

	  else if (strcmp(argv[loc],"-region") == 0) {
	    // allow user to specif elements via a region

	    if (argc < loc+2) {
	      opserr << "WARNING recorder Element .. -region tag?  .. - no region specified\n";
	      return TCL_ERROR;
	    }
	    int tag;
	    if (Tcl_GetInt(interp, argv[loc+1], &tag) != TCL_OK) {
	      opserr << "WARNING recorder Element -region tag? - invalid tag " << argv[loc+1] << endln;
	      return TCL_ERROR;
	    }      
	    MeshRegion *theRegion = theDomain.getRegion(tag);
	    if (theRegion == 0) {
	      opserr << "WARNING recorder Element -region " << tag << " - region does not exist" << endln;
	      return TCL_OK;
	    }      
	    const ID &eleRegion = theRegion->getElements();
	    for (int i=0; i<eleRegion.Size(); i++)
	      eleIDs[numEle++] = eleRegion(i);

	    loc += 2;
	  } 

	  else if ((strcmp(argv[loc],"-time") == 0) || (strcmp(argv[loc],"-load") == 0)) { 
	    // allow user to specify const load
	    echoTime = true;
	    loc++;
	  } 

	  else if (strcmp(argv[loc],"-dT") == 0) {
	    // allow user to specify time step size for recording
	    loc++;
	    if (Tcl_GetDouble(interp, argv[loc], &dT) != TCL_OK)	
	      return TCL_ERROR;	
	    loc++;
	  } 

	  else if (strcmp(argv[loc],"-file") == 0) {
	    fileName = argv[loc+1];
	    loc += 2;
	    if (strcmp(argv[loc],"-xml") == 0) {
	      eMode = XML_FILE;
	      loc+=1;
	    } else if (strcmp(argv[loc],"-headings") == 0) {
	      eMode = DATA_FILE;
	      loc +=1;
	    }
	  }
	  
	  else if (strcmp(argv[loc],"-database") == 0) {
	    theRecorderDatabase = theDatabase;
	    if (theRecorderDatabase != 0)
	      tableName = argv[loc+1];
	    else {
	      opserr << "WARNING recorder Node .. -database &lt;fileName&gt; - NO CURRENT DATABASE, results to File instead\n";
	      fileName = argv[loc+1];
	    }
	    
	    loc += 2;
	  }

	  else if (strcmp(argv[loc],"-nees") == 0) {
	    // allow user to specify load pattern other than current
	    fileName = argv[loc+1];
	    eMode = XML_FILE;
	    loc += 2;
	  }	    

	  else {
	    // first unknown string then is assumed to start 
	    // element response request starts
	    eleData = loc;
	    flags = 1;
	  }
	}

	// if user has specified no element tags lets assume he wants them all
	if (numEle == 0) {
	  ElementIter &theEleIter = theDomain.getElements();
	  Element *theEle;
	  while ((theEle = theEleIter()) != 0)
	    eleIDs[numEle++] = theEle->getTag();
	}

	if (eleData >= argc) {
	  opserr << "ERROR: No response type specified for element recorder. " << endln;
	  return TCL_ERROR;
	}

	const char **data = new const char *[argc-eleData];

	int i,j;
	for (i=eleData, j=0; i<argc; i++, j++)
	  data[j] = argv[i];


      // construct the DataHandler
      if (fileName != 0) {
	theDataOutputHandler = new DataOutputFileHandler(fileName, eMode);
      } else if (tableName != 0) {
	theDataOutputHandler = new DataOutputDatabaseHandler(theDatabase, tableName);
      } else
	theDataOutputHandler = new DataOutputStreamHandler();


      if (strcmp(argv[1],"Element") == 0) 
	(*theRecorder) = new ElementRecorder(eleIDs, 
					     data, 
					     argc-eleData, 
					     echoTime, 
					     theDomain, 
					     *theDataOutputHandler,
					     dT);
      else
	(*theRecorder) = new EnvelopeElementRecorder(eleIDs, 
						     data, 
						     argc-eleData, 
						     theDomain, 
						     *theDataOutputHandler,
						     dT);
      
      delete [] data;
    }
    

    //////////////////////Begining of ElementDamage recorder//////////////////////
    ///////////////////////////////By Arash Altoontash////////////////////////////

    else if ( (strcmp(argv[1],"Damage") == 0) || (strcmp(argv[1],"ElementDamage") == 0) ||
		(strcmp(argv[1],"damage") == 0) || (strcmp(argv[1],"elementDamage") == 0) ) {
		

      if (argc < 7) {
	opserr << "WARNING recorder ElementDamage eleID? <-time> "
	       << "<-file fileName?> <-section secID1? secID2? ...> <-dof dofID?> <-damage dmgID?>";
	return TCL_ERROR;
      }    
      
      double dT = 0.0;
      bool echoTime = false;
      int loc = 2;
      int eleID;
      
      if ( Tcl_GetInt(interp, argv[loc], &eleID) != TCL_OK) {
	opserr << "WARNING recorder ElementDamage: No element tag specified ";
	return TCL_ERROR;
      }
      loc++;
      
      if ( (strcmp(argv[loc],"-time") == 0) || (strcmp(argv[loc],"-load") == 0) ) { 
	// allow user to specify const load
	echoTime = true;
	loc++;
      }		
      else if ( strcmp(argv[loc],"-dT" ) == 0 ) {
	// allow user to specify time step size for recording
	loc++;
	if (Tcl_GetDouble(interp, argv[loc], &dT) != TCL_OK)	
	  return TCL_ERROR;	
	loc++;
      } 
      
      
      if (strcmp(argv[loc],"-file") == 0) {
	// allow user to specify load pattern other than current
	loc++;
	fileName = argv[loc];
	loc++;
      }
      
      if ( strcmp(argv[loc],"-section")!=0 && strcmp(argv[loc],"section")!=0 ) {
	opserr << "WARNING recorder ElementDamage: Section keyword not specified ";
	return TCL_ERROR;
      }
      loc++;
      
      int secID;		
      int endSecIDs = loc;
      int numSec = 0;
      while ( Tcl_GetInt(interp, argv[endSecIDs], &secID ) == TCL_OK) {
	endSecIDs++;
      }
      
      numSec = endSecIDs - loc ;
      // create an ID to hold section/material tags
      ID secIDs(numSec);
      
      // read in the sec tags to the ID
      for (int i=loc; i<endSecIDs; i++) {
	if (Tcl_GetInt(interp, argv[i], &secID) != TCL_OK)	
	  return TCL_ERROR;	
	secIDs[loc-i] = secID;	  
      }
      
      loc = endSecIDs;
      
      int dofID = 0;
      if ( strcmp(argv[loc],"-dof")==0 || strcmp(argv[loc],"dof")==0 ||
	   strcmp(argv[loc],"-DOF")==0 || strcmp(argv[loc],"DOF")==0 ) {
	loc++;
	if (Tcl_GetInt(interp, argv[loc], &dofID) != TCL_OK) {
	  opserr << "WARNING recorder ElementDamage: No dof tag specified ";
	  return TCL_ERROR;
	}
	loc++;
      }
      
      if ( strcmp(argv[loc],"-damage")!=0 && strcmp(argv[loc],"damage")!=0 ) {
	opserr << "WARNING recorder ElementDamage: No damege tag specified ";
	return TCL_ERROR;
      }
      loc++;
      
      int dmgID;
      if (Tcl_GetInt(interp, argv[loc], &dmgID) != TCL_OK) {
	opserr << "WARNING recorder ElementDamage: No damege tag specified ";
	return TCL_ERROR;
      }
      
      DamageModel *dmgPTR;
      dmgPTR = theDamageTclModelBuilder->getDamageModel(dmgID);
      
      if ( dmgPTR == NULL )
	{
	  opserr << "WARNING recorder ElementDamage: specified damage model not found\n";
	  return TCL_ERROR;
	}
      
      
      //	const char **data = new const char *[argc-eleData];
      
      
      // now construct the recorder
      (*theRecorder) = new DamageRecorder( eleID , secIDs, dofID , dmgPTR , theDomain, echoTime, dT , fileName );
      
      
    }
	//////////////////////End of ElementDamage recorder////////////////////////////
    
    
    // create a recorder to write nodal displacement quantities to a file
    else if ((strcmp(argv[1],"Node") == 0) || (strcmp(argv[1],"EnvelopeNode") == 0) 
	     || (strcmp(argv[1],"NodeEnvelope") == 0)) {	

      if (argc < 7) {
	    opserr << "WARNING recorder Node ";
	    opserr << "-node <list nodes> -dof <doflist> -file <fileName> -dT <dT> reponse";
	    return TCL_ERROR;
	}    

// AddingSensitivity:BEGIN ///////////////////////////////////
      int sensitivity = 0;
// AddingSensitivity:END /////////////////////////////////////

      TCL_Char *responseID = 0;
      echoMode eMode = NONE;       // enum found in DataOutputFileHandler.h

      int pos = 2;

      /* KEEP - FOR LEGACY REASONS NEED TO KEEP THE FOLLOWING UGLY STUFF */
      if ((strcmp(argv[pos],"-time") != 0) && (strcmp(argv[pos],"-load") != 0) &&
	  (strcmp(argv[pos],"-dT") !=  0) && (strcmp(argv[pos],"-node") != 0) &&
	  (strcmp(argv[pos],"-nees") !=  0) && (strcmp(argv[pos],"-database") != 0) &&
	  (strcmp(argv[pos],"-dof") != 0) && (strcmp(argv[pos],"-file") != 0)) {
	pos = 4;
	responseID = argv[3];
	fileName = argv[2];
      } 
      /* ********************** END OF KEEP ***************************  */

      bool echoTimeFlag = false;
      int flags = 0;
      double dT = 0.0;
      int numNodes = 0;

      // create ID's to contain the node tags & the dofs
      ID theNodes(0,16);
      ID theDofs(0, 6);

      while (flags == 0 && pos < argc) {

	if (strcmp(argv[pos],"-time") == 0) {
	  echoTimeFlag = true;
	  pos++;
	}
	
	else if (strcmp(argv[pos],"-load") == 0) {
	  echoTimeFlag = true;      
	  pos++;
	}

	else if (strcmp(argv[pos],"-file") == 0) {
	  fileName = argv[pos+1];
	  pos += 2;
	  if (strcmp(argv[pos],"-xml") == 0) {
	    eMode = XML_FILE;
	    pos+=1;
	  } else if (strcmp(argv[pos],"-headings") == 0) {
	    eMode = DATA_FILE;
	    pos +=1;
	  }
	}

	else if (strcmp(argv[pos],"-database") == 0) {
	  theRecorderDatabase = theDatabase;
	  if (theRecorderDatabase != 0)
	    tableName = argv[pos+1];
	  else {
	    opserr << "WARNING recorder Node .. -database &lt;fileName&gt; - NO CURRENT DATABASE, results to File instead\n";
	    fileName = argv[pos+1];
	  }

	  pos += 2;
	}
	else if (strcmp(argv[pos],"-nees") == 0) {
	  // allow user to specify load pattern other than current
	  fileName = argv[pos+1];
	  eMode = XML_FILE;
	  pos += 2;
	}	    

	else if (strcmp(argv[pos],"-dT") == 0) {
	  pos ++;
	  if (Tcl_GetDouble(interp, argv[pos], &dT) != TCL_OK)	
	    return TCL_ERROR;		  
	  pos++;
	}

	else if ((strcmp(argv[pos],"-node") == 0) || 
		 (strcmp(argv[pos],"-nodes") == 0)) {
	  pos++;
	  
	  // read in the node tags or 'all' can be used
	  if (strcmp(argv[pos],"all") == 0) {
	    numNodes = theDomain.getNumNodes();
	    
	    NodeIter &theNodeIter = theDomain.getNodes();
	    Node *theNode;
	    int loc=0;
	    while ((theNode= theNodeIter()) != 0) {
	      int tag = theNode->getTag();
	      theNodes[loc++] = tag;
	    }
	    pos++;
	  } else {
	    int node;
	    for (int j=pos; j< argc; j++) 
	      if (Tcl_GetInt(interp, argv[pos], &node) != TCL_OK) {
		j = argc;
		Tcl_ResetResult(interp);
	      } else {
		theNodes[numNodes] = node;
		numNodes++;
		pos++;
	      }
	  }
	} 

	else if (strcmp(argv[pos],"-nodeRange") == 0) {
	    
	  // ensure no segmentation fault if user messes up
	  if (argc < pos+3) {
	    opserr << "WARNING recorder Node .. -nodeRange start? end?  .. - no ele tags specified\n";
	    return TCL_ERROR;
	  }
	  
	  //
	  // read in start and end tags of two elements & add set [start,end]
	  //
	    
	  int start, end;
	  if (Tcl_GetInt(interp, argv[pos+1], &start) != TCL_OK) {
	    opserr << "WARNING recorder Node -nodeRange start? end? - invalid start " << argv[pos+1] << endln;
	    return TCL_ERROR;
	  }      
	  if (Tcl_GetInt(interp, argv[pos+2], &end) != TCL_OK) {
	    opserr << "WARNING recorder Node -nodeRange start? end? - invalid end " << argv[pos+2] << endln;
	    return TCL_ERROR;
	  }      
	  if (start > end) {
	    int swap = end;
	    end = start;
	    start = swap;
	  }
	  
	  for (int i=start; i<=end; i++)
	    theNodes[numNodes++] = i;	    
	  pos += 3;
	}

	else if (strcmp(argv[pos],"-region") == 0) {
	  // allow user to specif elements via a region
	  
	  if (argc < pos+2) {
	    opserr << "WARNING recorder Node .. -region tag?  .. - no region specified\n";
	    return TCL_ERROR;
	  }
	  int tag;
	  if (Tcl_GetInt(interp, argv[pos+1], &tag) != TCL_OK) {
	    opserr << "WARNING recorder Node -region tag? - invalid tag " << argv[pos+1] << endln;
	    return TCL_ERROR;
	  }      
	  MeshRegion *theRegion = theDomain.getRegion(tag);
	  if (theRegion == 0) {
	    opserr << "WARNING recorder Node -region " << tag << " - region does not exist" << endln;
	    return TCL_OK;
	  }      
	  const ID &nodeRegion = theRegion->getNodes();
	  for (int i=0; i<nodeRegion.Size(); i++)
	    theNodes[numNodes++] = nodeRegion(i);
	  
	  pos += 2;
	} 
	
	else if (strcmp(argv[pos],"-dof") == 0) {
	  pos++;
	  int numDOF = 0;
	  int dof;
	  for (int j=pos; j< argc; j++) 
	    if (Tcl_GetInt(interp, argv[pos], &dof) != TCL_OK) {
	      j = argc;
	      Tcl_ResetResult(interp);
	    } else {
	      theDofs[numDOF] = dof-1;  // -1 for c indexing of the dof's
	      numDOF++;
	      pos++;
	    }
	}
// AddingSensitivity:BEGIN //////////////////////////////////////
	else if (strcmp(argv[pos],"-sensitivity") == 0) {
		pos++;
		if (Tcl_GetInt(interp, argv[pos], &sensitivity) != TCL_OK) {
			opserr << "ERROR: Invalid gradient number to node recorder." << endln;
			return TCL_ERROR;
		}
		pos++;
	}
// AddingSensitivity:END ////////////////////////////////////////
	else	 
	  flags = 1;
	}
      
      if (pos >= argc) {
	opserr << "WARNING: No response type specified for node recorder, will assume you meant -disp\n" << endln;
      }
      
      if (responseID == 0 && pos < argc) {
	responseID  = argv[pos];
      }

      if (numNodes == 0) {
	NodeIter &theNodeIter = theDomain.getNodes();
	Node *theNode;
	while ((theNode= theNodeIter()) != 0) {
	  int tag = theNode->getTag();
	  theNodes[numNodes++] = tag;
	}
      }

      // construct the DataHandler
      if (fileName != 0) {
	theDataOutputHandler = new DataOutputFileHandler(fileName, eMode);
      } else if (tableName != 0) {
	theDataOutputHandler = new DataOutputDatabaseHandler(theDatabase, tableName);
      } else
	theDataOutputHandler = new DataOutputStreamHandler();


      if (strcmp(argv[1],"Node") == 0) {
	(*theRecorder) = new NodeRecorder(theDofs, 
					  theNodes, 
					  sensitivity,
					  responseID, 
					  theDomain, 
					  *theDataOutputHandler, 
					  dT, 
					  echoTimeFlag);
      } else
	  
	(*theRecorder) = new EnvelopeNodeRecorder(theDofs, 
						  theNodes, 
						  responseID, 
						  theDomain,
						  *theDataOutputHandler,
						  dT);
    }

    else if (strcmp(argv[1],"Pattern") == 0) {
      if (argc < 4) {
	opserr << "WARNING recorder Pattern filename? <startFlag> patternTag?";
	return TCL_ERROR;
      }
      
      int flag = 0;
      if (strcmp(argv[3],"-time") == 0)
	flag = 1;
      if (strcmp(argv[3],"-load") == 0)
	flag = 2;
      
      int pos = 3;
      if (flag != 0) pos = 4;
      
      int patternTag;
      
      if (Tcl_GetInt(interp, argv[pos++], &patternTag) != TCL_OK)
	return TCL_ERROR;
      
      (*theRecorder) = new PatternRecorder(patternTag, theDomain, argv[2], 0.0
					   , flag);
    }
    
    // Create a recorder to write nodal drifts to a file
    else if (strcmp(argv[1],"Drift") == 0) {

      if (argc < 7) {
	opserr << "WARNING recorder Drift filename? <startFlag> ";
	opserr << "node1? node2? dof? perpDirn?\n";
	return TCL_ERROR;
      }    

      if ((strcmp(argv[2],"-file") != 0) && (strcmp(argv[2],"-iNode") != 0) &&
	  (strcmp(argv[2],"-jNode") != 0) && (strcmp(argv[2],"-dirn") != 0) &&
	  (strcmp(argv[2],"-perpDirn") != 0) && (strcmp(argv[2],"-dof") != 0)) {


	//
	// for legacy reasons we enter the first branch of the if
	//

	echoMode eMode = NONE;       // enum found in DataOutputFileHandler.h
	fileName = argv[2];
	
	int flag = 0;
	if (strcmp(argv[3],"-time") == 0) 
	  flag = 1;
	if (strcmp(argv[3],"-load") == 0)
	  flag = 2;      
	
	int pos = 3;
	if (flag != 0) pos = 4;
	
	int node1, node2, dof, perpDirn;
	
	if (Tcl_GetInt(interp, argv[pos++], &node1) != TCL_OK)	
	  return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[pos++], &node2) != TCL_OK)	
	  return TCL_ERROR;	      
	if (Tcl_GetInt(interp, argv[pos++], &dof) != TCL_OK)	
	  return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[pos++], &perpDirn) != TCL_OK)	
	  return TCL_ERROR;

	theDataOutputHandler = new DataOutputFileHandler(fileName, eMode);


	// Subtract one from dof and perpDirn for C indexing
	(*theRecorder) = new DriftRecorder(node1, node2, dof-1, perpDirn-1,
					   theDomain, *theDataOutputHandler);
      } else {

	echoMode eMode = NONE;       // enum found in DataOutputFileHandler.h

	bool echoTimeFlag = false;
	ID iNodes(0,16);
	ID jNodes(0,16);
	int dof = 1;
	int perpDirn = 2;
	int pos = 2;
	while (pos < argc) {

	  if (strcmp(argv[pos],"-file") == 0) {
	    fileName = argv[pos+1];
	    pos += 2;
	    if (strcmp(argv[pos],"-xml") == 0) {
	      eMode = XML_FILE;
	      pos+=1;
	    } else if (strcmp(argv[pos],"-headings") == 0) {
	      eMode = DATA_FILE;
	      pos +=1;
	    }
	  }

	  else if (strcmp(argv[pos],"-database") == 0) {
	    theRecorderDatabase = theDatabase;
	    if (theRecorderDatabase != 0)
	      tableName = argv[pos+1];
	    else {
	      opserr << "WARNING recorder Node .. -database &lt;fileName&gt; - NO CURRENT DATABASE, results to File instead\n";
	      fileName = argv[pos+1];
	    }
	    
	    pos += 2;
	  }
	  else if (strcmp(argv[pos],"-nees") == 0) {
	    // allow user to specify load pattern other than current
	    fileName = argv[pos+1];
	    eMode = XML_FILE;
	    pos += 2;
	  }	    

	  else if ((strcmp(argv[pos],"-iNode") == 0) || 
		 (strcmp(argv[pos],"-iNodes") == 0)) {
	    pos++;
	    
	    int node;
	    int numNodes = 0;
	    for (int j=pos; j< argc; j++) 
	      if (Tcl_GetInt(interp, argv[pos], &node) != TCL_OK) {
		j = argc;
		Tcl_ResetResult(interp);
	      } else {
		iNodes[numNodes] = node;
		numNodes++;
		pos++;
	      }
	  }

	  else if ((strcmp(argv[pos],"-jNode") == 0) || 
		 (strcmp(argv[pos],"-jNodes") == 0)) {
	    pos++;
	    int node;
	    int numNodes = 0;
	    for (int j=pos; j< argc; j++) 
	      if (Tcl_GetInt(interp, argv[pos], &node) != TCL_OK) {
		j = argc;
		Tcl_ResetResult(interp);
	      } else {
		jNodes[numNodes] = node;
		numNodes++;
		pos++;
	      }
	  }

	  else if (strcmp(argv[pos],"-dof") == 0) {
	    if (Tcl_GetInt(interp, argv[pos+1], &dof) != TCL_OK) {
	      pos = argc;
	    } 
	    pos+=2;
	  }
	  
	  else if (strcmp(argv[pos],"-perpDirn") == 0) {
	    if (Tcl_GetInt(interp, argv[pos+1], &perpDirn) != TCL_OK) {
	      pos = argc;
	    } 
	    pos+=2;
	  } 
	  
	  else if (strcmp(argv[pos],"-time") == 0) {
	    echoTimeFlag = true;
	    pos+=1;

	  } else 
	    pos++;
	}
	  
	if (iNodes.Size() != jNodes.Size()) {
	  opserr << "WARNING recorder Drift - the number of iNodes and jNodes must be the same " << iNodes << " " << jNodes << endln;
	  return TCL_ERROR;
	}

	// construct the DataHandler
	if (fileName != 0) {
	  theDataOutputHandler = new DataOutputFileHandler(fileName, eMode);
	} else if (tableName != 0) {
	  theDataOutputHandler = new DataOutputDatabaseHandler(theDatabase, tableName);
	} else
	  theDataOutputHandler = new DataOutputStreamHandler();
	
	// Subtract one from dof and perpDirn for C indexing
	(*theRecorder) = new DriftRecorder(iNodes, jNodes, dof-1, perpDirn-1,
					   theDomain, *theDataOutputHandler, echoTimeFlag);
      }
    }
    
    // a recorder for the graphical display of the domain
    else if (strcmp(argv[1],"display") == 0) {

	int xLoc, yLoc, width, height;

	if (argc < 7) {
	    opserr << "WARNING recorder display title xLoc yLoc pixelsX pixelsY <-file fileName?>";
	    return TCL_ERROR;
	}    
	if (Tcl_GetInt(interp, argv[3], &xLoc) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[4], &yLoc) != TCL_OK)	
	    return TCL_ERROR;	      
	if (Tcl_GetInt(interp, argv[5], &width) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[6], &height) != TCL_OK)	
	    return TCL_ERROR;
	
	// check if we are to wipe image on each redraw
	int wipeFlag = 0;
	if (argc == 8) 
	  if (strcmp(argv[7],"-wipe") == 0) 
	    wipeFlag = 1;


#ifdef _NOGRAPHICS
      return TCL_OK;
#else
	if (argc == 7 || argc == 8)
	  (*theRecorder) = new TclFeViewer(argv[2], xLoc, yLoc, width, height, theDomain, wipeFlag, interp);
	else if (argc == 9)
	  (*theRecorder) = new TclFeViewer(argv[2], xLoc, yLoc, width, height, argv[8], theDomain, interp);
#endif
    }

    else if (strcmp(argv[1],"plot") == 0) {

	int xLoc, yLoc, width, height;
	if (argc < 9) {
	    opserr << "WARNING recorder display fileName? windowTitle? xLoc yLoc pixelsX pixelsY -columns colX1 colY1 -columns colX2 ...";
	    return TCL_ERROR;
	}    
	if (Tcl_GetInt(interp, argv[4], &xLoc) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[5], &yLoc) != TCL_OK)	
	    return TCL_ERROR;	      
	if (Tcl_GetInt(interp, argv[6], &width) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[7], &height) != TCL_OK)	
	    return TCL_ERROR;	      

	int loc = 8;

	double dT = 0.0;
	loc = 0;
	ID cols(0,16);
	int numCols = 0;
	while (loc < argc) {
	  if ((strcmp(argv[loc],"-columns") == 0) ||
	      (strcmp(argv[loc],"-cols") == 0) ||
	      (strcmp(argv[loc],"-col") == 0)) {
	    if (argc < loc+2)
	      return TCL_ERROR;

	    int colX, colY;
	    if (Tcl_GetInt(interp, argv[loc+1], &colX) != TCL_OK)	
	      return TCL_ERROR;	

	    if (Tcl_GetInt(interp, argv[loc+2], &colY) != TCL_OK)	
	      return TCL_ERROR;	

	    cols[numCols++] = colX;
	    cols[numCols++] = colY;
	    loc += 3;
	  } 
	  else if (strcmp(argv[loc],"-dT") == 0) {

	    if (Tcl_GetDouble(interp, argv[loc+1], &dT) != TCL_OK)	
	      return TCL_ERROR;	
	    loc += 2;	    
	  }
	  else
	    loc++;
	}

#ifdef _NOGRAPHICS
	return TCL_OK;
#else
	FilePlotter *thePlotter = new FilePlotter(argv[2], argv[3], xLoc, yLoc, width, height, dT);
	(*theRecorder) = thePlotter;    
	thePlotter->setCol(cols);
#endif
    }

    else if (strcmp(argv[1],"plotDifferent") == 0) {

	int xLoc, yLoc, width, height;
	if (argc < 10) {
	    opserr << "WARNING recorder display fileName? windowTitle? xLoc yLoc pixelsX pixelsY -columns colX1 colY1 -columns colX2 ...";
	    return TCL_ERROR;
	}    
	if (Tcl_GetInt(interp, argv[5], &xLoc) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[6], &yLoc) != TCL_OK)	
	    return TCL_ERROR;	      
	if (Tcl_GetInt(interp, argv[7], &width) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[8], &height) != TCL_OK)	
	    return TCL_ERROR;	      

	int loc = 9;

	double dT = 0.0;
	loc = 0;
	ID cols(0,16);
	int numCols = 0;
	while (loc < argc) {
	  if ((strcmp(argv[loc],"-columns") == 0) ||
	      (strcmp(argv[loc],"-cols") == 0) ||
	      (strcmp(argv[loc],"-col") == 0)) {
	    if (argc < loc+2)
	      return TCL_ERROR;

	    int colX, colY;
	    if (Tcl_GetInt(interp, argv[loc+1], &colX) != TCL_OK)	
	      return TCL_ERROR;	

	    if (Tcl_GetInt(interp, argv[loc+2], &colY) != TCL_OK)	
	      return TCL_ERROR;	

	    cols[numCols++] = colX;
	    cols[numCols++] = colY;
	    loc += 3;
	  } 
	  else if (strcmp(argv[loc],"-dT") == 0) {

	    if (Tcl_GetDouble(interp, argv[loc+1], &dT) != TCL_OK)	
	      return TCL_ERROR;	
	    loc += 2;	    
	  }
	  else
	    loc++;
	}

#ifdef _NOGRAPHICS
	return TCL_OK;
#else
	FilePlotter *thePlotter = new FilePlotter(argv[2], argv[3], argv[4], xLoc, yLoc, width, height, dT);
	(*theRecorder) = thePlotter;    
	thePlotter->setCol(cols);
#endif
    }


    else if (strcmp(argv[1],"increments") == 0) {

	int xLoc, yLoc, width, height;
	
	if (theAlgorithm == 0) {
	    opserr << "WARNING recorder increments - only allowed as algorithmRecorder";
	    return TCL_ERROR;
	}
	if (argc < 7) {
	    opserr << "WARNING recorder display windowTitle? xLoc yLoc pixelsX pixelsY ";
	    return TCL_ERROR;
	}    
	if (Tcl_GetInt(interp, argv[3], &xLoc) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[4], &yLoc) != TCL_OK)	
	    return TCL_ERROR;	      
	if (Tcl_GetInt(interp, argv[5], &width) != TCL_OK)	
	    return TCL_ERROR;	
	if (Tcl_GetInt(interp, argv[6], &height) != TCL_OK)	
	    return TCL_ERROR;	      

	TCL_Char *fileName = 0;
	bool displayRecord = false;
	int loc = 7;
	while (loc < argc) {
	  if ((strcmp(argv[loc],"-file") == 0) ||
	      (strcmp(argv[loc],"-file") == 0)) {
	    
	      if (argc < loc+1)
		return TCL_ERROR;
	      loc++;
	      fileName = argv[loc];
	      loc++;
	  } else if (strcmp(argv[loc],"-display") == 0) {
	    displayRecord = true;
	    loc++;				   
	  } 
	  else
	    loc++;
	}

#ifdef _NOGRAPHICS
	return TCL_OK;
#else
	AlgorithmIncrements *thePlotter =  new AlgorithmIncrements(theAlgorithm, 
								   argv[2], xLoc, yLoc, width, height, 
								   displayRecord, fileName);
	(*theRecorder) = thePlotter;
#endif
    }

    else if (strcmp(argv[1],"GSA") == 0) {
	if (argc < 3) {
	  opserr << argc;
	  opserr << "WARNING recorder GSA -file fileName? -dT deltaT? - not enough arguments\n";
	  return TCL_ERROR;
	}    
	TCL_Char *fileName = 0;
	TCL_Char *title1 =0;
	TCL_Char *title2 =0;
	TCL_Char *title3 =0;
	TCL_Char *jobno =0;
	TCL_Char *initials =0;
	TCL_Char *spec =0;
	TCL_Char *currency =0;
	TCL_Char *length =0;
	TCL_Char *force =0;
	TCL_Char *temp =0;
	double dT = 0.0;
	int loc = 2;

	while (loc < argc) {
	  if ((strcmp(argv[loc],"-file") == 0) ||
	      (strcmp(argv[loc],"-file") == 0)) {
	    fileName = argv[loc+1];
	    loc += 2;
	  } else if ((strcmp(argv[loc],"-title1") == 0) ||
	      (strcmp(argv[loc],"-Title1e") == 0)) {
	    title1 = argv[loc+1];
	    loc += 2;
	  } else if ((strcmp(argv[loc],"-title2") == 0) ||
	      (strcmp(argv[loc],"-Title2e") == 0)) {
	    title2 = argv[loc+1];
	    loc += 2;
	  } else if ((strcmp(argv[loc],"-title3") == 0) ||
	      (strcmp(argv[loc],"-Title3e") == 0)) {
	    title3 = argv[loc+1];
	    loc += 2;
	  } else if ((strcmp(argv[loc],"-jobno") == 0) ||
	      (strcmp(argv[loc],"-JobNo") == 0)) {
	    jobno = argv[loc+1];
	    loc += 2;
	  } else if ((strcmp(argv[loc],"-initials") == 0) ||
	      (strcmp(argv[loc],"-Initials") == 0)) {
	    initials = argv[loc+1];
	    loc += 2;
	  } else if ((strcmp(argv[loc],"-spec") == 0) ||
	      (strcmp(argv[loc],"-Spec") == 0)) {
	    spec = argv[loc+1];
	    loc += 2;
	  } else if ((strcmp(argv[loc],"-currency") == 0) ||
	      (strcmp(argv[loc],"-Currency") == 0)) {
	    currency = argv[loc+1];
	    loc += 2;
	  } else if ((strcmp(argv[loc],"-length") == 0) ||
	      (strcmp(argv[loc],"-Length") == 0)) {
	    length = argv[loc+1];
	    loc += 2;
	  } else if ((strcmp(argv[loc],"-force") == 0) ||
	      (strcmp(argv[loc],"-Force") == 0)) {
	    force = argv[loc+1];
	    loc += 2;
	  } else if ((strcmp(argv[loc],"-temp") == 0) ||
	      (strcmp(argv[loc],"-Temp") == 0)) {
	    temp = argv[loc+1];
	    loc += 2;
	  }
	  else if (strcmp(argv[loc],"-dT") == 0) {
	    if (Tcl_GetDouble(interp, argv[loc+1], &dT) != TCL_OK)	
	      return TCL_ERROR;	      
	    loc += 2;
	  }
	  else
	    loc++;
	}

	GSA_Recorder *theR = new GSA_Recorder(theDomain, fileName, title1, title2, title3,
					      jobno, initials, spec, currency, length, force,
					      temp, dT);
	(*theRecorder) = theR;
    }
    
    else if (strcmp(argv[1],"YsVisual") == 0)
      { //!!
	
	int eleTag, xLoc, yLoc, width, height;
	double scale;
	
	if (argc < 7)
	  {
	    opserr << "WARNING recorder YsVisual eleTag title scale, xLoc yLoc pixelsX pixelsY \n";
	    return TCL_ERROR;
	  }
	
	if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK)
	  return TCL_ERROR;
	if (Tcl_GetDouble(interp, argv[4], &scale) != TCL_OK)
	  return TCL_ERROR;
	if (Tcl_GetInt(interp, argv[5], &xLoc) != TCL_OK)
	  return TCL_ERROR;
	if (Tcl_GetInt(interp, argv[6], &yLoc) != TCL_OK)
	  return TCL_ERROR;
	if (Tcl_GetInt(interp, argv[7], &width) != TCL_OK)
	  return TCL_ERROR;
	if (Tcl_GetInt(interp, argv[8], &height) != TCL_OK)
	  return TCL_ERROR;
	
	Element *theEle = theDomain.getElement(eleTag);
	if(theEle==0)
	  {
	    opserr << "WARNING no element of tag " << eleTag << "\n";
	    return TCL_ERROR;
	  }
	//void createView(char *title, int x, int y, int cx, int cy, char displaytype = 'l');
#ifdef _NOGRAPHICS
	return TCL_OK;
#else
	YsVisual *theVisual = new YsVisual(theEle, argv[3], scale, xLoc, yLoc, width, height);
	(*theRecorder) = theVisual;
	return TCL_OK;
#endif
      }
    
    // no recorder type specified yet exists
    else {
      opserr << "WARNING No recorder type exists ";
      opserr << "for recorder of type:" << argv[1];
      
      return TCL_ERROR;
    }    
    
    // check we instantiated a recorder .. if not ran out of memory
    if ((*theRecorder) == 0) {
	opserr << "WARNING ran out of memory - recorder " << argv[1]<< endln;
	return TCL_ERROR;
    } 

    // operation successfull
    return TCL_OK;
}


int 
TclAddRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
	       TCL_Char **argv, Domain &theDomain)
{
	Recorder *theRecorder;
	TclCreateRecorder(clientData, interp, argc, argv, theDomain, &theRecorder);
	
	if ((theRecorder == 0) || (theDomain.addRecorder(*theRecorder)) < 0) {
		opserr << "WARNING could not add to domain - recorder " << argv[1]<< endln;
		if (theRecorder == 0) 
			opserr << "could not create recorder\n";
		else
			delete theRecorder;
		return TCL_ERROR;
	} 
	return TCL_OK;
	
}


int 
TclAddAlgorithmRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
			TCL_Char **argv, Domain &theDomain, EquiSolnAlgo *theAlgo)
{
	Recorder *theRecorder = 0;
	theAlgorithm = theAlgo;
	if (TclCreateRecorder(clientData, interp, argc, argv, theDomain,
			&theRecorder) == TCL_ERROR) {
		return TCL_ERROR;
	} else {
		// add the recorder to the domain, 
		// NOTE: will not be called with theALgo == 0
		// see ~/g3/SRC/tcl/commands.C file
		if (theRecorder == 0 || theAlgo->addRecorder(*theRecorder) < 0) {
			opserr << "WARNING could not add to algorithm - recorder " << argv[1]<< endln;
			delete theRecorder;
			return TCL_ERROR;
		} 
		return TCL_OK;
	}
}


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

// $Revision$
// $Date$
// $URL$

// Written: fmk 
// Created: 11/00
// Revision: A
//
// Description: This file contains the function invoked when the user invokes
// the Pattern command in the interpreter. It is invoked by the 
// TclModelBuilder_addPattern function in the TclModelBuilder.C file. Current 
// valid Pattern types are:

#include <tcl.h>
#include <Domain.h>
#include <LinearSeries.h>
#include <ConstantSeries.h>
#include <RectangularSeries.h>
#include <TrigSeries.h>
#include <PulseSeries.h>
#include <TriangleSeries.h>
#include <PathTimeSeries.h>
#include <PathSeries.h>
#include <PeerMotion.h>
#include <PeerNGAMotion.h>
#include <string.h>

#ifdef _RELIABILITY
#include <DiscretizedRandomProcessSeries.h>
#include <SimulatedRandomProcessSeries.h>
#include <Spectrum.h>
#include <RandomNumberGenerator.h>
#include <ReliabilityDomain.h>
//#include <NewDiscretizedRandomProcessSeries.h>

extern ReliabilityDomain *theReliabilityDomain;
extern RandomNumberGenerator *theRandomNumberGenerator;
#endif

#include <SimulationInformation.h>
extern SimulationInformation simulationInfo;
//extern const char * getInterpPWD(Tcl_Interp *interp);  // commands.cpp

// little function to free memory after invoke Tcl_SplitList
//   note Tcl_Split list stores the array of pointers and the strings in 
//   one array, which is why Tcl_Free needs only be called on the array.
static void cleanup(TCL_Char **argv) {
  Tcl_Free((char *) argv);
}

extern void *OPS_ConstantSeries(void);
extern void *OPS_LinearSeries(void);
extern void *OPS_TriangleSeries(void);
extern void *OPS_TrigSeries(void);
extern void *OPS_RectangularSeries(void);
extern void *OPS_PulseSeries(void);
extern void *OPS_PeerMotion(void);
extern void *OPS_PeerNGAMotion(void);

#include <elementAPI.h>
#include <TclModelBuilder.h>

TimeSeries *
TclTimeSeriesCommand(ClientData clientData, 
		     Tcl_Interp *interp,
		     int argc, 
		     TCL_Char **argv,
		     Domain *theDomain)
{
  // note the 1 instead of usual 2
    OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, theDomain);
			    
  TimeSeries *theSeries = 0;

  if ((strcmp(argv[0],"Constant") == 0) || (strcmp(argv[0],"ConstantSeries") == 0)) {

    void *theResult = OPS_ConstantSeries();
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;
    
  } else if ((strcmp(argv[0],"Trig") == 0) || (strcmp(argv[0],"TrigSeries") == 0) || 
	     (strcmp(argv[0],"Sine") == 0) || (strcmp(argv[0],"SineSeries") == 0)) {

    void *theResult = OPS_TrigSeries();
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

  }	

  else if ((strcmp(argv[0],"Linear") == 0) || (strcmp(argv[0],"LinearSeries") == 0)) {

    void *theResult = OPS_LinearSeries();
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

  }

  else if (strcmp(argv[0],"Rectangular") == 0) {

    void *theResult = OPS_RectangularSeries();
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

  }

  else if ((strcmp(argv[0],"Pulse") == 0) || (strcmp(argv[0],"PulseSeries") == 0))  {

    void *theResult = OPS_PulseSeries();
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

    }

  else if ((strcmp(argv[0],"Triangle") == 0) || (strcmp(argv[0],"TriangleSeries") == 0))  {

    void *theResult = OPS_TriangleSeries();
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

  }
  
  else if ((strcmp(argv[0],"Series") == 0) ||
	   (strcmp(argv[0],"Path") == 0)) {

    // LoadPattern and PathSeries - read args and create RectangularSeries object
    double cFactor = 1.0;

    if (argc < 3) {
      opserr << "WARNING not enough args - ";
      opserr << " Series -dt timeIncr -values {list of points }\n"; 
      return 0;	
    }

    int tag = 0;
    double timeIncr = 0.0;
    int endMarker =  1;
    bool done = false;
    int fileName = 0;
    int fileTimeName = 0;
    int filePathName = 0;
    Vector *dataPath = 0;
    Vector *dataTime = 0;
    bool useLast = false;
    bool prependZero = false;
    double startTime = 0.0;

    if (Tcl_GetInt(interp, argv[endMarker], &tag) == TCL_OK) {
      endMarker++;
    }

    while (endMarker < argc && done == false) {
	
      if (strcmp(argv[endMarker],"-dt") == 0) {
	// allow user to specify the time increment
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetDouble(interp, argv[endMarker], &timeIncr) != TCL_OK) {
	  
	  opserr << "WARNING invalid dt " << argv[endMarker] << " - ";
	  opserr << " Series -dt dt ... \n";
	  return 0;
	}
      } 

      else if (strcmp(argv[endMarker],"-tag") == 0) {
	// allow user to specify the tag
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetInt(interp, argv[endMarker], &tag) != TCL_OK) {
	  
	  opserr << "WARNING invalid tag " << argv[endMarker] << " - ";
	  return 0;
	}
      } 

      else if (strcmp(argv[endMarker],"-factor") == 0) {
	// allow user to specify the factor
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	  
	  opserr << "WARNING invalid cFactor " << argv[endMarker] << " - ";
	  opserr << " Series -factor ... \n";
	  return 0;
	}
      } 

      else if (strcmp(argv[endMarker],"-file") == 0) {
	// allow user to specify the file name containing time and data points
	endMarker++;
	if (endMarker != argc) {
	  fileName = endMarker; // argv[endMarker];
	}
      }

      else if (strcmp(argv[endMarker],"-filePath") == 0) {
	// allow user to specify the file name containing the data points
	endMarker++;
	if (endMarker != argc) {
	  filePathName = endMarker; // argv[endMarker];
	}
      }

      else if (strcmp(argv[endMarker],"-fileTime") == 0) {
	// allow user to specify the file name containing the data points
	endMarker++;
	if (endMarker != argc) {
	  fileTimeName = endMarker; // argv[endMarker];
	}
      }

      else if (strcmp(argv[endMarker],"-values") == 0) {
	// allow user to specify the data points in tcl list
	endMarker++;
	if (endMarker != argc) {
	  int pathSize;
	  TCL_Char **pathStrings;
	  
	  if (Tcl_SplitList(interp, argv[endMarker], 
			    &pathSize, &pathStrings) != TCL_OK) {
		      
	    opserr << "WARNING problem splitting path list " << argv[endMarker] << " - ";
	    opserr << " Series -values {path} ... \n";
	    return 0;
	  }
	  
	  dataPath = new Vector(pathSize);
	  for (int i=0; i<pathSize; i++) {
	    double value;
	    if ( Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
	      opserr << "WARNING problem reading path data value " << pathStrings[i] << " - ";
	      opserr << " Series -values {path} ... \n";
	      cleanup(pathStrings);
	      return 0;
	    }
	    (*dataPath)(i) = value;
	  }
	  // free up the array of pathsStrings .. see tcl man pages as to why
	  cleanup(pathStrings);
	}
      }

      else if (strcmp(argv[endMarker],"-time") == 0) {
	// allow user to specify the data points in tcl list
	endMarker++;
	if (endMarker != argc) {
	  int pathSize;
	  TCL_Char **pathStrings;
	  
	  if (Tcl_SplitList(interp, argv[endMarker], 
			    &pathSize, &pathStrings) != TCL_OK) {
			
	    opserr << "WARNING problem spltting time path " << argv[endMarker] << " - ";
	    opserr << " Series -time {times} ... \n";
	    return 0;
	  }
	  
	  dataTime = new Vector(pathSize);
	  for (int i=0; i<pathSize; i++) {
	    double value;
	    if ( Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
	      opserr << "WARNING problem reading time path value " << pathStrings[i] << " - ";
	      opserr << " Series -values {path} ... \n";

	      cleanup(pathStrings);
	      return 0;
	    }
	    (*dataTime)(i) = value;
	  }
	  // free up the array of pathsStrings .. see tcl man pages as to why
	  cleanup(pathStrings);
	}
      }
      
      else if (strcmp(argv[endMarker],"-useLast") == 0) {
	    useLast = true;
      }

      else if (strcmp(argv[endMarker],"-prependZero") == 0) {
	    prependZero = true;
      }

      else if (strcmp(argv[endMarker],"-startTime") == 0 ||
          strcmp(argv[endMarker],"-tStart") == 0) {
	// allow user to specify the start time
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetDouble(interp, argv[endMarker], &startTime) != TCL_OK) {
	  
	  opserr << "WARNING invalid tStart " << argv[endMarker] << " - ";
	  opserr << " Series -startTime tStart ... \n";
	  return 0;
	}
      } 

      endMarker++;
    }

    if (filePathName != 0 && fileTimeName == 0 && timeIncr != 0.0) {
      //      const char *pwd = getInterpPWD(interp);
      //      simulationInfo.addInputFile(argv[filePathName], pwd);  
      theSeries = new PathSeries(tag, argv[filePathName], timeIncr, cFactor,
                                 useLast, prependZero, startTime);
    }

    else if (fileName != 0) {
      //const char *pwd = getInterpPWD(interp);
      //      simulationInfo.addInputFile(argv[fileName], pwd);  
      theSeries = new PathTimeSeries(tag, argv[fileName], cFactor, useLast);

    } else if (filePathName != 0 && fileTimeName != 0) {
      //      const char *pwd = getInterpPWD(interp);
      //      simulationInfo.addInputFile(argv[filePathName], pwd);  
      //      simulationInfo.addInputFile(argv[fileTimeName], pwd);  
      theSeries = new PathTimeSeries(tag, argv[filePathName], argv[fileTimeName],
                                     cFactor, useLast); 

    } else if (dataPath != 0 && dataTime == 0 && timeIncr != 0.0) {
      theSeries = new PathSeries(tag, *dataPath, timeIncr, cFactor,
                                 useLast, prependZero, startTime);
      delete dataPath;

    } else if (dataPath != 0 && dataTime != 0) {
      theSeries = new PathTimeSeries(tag, *dataPath, *dataTime, cFactor, useLast);  
      delete dataPath;
      delete dataTime;

    } else {
      opserr << "WARNING choice of options for Path Series invalid - valid options for ";
      opserr << " Path are\n";
      opserr << " \t -fileT fileTimeName -fileP filePathName \n";
      opserr << " \t -dt constTimeIncr -file filePathName\n";
      opserr << " \t -dt constTimeIncr -values {list of points on path}\n";
      opserr << " \t -time {list of time points} -values {list of points on path}\n";
      return 0;
    }      
	
  } 

  else if ((strcmp(argv[0],"PeerDatabase") == 0) || (strcmp(argv[0],"PeerMotion") == 0)) {

    void *theResult = OPS_PeerMotion();
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;
    
    PeerMotion *thePeerMotion = (PeerMotion *)theSeries;       	

    if (argc > 4 && theSeries != 0) {
      int argCount = 4;

      while (argCount+1 < argc) {
	if ((strcmp(argv[argCount],"-dT") == 0) || (strcmp(argv[argCount],"-dt") == 0) || 
	    (strcmp(argv[argCount],"-DT") == 0)) {
	  const char *variableName = argv[argCount+1];
	  double dT = thePeerMotion->getDt();
	  char string[30];
	  sprintf(string,"set %s %.18e", variableName, dT);   Tcl_Eval(interp, string);
	  argCount+=2;
	} else if ((strcmp(argv[argCount],"-nPts") == 0) || (strcmp(argv[argCount],"-NPTS") == 0)) {
	  const char *variableName = argv[argCount+1];
	  int nPts = thePeerMotion->getNPts();
	  char string[30];
	  sprintf(string,"set %s %d", variableName, nPts);   Tcl_Eval(interp, string);
	  argCount+=2;
	} else
	  argCount++;
      }
    }
  }

  else if ((strcmp(argv[0],"PeerNGADatabase") == 0) || (strcmp(argv[0],"PeerNGAMotion") == 0)) {


    void *theResult = OPS_PeerNGAMotion();
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;
    
    PeerNGAMotion *thePeerMotion = (PeerNGAMotion*)(theSeries);       	


    if (argc > 3 && theSeries != 0) {
      int argCount = 3;

      while (argCount+1 < argc) {
	if ((strcmp(argv[argCount],"-dT") == 0) || (strcmp(argv[argCount],"-dt") == 0) || 
	    (strcmp(argv[argCount],"-DT") == 0)) {
	  const char *variableName = argv[argCount+1];
	  double dT = thePeerMotion->getDt();
	  char string[30];
	  sprintf(string,"set %s %.18e", variableName, dT);   Tcl_Eval(interp, string);
	  argCount+=2;

	} else if ((strcmp(argv[argCount],"-nPts") == 0) || (strcmp(argv[argCount],"-NPTS") == 0)) {
	  const char *variableName = argv[argCount+1];
	  int nPts = thePeerMotion->getNPts();
	  char string[30];
	  sprintf(string,"set %s %d", variableName, nPts);   Tcl_Eval(interp, string);

	  argCount+=2;
	} else
	  argCount++;
      }
    }
  }

#ifdef _RELIABILITY

  else if (strcmp(argv[0],"DiscretizedRandomProcess") == 0) {
    
    double mean, maxStdv;
    ModulatingFunction *theModFunc;
    
    if (Tcl_GetDouble(interp, argv[1], &mean) != TCL_OK) {
      opserr << "WARNING invalid input: random process mean \n";
      return 0;
    }
    
    if (Tcl_GetDouble(interp, argv[2], &maxStdv) != TCL_OK) {
      opserr << "WARNING invalid input: random process max stdv \n";
      return 0;
    }
    
    // Number of modulating functions
    int argsBeforeModList = 3;
    int numModFuncs = argc-argsBeforeModList;
    
    // Create an array to hold pointers to modulating functions
    ModulatingFunction **theModFUNCS = new ModulatingFunction *[numModFuncs];
    
    // For each modulating function, get the tag and ensure it exists
    int tagI;
    for (int i=0; i<numModFuncs; i++) {
      if (Tcl_GetInt(interp, argv[i+argsBeforeModList], &tagI) != TCL_OK) {
	opserr << "WARNING invalid modulating function tag. " << endln;
	return 0;
      }
      
      theModFunc = 0;
      theModFunc = theReliabilityDomain->getModulatingFunction(tagI);
      
      if (theModFunc == 0) {
	opserr << "WARNING modulating function number "<< argv[i+argsBeforeModList] << "does not exist...\n";
	delete [] theModFUNCS;
	return 0;
      }
      else {
	theModFUNCS[i] = theModFunc;
      }
    }	
    
    // Parsing was successful, create the random process series object
    theSeries = new DiscretizedRandomProcessSeries(0,numModFuncs,theModFUNCS,mean,maxStdv);       	
  }
  
   ///// added by K Fujimura /////
   /*FMK RELIABILITY
  else if (strcmp(argv[0],"NewDiscretizedRandomProcess") == 0) {
    
    double mean, maxStdv;
    ModulatingFunction *theModFunc;
    
    if (Tcl_GetDouble(interp, argv[1], &mean) != TCL_OK) {
      opserr << "WARNING invalid input: random process mean \n";
      return 0;
    }
    
    if (Tcl_GetDouble(interp, argv[2], &maxStdv) != TCL_OK) {
      opserr << "WARNING invalid input: random process max stdv \n";
      return 0;
    }
    
    // Number of modulating functions
    int argsBeforeModList = 3;
    int numModFuncs = argc-argsBeforeModList;
    
    // Create an array to hold pointers to modulating functions
    ModulatingFunction **theModFUNCS = new ModulatingFunction *[numModFuncs];
    
    // For each modulating function, get the tag and ensure it exists
    int tagI;
    for (int i=0; i<numModFuncs; i++) {
      if (Tcl_GetInt(interp, argv[i+argsBeforeModList], &tagI) != TCL_OK) {
	opserr << "WARNING invalid modulating function tag. " << endln;
	return 0;
      }
      
      theModFunc = 0;
      theModFunc = theReliabilityDomain->getModulatingFunction(tagI);
      
      if (theModFunc == 0) {
	opserr << "WARNING modulating function number "<< argv[i+argsBeforeModList] << "does not exist...\n";
	delete [] theModFUNCS;
	return 0;
      }
      else {
	theModFUNCS[i] = theModFunc;
      }
    }	
    
    // Parsing was successful, create the random process series object
  theSeries = new NewDiscretizedRandomProcessSeries(numModFuncs,theModFUNCS,mean,maxStdv);       	
  }
  */
  else if (strcmp(argv[0],"SimulatedRandomProcess") == 0) {
    
    int spectrumTag, numFreqIntervals;
    double mean;
    
    if (Tcl_GetInt(interp, argv[1], &spectrumTag) != TCL_OK) {
      opserr << "WARNING invalid input to SimulatedRandomProcess: spectrumTag" << endln;
      return 0;
    }
    
    if (Tcl_GetDouble(interp, argv[2], &mean) != TCL_OK) {
      opserr << "WARNING invalid input to SimulatedRandomProcess: mean" << endln;
      return 0;
    }
    
    if (Tcl_GetInt(interp, argv[3], &numFreqIntervals) != TCL_OK) {
      opserr << "WARNING invalid input to SimulatedRandomProcess: numFreqIntervals" << endln;
      return 0;
    }
    
    // Check that the random number generator exists
    if (theRandomNumberGenerator == 0) {
      opserr << "WARNING: A random number generator must be instantiated before SimulatedRandomProcess." << endln;
      return 0;
    }
    
    // Check that the spectrum exists
    Spectrum *theSpectrum = 0;
    theSpectrum = theReliabilityDomain->getSpectrum(spectrumTag);
    if (theSpectrum == 0) {
      opserr << "WARNING: Could not find the spectrum for the SimulatedRandomProcess." << endln;
      return 0;
    }
    
    // Parsing was successful, create the random process series object
    theSeries = new SimulatedRandomProcessSeries(0,theRandomNumberGenerator,theSpectrum,numFreqIntervals,mean);
  }

#endif

  else {
    for (int i = 0; i < argc; i++)
      opserr << argv[i] << ' ';
    opserr << endln;
    // type of load pattern type unknown
    opserr << "WARNING unknown Series type " << argv[0] << " - ";
    opserr << " valid types: Linear, Rectangular, Path, Constant, Trig, Sine\n";
    return 0;
  }

  return theSeries;
}

TimeSeries *
TclSeriesCommand(ClientData clientData, Tcl_Interp *interp, TCL_Char *arg)
{
  int argc;
  TCL_Char **argv;

  int timeSeriesTag = 0;
  if (Tcl_GetInt(interp, arg, &timeSeriesTag) == TCL_OK) {
    return OPS_getTimeSeries(timeSeriesTag);
  }

  // split the list
  if (Tcl_SplitList(interp, arg, &argc, &argv) != TCL_OK) {
    opserr << "WARNING could not split series list " << arg << endln;
    return 0;
  }

  TimeSeries *theSeries = TclTimeSeriesCommand(clientData, interp, argc, argv, 0);

  // clean up after ourselves and return the series
  cleanup(argv);
  return theSeries;
}

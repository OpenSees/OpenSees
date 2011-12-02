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
// $Date: 2000-12-14 08:43:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TclSeriesCommand.cpp,v $

// Written: fmk 
// Created: 11/00
// Revision: A
//
// Description: This file contains the function invoked when the user invokes
// the Pattern command in the interpreter. It is invoked by the 
// TclModelBuilder_addPattern function in the TclModelBuilder.C file. Current 
// valid Pattern types are:

// What: "@(#) TclPatternCommand.C, revA"

#include <tcl.h>
#include <Domain.h>
#include <LinearSeries.h>
#include <ConstantSeries.h>
#include <RectangularSeries.h>
#include <TrigSeries.h>
#include <PathTimeSeries.h>
#include <PathSeries.h>
#include <string.h>


// little function to free memory after invoke Tcl_SplitList
//   note Tcl_Split list stores the array of pointers and the strings in 
//   one array, which is why Tcl_Free needs only be called on the array.
static void cleanup(char **argv) {
#ifdef TCL_Free	    
	  Tcl_Free((char *) argv);
#endif	    
}


TimeSeries *
TclSeriesCommand(ClientData clientData, Tcl_Interp *interp, char *arg)
{

  int argc;
  char **argv;

  // split the list
  if (Tcl_SplitList(interp, arg, &argc, &argv) != TCL_OK) {
    cerr << "WARNING could not split series list " << arg << endl;
    return 0;
  }
			    
  TimeSeries *theSeries = 0;

  if (strcmp(argv[0],"Constant") == 0) {
    // LoadPattern and ConstantSeries - read args & create LinearSeries object
    double cFactor = 1.0;
	
    int endMarker = 1;
    if ((endMarker != argc) && (strcmp(argv[endMarker],"-factor") == 0)) {
      // allow user to specify the factor
      endMarker++;
      if (endMarker == argc || 
	  Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	    
	cerr << "WARNING invalid cFactor " << argv[endMarker] << " - ";
	cerr << " Constant -factor cFactor\n";
	cleanup(argv);
	return 0;
      }
      endMarker++;
    }

    theSeries = new ConstantSeries(cFactor);       	
    
  } else if (strcmp(argv[0],"Trig") == 0 || 
	     strcmp(argv[0],"Sine") == 0) {
    // LoadPattern and TrigSeries - read args & create TrigSeries object
    double cFactor = 1.0;
    double tStart, tFinish, period;
    double shift = 0.0;
      
    if (argc < 4) {
      cerr << "WARNING not enough TimeSeries args - ";
      cerr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
      cleanup(argv);
      return 0;	
    }	
    if (Tcl_GetDouble(interp, argv[1], &tStart) != TCL_OK) {
      cerr << "WARNING invalid tStart " << argv[1] << " - ";
      cerr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
      cleanup(argv);
      return 0;				
    }
    if (Tcl_GetDouble(interp, argv[2], &tFinish) != TCL_OK) {
      cerr << "WARNING invalid tFinish " << argv[2] << " - ";
      cerr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
      cleanup(argv);
      return 0;	
    }     
    if (Tcl_GetDouble(interp, argv[3], &period) != TCL_OK) {
      cerr << "WARNING invalid period " << argv[3] << " - ";
      cerr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
      cleanup(argv);
      return 0;	
    }     
    
    int endMarker = 4;
    
    while (endMarker < argc && endMarker < argc) {
      if (strcmp(argv[endMarker],"-factor") == 0) {
	// allow user to specify the factor
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	  
	  cerr << "WARNING invalid cFactor " << argv[endMarker] << " -";
	  cerr << " Trig  tStart tFinish period -factor cFactor\n";
	  cleanup(argv);
	  return 0;
	}
      }

      else if (strcmp(argv[endMarker],"-shift") == 0) {
	// allow user to specify phase shift
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetDouble(interp, argv[endMarker], &shift) != TCL_OK) {
	    
	  cerr << "WARNING invalid phase shift " << argv[endMarker] << " - ";
	  cerr << " Trig tStart tFinish period -shift shift\n";
	  cleanup(argv);
	  return 0;
	}
      }
      endMarker++;
    }

    theSeries = new TrigSeries(tStart, tFinish, period, shift, cFactor);
	
  }	

  else if (strcmp(argv[0],"Linear") == 0) {
    // LoadPattern and LinearSeries - read args & create LinearSeries object
    double cFactor = 1.0;
      
    int endMarker = 1;
    
    if ((endMarker < argc) && (strcmp(argv[endMarker],"-factor") == 0)) {
      // allow user to specify the factor
      endMarker++;
      if (endMarker == argc || 
	  Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	
	cerr << "WARNING invalid cFactor " << argv[endMarker] << " - ";
	cerr << " Linear  -factor cFactor\n";
	cleanup(argv);
	return 0;
      }	
      endMarker++;
    }

    theSeries = new LinearSeries(cFactor);       	
  }

  else if (strcmp(argv[0],"Rectangular") == 0) {
    // LoadPattern and RectangularSeries - read args and create RectangularSeries object
    double tStart, tFinish;
    double cFactor = 1.0;
    if (argc < 3) {
      cerr << "WARNING not enogh args - ";
      cerr << " Rectangular tStart tFinish <-factor cFactor>\n";
      cleanup(argv);
      return 0;	
    }	
    if (Tcl_GetDouble(interp, argv[1], &tStart) != TCL_OK) {
      cerr << "WARNING invalid tStart " << argv[1] << " - ";
      cerr << " Rectangular tStart tFinish <-factor factor>\n";
      cleanup(argv);
      return 0;
    }
    if (Tcl_GetDouble(interp, argv[2], &tFinish) != TCL_OK) {
      cerr << "WARNING invalid tStart " << argv[2] << " - ";
      cerr << " Rectangular tStart tFinish <-factor fcator>\n";
      cleanup(argv);
      return 0;	
    }     

    int endMarker =  3;
    
    if ((endMarker != argc) && (strcmp(argv[endMarker],"-factor") == 0)) {
      // allow user to specify the factor
      endMarker++;
      if (endMarker == argc || 
	  Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	
	cerr << "WARNING invalid cFactor " << argv[endMarker] << " - ";
	cerr << " Rectangular tStart tFinish -factor cFactor\n";
	cleanup(argv);
	return 0;
      }
      endMarker++;
    }  

    theSeries = new RectangularSeries(tStart, tFinish, cFactor); 
      
  }
    
  else if ((strcmp(argv[0],"Series") == 0) ||
	   (strcmp(argv[0],"Path") == 0)) {

    // LoadPattern and PathSeries - read args and create RectangularSeries object
    double cFactor = 1.0;
    if (argc < 3) {

      cerr << "WARNING not enough args - ";
      cerr << " Series -dt timeIncr -values {list of points }\n"; 
      cleanup(argv);
      return 0;	
    }

    double timeIncr = 0.0;
    int endMarker =  1;
    bool done = false;
    char *fileTimeName = 0;
    char *filePathName = 0;
    Vector *dataPath = 0;
    Vector *dataTime = 0;

    while (endMarker < argc && done == false) {
	
      if (strcmp(argv[endMarker],"-dt") == 0) {
	// allow user to specify the time increment
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetDouble(interp, argv[endMarker], &timeIncr) != TCL_OK) {
	  
	  cerr << "WARNING invalid dt " << argv[endMarker] << " - ";
	  cerr << " Series -dt dt ... \n";
	  cleanup(argv);
	  return 0;
	}
      } 

      else if (strcmp(argv[endMarker],"-factor") == 0) {
	// allow user to specify the factor
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	  
	  cerr << "WARNING invalid cFactor " << argv[endMarker] << " - ";
	  cerr << " Series -factor ... \n";
	  cleanup(argv);
	  return 0;
	}
      } 

      else if (strcmp(argv[endMarker],"-filePath") == 0) {
	// allow user to specify the file name containg the data points
	endMarker++;
	if (endMarker != argc) {
	  filePathName = argv[endMarker];
	}
      }

      else if (strcmp(argv[endMarker],"-fileTime") == 0) {
	// allow user to specify the file name containg the data points
	endMarker++;
	if (endMarker != argc) {
	  fileTimeName = argv[endMarker];
	}
      }

      else if (strcmp(argv[endMarker],"-values") == 0) {
	// allow user to specify the data points in tcl list
	endMarker++;
	if (endMarker != argc) {
	  int pathSize;
	  char **pathStrings;
	  
	  if (Tcl_SplitList(interp, argv[endMarker], 
			    &pathSize, &pathStrings) != TCL_OK) {
		      
	    cerr << "WARNING problem splitting path list " << argv[endMarker] << " - ";
	    cerr << " Series -values {path} ... \n";
	    cleanup(argv);
	    return 0;
	  }
	  
	  dataPath = new Vector(pathSize);
	  for (int i=0; i<pathSize; i++) {
	    double value;
	    if ( Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
	      cerr << "WARNING problem reading path data value " << pathStrings[i] << " - ";
	      cerr << " Series -values {path} ... \n";
	      cleanup(argv);
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
	  char **pathStrings;
	  
	  if (Tcl_SplitList(interp, argv[endMarker], 
			    &pathSize, &pathStrings) != TCL_OK) {
			
	    cerr << "WARNING problem spltting time path " << argv[endMarker] << " - ";
	    cerr << " Series -time {times} ... \n";
	    cleanup(argv);
	    return 0;
	  }
	  
	  dataTime = new Vector(pathSize);
	  for (int i=0; i<pathSize; i++) {
	    double value;
	    if ( Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
	      cerr << "WARNING problem reading time path value " << pathStrings[i] << " - ";
	      cerr << " Series -values {path} ... \n";

	      cleanup(argv);
	      cleanup(pathStrings);
	      return 0;
	    }
	    (*dataTime)(i) = value;
	  }
	  // free up the array of pathsStrings .. see tcl man pages as to why
	  cleanup(pathStrings);
	}
      }
      
      endMarker++;
    }
    

    if (filePathName != 0 && fileTimeName == 0 && timeIncr != 0.0)
      theSeries = new PathSeries(filePathName, timeIncr, cFactor);
    else if (filePathName != 0 && fileTimeName != 0)
      theSeries = new PathTimeSeries(filePathName, fileTimeName, cFactor); 
    else if (dataPath != 0 && dataTime == 0 && timeIncr != 0.0) {
      theSeries = new PathSeries(*dataPath, timeIncr, cFactor); 
      delete dataPath;
    } else if (dataPath != 0 && dataTime != 0) {
      theSeries = new PathTimeSeries(*dataPath, *dataTime, cFactor);  
      delete dataPath;      
      delete dataTime;      
    } else {
      cerr << "WARNING choice of options for Path Series invalid - valid options for ";
      cerr << " Path are\n";
      cerr << " \t -fileT fileTimeName -fileP filePathName \n";
      cerr << " \t -dt constTimeIncr -file filePathName\n";
      cerr << " \t -dt constTimeIncr -values {list of points on path}\n";
      cerr << " \t -time {list of time points} -values {list of points on path}\n";
      cleanup(argv);
      return 0;
    }      
	
  } 
	
  else {
    // type of load pattern type unknown
    cerr << "WARNING unknown Series type " << argv[0] << " - ";
    cerr << " valid types: Linear, Rectangular, Path, Constant, Trig, Sine\n";
    cleanup(argv);
    return 0;
  }

  // clean up after ourselves and return the series
  cleanup(argv);
  return theSeries;
}





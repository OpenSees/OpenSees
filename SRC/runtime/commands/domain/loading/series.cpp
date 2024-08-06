//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the function invoked when the user invokes
// the pattern command in the interpreter. It is invoked by the
// TclBasicBuilder_addPattern function.
//
// Written: fmk
// Created: 11/00
//
#include <sys/types.h>
#include <sys/stat.h>

#include <tcl.h>
#include <string.h>
#include <runtimeAPI.h>
#include <Logging.h>
#include <Parsing.h>
#include <BasicModelBuilder.h>
// #include <elementAPI.h>

#include <Domain.h>
#include <LinearSeries.h>
#include <ConstantSeries.h>
#include <PathTimeSeries.h>
#include <PathSeries.h>
#include <TrigSeries.h>
#include <RectangularSeries.h>
#include <PulseSeries.h>
#include <TriangleSeries.h>
// #include <PeerMotion.h>
// #include <PeerNGAMotion.h>



// extern OPS_Routine OPS_ConstantSeries;
// extern OPS_Routine OPS_LinearSeries;
// extern OPS_Routine OPS_TrigSeries;
// extern OPS_Routine OPS_PulseSeries;
// extern OPS_Routine OPS_PeerMotion;
// extern OPS_Routine OPS_PeerNGAMotion;
// extern OPS_Routine OPS_TriangleSeries;
// extern OPS_Routine OPS_RectangularSeries;

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp,
                                       int cArg, int mArg, TCL_Char ** const argv,
                                       Domain *domain);


static void *
TclDispatch_newLinearSeries(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
  int numRemainingArgs = argc;

  int tag = 0;
  double cFactor = 1.0;

  if (numRemainingArgs != 0) {

    if (numRemainingArgs == 1 || numRemainingArgs == 3) {
      if (Tcl_GetInt(interp, argv[0], &tag) != 0) {
        opserr << G3_ERROR_PROMPT << "invalid series tag in LinearSeries tag? <-factor "
                  "factor?>"
               << "\n";
        return nullptr;
      }
      numRemainingArgs--;
    }

    if (numRemainingArgs > 1) {
      const char *argvS = argv[1];
      if (argvS == 0) {
        opserr << G3_ERROR_PROMPT << "string error in LinearSeries with tag: " << tag
               << "\n";
        return nullptr;
      }
      if (Tcl_GetDouble(interp, argv[2], &cFactor) != 0) {
        opserr << G3_ERROR_PROMPT << "invalid factor in LinearSeries with tag: " << tag
               << "\n";
        return nullptr;
      }
    }
  }

  return new LinearSeries(tag, cFactor);
}

static TimeSeries *
TclDispatch_newTimeSeries(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{

  // note the 1 instead of usual 2
  OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, nullptr);

  TimeSeries *theSeries = nullptr;

  if ((strcmp(argv[0], "Constant") == 0) ||
      (strcmp(argv[0], "ConstantSeries") == 0)) {
    // LoadPattern and ConstantSeries - read args & create LinearSeries object
    double cFactor = 1.0;
        
    int endMarker = 1;
    if ((endMarker != argc) && (strcmp(argv[endMarker],"-factor") == 0)) {
      // allow user to specify the factor
      endMarker++;
      if (endMarker == argc || 
          Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
            
        opserr << "WARNING invalid cFactor " << argv[endMarker] << " - ";
        opserr << " Constant -factor cFactor\n";
        Tcl_Free((char*)argv);
        return 0;
      }
      endMarker++;
    }

    theSeries = new ConstantSeries(cFactor);            
    

//  void *theResult = OPS_ConstantSeries(rt, argc, argv);
//  if (theResult != nullptr)
//    theSeries = (TimeSeries *)theResult;
  }

  else if (strcmp(argv[0],"Trig") == 0 || 
             strcmp(argv[0],"Sine") == 0) {
     // LoadPattern and TrigSeries - read args & create TrigSeries object
     int tag = 0;
     double cFactor = 1.0;
     double tStart, tFinish, period;
     double shift = 0.0;
       
     if (argc < 3) {
       opserr << "WARNING not enough TimeSeries args - ";
       opserr << " Trig <tag?> tStart tFinish period <-shift shift> <-factor cFactor>\n";
       return nullptr;
     }
     int argi = 1;

     if (argc == 5 || argc == 7 || argc == 9 || argc == 11) {
      if (Tcl_GetInt(interp, argv[argi++], &tag) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid series tag in Trig tag?" << "\n";
        return nullptr;
      }
     }

     if (Tcl_GetDouble(interp, argv[argi++], &tStart) != TCL_OK) {
       opserr << "WARNING invalid tStart " << argv[argi-1] << " - ";
       opserr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
       return nullptr;
     }

     if (Tcl_GetDouble(interp, argv[argi++], &tFinish) != TCL_OK) {
       opserr << "WARNING invalid tFinish " << argv[argi-1] << " - ";
       opserr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
       return nullptr; 
     }

     if (Tcl_GetDouble(interp, argv[argi++], &period) != TCL_OK) {
       opserr << "WARNING invalid period " << argv[argi-1] << " - ";
       opserr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
       return nullptr; 

     } else if (period == 0.0) {
       opserr << G3_WARN_PROMPT << "Period for '" << argv[0] << "' is zero.\n";
     }
      
     while (argi < argc) {
       if (strcmp(argv[argi], "-factor") == 0) {
         // allow user to specify a scaling factor
         argi++;
         if (argi == argc || 
             Tcl_GetDouble(interp, argv[argi], &cFactor) != TCL_OK) {
           
           opserr << "WARNING invalid cFactor " << argv[argi] << " -";
           opserr << " Trig  tStart tFinish period -factor cFactor\n";
           return nullptr;
         }
       }
 
       else if (strcmp(argv[argi],"-shift") == 0) {
         // allow user to specify phase shift
         argi++;
         if (argi == argc || 
             Tcl_GetDouble(interp, argv[argi], &shift) != TCL_OK) {
             
           opserr << "WARNING invalid phase shift " << argv[argi] << " - ";
           opserr << " Trig tStart tFinish period -shift shift\n";
           return nullptr;
         }
       }
       argi++;
     }
 
     theSeries = new TrigSeries(tag, tStart, tFinish, period, shift, cFactor);
         
   }

#if 0
   else if ((strcmp(argv[0], "Trig") == 0) ||
             (strcmp(argv[0], "TrigSeries") == 0) ||
             (strcmp(argv[0], "Sine") == 0) ||
             (strcmp(argv[0], "SineSeries") == 0)) {

    void *theResult = OPS_TrigSeries(rt, argc, argv);
    if (theResult != nullptr)
      theSeries = (TimeSeries *)theResult;

  }
#endif

  else if ((strcmp(argv[0], "Linear") == 0) ||
           (strcmp(argv[0], "LinearSeries") == 0)) {

    void *theResult = TclDispatch_newLinearSeries(clientData, interp, argc - 1, &argv[1]);

    if (theResult != nullptr)
      theSeries = (TimeSeries *)theResult;
    else
      opserr << "ERROR\n";

  }

  else if (strcmp(argv[0],"Pulse") == 0)  {
    // LoadPattern and PulseSeries - read args & create PulseSeries object
    double cFactor = 1.0;
    double tStart, tFinish, period;
    double width = 0.5;
    double shift = 0.0;
      
    if (argc < 4) {
      opserr << "WARNING not enough PulseSeries args - ";
      opserr << " Pulse tStart tFinish period <-width pulseWidth> <-shift shift> <-factor cFactor>\n";
      Tcl_Free((char*)argv);
      return 0; 
    }   
    if (Tcl_GetDouble(interp, argv[1], &tStart) != TCL_OK) {
      opserr << "WARNING invalid tStart " << argv[1] << " - ";
      opserr << " Pulse tStart tFinish period <-width pulseWidth> <-shift shift> <-factor cFactor>\n";
      Tcl_Free((char*)argv);
      return 0;                         
    }
    if (Tcl_GetDouble(interp, argv[2], &tFinish) != TCL_OK) {
      opserr << "WARNING invalid tFinish " << argv[2] << " - ";
      opserr << " Pulse tStart tFinish period <-width pulseWidth> <-shift shift> <-factor cFactor>\n";
      Tcl_Free((char*)argv);
      return 0; 
    }     
    if (Tcl_GetDouble(interp, argv[3], &period) != TCL_OK) {
      opserr << "WARNING invalid period " << argv[3] << " - ";
      opserr << " Pulse tStart tFinish period <-width pulseWidth> <-shift shift> <-factor cFactor>\n";
      Tcl_Free((char*)argv);
      return 0; 
    }     
    
    int endMarker = 4;
    
    while (endMarker < argc && endMarker < argc) {
      if (strcmp(argv[endMarker],"-factor") == 0) {
        // allow user to specify the factor
        endMarker++;
        if (endMarker == argc || 
            Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
          
          opserr << "WARNING invalid cFactor " << argv[endMarker] << " -";
          opserr << " Pulse tStart tFinish period <-width pulseWidth> <-shift shift> <-factor cFactor>\n";
          Tcl_Free((char*)argv);
          return 0;
        }
      }

      else if (strcmp(argv[endMarker],"-width") == 0) {
        // allow user to specify pulse width
        endMarker++;
        if (endMarker == argc || 
            Tcl_GetDouble(interp, argv[endMarker], &width) != TCL_OK) {
            
          opserr << "WARNING invalid pulse width " << argv[endMarker] << " - ";
          opserr << " Pulse tStart tFinish period <-width pulseWidth> <-shift shift> <-factor cFactor>\n";
          Tcl_Free((char*)argv);
          return 0;
        }
      }

      else if (strcmp(argv[endMarker],"-shift") == 0) {
        // allow user to specify phase shift
        endMarker++;
        if (endMarker == argc || 
            Tcl_GetDouble(interp, argv[endMarker], &shift) != TCL_OK) {
            
          opserr << "WARNING invalid phase shift " << argv[endMarker] << " - ";
          opserr << " Pulse tStart tFinish period <-width pulseWidth> <-shift shift> <-factor cFactor>\n";
          Tcl_Free((char*)argv);
          return 0;
        }
      }
      endMarker++;
    }

    theSeries = new PulseSeries(tStart, tFinish, period, width, shift, cFactor);
  }     

  else if (strcmp(argv[0],"Triangle") == 0)  {
    // LoadPattern and TriangleSeries - read args & create TriangleSeries object
    double cFactor = 1.0;
    double tStart, tFinish, period;
    double shift = 0.0;
      
    if (argc < 4) {
      opserr << "WARNING not enough TriangleSeries args - ";
      opserr << " Triangle tStart tFinish period <-shift shift> <-factor cFactor>\n";
      Tcl_Free((char*)argv);
      return 0; 
    }   
    if (Tcl_GetDouble(interp, argv[1], &tStart) != TCL_OK) {
      opserr << "WARNING invalid tStart " << argv[1] << " - ";
      opserr << " Triangle tStart tFinish period <-shift shift> <-factor cFactor>\n";
      Tcl_Free((char*)argv);
      return 0;                         
    }
    if (Tcl_GetDouble(interp, argv[2], &tFinish) != TCL_OK) {
      opserr << "WARNING invalid tFinish " << argv[2] << " - ";
      opserr << " Triangle tStart tFinish period <-shift shift> <-factor cFactor>\n";
      Tcl_Free((char*)argv);
      return 0; 
    }     
    if (Tcl_GetDouble(interp, argv[3], &period) != TCL_OK) {
      opserr << "WARNING invalid period " << argv[3] << " - ";
      opserr << " Triangle tStart tFinish period <-shift shift> <-factor cFactor>\n";
      Tcl_Free((char*)argv);
      return 0; 
    }     
    
    int endMarker = 4;
    
    while (endMarker < argc && endMarker < argc) {
      if (strcmp(argv[endMarker],"-factor") == 0) {
        // allow user to specify the factor
        endMarker++;
        if (endMarker == argc || 
            Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
          
          opserr << "WARNING invalid cFactor " << argv[endMarker] << " -";
          opserr << " Triangle tStart tFinish period <-shift shift> <-factor cFactor>\n";
          Tcl_Free((char*)argv);
          return 0;
        }
      }

      else if (strcmp(argv[endMarker],"-shift") == 0) {
        // allow user to specify phase shift
        endMarker++;
        if (endMarker == argc || 
            Tcl_GetDouble(interp, argv[endMarker], &shift) != TCL_OK) {
            
          opserr << "WARNING invalid phase shift " << argv[endMarker] << " - ";
          opserr << " Triangle tStart tFinish period <-shift shift> <-factor cFactor>\n";
          Tcl_Free((char*)argv);
          return 0;
        }
      }
      endMarker++;
    }

    theSeries = new TriangleSeries(tStart, tFinish, period, shift, cFactor);    
  }

  else if (strcmp(argv[0],"Rectangular") == 0) {
    // LoadPattern and RectangularSeries - read args and create RectangularSeries object
    double tStart, tFinish;
    double cFactor = 1.0;
    if (argc < 3) {
      opserr << "WARNING not enogh args - ";
      opserr << " Rectangular tStart tFinish <-factor cFactor>\n";
      Tcl_Free((char*)argv);
      return 0; 
    }   
    if (Tcl_GetDouble(interp, argv[1], &tStart) != TCL_OK) {
      opserr << "WARNING invalid tStart " << argv[1] << " - ";
      opserr << " Rectangular tStart tFinish <-factor factor>\n";
      Tcl_Free((char*)argv);
      return 0;
    }
    if (Tcl_GetDouble(interp, argv[2], &tFinish) != TCL_OK) {
      opserr << "WARNING invalid tStart " << argv[2] << " - ";
      opserr << " Rectangular tStart tFinish <-factor fcator>\n";
      Tcl_Free((char*)argv);
      return 0; 
    }     

    int endMarker =  3;
    
    if ((endMarker != argc) && (strcmp(argv[endMarker],"-factor") == 0)) {
      // allow user to specify the factor
      endMarker++;
      if (endMarker == argc || 
          Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
        
        opserr << "WARNING invalid cFactor " << argv[endMarker] << " - ";
        opserr << " Rectangular tStart tFinish -factor cFactor\n";
        Tcl_Free((char*)argv);
        return 0;
      }
      endMarker++;
    }  

    theSeries = new RectangularSeries(tStart, tFinish, cFactor); 
  }

#if 0
  else if (strcmp(argv[0], "Rectangular") == 0) {

    void *theResult = OPS_RectangularSeries(rt, argc, argv);
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

  }
  else if ((strcmp(argv[0], "Pulse") == 0) ||
           (strcmp(argv[0], "PulseSeries") == 0)) {

    void *theResult = OPS_PulseSeries(rt, argc, argv);
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

  }
  else if ((strcmp(argv[0], "Triangle") == 0) ||
           (strcmp(argv[0], "TriangleSeries") == 0)) {

    void *theResult = OPS_TriangleSeries(rt, argc, argv);
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

  }
#endif

  else if ((strcmp(argv[0], "Series") == 0) || (strcmp(argv[0], "Path") == 0)) {

    double cFactor = 1.0;

    if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "not enough args - ";
      opserr << " Series -dt timeIncr -values {list of points }\n";
      return 0;
    }

    int tag = 0;
    double timeIncr = 0.0;
    int endMarker = 1;
    bool done = false;
    int fileName = 0;
    int fileTimeName = 0;
    int filePathName = 0;
    Vector *dataPath = nullptr;
    Vector *dataTime = nullptr;
    bool useLast = false;
    bool prependZero = false;
    double startTime = 0.0;

    struct stat fileInfo;

    if (Tcl_GetInt(interp, argv[endMarker], &tag) == TCL_OK) {
      endMarker++;
    }

    while (endMarker < argc && done == false) {

      if (strcmp(argv[endMarker], "-dt") == 0) {
        // allow user to specify the time increment
        endMarker++;
        if (endMarker == argc ||
            Tcl_GetDouble(interp, argv[endMarker], &timeIncr) != TCL_OK) {

          opserr << G3_ERROR_PROMPT << "invalid dt " << argv[endMarker] << " - ";
          opserr << " Series -dt dt ... \n";
          return 0;
        }
      }

      else if (strcmp(argv[endMarker], "-tag") == 0) {
        // allow user to specify the tag
        endMarker++;
        if (endMarker == argc ||
            Tcl_GetInt(interp, argv[endMarker], &tag) != TCL_OK) {

          opserr << G3_ERROR_PROMPT << "invalid tag " << argv[endMarker] << " - ";
          return 0;
        }
      }

      else if (strcmp(argv[endMarker], "-factor") == 0) {
        // allow user to specify the factor
        endMarker++;
        if (endMarker == argc ||
            Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {

          opserr << G3_ERROR_PROMPT << "invalid cFactor " << argv[endMarker] << " - ";
          opserr << " Series -factor ... \n";
          return 0;
        }
      }

      else if (strcmp(argv[endMarker], "-file") == 0) {
        // allow user to specify the file name containing time and data points
        endMarker++;
        if (endMarker != argc) {
          fileName = endMarker; // argv[endMarker];
          if (stat(argv[endMarker], &fileInfo ) != 0) {
            opserr << G3_ERROR_PROMPT << "Cannot open file "
                   << argv[endMarker] << "\n";
            return nullptr;
          }
        }
      }

      else if (strcmp(argv[endMarker], "-filePath") == 0) {
        // allow user to specify the file name containing the data points
        endMarker++;
        if (endMarker != argc) {
          filePathName = endMarker; // argv[endMarker];
          if (stat(argv[endMarker], &fileInfo ) != 0) {
            opserr << G3_ERROR_PROMPT << "Cannot open file "
                   << argv[endMarker] << "\n";
            return nullptr;
          }
        }
      }

      else if (strcmp(argv[endMarker], "-fileTime") == 0) {
        // allow user to specify the file name containing the data points
        endMarker++;
        if (endMarker != argc) {
          fileTimeName = endMarker; // argv[endMarker];
          if (stat(argv[endMarker], &fileInfo ) != 0) {
            opserr << G3_ERROR_PROMPT << "Cannot open file "
                   << argv[endMarker] << "\n";
            return nullptr;
          }
        }
      }

      else if (strcmp(argv[endMarker], "-values") == 0) {
        // allow user to specify the data points in tcl list
        endMarker++;
        if (endMarker != argc) {
          int pathSize;
          TCL_Char **pathStrings;

          if (Tcl_SplitList(interp, argv[endMarker], &pathSize, &pathStrings) !=
              TCL_OK) {

            opserr << G3_ERROR_PROMPT << "problem splitting path list " << argv[endMarker]
                   << " - ";
            opserr << " Series -values {path} ... \n";
            return nullptr;
          }

          dataPath = new Vector(pathSize);
          for (int i = 0; i < pathSize; ++i) {
            double value;
            if (Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
              opserr << G3_ERROR_PROMPT << "problem reading path data value "
                     << pathStrings[i] << " - ";
              opserr << " Series -values {path} ... \n";
              Tcl_Free((char *)pathStrings);
              return nullptr;
            }
            (*dataPath)(i) = value;
          }
          // free up the array of pathsStrings .. see tcl man pages as to why
          Tcl_Free((char *)pathStrings);
        }
      }

      else if (strcmp(argv[endMarker], "-time") == 0) {
        // allow user to specify the data points in tcl list
        endMarker++;
        if (endMarker != argc) {
          int pathSize;
          TCL_Char **pathStrings;

          if (Tcl_SplitList(interp, argv[endMarker], &pathSize, &pathStrings) !=
              TCL_OK) {

            opserr << G3_ERROR_PROMPT << "problem spltting time path " << argv[endMarker]
                   << " - ";
            opserr << " Series -time {times} ... \n";
            return 0;
          }

          dataTime = new Vector(pathSize);
          for (int i = 0; i < pathSize; ++i) {
            double value;
            if (Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
              opserr << G3_ERROR_PROMPT << "problem reading time path value "
                     << pathStrings[i] << " - ";
              opserr << " Series -values {path} ... \n";

              Tcl_Free((char *)pathStrings);
              return 0;
            }
            (*dataTime)(i) = value;
          }
          // free up the array of pathsStrings .. see tcl man pages as to why
          Tcl_Free((char *)pathStrings);
        }
      }

      else if (strcmp(argv[endMarker], "-useLast") == 0) {
        useLast = true;
      }

      else if (strcmp(argv[endMarker], "-prependZero") == 0) {
        prependZero = true;
      }

      else if (strcmp(argv[endMarker], "-startTime") == 0 ||
               strcmp(argv[endMarker], "-tStart") == 0) {
        // allow user to specify the start time
        endMarker++;
        if (endMarker == argc ||
            Tcl_GetDouble(interp, argv[endMarker], &startTime) != TCL_OK) {

          opserr << G3_ERROR_PROMPT << "invalid tStart " << argv[endMarker] << " - ";
          opserr << " Series -startTime tStart ... \n";
          return 0;
        }
      }

      endMarker++;
    }

    if (filePathName != 0 && fileTimeName == 0 && timeIncr != 0.0) {
      theSeries = new PathSeries(tag, argv[filePathName], timeIncr, cFactor,
                                 useLast, prependZero, startTime);
    }

    else if (fileName != 0) {
      theSeries = new PathTimeSeries(tag, argv[fileName], cFactor, useLast);

    } else if (filePathName != 0 && fileTimeName != 0) {
      theSeries = new PathTimeSeries(tag, argv[filePathName],
                                     argv[fileTimeName], cFactor, useLast);

    } else if (dataPath != 0 && dataTime == 0 && timeIncr != 0.0) {
      theSeries = new PathSeries(tag, *dataPath, timeIncr, cFactor, useLast,
                                 prependZero, startTime);
      delete dataPath;

    } else if (dataPath != 0 && dataTime != 0) {
      if (dataTime->Size() != dataPath->Size()) {
        opserr << G3_ERROR_PROMPT << "size of time vector (" << dataTime->Size()
               << ") must be equal to size of values (" << dataPath->Size() << ")\n";
        return nullptr;
      }
      theSeries =
          new PathTimeSeries(tag, *dataPath, *dataTime, cFactor, useLast);
      delete dataPath;
      delete dataTime;

    } else {
      opserr << G3_ERROR_PROMPT << "choice of options for Path Series invalid - valid "
                "options for ";
      opserr << " Path are\n";
      opserr << " \t -fileT fileTimeName -fileP filePathName \n";
      opserr << " \t -dt constTimeIncr -file filePathName\n";
      opserr << " \t -dt constTimeIncr -values {list of points on path}\n";
      opserr << " \t -time {list of time points} -values {list of points on "
                "path}\n";
      return 0;
    }

  }

#if 0
  else if ((strcmp(argv[0], "PeerDatabase") == 0) ||
           (strcmp(argv[0], "PeerMotion") == 0)) {

    void *theResult = OPS_PeerMotion(rt, argc, argv);
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

    PeerMotion *thePeerMotion = (PeerMotion *)theSeries;

    if (argc > 4 && theSeries != 0) {
      int argCount = 4;

      while (argCount + 1 < argc) {
        if ((strcmp(argv[argCount], "-dT") == 0) ||
            (strcmp(argv[argCount], "-dt") == 0) ||
            (strcmp(argv[argCount], "-DT") == 0)) {
          const char *variableName = argv[argCount + 1];
          double dT = thePeerMotion->getDt();
          char string[30];
          sprintf(string, "set %s %.18e", variableName, dT);
          if (Tcl_Eval(interp, string) != TCL_OK) {
            opserr << G3_WARN_PROMPT << Tcl_GetStringResult(interp);
            Tcl_Exit(TCL_ERROR);
          }
          argCount += 2;
        } else if ((strcmp(argv[argCount], "-nPts") == 0) ||
                   (strcmp(argv[argCount], "-NPTS") == 0)) {
          const char *variableName = argv[argCount + 1];
          int nPts = thePeerMotion->getNPts();
          char string[30];
          sprintf(string, "set %s %d", variableName, nPts);
          if (Tcl_Eval(interp, string) != TCL_OK) {
            opserr << G3_WARN_PROMPT << Tcl_GetStringResult(interp);
            Tcl_Exit(TCL_ERROR);
          }
          argCount += 2;
        } else
          argCount++;
      }
    }
  }
#endif


#if 0
  else if ((strcmp(argv[0], "PeerNGADatabase") == 0) ||
           (strcmp(argv[0], "PeerNGAMotion") == 0)) {

    void *theResult = OPS_PeerNGAMotion(rt, argc, argv);
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

    PeerNGAMotion *thePeerMotion = (PeerNGAMotion *)(theSeries);

    if (argc > 3 && theSeries != 0) {
      int argCount = 3;

      while (argCount + 1 < argc) {
        if ((strcmp(argv[argCount], "-dT") == 0) ||
            (strcmp(argv[argCount], "-dt") == 0) ||
            (strcmp(argv[argCount], "-DT") == 0)) {
          const char *variableName = argv[argCount + 1];
          double dT = thePeerMotion->getDt();
          char string[30];
          sprintf(string, "set %s %.18e", variableName, dT);
          if (Tcl_Eval(interp, string) != TCL_OK) {
            opserr << G3_WARN_PROMPT << Tcl_GetStringResult(interp);
            Tcl_Exit(TCL_ERROR);
          }
          argCount += 2;

        } else if ((strcmp(argv[argCount], "-nPts") == 0) ||
                   (strcmp(argv[argCount], "-NPTS") == 0)) {
          const char *variableName = argv[argCount + 1];
          int nPts = thePeerMotion->getNPts();
          char string[30];
          sprintf(string, "set %s %d", variableName, nPts);
          if (Tcl_Eval(interp, string) != TCL_OK) {
            opserr << G3_WARN_PROMPT << Tcl_GetStringResult(interp);
            Tcl_Exit(TCL_ERROR);
          }
          argCount += 2;
        } else
          argCount++;
      }
    }
  }
#endif

#ifdef _RELIABILITY

  else if (strcmp(argv[0], "DiscretizedRandomProcess") == 0) {

    double mean, maxStdv;
    ModulatingFunction *theModFunc;

    if (Tcl_GetDouble(interp, argv[1], &mean) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid input: random process mean \n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[2], &maxStdv) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid input: random process max stdv \n";
      return 0;
    }

    // Number of modulating functions
    int argsBeforeModList = 3;
    int numModFuncs = argc - argsBeforeModList;

    // Create an array to hold pointers to modulating functions
    ModulatingFunction **theModFUNCS = new ModulatingFunction *[numModFuncs];

    // For each modulating function, get the tag and ensure it exists
    int tagI;
    for (int i = 0; i < numModFuncs; ++i) {
      if (Tcl_GetInt(interp, argv[i + argsBeforeModList], &tagI) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid modulating function tag. " << "\n";
        return 0;
      }

      theModFunc = 0;
      theModFunc = theReliabilityDomain->getModulatingFunction(tagI);

      if (theModFunc == 0) {
        opserr << G3_ERROR_PROMPT << "modulating function number "
               << argv[i + argsBeforeModList] << "does not exist...\n";
        delete[] theModFUNCS;
        return 0;
      } else {
        theModFUNCS[i] = theModFunc;
      }
    }

    // Parsing was successful, create the random process series object
    theSeries = new DiscretizedRandomProcessSeries(0, numModFuncs, theModFUNCS,
                                                   mean, maxStdv);
  }

  else if (strcmp(argv[0], "SimulatedRandomProcess") == 0) {

    int spectrumTag, numFreqIntervals;
    double mean;

    if (Tcl_GetInt(interp, argv[1], &spectrumTag) != TCL_OK) {
      opserr << "WARNING invalid input to SimulatedRandomProcess: spectrumTag"
             << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[2], &mean) != TCL_OK) {
      opserr << "WARNING invalid input to SimulatedRandomProcess: mean"
             << "\n";
      return 0;
    }

    if (Tcl_GetInt(interp, argv[3], &numFreqIntervals) != TCL_OK) {
      opserr
          << "WARNING invalid input to SimulatedRandomProcess: numFreqIntervals"
          << "\n";
      return 0;
    }

    // Check that the random number generator exists
    if (theRandomNumberGenerator == 0) {
      opserr << "WARNING: A random number generator must be instantiated "
                "before SimulatedRandomProcess."
             << "\n";
      return 0;
    }

    // Check that the spectrum exists
    Spectrum *theSpectrum = 0;
    theSpectrum = theReliabilityDomain->getSpectrum(spectrumTag);
    if (theSpectrum == 0) {
      opserr << "WARNING: Could not find the spectrum for the "
                "SimulatedRandomProcess."
             << "\n";
      return 0;
    }

    // Parsing was successful, create the random process series object
    theSeries = new SimulatedRandomProcessSeries(
        0, theRandomNumberGenerator, theSpectrum, numFreqIntervals, mean);
  }

#endif

  else {
    for (int i = 0; i < argc; ++i)
      opserr << argv[i] << ' ';
    opserr << "\n";
    // type of load pattern type unknown
    opserr << "WARNING unknown Series type " << argv[0] << " - ";
    opserr << " valid types: Linear, Rectangular, Path, Constant, Trig, Sine\n";
    return 0;
  }

  return theSeries;
}

TimeSeries *
TclSeriesCommand(ClientData clientData, Tcl_Interp *interp, TCL_Char * const arg)
{
  int argc;
  TCL_Char ** argv;
  TimeSeries *series;
  int timeSeriesTag = 0;

  if (Tcl_GetInt(interp, arg, &timeSeriesTag) == TCL_OK) {
    if (clientData && (series = ((BasicModelBuilder*)clientData)->getTypedObject<TimeSeries>(timeSeriesTag)))
      return series->getCopy();
  }

  // split the list
  if (Tcl_SplitList(interp, arg, &argc, &argv) != TCL_OK) {
    opserr << "WARNING could not split series list " << arg << "\n";
    return nullptr;
  }

  TimeSeries *theSeries = TclDispatch_newTimeSeries(clientData, interp, argc, argv);

  // clean up after ourselves and return the series
  Tcl_Free((char *)argv);
  return theSeries;
}

int
TclCommand_addTimeSeries(ClientData clientData, Tcl_Interp *interp, int argc,
                         TCL_Char ** const argv)
{
  TimeSeries *theSeries = TclDispatch_newTimeSeries(clientData, interp, argc - 1, &argv[1]);

  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  if (theSeries != nullptr) {
    int tag;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "failed to read tag \"" << argv[2] << "\"\n";
      return TCL_ERROR;
    }
    if (builder->addTypedObject<TimeSeries>(tag, theSeries) < 0)
      return TCL_ERROR;
    else
      return TCL_OK;
  }
  return TCL_ERROR;
}


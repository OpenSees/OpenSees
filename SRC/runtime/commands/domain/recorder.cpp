//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the function that is invoked
// by the interpreter when the comand 'record' is invoked by the
// user.
//
// Written: fmk, cmp
// Created: 04/98
//
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#define strcmp strcasecmp

#include <tcl.h>
#include <G3_Logging.h>
#include <Domain.h>
#include <NodeIter.h>
#include <NodeData.h>
#include <ElementIter.h>
#include <Node.h>
#include <Element.h>
#include <Parameter.h>
#include <DamageModel.h>
#include <packages.h>

// Streams
#include <StandardStream.h>
#include <DataFileStream.h>
#include <DataFileStreamAdd.h>
#include <XmlFileStream.h>
#include <BinaryFileStream.h>
#include <DatabaseStream.h>
#include <DummyStream.h>
#include <TCP_Stream.h>

// Recorders
#include <NodeRecorder.h>
#include <NodeRecorderRMS.h>
#include <EnvelopeNodeRecorder.h>
#include <PatternRecorder.h>
#include <DriftRecorder.h>
#include <EnvelopeDriftRecorder.h>
#include <ElementRecorder.h>
#include <EnvelopeElementRecorder.h>
#include <NormElementRecorder.h>
#include <NormEnvelopeElementRecorder.h>
#include <DamageRecorder.h>
#include <MeshRegion.h>
#include <RemoveRecorder.h>

#define MAX_NDF 6
extern FE_Datastore    *theDatabase;
extern FEM_ObjectBroker theBroker;

OPS_Routine OPS_PVDRecorder;
OPS_Routine OPS_GmshRecorder;
OPS_Routine OPS_MPCORecorder;
OPS_Routine OPS_VTK_Recorder;
OPS_Routine OPS_ElementRecorderRMS;

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp,
                        int cArg, int mArg, TCL_Char ** const argv,
                        Domain *domain);

extern TimeSeries *TclSeriesCommand(ClientData clientData, Tcl_Interp *interp,
                                    TCL_Char *arg);

//
// Dynamically loaded recorders
typedef struct externalRecorderCommand {
  char *funcName;
  void *(*funcPtr)();
  struct externalRecorderCommand *next;
} ExternalRecorderCommand;

static ExternalRecorderCommand *theExternalRecorderCommands = NULL;

//
// Small structure to hold common config options
struct OutputOptions {
  TCL_Char   *filename  = nullptr;
  TCL_Char   *tableName = nullptr;
  const char *inetAddr  = nullptr;
  int         inetPort  = 0;
  int  precision        = 6;
  int writeBufferSize   = 0;
  bool doScientific     = false;
  bool closeOnWrite     = false;

  FE_Datastore *theDatabase = nullptr;

  enum Mode {
    STANDARD_STREAM,
    DATA_STREAM,
    XML_STREAM,
    DATABASE_STREAM,
    BINARY_STREAM,
    DATA_STREAM_CSV,
    TCP_STREAM,
    DATA_STREAM_ADD,
    MODE_UNSPECIFIED
  } eMode = STANDARD_STREAM;
};


static int
createNodeRecorder(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char ** const argv, Recorder **theRecorder);

static OPS_Stream *
createOutputStream(OutputOptions &options)
{
  OPS_Stream *theOutputStream = nullptr;

  // construct the DataHandler
  if (options.filename != nullptr) {
    if (options.eMode == OutputOptions::DATA_STREAM) {
      theOutputStream = new DataFileStream(
          options.filename, 
          openMode::OVERWRITE, 2, 0, 
          options.closeOnWrite, 
          options.precision, 
          options.doScientific);

    } else if (options.eMode == OutputOptions::DATA_STREAM_ADD) {
      theOutputStream = new DataFileStreamAdd(
          options.filename, 
          openMode::OVERWRITE, 2, 0, 
          options.closeOnWrite, 
          options.precision, 
          options.doScientific);

    } else if (options.eMode == OutputOptions::DATA_STREAM_CSV) {
      theOutputStream = new DataFileStream(
          options.filename, 
          openMode::OVERWRITE, 2, 1, 
          options.closeOnWrite, 
          options.precision, 
          options.doScientific);

    } else if (options.eMode == OutputOptions::XML_STREAM) {
      theOutputStream = new XmlFileStream(options.filename);

    } else if (options.eMode == OutputOptions::BINARY_STREAM) {
      theOutputStream = new BinaryFileStream(options.filename);
    }

  } else if (options.eMode == OutputOptions::TCP_STREAM && options.inetAddr != 0) {
    theOutputStream = new TCP_Stream(options.inetPort, options.inetAddr);

  } else if (options.eMode == OutputOptions::DATABASE_STREAM && options.tableName != 0) {
    theOutputStream = new DatabaseStream(options.theDatabase, options.tableName);

  } else {
    theOutputStream = new StandardStream();
  }

  theOutputStream->setPrecision(options.precision);

  return theOutputStream;
}

static NodeData
getNodeDataFlag(const char *dataToStore, Domain& theDomain, int* dataIndex)
{
  NodeData dataFlag = NodeData::DisplTrial; // 10;

  if (dataToStore == nullptr)
    dataFlag = NodeData::DisplTrial; // 0;
  else if ((strcmp(dataToStore, "disp") == 0)) {
    dataFlag = NodeData::DisplTrial; // 0
  } else if ((strcmp(dataToStore, "vel") == 0)) {
    dataFlag = NodeData::VelocTrial; // 1
  } else if ((strcmp(dataToStore, "accel") == 0)) {
    dataFlag = NodeData::AccelTrial; // 2
  } else if ((strcmp(dataToStore, "incrDisp") == 0)) {
    dataFlag = NodeData::IncrDisp;  // 3;
  } else if ((strcmp(dataToStore, "incrDeltaDisp") == 0)) {
    dataFlag = NodeData::IncrDeltaDisp; // 4;
  } else if ((strcmp(dataToStore, "unbalance") == 0)) {
    dataFlag = NodeData::UnbalancedLoad; // 5
  } else if ((strcmp(dataToStore, "unbalanceInclInertia") == 0) ||
	     (strcmp(dataToStore, "unbalanceIncInertia") == 0) ||
	     (strcmp(dataToStore, "unbalanceIncludingInertia") == 0))  {
    dataFlag = NodeData::UnbalanceInclInertia; // 6
  } else if ((strcmp(dataToStore, "reaction") == 0)) {
    dataFlag = NodeData::Reaction;   // 7
  } else if (((strcmp(dataToStore, "reactionIncInertia") == 0))
	     || ((strcmp(dataToStore, "reactionInclInertia") == 0))
	     || ((strcmp(dataToStore, "reactionIncludingInertia") == 0))) {
    dataFlag = NodeData::ReactionInclInertia; // 8;
  } else if (((strcmp(dataToStore, "rayleighForces") == 0))
	     || ((strcmp(dataToStore, "rayleighDampingForces") == 0))) {
    dataFlag = NodeData::ReactionInclRayleigh; // 9;

  } else if ((strcmp(dataToStore, "dispNorm") == 0)) {
    dataFlag = NodeData::DisplNorm;

  } else if (((strcmp(dataToStore, "nodalRayleighForces") == 0))) {
    dataFlag = NodeData::RayleighForces; // 10001;

  } else if (((strcmp(dataToStore, "pressure") == 0))) {
    dataFlag = NodeData::Pressure;

  } else if ((strncmp(dataToStore, "eigen",5) == 0)) {
    int mode = atoi(&(dataToStore[5]));
    *dataIndex = mode;
    if (mode > 0)
      dataFlag = NodeData::EigenVector; // 10 + mode;
    else
      dataFlag = NodeData::Empty; // 10;

  } else if ((strncmp(dataToStore, "sensitivity",11) == 0)) {
    int paramTag = atoi(&(dataToStore[11]));
    Parameter *theParameter = theDomain.getParameter(paramTag);
    int grad = -1;
    if (theParameter != nullptr)
      grad = theParameter->getGradIndex();
    *dataIndex = grad;
    if (grad > 0)
      dataFlag = NodeData::DisplSensitivity; // 1000 + grad;
    else
      dataFlag = NodeData::Empty; // 10;

  } else if ((strncmp(dataToStore, "velSensitivity",14) == 0)) {
    int paramTag = atoi(&(dataToStore[14]));
    Parameter *theParameter = theDomain.getParameter(paramTag);
    int grad = -1;
    if (theParameter != nullptr)
      grad = theParameter->getGradIndex();
    
    *dataIndex = grad;
    if (grad > 0)
      dataFlag = NodeData::VelocSensitivity; // 2000 + grad;
    else
      dataFlag = NodeData::Empty; // 10;

  } else if ((strncmp(dataToStore, "accSensitivity",14) == 0)) {
    int paramTag = atoi(&(dataToStore[14]));
    Parameter *theParameter = theDomain.getParameter(paramTag);
    int grad = -1;
    if (theParameter != nullptr)
      grad = theParameter->getGradIndex();

    *dataIndex = grad;
    if (grad > 0)
      dataFlag = NodeData::AccelSensitivity; // 3000 + grad;
    else
      dataFlag = NodeData::Empty; // 10;

  } else {
    // TODO
    dataFlag = NodeData::Empty; // 10;
    opserr << "NodeRecorder::NodeRecorder - dataToStore '" << dataToStore;
    opserr << "' not recognized (disp, vel, accel, incrDisp, incrDeltaDisp)\n";
  }
  return dataFlag;
}


static int
parseOutputOption(OutputOptions *options, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
    OutputOptions::Mode eMode = OutputOptions::MODE_UNSPECIFIED;

    int loc = 0;
    if (strcmp(argv[loc], "-precision") == 0) {
      if (++loc >= argc || Tcl_GetInt(interp, argv[loc], &options->precision) != TCL_OK)
        return -1;
      loc++;
    }
 
    else if (strcmp(argv[loc], "-scientific") == 0) {
      options->doScientific = true;
      loc++;
    }

    else if (strcmp(argv[loc], "-closeOnWrite") == 0) {
      options->closeOnWrite = true;
      loc++;
    }

    else if (strcmp(argv[loc], "-buffer") == 0 ||
             strcmp(argv[loc], "-bufferSize") == 0) {
      loc++;
      if (loc >= argc || Tcl_GetInt(interp, argv[loc], &options->writeBufferSize) != TCL_OK)
        return -1;
      loc++;
    }
 
    else {
      // pick out filename
      if ((strcmp(argv[loc], "-file") == 0) ||
          (strcmp(argv[loc], "-txt")  == 0))
        eMode = OutputOptions::DATA_STREAM;

      else if (strcmp(argv[loc], "-fileAdd") == 0)
        eMode = OutputOptions::DATA_STREAM_ADD;

      else if ((strcmp(argv[loc], "-fileCSV") == 0) ||
               (strcmp(argv[loc], "-csv") == 0))
        eMode = OutputOptions::DATA_STREAM_CSV;

      else if ((strcmp(argv[loc], "-nees") == 0) ||
               (strcmp(argv[loc], "-xml")  == 0)) {
        eMode = OutputOptions::XML_STREAM;
      } 
      else if ((strcmp(argv[loc], "-binary") == 0)) {
        eMode = OutputOptions::BINARY_STREAM;
      }
      else if ((strcmp(argv[loc], "-TCP") == 0) ||
               (strcmp(argv[loc], "-tcp") == 0)) {
        options->inetAddr = argv[loc + 1];
        if (Tcl_GetInt(interp, argv[loc + 2], &options->inetPort) != TCL_OK) {
          ;
        }
        eMode = OutputOptions::TCP_STREAM;
        loc += 3;
      } 
      else if (strcmp(argv[loc], "-database") == 0) {
        // TODO:
        // theRecorderDatabase = theDatabase;
        // if (theRecorderDatabase != 0) {
        //   tableName = argv[loc + 1];
        //   eMode = OutputOptions::DATABASE_STREAM;
        // } else {
        //   opserr << "WARNING recorder Node .. -database " << filename
        //          << "  - NO CURRENT DATABASE, results to File instead\n";
        //   filename = argv[loc + 1];
        // }

        // loc += 2;
      }

      // If one of these matched, move on to filename
      if (eMode != OutputOptions::MODE_UNSPECIFIED) {
        if (loc + 1 < argc) {
          options->filename = argv[loc + 1];
          options->eMode = eMode;
        } else {
          opserr << G3_ERROR_PROMPT
                 << "expected file name after flag '" << argv[loc] << "\n";
          return -1;
        }
        loc += 2;
      }
    }
    return loc;
}

static int
TclCreateRecorder(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char ** const argv, Domain &theDomain, Recorder **theRecorder)
{
  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData;
  G3_Runtime *rt = G3_getRuntime(interp);
  (*theRecorder) = nullptr;

  // make sure at least one other argument to contain integrator
  if (argc < 2) {
    opserr << "WARNING need to specify a Recorder type\n";
    return TCL_ERROR;
  }

  //
  // check argv[1] for type of Recorder, parse in rest of arguments
  // needed for the type of Recorder, create the object and add to Domain
  //
  OPS_Stream   *theOutputStream     = nullptr;

  // an Element Recorder or ElementEnvelope Recorder
  if ((strcmp(argv[1], "Element") == 0) ||
      (strcmp(argv[1], "EnvelopeElement") == 0) ||
      (strcmp(argv[1], "NormElement") == 0) ||
      (strcmp(argv[1], "NormEnvelopeElement") == 0)) {

    OutputOptions options;

    int numEle = 0;
    int endEleIDs = 2;
    double dT   = 0.0;
    double rTolDt = 1e-5;
    bool echoTime = false;
    int loc = endEleIDs;
    int flags   = 0;
    int eleData = 0;
    ID *eleIDs  = 0;

    ID *specificIndices = nullptr;

    while (flags == 0 && loc < argc) {
      int consumed;
      if ((consumed = parseOutputOption(&options, interp, argc-loc, &argv[loc])) != 0) {
        if (consumed > 0)
          loc += consumed;
        else
          return TCL_ERROR;
      }

      else if (strcmp(argv[loc], "-rTolDt") == 0) {
        loc++;
        if (Tcl_GetDouble(interp, argv[loc], &rTolDt) != TCL_OK)
          return TCL_ERROR;
        loc++;
      }

      else if ((strcmp(argv[loc], "-ele") == 0) ||
               (strcmp(argv[loc], "-tag") == 0) ||
               (strcmp(argv[loc], "-eles") == 0) ||
               (strcmp(argv[loc], "-elem") == 0) ||
               (strcmp(argv[loc], "-element") == 0)) {

        // ensure no segmentation fault if user messes up
        if (argc < loc + 2) {
          opserr << "WARNING recorder Element .. -ele tag1? .. - no ele tags "
                    "specified\n";
          return TCL_ERROR;
        }

        // read in a list of ele until end of command or other flag
        loc++;
        int eleTag;
        eleIDs = new ID(0, 32);
        while (loc < argc && Tcl_GetInt(interp, argv[loc], &eleTag) == TCL_OK) {
          (*eleIDs)[numEle] = eleTag;
          numEle++;
          loc++;
        }
        Tcl_ResetResult(interp);

        if (loc == argc) {
          opserr << "ERROR: No response type specified for element recorder. "
                 << endln;
          delete eleIDs;
          return TCL_ERROR;
        }

        if (strcmp(argv[loc], "all") == 0) {
          eleIDs = nullptr;
          loc++;
        }

      } else if ((strcmp(argv[loc], "-eleRange") == 0) ||
                 (strcmp(argv[loc], "-range") == 0)) {

        // ensure no segmentation fault if user messes up
        if (argc < loc + 3) {
          opserr << "WARNING recorder Element .. -eleRange start? end?  .. - "
                    "no ele tags specified\n";
          return TCL_ERROR;
        }

        // read in start and end tags of two elements & add set [start,end]
        int start, end;
        if (Tcl_GetInt(interp, argv[loc + 1], &start) != TCL_OK) {
          opserr << "WARNING recorder Element -eleRange start? end? - invalid "
                    "start "
                 << argv[loc + 1] << endln;
          return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[loc + 2], &end) != TCL_OK) {
          opserr
              << "WARNING recorder Element -eleRange start? end? - invalid end "
              << argv[loc + 2] << endln;
          return TCL_ERROR;
        }

        if (start > end) {
          int swap = end;
          end = start;
          start = swap;
        }

        eleIDs = new ID(end - start);
        for (int i = start; i <= end; ++i)
          (*eleIDs)[numEle++] = i;

        loc += 3;
      }

      else if (strcmp(argv[loc], "-region") == 0) {
        // allow user to specif elements via a region
        if (argc < loc + 2) {
          opserr << "WARNING recorder Element .. -region tag?  .. - no region "
                    "specified\n";
          return TCL_ERROR;
        }
        int tag;
        if (Tcl_GetInt(interp, argv[loc + 1], &tag) != TCL_OK) {
          opserr << "WARNING recorder Element -region tag? - invalid tag "
                 << argv[loc + 1] << endln;
          return TCL_ERROR;
        }
        MeshRegion *theRegion = domain->getRegion(tag);
        if (theRegion == nullptr) {
          opserr << "WARNING recorder Element -region " << tag
                 << " - region does not exist" << endln;
          return TCL_ERROR; // was TCL_OK
        }
        const ID &eleRegion = theRegion->getElements();

        eleIDs = new ID(eleRegion.Size());
        for (int i = 0; i < eleRegion.Size(); ++i)
          (*eleIDs)[numEle++] = eleRegion(i);

        loc += 2;
      }

      else if ((strcmp(argv[loc], "-dof") == 0) ||
               (strcmp(argv[loc], "-dofs") == 0)) {

        // ensure no segmentation fault if user messes up
        if (argc < loc + 2) {
          opserr << "WARNING recorder Element .. -ele tag1? .. - no ele tags "
                    "specified\n";
          return TCL_ERROR;
        }

        // read in a list of indices until end of command or other flag
        loc++;
        int index;
        int numIndex = 0;
        specificIndices = new ID(0, MAX_NDF);
        while (loc < argc && Tcl_GetInt(interp, argv[loc], &index) == TCL_OK) {
          (*specificIndices)[numIndex] = index - 1; // opensees to c indexing
          numIndex++;
          loc++;
        }
        Tcl_ResetResult(interp);

      } else if ((strcmp(argv[loc], "-time") == 0) ||
                 (strcmp(argv[loc], "-load") == 0)) {
        echoTime = true;
        loc++;
      }

      else if (strcmp(argv[loc], "-dT") == 0) {
        // allow user to specify time step size for recording
        loc++;
        if (Tcl_GetDouble(interp, argv[loc], &dT) != TCL_OK)
          return TCL_ERROR;
        loc++;
      }

      else {
        // TODO: handle the same as Node recorder; see Example1.1.py
        // first unknown string then is assumed to start
        // element response request starts
        eleData = loc;
        flags = 1;
      }
    }

    if (eleData >= argc) {
      opserr << "ERROR: No response type specified for element recorder. "
             << endln;
      return TCL_ERROR;
    }

    const char **data = new const char *[argc - eleData];
    for (int i = eleData, j = 0; i < argc; i++, j++)
      data[j] = argv[i];

    // construct the DataHandler
    theOutputStream = createOutputStream(options);

    if (strcmp(argv[1], "Element") == 0)
      (*theRecorder) = new ElementRecorder(eleIDs, data, argc - eleData, echoTime, *domain,
                                           *theOutputStream, dT, rTolDt, specificIndices);

    else if (strcmp(argv[1], "EnvelopeElement") == 0)
      (*theRecorder) = new EnvelopeElementRecorder(eleIDs, data, argc - eleData,
                                                   *domain, *theOutputStream, dT, rTolDt,
                                                   echoTime, specificIndices);

    else if (strcmp(argv[1], "NormElement") == 0)
      (*theRecorder) = new NormElementRecorder(eleIDs, data, argc - eleData, 
                                               echoTime, *domain, *theOutputStream, 
                                               dT, rTolDt, specificIndices);

    else
      (*theRecorder) = new NormEnvelopeElementRecorder(eleIDs, data, argc - eleData,
                                                       *domain, *theOutputStream, dT, 1e-6,
                                                       echoTime, specificIndices);

    if (eleIDs != nullptr)
      delete eleIDs;

    delete[] data;
  }

  else if ((strcmp(argv[1], "Damage") == 0) ||
           (strcmp(argv[1], "ElementDamage") == 0) ||
           (strcmp(argv[1], "damage") == 0) ||
           (strcmp(argv[1], "elementDamage") == 0)) {
    //////////  By Arash Altoontash /////////////////
    TCL_Char *filename  = nullptr;

    if (argc < 7) {
      opserr << "WARNING recorder ElementDamage eleID? <-time> "
             << "<-file filename?> <-section secID1? secID2? ...> <-dof "
                "dofID?> <-damage dmgID?>";
      return TCL_ERROR;
    }

    double dT = 0.0;
    bool echoTime = false;
    int loc = 2;
    int eleID;

    if (Tcl_GetInt(interp, argv[loc], &eleID) != TCL_OK) {
      opserr << "WARNING recorder ElementDamage: No element tag specified ";
      return TCL_ERROR;
    }
    loc++;

    if ((strcmp(argv[loc], "-time") == 0) ||
        (strcmp(argv[loc], "-load") == 0)) {
      // allow user to specify const load
      echoTime = true;
      loc++;
    } 
    else if (strcmp(argv[loc], "-dT") == 0) {
      // allow user to specify time step size for recording
      loc++;
      if (Tcl_GetDouble(interp, argv[loc], &dT) != TCL_OK)
        return TCL_ERROR;
      loc++;
    }

    if (strcmp(argv[loc], "-file") == 0) {
      // allow user to specify load pattern other than current
      loc++;
      filename = argv[loc];
      loc++;
    }

    if (strcmp(argv[loc], "-section") != 0 &&
        strcmp(argv[loc], "section") != 0) {
      opserr
          << "WARNING recorder ElementDamage: Section keyword not specified ";
      return TCL_ERROR;
    }
    loc++;

    int secID;
    int endSecIDs = loc;
    int numSec = 0;
    while (Tcl_GetInt(interp, argv[endSecIDs], &secID) == TCL_OK) {
      endSecIDs++;
    }

    numSec = endSecIDs - loc;
    // create an ID to hold section/material tags
    ID secIDs(numSec);

    // read in the sec tags to the ID
    for (int i = loc; i < endSecIDs; ++i) {
      if (Tcl_GetInt(interp, argv[i], &secID) != TCL_OK)
        return TCL_ERROR;
      secIDs[loc - i] = secID;
    }

    loc = endSecIDs;

    int dofID = 0;
    if (strcmp(argv[loc], "-dof") == 0 || strcmp(argv[loc], "dof") == 0 ||
        strcmp(argv[loc], "-DOF") == 0 || strcmp(argv[loc], "DOF") == 0) {
      loc++;
      if (Tcl_GetInt(interp, argv[loc], &dofID) != TCL_OK) {
        opserr << "WARNING recorder ElementDamage: No dof tag specified ";
        return TCL_ERROR;
      }
      loc++;
    }

    if (strcmp(argv[loc], "-damage") != 0 && strcmp(argv[loc], "damage") != 0) {
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
    dmgPTR = OPS_getDamageModel(dmgID);

    if (dmgPTR == NULL) {
      opserr << "WARNING recorder ElementDamage: specified damage model not "
                "found\n";
      return TCL_ERROR;
    }


    OPS_Stream *theOutput = new DataFileStream(filename);

    // now construct the recorder
    (*theRecorder) = new DamageRecorder(eleID, secIDs, dofID, dmgPTR, *domain,
                                        echoTime, dT, 1e-6, *theOutput);

  }

  else if (/* (strcmp(argv[1], "Remove") == 0) || */
           (strcmp(argv[1], "ElementRemoval") == 0) ||
           (strcmp(argv[1], "NodeRemoval") == 0) ||
           (strcmp(argv[1], "Collapse") == 0)) {

    if (argc < 4) {
      opserr << "WARNING recorder Collapse -ele eleID <eleID2? ...>  -node "
                "nodeID <-time> <-file filename?> ? "
             << "\n or recorder Collapse -ele eleID1 <eleID2? ...>? <-sec "
                "secID1? secID2? ...> -crit crit1? value1?"
             << " <-crit crit2? value2?> <-time> <-file filename?> <-mass "
                "mass1? mass2? ...> <-g gAcc gDir? gPat?>?"
             << endln;
      return TCL_ERROR;
    }

    // current maximum number of either-or criteria for an element
    int maxNumCriteria = 2; 

    double dT = 0.0;
    bool echoTime = false;
    int endEleIDs = 2;
    int numEle = endEleIDs - 2;
    int loc = endEleIDs;
    int flags    = 0;
    int nodeTag  = 0;
    int nTagbotn = 0;
    int nTagmidn = 0;
    int nTagtopn = 0;
    int globgrav = 0;
    const char *filenameinf = 0;

    TCL_Char *filename  = nullptr;

    //  end of new
    int numSecondaryEle = 0;

    // create an ID to hold ele tags
    // ID eleIDs(numEle, numEle+1);
    ID eleIDs;
    eleIDs = ID(1);
    ID secondaryEleIDs = ID(1);
    secondaryEleIDs[0] = 0;
    bool secondaryFlag = false;
    ID secIDs = 0;

    // optional mass and weight definition
    Vector eleMass(1);
    eleMass.Zero();
    double gAcc = 0;
    int gDir = 0, gPat = 0;

    Vector remCriteria(2 * maxNumCriteria);
    remCriteria.Zero();
    int numCrit = 0;

    while (flags == 0 && loc < argc) {

      if ((strcmp(argv[loc], "-node") == 0) ||
          (strcmp(argv[loc], "-tag") == 0)) {

        if (Tcl_GetInt(interp, argv[loc + 1], &nodeTag) != TCL_OK) {
          opserr << "WARNING recorder Collapse -node - invalid node tag "
                 << argv[loc + 1] << endln;
          return TCL_ERROR;
        }

        Node *theNode = domain->getNode(nodeTag);
        if (theNode == nullptr) {
          opserr << "WARNING recorder Collapse -node - invalid node "
                 << argv[loc + 1] << endln;
          return TCL_ERROR;
        }
        loc += 2;
      }
      // new

      else if (strcmp(argv[loc], "-file_infill") == 0) {
        filenameinf = argv[loc + 1];
        loc += 2;
      }

      else if (strcmp(argv[loc], "-checknodes") == 0) {
        if (Tcl_GetInt(interp, argv[loc + 1], &nTagbotn) != TCL_OK) {
          opserr << "WARNING recorder Collapse -node - invalid node tag "
                 << argv[loc + 1] << endln;
          return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[loc + 2], &nTagmidn) != TCL_OK) {
          opserr << "WARNING recorder Collapse -node - invalid node tag "
                 << argv[loc + 1] << endln;
          return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[loc + 3], &nTagtopn) != TCL_OK) {
          opserr << "WARNING recorder Collapse -node - invalid node tag "
                 << argv[loc + 1] << endln;
          return TCL_ERROR;
        }
        loc += 4;
      }

      else if (strcmp(argv[loc], "-global_gravaxis") == 0) {
        if (Tcl_GetInt(interp, argv[loc + 1], &globgrav) != TCL_OK) {
          opserr << "WARNING recorder Collapse -global_gravaxis - invalid "
                    "global axis for gravity "
                 << argv[loc + 1] << endln;
          return TCL_ERROR;
        }
        loc += 2;
      }

      //    end of new

      else if ((strcmp(argv[loc], "-slave") == 0) ||
               (strcmp(argv[loc], "-secondary") == 0)) {
        secondaryFlag = true;
        loc++;
      }

      else if ((strcmp(argv[loc], "-ele") == 0) ||
               (strcmp(argv[loc], "-eles") == 0) ||
               (strcmp(argv[loc], "-element") == 0)) {

        // ensure no segmentation fault if user messes up
        if (argc < loc + 2) {
          opserr << "WARNING recorder Collapse .. -ele tag1? .. - no ele tags "
                    "specified\n";
          return TCL_ERROR;
        }

        //
        // read in a list of ele until end of command or other flag
        //
        loc++;
        int eleTag;
        while (loc < argc && Tcl_GetInt(interp, argv[loc], &eleTag) == TCL_OK) {
          if (secondaryFlag == false)
            eleIDs[numEle++] = eleTag;
          else
            secondaryEleIDs[numSecondaryEle++] = eleTag;
          secondaryFlag = false;
          loc++;
        }

        Tcl_ResetResult(interp);

        if (strcmp(argv[loc], "all") == 0) {
          ElementIter &theEleIter = domain->getElements();
          Element *theEle;
          while ((theEle = theEleIter()) != 0)
            eleIDs[numEle++] = theEle->getTag();
          loc++;
        }

      }

      else if (strcmp(argv[loc], "-eleRange") == 0) {

        // ensure no segmentation fault if user messes up
        if (argc < loc + 3) {
          opserr << "WARNING recorder Element .. -eleRange start? end?  .. - "
                    "no ele tags specified\n";
          return TCL_ERROR;
        }

        //
        // read in start and end tags of two elements & add set [start,end]
        //
        int start, end;
        if (Tcl_GetInt(interp, argv[loc + 1], &start) != TCL_OK) {
          opserr << "WARNING recorder Element -eleRange start? end? - invalid "
                    "start "
                 << argv[loc + 1] << endln;
          return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[loc + 2], &end) != TCL_OK) {
          opserr
              << "WARNING recorder Element -eleRange start? end? - invalid end "
              << argv[loc + 2] << endln;
          return TCL_ERROR;
        }
        if (start > end) {
          int swap = end;
          end = start;
          start = swap;
        }

        for (int i = start; i <= end; ++i)
          if (secondaryFlag == false)
            eleIDs[numEle++] = i;
          else
            secondaryEleIDs[numSecondaryEle++] = i;

        secondaryFlag = false;
        loc += 3;
      }

      else if (strcmp(argv[loc], "-region") == 0) {
        // allow user to specif elements via a region
        if (argc < loc + 2) {
          opserr << "WARNING recorder Element .. -region tag?  .. - no region "
                    "specified\n";
          return TCL_ERROR;
        }
        int tag;
        if (Tcl_GetInt(interp, argv[loc + 1], &tag) != TCL_OK) {
          opserr << "WARNING recorder Element -region tag? - invalid tag "
                 << argv[loc + 1] << endln;
          return TCL_ERROR;
        }
        MeshRegion *theRegion = domain->getRegion(tag);
        if (theRegion == 0) {
          opserr << "WARNING recorder Element -region " << tag
                 << " - region does not exist" << endln;
          return TCL_OK;
        }
        const ID &eleRegion = theRegion->getElements();
        for (int i = 0; i < eleRegion.Size(); ++i)
          if (secondaryFlag == false)
            eleIDs[numEle++] = eleRegion(i);
          else
            secondaryEleIDs[numSecondaryEle++] = eleRegion(i);

        secondaryFlag = false;
        loc += 2;
      }

      else if ((strcmp(argv[loc], "-time") == 0) ||
               (strcmp(argv[loc], "-load") == 0)) {
        // allow user to specify const load
        echoTime = true;
        loc++;
      }

      else if (strcmp(argv[loc], "-dT") == 0) {
        // allow user to specify time step size for recording
        loc++;
        if (Tcl_GetDouble(interp, argv[loc], &dT) != TCL_OK)
          return TCL_ERROR;
        loc++;
      }

      else if (strcmp(argv[loc], "-file") == 0) {
        filename = argv[loc + 1];
        loc += 2;
      }

      else if ((strcmp(argv[loc], "-mass") == 0)) {
        eleMass.resize(numEle);
        eleMass.Zero();
        loc++;
        double eleM = 0;

        if (loc < argc && Tcl_GetDouble(interp, argv[loc + 1], &eleM) != TCL_OK) {
          Tcl_GetDouble(interp, argv[loc], &eleM);
          for (int i = 0; i < numEle; ++i)
            eleMass(i) = eleM;
          loc++;
        } else {
          int i = 0;
          while (loc < argc && Tcl_GetDouble(interp, argv[loc], &eleM) == TCL_OK) {
            eleMass(i) = eleM;
            loc++;
            i++;
          }
        }

        // Tcl_ResetResult(interp);
      }

      else if ((strcmp(argv[loc], "-g") == 0)) {
        loc++;

        if (Tcl_GetDouble(interp, argv[loc], &gAcc) != TCL_OK) {
          opserr << "WARNING recorder Remove -g gValue? gDir? gPat?... invalid "
                    "gValue ";
          opserr << argv[loc] << endln;
          return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[loc + 1], &gDir) != TCL_OK) {
          opserr << "WARNING recorder Remove -g gValue? gDir? gPat?... invalid "
                    "gDir ";
          opserr << argv[loc + 1] << endln;
          return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[loc + 2], &gPat) != TCL_OK) {
          opserr << "WARNING recorder Remove -g gValue? gDir? gPat?... invalid "
                    "gPat ";
          opserr << argv[loc + 1] << endln;
          return TCL_ERROR;
        }
        loc += 3;
      }

      else if ((strcmp(argv[loc], "-section") == 0) ||
               (strcmp(argv[loc], "-sec") == 0) ||
               (strcmp(argv[loc], "-comp") == 0)) {

        loc++;

        int secID;
        int endSecIDs = loc;
        int numSec = 0;
        while (Tcl_GetInt(interp, argv[endSecIDs], &secID) == TCL_OK) {
          endSecIDs++;
        }

        numSec = endSecIDs - loc;
        // create an ID to hold section/material tags
        secIDs = ID(numSec);

        // read in the sec tags to the ID
        for (int i = loc; i < endSecIDs; ++i) {
          if (Tcl_GetInt(interp, argv[i], &secID) != TCL_OK)
            return TCL_ERROR;
          secIDs[i - loc] = secID;
        }

        loc = endSecIDs;
      }

      else if (strcmp(argv[loc], "-criteria") == 0 ||
               strcmp(argv[loc], "-crit") == 0) {

        int critTag = 0;
        double critValue = 0;

        if (strcmp(argv[loc + 1], "minStrain") == 0 ||
            strcmp(argv[loc + 1], "1") == 0)
          critTag = 1;
        else if (strcmp(argv[loc + 1], "maxStrain") == 0 ||
                 strcmp(argv[loc + 1], "2") == 0)
          critTag = 2;
        else if (strcmp(argv[loc + 1], "axialDI") == 0 ||
                 strcmp(argv[loc + 1], "3") == 0)
          critTag = 3;
        else if (strcmp(argv[loc + 1], "flexureDI") == 0 ||
                 strcmp(argv[loc + 1], "4") == 0)
          critTag = 4;
        else if (strcmp(argv[loc + 1], "axialLS") == 0 ||
                 strcmp(argv[loc + 1], "5") == 0)
          critTag = 5;
        else if (strcmp(argv[loc + 1], "shearLS") == 0 ||
                 strcmp(argv[loc + 1], "6") == 0)
          critTag = 6;
        else if (strcmp(argv[loc + 1], "INFILLWALL") == 0 ||
                 strcmp(argv[loc + 1], "7") == 0)
          critTag = 7;
        else {
          opserr << "Error: RemoveRecorder - Removal Criteria " << argv[loc + 1]
                 << " not recognized" << endln;
          return TCL_ERROR;
        }
        if (critTag != 7) {
          if (Tcl_GetDouble(interp, argv[loc + 2], &critValue) != TCL_OK) {
            opserr << "WARNING recorder Remove -crit critTag? critValue?... "
                      "invalid critValue ";
            opserr << argv[loc + 1] << endln;
            return TCL_ERROR;
          }
        }

        remCriteria[2 * numCrit] = critTag;
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
          secIDs[1] = 1; // this is not used directly, gets rid of the "-sec"
                         // for infillwall
        }
      }

      else {
        // first unknown string then is assumed to start
        // element response request starts
        // eleData = loc;
        flags = 1;
      }
    }

    // if user has specified no element tags lets assume he wants them all
    if (numEle == 0) {
      ElementIter &theEleIter = domain->getElements();
      Element *theEle;
      while ((theEle = theEleIter()) != nullptr)
        eleIDs[numEle++] = theEle->getTag();
    }

    theOutputStream = new DummyStream();
    // now construct the recorder
    (*theRecorder) = new RemoveRecorder(
        nodeTag, eleIDs, secIDs, secondaryEleIDs, remCriteria, *domain,
        *theOutputStream, echoTime, dT, 1e-6, filename, eleMass, gAcc, gDir, gPat,
        nTagbotn, nTagmidn, nTagtopn, globgrav, filenameinf);

  }

  else if ((strcasecmp(argv[1], "Node") == 0) ||
           (strcasecmp(argv[1], "NodeRMS") == 0) ||
           (strcasecmp(argv[1], "EnvelopeNode") == 0) ||
           (strcasecmp(argv[1], "NodeEnvelope") == 0)) {
    return createNodeRecorder(clientData, interp, argc, argv, theRecorder);
  }

  else if (strcmp(argv[1], "Pattern") == 0) {
    if (argc < 4) {
      opserr << "WARNING recorder Pattern filename? <startFlag> patternTag?";
      return TCL_ERROR;
    }

    int flag = 0;
    if (strcmp(argv[3], "-time") == 0)
      flag = 1;
    if (strcmp(argv[3], "-load") == 0)
      flag = 2;

    int pos = 3;
    if (flag != 0)
      pos = 4;

    int patternTag;

    if (Tcl_GetInt(interp, argv[pos++], &patternTag) != TCL_OK)
      return TCL_ERROR;

    (*theRecorder) =
        new PatternRecorder(patternTag, *domain, argv[2], 0.0, flag);
  }

  // Create a recorder to write nodal drifts to a file
  else if ((strcmp(argv[1], "Drift") == 0) ||
           (strcmp(argv[1], "EnvelopeDrift") == 0)) {

    OutputOptions options;

    bool echoTimeFlag = false;
    ID iNodes(0, 16);
    ID jNodes(0, 16);
    int dof = 1;
    int perpDirn = 2;
    double dT = 0.0;

    int pos = 2;
    while (pos < argc) {
      int consumed;
      if ((consumed = parseOutputOption(&options, interp, argc-pos, &argv[pos])) != 0) {
        if (consumed > 0)
          pos += consumed;
        else
          return TCL_ERROR;
      }


      else if (strcmp(argv[pos], "-time") == 0) {
        echoTimeFlag = true;
        pos += 1;

      }

      else if (strcmp(argv[pos], "-dT") == 0) {
        // allow user to specify time step size for recording
        pos++;
        if (Tcl_GetDouble(interp, argv[pos], &dT) != TCL_OK)
          return TCL_ERROR;
        pos++;
      } 

      else if ((strcmp(argv[pos], "-iNode") == 0) ||
               (strcmp(argv[pos], "-iNodes") == 0)) {
        pos++;

        int node;
        int numNodes = 0;
        for (int j = pos; j < argc; j++)
          if (Tcl_GetInt(interp, argv[pos], &node) != TCL_OK) {
            j = argc;
            Tcl_ResetResult(interp);
          } else {
            iNodes[numNodes] = node;
            numNodes++;
            pos++;
          }
      }

      else if ((strcmp(argv[pos], "-jNode") == 0) ||
               (strcmp(argv[pos], "-jNodes") == 0)) {
        pos++;
        int node;
        int numNodes = 0;
        for (int j = pos; j < argc; j++)
          if (Tcl_GetInt(interp, argv[pos], &node) != TCL_OK) {
            j = argc;
            Tcl_ResetResult(interp);
          } else {
            jNodes[numNodes] = node;
            numNodes++;
            pos++;
          }
      }

      else if (strcmp(argv[pos], "-dof") == 0) {
        if (Tcl_GetInt(interp, argv[pos + 1], &dof) != TCL_OK) {
          pos = argc;
        }
        pos += 2;
      }

      else if (strcmp(argv[pos], "-perpDirn") == 0) {
        if (Tcl_GetInt(interp, argv[pos + 1], &perpDirn) != TCL_OK) {
          pos = argc;
        }
        pos += 2;
      }

      else
        pos++;
    }

    if (iNodes.Size() != jNodes.Size()) {
      opserr << "WARNING recorder Drift - the number of iNodes and jNodes must "
                "be the same "
             << iNodes << " " << jNodes << endln;
      return TCL_ERROR;
    }

    // construct the DataHandler
    theOutputStream = createOutputStream(options);


    // Subtract one from dof and perpDirn for C indexing
    if (strcmp(argv[1], "Drift") == 0)
      (*theRecorder) =
          new DriftRecorder(iNodes, jNodes, dof - 1, perpDirn - 1, *domain,
                            *theOutputStream, echoTimeFlag, dT);
    else
      (*theRecorder) =
          new EnvelopeDriftRecorder(iNodes, jNodes, dof - 1, perpDirn - 1,
                                    *domain, *theOutputStream, echoTimeFlag);

  }

  else if (strcmp(argv[1], "display") == 0) {
    // TODO: Handle "recorder display"
    return TCL_OK;
  }

  else if (strcmp(argv[1], "plot") == 0) {
    // TODO: handle "recorder plot"
    return TCL_OK;
  }

  else if (strcmp(argv[1], "vtk") == 0 || strcmp(argv[1], "VTK") == 0) {
    OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, nullptr);
    (*theRecorder) = (Recorder *)OPS_VTK_Recorder(rt, argc, argv);
  } 
  else if (strcmp(argv[1], "ElementRMS") == 0) {
    OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, nullptr);
    (*theRecorder) = (Recorder *)OPS_ElementRecorderRMS(rt, argc, argv);
  }

  else if (strcmp(argv[1], "gmsh") == 0 || strcmp(argv[1], "GMSH") == 0) {
    OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, nullptr);
    (*theRecorder) = (Recorder *)OPS_GmshRecorder(rt, argc, argv);
  }

  // else if (strcmp(argv[1],"gmshparallel") == 0 ||
  // strcmp(argv[1],"GMSHPARALLEL") == 0) {
  //  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
  //  (*theRecorder) = (Recorder*) OPS_GmshRecorderParallel();
  //  }

#if 0
  else if (strcmp(argv[1], "mpco") == 0) {
      OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
    (*theRecorder) = (Recorder*)OPS_MPCORecorder(rt, argc, argv);
    if (theRecorder == 0) {
      return TCL_ERROR;
    }
  }
#endif

#if 0
  else if (strcmp(argv[1],"GSA") == 0) {
      if (argc < 3) {
        opserr << argc;
        opserr << "WARNING recorder GSA -file filename? -dT deltaT? - not enough arguments\n"; 
        return TCL_ERROR;
      }
      TCL_Char *filename = 0;
      TCL_Char *title1 =0;
      TCL_Char *title2 =0;
      TCL_Char *title3 =0;
      TCL_Char *jobno  =0;
      TCL_Char *initials =0;
      TCL_Char *spec     =0;
      TCL_Char *currency =0;
      TCL_Char *length =0;
      TCL_Char *force  =0;
      TCL_Char *temp   =0;
      double dT = 0.0;
      int loc = 2;

      while (loc < argc) {
        if ((strcmp(argv[loc],"-file") == 0) ||
            (strcmp(argv[loc],"-file") == 0)) {
          filename = argv[loc+1];
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

      GSA_Recorder *theR = new GSA_Recorder(theDomain, filename, title1, title2,
              title3, jobno, initials, spec, currency, length, force, temp, dT);
      (*theRecorder) = theR;
  }
#endif

  else {

    // try existing loaded packages
    ExternalRecorderCommand *recorderCommands = theExternalRecorderCommands;
    bool found = false;

    while (recorderCommands != NULL && found == false) {
      if (strcmp(argv[1], recorderCommands->funcName) == 0) {

        OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, domain);
        void *theRes = (*(recorderCommands->funcPtr))();
        if (theRes != 0) {
          *theRecorder = (Recorder *)theRes;
          found = true;
        }
      } else
        recorderCommands = recorderCommands->next;
    }

    // if not there try loading package
  }

  if (*theRecorder == nullptr) {
    opserr << G3_ERROR_PROMPT << "No recorder exists "
           << "with type '" << argv[1] << "'\n";

    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclAddRecorder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{

  Domain* domain = (Domain*)clientData;

  Recorder *theRecorder = nullptr;

  if (TclCreateRecorder(clientData, interp, argc, argv, *domain, &theRecorder) != TCL_OK)
    return TCL_ERROR;

  else if (theRecorder == nullptr) {
    Tcl_SetObjResult(interp, Tcl_NewIntObj(-1));
    return TCL_OK;
  }

  else if ((domain->addRecorder(*theRecorder)) < 0) {
    opserr << G3_ERROR_PROMPT << "Failed to add recorder to domain" << endln;
    delete theRecorder;
    Tcl_SetObjResult(interp, Tcl_NewIntObj(-1));
    return TCL_ERROR;
  }

  int recorderTag = theRecorder->getTag();
  Tcl_SetObjResult(interp, Tcl_NewIntObj(recorderTag));
  return TCL_OK;
}

int
TclCommand_record(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ((Domain*)clientData)->record(false);
  return TCL_OK;
}


// by SAJalali
int
OPS_recorderValue(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char ** const argv)
{
  Domain *domain = (Domain*)clientData;

  // clmnID starts from 1
  if (argc < 3) {
    opserr << "WARNING want - recorderValue recorderTag clmnID <rowOffset> "
              "<-reset>\n";
    return TCL_ERROR;
  }

  int tag, rowOffset;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING recorderValue recorderTag? clmnID <rowOffset> <-reset> "
              "could not read recorderTag\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING recorderValue recorderTag? clmnID - could not read "
              "clmnID \n";
    return TCL_ERROR;
  }
  dof--;
  rowOffset = 0;
  int curArg = 3;
  if (argc > curArg) {
    if (Tcl_GetInt(interp, argv[curArg], &rowOffset) != TCL_OK) {
      opserr << "WARNING recorderValue recorderTag? clmnID <rowOffset> "
                "<-reset> could not read rowOffset \n";
      return TCL_ERROR;
    }
    curArg++;
  }
  bool reset = false;
  if (argc > curArg) {
    if (strcmp(argv[curArg], "-reset") == 0)
      reset = true;
    curArg++;
  }
  Recorder *theRecorder = domain->getRecorder(tag);
  double res = theRecorder->getRecordedValue(dof, rowOffset, reset);

  // now we copy the value to the tcl string that is returned
  char buffer[40];
  sprintf(buffer, "%35.8f", res);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);
  return TCL_OK;
}

static int
createNodeRecorder(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char ** const argv, Recorder **theRecorder)
{
  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData;
  G3_Runtime *rt = G3_getRuntime(interp);

  // make sure at least one other argument to contain integrator
  if (argc < 2) {
    opserr << "WARNING need to specify a Recorder type\n";
    return TCL_ERROR;
  }

  //
  // check argv[1] for type of Recorder, parse in rest of arguments
  // needed for the type of Recorder, create the object and add to Domain
  //
  OPS_Stream   *theOutputStream = nullptr;
  TCL_Char     *responseID      = nullptr;
  bool         echoTimeFlag     = false;
  double       dT               = 0.0;
  double       rTolDt           = 1e-5;
  int          numNodes         = 0;
  int          gradIndex        = -1;
  int          dataIndex        = -1;
  NodeData     dataFlag;

  // create ID's to contain the node tags and dofs
  ID *theNodes = nullptr;
  ID theDofs(0, MAX_NDF);
  ID theTimeSeriesID(0, MAX_NDF);

  TimeSeries **theTimeSeries = nullptr;

  OutputOptions options;

  if (argc < 7) {
    opserr << G3_ERROR_PROMPT << "recorder Node ";
    opserr << "-node <list nodes> -dof <doflist> -file <filename> -dT <dT> "
              "<reponse>";
    return TCL_ERROR;
  }


  int pos = 2;
  while (pos < argc) {
    int consumed;
    if ((consumed = parseOutputOption(&options, interp, argc-pos, &argv[pos])) != 0) {
      if (consumed > 0)
        pos += consumed;
      else
        return TCL_ERROR;
    }

    else if ((strcmp(argv[pos], "-time") == 0) ||
        (strcmp(argv[pos], "-load") == 0)) {
      echoTimeFlag = true;
      pos++;
    }

    else if (strcmp(argv[pos], "-dT") == 0) {
      // allow user to specify time step size for recording
      pos++;
      if (Tcl_GetDouble(interp, argv[pos], &dT) != TCL_OK)
        return TCL_ERROR;
      pos++;
    }
    else if (strcmp(argv[pos], "-rTolDt") == 0) {
      pos++;
      if (Tcl_GetDouble(interp, argv[pos], &rTolDt) != TCL_OK)
        return TCL_ERROR;
      pos++;
    }
    else if (strcmp(argv[pos], "-timeSeries") == 0) {

      pos++;
      int numTimeSeries = 0;
      int dof;

      for (int j = pos; j < argc; j++) {
        if (Tcl_GetInt(interp, argv[pos], &dof) != TCL_OK) {
          j = argc;
          Tcl_ResetResult(interp);
        } else {
          theTimeSeriesID[numTimeSeries] = dof; // -1 for c indexing of the dof's
          numTimeSeries++;
          pos++;
        }
      }

      theTimeSeries = new TimeSeries *[numTimeSeries];
      for (int j = 0; j < numTimeSeries; j++) {
        int timeSeriesTag = theTimeSeriesID(j);
        if (timeSeriesTag != 0 && timeSeriesTag != -1) {
          theTimeSeries[j] = G3_getTimeSeries(rt, theTimeSeriesID(j));
        } else {
          theTimeSeries[j] = 0;
        }
      }
    }

    else if ((strcmp(argv[pos], "-node") == 0) ||
             (strcmp(argv[pos], "-tag") == 0)  ||
             (strcmp(argv[pos], "-nodes") == 0)) {
      pos++;

      // read in the node tags or 'all' can be used
      if (strcmp(argv[pos], "all") == 0) {
        opserr << "recoder Node - error -all option has been removed, use "
                  "-nodeRange instaed\n";
        return TCL_ERROR;

      } else {
        theNodes = new ID(0, 16);
        int node;
        for (int j = pos; j < argc; j++)
          if (Tcl_GetInt(interp, argv[pos], &node) != TCL_OK) {
            j = argc;
            Tcl_ResetResult(interp);

          } else if (domain->getNode(node) == nullptr) {
            delete theNodes;
            theNodes = nullptr;
            opserr << G3_ERROR_PROMPT << "cannot find node with tag " << node << "\n";
            return TCL_ERROR;

          } else {
            (*theNodes)[numNodes] = node;
            numNodes++;
            pos++;
          }
      }
    }

    else if ((strcmp(argv[pos], "-nodeRange") == 0) || 
             (strcmp(argv[pos], "-range")==0)) {
      // ensure no segmentation fault if user messes up
      if (argc < pos + 3) {
        opserr << G3_ERROR_PROMPT << "recorder " << argv[1] 
               << " .. -range start? end?  .. - missing start/end tags\n";
        return TCL_ERROR;
      }

      // read in start and end tags of two elements & add set [start,end]
      int start, end;
      if (Tcl_GetInt(interp, argv[pos + 1], &start) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "recorder " << argv[1] 
               << " -range start? end? - invalid start "
               << argv[pos + 1] << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[pos + 2], &end) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "recorder " << argv[1] 
               << " -range start? end? - invalid end "
               << argv[pos + 2] << endln;
        return TCL_ERROR;
      }

      if (start > end) {
        int swap = end;
        end = start;
        start = swap;
      }

      theNodes = new ID(end - start + 1);
      for (int i = start; i <= end; ++i)
        (*theNodes)[numNodes++] = i;

      pos += 3;
    }

    else if (strcmp(argv[pos], "-region") == 0) {
      // allow user to specif elements via a region
      if (argc < pos + 2) {
        opserr << "WARNING recorder Node .. -region tag?  .. - no region "
                  "specified\n";
        return TCL_ERROR;
      }
      int tag;
      if (Tcl_GetInt(interp, argv[pos + 1], &tag) != TCL_OK) {
        opserr << "WARNING recorder Node -region tag? - invalid tag "
               << argv[pos + 1] << endln;
        return TCL_ERROR;
      }
      MeshRegion *theRegion = domain->getRegion(tag);
      if (theRegion == nullptr) {
        opserr << "WARNING recorder Node -region " << tag
               << " - region does not exist" << endln;
        return TCL_OK;
      }

      const ID &nodeRegion = theRegion->getNodes();
      theNodes = new ID(nodeRegion);

      pos += 2;
    }

    else if (strcmp(argv[pos], "-dof") == 0) {
      pos++;
      int numDOF = 0;
      int dof;
      for (int j = pos; j < argc; j++)
        if (Tcl_GetInt(interp, argv[pos], &dof) != TCL_OK) {
          j = argc;
          Tcl_ResetResult(interp);
        } else {
          theDofs[numDOF] = dof - 1; // -1 for c indexing of the dof's
          numDOF++;
          pos++;
        }
    }
    // AddingSensitivity:BEGIN //////////////////////////////////////
    else if (strcmp(argv[pos], "-sensitivity") == 0) {
      pos++;
      int paramTag;
      if (Tcl_GetInt(interp, argv[pos], &paramTag) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid parameter tag to node recorder." << endln;
        return TCL_ERROR;
      }
      pos++;

      // Now get gradIndex from parameter tag
      Parameter *theParameter = domain->getParameter(paramTag);
      if (theParameter == nullptr) {
        opserr << G3_ERROR_PROMPT << "parameter " << paramTag << " not found"
               << endln;
        return TCL_ERROR;
      }
      gradIndex = theParameter->getGradIndex();
    }
    // AddingSensitivity:END ////////////////////////////////////////
    else if (responseID == nullptr && pos < argc) {
      responseID = argv[pos];
      pos++;
    }
    else if (pos < argc) {
      opserr << "WARNING Unknown argument " << argv[pos] << "\n";
    }
  } // while (pos < argc)

  if (responseID == nullptr) { // pos >= argc) {
    opserr << "WARNING: No response type specified for node recorder, will "
              "assume you meant -disp\n";
  }

  theOutputStream = createOutputStream(options);

  if (theTimeSeries != nullptr && theTimeSeriesID.Size() < theDofs.Size()) {
    opserr << G3_ERROR_PROMPT << "recorder Node/EnvelopNode # TimeSeries must equal # "
              "dof - IGNORING TimeSeries OPTION\n";
    for (int i = 0; i < theTimeSeriesID.Size(); ++i) {
      if (theTimeSeries[i] != nullptr)
        delete theTimeSeries[i];
      delete[] theTimeSeries;
      theTimeSeries = nullptr;
    }
  }


  if ((dataFlag = getNodeDataFlag(responseID, *domain, &dataIndex)) == NodeData::Unknown) {
    opserr << G3_ERROR_PROMPT << "invalid response ID '" << responseID << "'\n";
    return TCL_ERROR;
  }

  // Create and set the recorder
  if (strcasecmp(argv[1], "Node") == 0) {
    (*theRecorder) =
        new NodeRecorder(theDofs, theNodes, gradIndex, dataFlag, dataIndex, *domain,
                         *theOutputStream, dT, rTolDt, echoTimeFlag, theTimeSeries);

  } else if (strcasecmp(argv[1], "NodeRMS") == 0) {
    (*theRecorder) =
        new NodeRecorderRMS(theDofs, theNodes, dataFlag, dataIndex, *domain,
                            *theOutputStream, dT, rTolDt, theTimeSeries);

  } else {
    (*theRecorder) = new EnvelopeNodeRecorder(theDofs, theNodes, dataFlag, dataIndex,
                                              *domain, *theOutputStream, dT, rTolDt,
                                              echoTimeFlag, theTimeSeries);
  }

  if (*theRecorder != nullptr) {
    opsdbg << G3_DEBUG_PROMPT << "Created recorder \n";
    (*theRecorder)->Print(opsdbg, 0);
  }

  if (theNodes != nullptr)
    delete theNodes;

  return TCL_OK;
}

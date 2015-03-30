// $Revision: 1.2 $
// $Date: 2007-10-16 00:11:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/TclUpdateMaterialCommand.cpp,v $

// fmk

#include <TclModelBuilder.h>
#include <Domain.h>
#include <MatParameter.h>

#include <string.h>

int
TclCommand_UpdateMaterialsCommand(ClientData clientData, 
				  Tcl_Interp *interp, 
				  int argc,
				  TCL_Char **argv, 
				  TclModelBuilder *theTclBuilder,
				  Domain *theDomain)
{

  if (argc < 5) {
    opserr << "WARNING insufficient number of UpdateMaterialStage arguments\n";
    opserr << "Want: UpdateMaterialStage material matTag? stage value?" << endln;
    return TCL_ERROR;
  }
  
  if (strcmp(argv[1],"-material") != 0) {
    opserr << "WARNING UpdateMaterialStage: Only accept parameter '-material' for now" << endln;
    return TCL_ERROR;		
  }		
  
  int materialTag, value;
  double valueD;
  
  if (Tcl_GetInt(interp, argv[2], &materialTag) != TCL_OK) {
    opserr << "WARNING MYSstage: invalid material tag" << endln;
    return TCL_ERROR;		
  }

  int parTag = theDomain->getNumParameters();
  parTag++;

  if (argc > 5) {
    if (strcmp(argv[5],"-parameter") == 0) {
      if (Tcl_GetInt(interp, argv[6], &parTag) != TCL_OK) {
	opserr << "WARNING UpdateMaterialStage: invalid parameter tag" << endln;
	return TCL_ERROR;		
      }
    }
  }	

  MatParameter *theParameter = new MatParameter(parTag, materialTag, argv[3]);

  if (theDomain->addParameter(theParameter) == false) {
    opserr << "WARNING could not add updateMaterialStage - MaterialStageParameter to domain" << endln;
    return TCL_ERROR;		
  }
 
  int res = 0;
  if (Tcl_GetInt(interp, argv[4], &value) != TCL_OK) {
 
	  if (Tcl_GetDouble(interp, argv[4], &valueD) != TCL_OK) {
      opserr << "WARNING UpdateMaterialStage: could not read value" << endln;
      return TCL_ERROR;		
    } else {
	
      res = theDomain->updateParameter(parTag, valueD);
	 
	  theDomain->removeParameter(parTag);
	 
    }
  } else {
	 
  res = theDomain->updateParameter(parTag, value);
 
  theDomain->removeParameter(parTag);
 
  }
  return res;
}



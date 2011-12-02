// @ rkaul@stanford.edu
// @ ggd@stanford.edu

#include <TclModelBuilder.h>
#include <string.h>
#include <Vector.h>

#include <YieldSurface_BC.h>
#include <YS_Section2D01.h>
#include <YS_Section2D02.h>

static void printCommand(int argc, TCL_Char **argv)
{
  opserr << "Input command: ";
  for (int i=0; i<argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
} 

SectionForceDeformation*
TclModelBuilderYS_SectionCommand(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theBuilder)
{
  if (argc < 3) {
    opserr << "WARNING insufficient number of arguments\n";
    printCommand(argc, argv);
    return 0;
  }
  
  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid section tag\n";
    printCommand(argc, argv);
    return 0;
  }
  
  SectionForceDeformation *theModel = 0;

  if (strcmp(argv[1],"YS_Section2D01") == 0 ||
      strcmp(argv[1],"YS_Section2d01") == 0) {

    if (argc < 7) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc,argv);
      opserr << "Want: section YS_Section2D01 tag? E? A? Iz? ysTag? <algo?>" << endln;
      return 0;
    }

    int algo, ysTag;
    double E, A, Iz;  
    int indx = 3;
  
    if (Tcl_GetDouble (interp, argv[indx++], &E) != TCL_OK) {
      opserr << "WARNING invalid E" << endln;
      opserr << " section: " << tag << endln;
      return 0;
    }
  
    if (Tcl_GetDouble (interp, argv[indx++], &A) != TCL_OK) {
      opserr << "WARNING invalid A" << endln;
      opserr << " section: " << tag << endln;
      return 0;
    }
  
    if (Tcl_GetDouble (interp, argv[indx++], &Iz) != TCL_OK) {
      opserr << "WARNING invalid Iz" << endln;
      opserr << " section: " << tag << endln;
      return 0;
    }

    if (Tcl_GetInt (interp, argv[indx++], &ysTag) != TCL_OK) {
      opserr << "WARNING invalid ysTag" << endln;
      opserr << " section: " << tag << endln;
      return 0;
    }
    
    YieldSurface_BC *ys = theBuilder->getYieldSurface_BC(ysTag);

    if (ys == 0) {
      opserr << "WARNING yield surface does not exist\n";
      opserr << "yieldSurface: " << ysTag; 
      opserr << "\nsection YieldSurface: " << tag << endln;
      return 0;
    }

    bool useKr = true;
    if(argc > indx) {
      if (Tcl_GetInt(interp, argv[indx++], &algo) != TCL_OK) {
	opserr << "WARNING invalid algo" << endln;
	opserr << " section: " << tag << endln;
	return 0;
      }
      if(algo == 0)
	useKr = false;
    }
    
    theModel = new YS_Section2D01(tag,E, A, Iz, ys, useKr);
  }

  else if (strcmp(argv[1],"YS_Section2D02") == 0 ||
	   strcmp(argv[1],"YS_Section2d02") == 0) {

    if (argc < 8) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc,argv);
      opserr << "Want: section YS_Section2D01 tag? E? A? Iz? maxPlastRot? ysTag? <algo?>" << endln;
      return 0;
    }

    int algo, ysTag;
    double E, A, Iz, maxPlstkRot;    
    int indx = 3;
    
    if (Tcl_GetDouble (interp, argv[indx++], &E) != TCL_OK) {
      opserr << "WARNING invalid E" << endln;
      opserr << " section: " << tag << endln;
      return 0;
    }
    
    if (Tcl_GetDouble (interp, argv[indx++], &A) != TCL_OK) {
      opserr << "WARNING invalid A" << endln;
      opserr << " section: " << tag << endln;
      return 0;
    }
    
    if (Tcl_GetDouble (interp, argv[indx++], &Iz) != TCL_OK) {
      opserr << "WARNING invalid Iz" << endln;
      opserr << " section: " << tag << endln;
      return 0;
    }
    
    if (Tcl_GetDouble (interp, argv[indx++], &maxPlstkRot) != TCL_OK) {
      opserr << "WARNING maxPlstkRot " << endln;
      opserr << " section: " << tag << endln;
      return 0;
    }

    if (Tcl_GetInt (interp, argv[indx++], &ysTag) != TCL_OK) {
      opserr << "WARNING invalid ysTag" << endln;
      opserr << " section: " << tag << endln;
      return 0;
    }

    YieldSurface_BC *ys = theBuilder->getYieldSurface_BC(ysTag);

    if (ys == 0) {
      opserr << "WARNING yield surface does not exist\n";
      opserr << "yieldSurface: " << ysTag; 
      opserr << "\nsection YieldSurface: " << tag << endln;
      return 0;
    }

    bool useKr = true;
    if(argc > indx) {
      if (Tcl_GetInt(interp, argv[indx++], &algo) != TCL_OK) {
	opserr << "WARNING invalid algo" << endln;
	opserr << " section: " << tag << endln;
	return 0;
      }
      if(algo == 0)
	useKr = false;
    }
     
    theModel = new YS_Section2D02(tag,E, A, Iz, maxPlstkRot, ys, useKr);    
  }

  return theModel;
}

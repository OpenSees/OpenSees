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
// $Date: 2006-03-17 22:47:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/TclLimitState.cpp,v $
                                                                        
                                                                        
// Written: fmk
// Created: 02/06


#include <ArrayOfTaggedObjects.h>

#include <LimitStateMaterial.h>
#include <AxialCurve.h>	// KJE
#include <ThreePointCurve.h> // KJE
#include <ShearCurve.h> // KJE


#include <string.h>

static TaggedObjectStorage *theLimitCurves = 0;
static Domain *theDomain = 0;

int
Tcl_AddLimitCurve(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);


//
// the following are invoked to add and remove the limitCurve command to the interpreter
//

int
Tcl_AddLimitCurveCommand(Tcl_Interp *theInterp, Domain *theTclDomain)
{
  Tcl_CreateCommand(theInterp, "limitCurve", 
		    Tcl_AddLimitCurve,
		    (ClientData)NULL, NULL);
  theDomain = theTclDomain;
  return 0;
}

int
Tcl_RemoveLimitCurveCommand(Tcl_Interp *theInterp)
{
  if (theLimitCurves != 0) {
      theLimitCurves->clearAll();
      delete theLimitCurves;
      theLimitCurves = 0;
  }
  
  return Tcl_DeleteCommand(theInterp, "limitCurve");
}


//
// the following are used to add and get copy to the LimitCurves
//

int 
AddLimitCurve(LimitCurve &theCurve)
{
  if (theLimitCurves == 0) {
    theLimitCurves = new ArrayOfTaggedObjects(32);
    if (theLimitCurves == 0) {
      opserr << "AddLimitCurve() - failed to add limit curve: " << theCurve;
    }
  }
  bool result = theLimitCurves->addComponent(&theCurve);
  if (result == true)
    return 0;
  else {
    opserr << "AddLimitCurve() - failed to add limit curve: " << theCurve;
    return -1;
  }
}


LimitCurve *
GetLimitCurve(int tag)
{
  if (theLimitCurves == 0)
    return 0;

  TaggedObject *mc = theLimitCurves->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  LimitCurve *result = (LimitCurve *)mc;
  return result;
}

static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
	opserr << argv[i] << " ";
    opserr << endln;
} 

// 
// the proceure invoked by the interpreter when limitCurve is invoked
//

int
Tcl_AddLimitCurve(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

  // Make sure there is a minimum number of arguments
  if (argc < 8) {
    opserr << "WARNING insufficient number of limit curve arguments\n";
    opserr << "Want: limitCurve type? tag? <specific curve args>" << endln;
    return TCL_ERROR;
  }
  
  // Pointer to a limit curve that will be added to the model builder
  LimitCurve *theCurve = 0;
  
  // Check argv[1] for limit curve type
  
  ///////////////////////
  // Axial Limit Curve
  ///////////////////////
  if (strcmp(argv[1],"Axial") == 0) {
    if (argc != 9 && argc != 12 && argc != 14 && argc != 15) { //SDK
      opserr << "WARNING invalid number of arguments\n";
      //	    printCommand(argc,argv); // Commented out by Terje
      opserr << "Want: limitCurve Axial tag? eleTag? Fsw? Kdeg? Fres? defType? forType?" << endln; //SDK
      opserr << "<ndI? ndJ? dof? perpDirn? delta? eleRemove?>" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    int eleTag;
    double Fsw;//SDK
    
    
    double Kdeg;
    double Fres;
    int defType, forType;
    int ndI = 0;
    int ndJ = 0;
    int dof = 0;
    int perpDirn = 0;
    int eleRemove = 0;
    double delta = 0.0;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid Axial LimitCurve tag" << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetInt(interp, argv[3], &eleTag) != TCL_OK) {
      opserr << "WARNING invalid element tag for associated beam-column element (eleTag)\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp,argv[4], &Fsw) != TCL_OK) { //SDK
      opserr << "WARNING invalid Fsw\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return TCL_ERROR;
    }
    
    
    
    if (Tcl_GetDouble(interp,argv[5], &Kdeg) != TCL_OK) {
      opserr << "WARNING invalid degrading slope Kdeg\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp,argv[6], &Fres) != TCL_OK) {
      opserr << "WARNING invalid residual capacity Fres\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return TCL_ERROR;
    }
    
    
    if (Tcl_GetInt(interp,argv[7], &defType) != TCL_OK) {
      opserr << "WARNING invalid deformation type defType\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetInt(interp,argv[8], &forType) != TCL_OK) {
      opserr << "WARNING invalid force type forType\n";
      opserr << "LimitCurve Axial: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (defType == 2)
      {
	if (Tcl_GetInt(interp,argv[9], &ndI) != TCL_OK) {
	  opserr << "WARNING invalid node I\n";
	  opserr << "LimitCurve Axial: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp,argv[10], &ndJ) != TCL_OK) {
	  opserr << "WARNING invalid node J\n";
	  opserr << "LimitCurve Axial: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp,argv[11], &dof) != TCL_OK) {
	  opserr << "WARNING invalid degree of freedom for drift\n";
	  opserr << "LimitCurve Axial: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp,argv[12], &perpDirn) != TCL_OK) {
	  opserr << "WARNING invalid direction for column length\n";
	  opserr << "LimitCurve Axial: " << tag << endln;
	  return TCL_ERROR;
	}
      }
    
    
    if (argc >= 14) {
      
      if (Tcl_GetDouble(interp,argv[13], &delta) != TCL_OK) {
	opserr << "WARNING invalid shift in drift surface (delta)\n";
	opserr << "LimitCurve Axial: " << tag << endln;
	return TCL_ERROR;
      }
    } 
    if (argc >= 15) {
      
      if (Tcl_GetInt(interp,argv[14], &eleRemove) != TCL_OK) {
	opserr << "WARNING invalid element removal option\n";
	opserr << "LimitCurve Axial: " << tag << endln;
	return TCL_ERROR;
      }
    }
    
    
    // Parsing was successful, allocate the limit curve
    // Subtract one from dof and perpDirn for C indexing
    theCurve = new AxialCurve(interp, tag, eleTag, theDomain, Fsw, //SDK
			      Kdeg, Fres, defType, forType, ndI, ndJ, dof-1, perpDirn-1, 
			      delta, eleRemove);  
  }
  
  //////////////////////////
  // Three Point LimitCurve
  //////////////////////////
  else if (strcmp(argv[1],"ThreePoint") == 0) {
    if (argc < 14 || argc > 18) {
      opserr << "WARNING insufficient arguments\n";
      //	    printCommand(argc,argv); // Commented out by Terje
      opserr << "Want: limitCurve ThreePoint tag? eleTag? x1? y1? x2? y2? x3? y3?";
      opserr << "Kdeg? Fres? defType? forType?" << endln;
      opserr << "<ndI? ndJ? dof? perpDirn?>" << endln;
      return TCL_ERROR;
    }
    
    int tag;
    int eleTag;
    double Kdeg;
    double Fres;
    int defType, forType;
    double x1, y1;
    double x2, y2;
    double x3, y3;
    int ndI = 0;
    int ndJ = 0;
    int dof = 0;
    int perpDirn = 0;
    
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid limitCurve ThreePoint tag" << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetInt(interp, argv[3], &eleTag) != TCL_OK) {
      opserr << "WARNING invalid element tag for associated beam-column element (eleTag)\n";
      opserr << "LimitCurve ThreePoint: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[4], &x1) != TCL_OK) {
      opserr << "WARNING invalid x1\n";
      opserr << "limitCurve ThreePoint: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[5], &y1) != TCL_OK) {
      opserr << "WARNING invalid y1\n";
      opserr << "limitCurve ThreePoint: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[6], &x2) != TCL_OK) {
      opserr << "WARNING invalid x2\n";
      opserr << "limitCurve ThreePoint: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[7], &y2) != TCL_OK) {
      opserr << "WARNING invalid y2\n";
      opserr << "limitCurve ThreePoint: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[8], &x3) != TCL_OK) {
      opserr << "WARNING invalid x3\n";
      opserr << "limitCurve ThreePoint: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[9], &y3) != TCL_OK) {
      opserr << "WARNING invalid y3\n";
      opserr << "limitCurve ThreePoint: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp,argv[10], &Kdeg) != TCL_OK) {
      opserr << "WARNING invalid degrading slope Kdeg\n";
      opserr << "LimitCurve ThreePoint: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp,argv[11], &Fres) != TCL_OK) {
      opserr << "WARNING invalid residual capacity Fres\n";
      opserr << "LimitCurve ThreePoint: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetInt(interp,argv[12], &defType) != TCL_OK) {
      opserr << "WARNING invalid deformation type defType\n";
      opserr << "LimitCurve ThreePoint: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetInt(interp,argv[13], &forType) != TCL_OK) {
      opserr << "WARNING invalid force type forType\n";
      opserr << "LimitCurve ThreePoint: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (defType == 2)
      {
	if (Tcl_GetInt(interp,argv[14], &ndI) != TCL_OK) {
	  opserr << "WARNING invalid node I\n";
	  opserr << "LimitCurve ThreePoint: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp,argv[15], &ndJ) != TCL_OK) {
	  opserr << "WARNING invalid node J\n";
	  opserr << "LimitCurve ThreePoint: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp,argv[16], &dof) != TCL_OK) {
	  opserr << "WARNING invalid degree of freedom for drift\n";
	  opserr << "LimitCurve ThreePoint: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp,argv[17], &perpDirn) != TCL_OK) {
	  opserr << "WARNING invalid direction for column length\n";
	  opserr << "LimitCurve ThreePoint: " << tag << endln;
	  return TCL_ERROR;
	}
	//		cout << ndI << endln;
	//		cout << ndJ << endln; // Commented out by Terje
      }
    
    
    // Parsing was successful, allocate the material
    // Subtract one from dof and perpDirn for C indexing
    theCurve = new ThreePointCurve(tag, eleTag, theDomain, 
				   x1, y1, x2, y2, x3, y3, Kdeg, Fres, defType, forType,
				   ndI, ndJ, dof-1, perpDirn-1);         
  }
  
  /////////////////////////////
  // Shear Limit Curve
  /////////////////////////////
  else if (strcmp(argv[1],"Shear") == 0) {
    if (argc < 14 || argc > 19) { //SDK
      opserr << "WARNING insufficient arguments\n";
      //	    printCommand(argc,argv); // Commented out by Terje
      opserr << "Want: limitCurve Shear tag? eleTag? rho? fc? b? h? d? Fsw? "; //SDK
      opserr << "Kdeg? Fres? defType? forType?" << endln;
      opserr << "<ndI? ndJ? dof? perpDirn? delta?>" << endln;
      return TCL_ERROR;
    }
    
    int tag;
    int eleTag;
    double Kdeg;
    double Fres;
    int defType, forType;
    double rho;
    double fc;
    double b, h, d;
    int ndI = 0;
    int ndJ = 0;
    int dof = 0;
    int perpDirn = 0;
    double Fsw = 0.0; //SDK
    
    double delta =0.0;
    
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid limitCurve Shear tag" << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetInt(interp, argv[3], &eleTag) != TCL_OK) {
      opserr << "WARNING invalid element tag for associated beam-column element (eleTag)\n";
      opserr << "LimitCurve Shear: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[4], &rho) != TCL_OK) {
      opserr << "WARNING invalid trans reinf ratio\n";
      opserr << "limitCurve Shear: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[5], &fc) != TCL_OK) {
      opserr << "WARNING invalid concrete strength\n";
      opserr << "limitCurve Shear: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[6], &b) != TCL_OK) {
      opserr << "WARNING invalid b\n";
      opserr << "limitCurve Shear: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[7], &h) != TCL_OK) {
      opserr << "WARNING invalid h\n";
      opserr << "limitCurve Shear: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[8], &d) != TCL_OK) {
      opserr << "WARNING invalid d\n";
      opserr << "limitCurve Shear: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[9], &Fsw) != TCL_OK) { //SDK
      opserr << "WARNING invalid Fsw\n";
      opserr << "limitCurve Shear: " << tag << endln;
      return TCL_ERROR;	
    }
    
    
    if (Tcl_GetDouble(interp,argv[10], &Kdeg) != TCL_OK) {
      opserr << "WARNING invalid degrading slope Kdeg\n";
      opserr << "LimitCurve Shear: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp,argv[11], &Fres) != TCL_OK) {
      opserr << "WARNING invalid residual capacity Fres\n";
      opserr << "LimitCurve Shear: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetInt(interp,argv[12], &defType) != TCL_OK) {
      opserr << "WARNING invalid deformation type defType\n";
      opserr << "LimitCurve Shear: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetInt(interp,argv[13], &forType) != TCL_OK) {
      opserr << "WARNING invalid force type forType\n";
      opserr << "LimitCurve Shear: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (defType == 2)
      {
	if (Tcl_GetInt(interp,argv[14], &ndI) != TCL_OK) {
	  opserr << "WARNING invalid node I\n";
	  opserr << "LimitCurve Shear: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp,argv[15], &ndJ) != TCL_OK) {
	  opserr << "WARNING invalid node J\n";
	  opserr << "LimitCurve Shear: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp,argv[16], &dof) != TCL_OK) {
	  opserr << "WARNING invalid degree of freedom for drift\n";
	  opserr << "LimitCurve Shear: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp,argv[17], &perpDirn) != TCL_OK) {
	  opserr << "WARNING invalid direction for column length\n";
	  opserr << "LimitCurve Shear: " << tag << endln;
	  return TCL_ERROR;
	}
	
      }
    
    if (argc == 19) 
      {
	if (Tcl_GetDouble(interp,argv[18], &delta) != TCL_OK) {
	  opserr << "WARNING invalid shift in drift surface (delta)\n";
	  opserr << "LimitCurve Shear: " << tag << endln;
	  return TCL_ERROR;
	}
      }
    
    
    
    
    // Parsing was successful, allocate the material
    // Subtract one from dof and perpDirn for C indexing
    theCurve = new ShearCurve(tag, eleTag, theDomain, 
			      rho, fc, b, h, d, Fsw, Kdeg, Fres, defType, forType, //SDK
			      ndI, ndJ, dof-1, perpDirn-1, delta);         
  }
  
  
  
  else  {
    opserr << "WARNING unknown type of limitCurve: " << argv[1];
    opserr << "\nValid types: Axial, ThreePoint, Shear" << endln;
    return TCL_ERROR;
  }
  
  // Ensure we have created the Material, out of memory if got here and no material
  if (theCurve == 0) {
    opserr << "WARNING ran out of memory creating limitCurve\n";
    opserr << argv[1] << endln;
    return TCL_ERROR;
  }
  
  // Now add the curve to the modelBuilder
  if (AddLimitCurve(*theCurve) < 0) {
    opserr << "WARNING could not add limitCurve to the domain\n";
    opserr << *theCurve << endln;
    delete theCurve; // invoke the curve objects destructor, otherwise mem leak
    return TCL_ERROR;
  }
  
  return TCL_OK;
}


// 
// the procedure invoked by the interpreter when a uniaxialMetrail limtState is found
//

UniaxialMaterial *
Tcl_AddLimitStateMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (strcmp(argv[1],"LimitState") != 0)
    return 0;

  UniaxialMaterial *theMaterial = 0;

  if (argc != 20 && argc != 19 && argc != 16 && argc != 15 && argc != 22 && argc != 23) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc,argv);
    opserr << "Want: uniaxialMaterial LimitState tag? mom1p? rot1p? mom2p? rot2p? mom3p? rot3p? "
	 << "\nmom1n? rot1n? mom2n? rot2n? mom3n? rot3n? pinchX? pinchY? damfc1? damfc2? beta? "
	 << "\n<curveTag? curveType? degrade?>";
    return theMaterial;
  }

  int tag;
  double mom1p, mom2p, mom3p;
  double rot1p, rot2p, rot3p;
  double mom1n, mom2n, mom3n;
  double rot1n, rot2n, rot3n;
  double pinchX, pinchY;
  double damfc1, damfc2;
  double beta = 0.0;
  int curveTag;
  int curveType;
  int degrade = 0;
  
  int i = 2;
  
  if (Tcl_GetInt(interp, argv[i++], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial LimitState tag" << endln;
    return theMaterial;
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &mom1p) != TCL_OK) {
    opserr << "WARNING invalid mom1p\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &rot1p) != TCL_OK) {
    opserr << "WARNING invalid rot1p\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &mom2p) != TCL_OK) {
    opserr << "WARNING invalid mom2p\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &rot2p) != TCL_OK) {
    opserr << "WARNING invalid rot2p\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (argc > 16) {
    if (Tcl_GetDouble(interp, argv[i++], &mom3p) != TCL_OK) {
      opserr << "WARNING invalid mom3p\n";
      opserr << "LimitState material: " << tag << endln;
      return theMaterial;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &rot3p) != TCL_OK) {
      opserr << "WARNING invalid rot3p\n";
      opserr << "LimitState material: " << tag << endln;
      return theMaterial;
    }
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &mom1n) != TCL_OK) {
    opserr << "WARNING invalid mom1n\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &rot1n) != TCL_OK) {
    opserr << "WARNING invalid rot1n\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &mom2n) != TCL_OK) {
    opserr << "WARNING invalid mom2n\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &rot2n) != TCL_OK) {
    opserr << "WARNING invalid rot2n\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (argc > 16) {
    if (Tcl_GetDouble(interp, argv[i++], &mom3n) != TCL_OK) {
      opserr << "WARNING invalid mom3n\n";
      opserr << "LimitState material: " << tag << endln;
      return theMaterial;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &rot3n) != TCL_OK) {
      opserr << "WARNING invalid rot3n\n";
      opserr << "LimitState material: " << tag << endln;
      return theMaterial;
    }
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &pinchX) != TCL_OK) {
    opserr << "WARNING invalid pinchX\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &pinchY) != TCL_OK) {
    opserr << "WARNING invalid pinchY\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &damfc1) != TCL_OK) {
    opserr << "WARNING invalid damfc1\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (Tcl_GetDouble(interp, argv[i++], &damfc2) != TCL_OK) {
    opserr << "WARNING invalid damfc2\n";
    opserr << "LimitState material: " << tag << endln;
    return theMaterial;
  }
  
  if (argc == 20 || argc == 16 || argc >= 22 ) {
    if (Tcl_GetDouble(interp, argv[i++], &beta) != TCL_OK) {
      opserr << "WARNING invalid beta\n";
      opserr << "LimitState material: " << tag << endln;
      return theMaterial;
    }
  }
  
  if (argc == 22 || argc == 23) {
    if (Tcl_GetInt(interp, argv[i++], &curveTag) != TCL_OK) {
      opserr << "WARNING invalid tag for LimitCurve (curveTag)\n";
      opserr << "LimitState material: " << tag << endln;
      return theMaterial;
    }
    
    LimitCurve *theCurve = GetLimitCurve(curveTag);
    
    if (theCurve == 0) {
      opserr << "WARNING limit curve does not exist\n";
      opserr << "limit curve: " << curveTag; 
      opserr << "\nLimitStateMaterial: " << tag << endln;
      return theMaterial;
    }
    
    if (Tcl_GetInt(interp, argv[i++], &curveType) != TCL_OK) {
      opserr << "WARNING invalid curveType\n";
      opserr << "LimitState material: " << tag << endln;
      return theMaterial;
    }
    
    if (argc == 23) {
      if (Tcl_GetInt(interp, argv[i++], &degrade) != TCL_OK) {
	opserr << "WARNING invalid degrade option\n";
	opserr << "LimitState material: " << tag << endln;
	return theMaterial;
      }
    } 
    
    theMaterial = new LimitStateMaterial (tag,
					  mom1p, rot1p, mom2p, rot2p, mom3p, rot3p,
					  mom1n, rot1n, mom2n, rot2n, mom3n, rot3n,
					  pinchX, pinchY, damfc1, damfc2, beta,
					  *theCurve, curveType, degrade);
  }
  
  // Parsing was successful, allocate the material
  if (argc == 20 || argc == 19)		
    theMaterial = new LimitStateMaterial (tag, 
					  mom1p, rot1p, mom2p, rot2p, mom3p, rot3p,
					  mom1n, rot1n, mom2n, rot2n, mom3n, rot3n,
					  pinchX, pinchY, damfc1, damfc2, beta);
  
  else if (argc == 16 || argc == 15)
    theMaterial = new LimitStateMaterial (tag,
					  mom1p, rot1p, mom2p, rot2p,
					  mom1n, rot1n, mom2n, rot2n,
					  pinchX, pinchY, damfc1, damfc2, beta);
  
  return theMaterial;
}


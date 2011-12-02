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
// $Date: 2006-01-19 19:19:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/TclReinforcingSteel.cpp,v $

/* ****************************************************************** **
** THIS FILE WAS DEVELOPED AT UC DAVIS                                **
**                                                                    **
** Programmed by: Jon Mohle (jfmohle@ucdavis.edu)                     **
** Supervisor: Sashi Kunnath (skkunnath@ucdavis.edu)                  **
**                                                                    **
********************************************************************* */
// Written: Jon Mohle
// Created: August 2005

                                                                       

#include <TclModelBuilder.h>
#include <ReinforcingSteel.h>   // Jon Mohle

#include <Vector.h>
#include <string.h>
#include <tcl.h>

int
TclCommand_ReinforcingSteel(ClientData clientData, Tcl_Interp *interp, int argc, 
			    TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
  UniaxialMaterial *theMaterial = 0;
  if (argc < 9) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial ReinforcingSteel tag? fy? fu? Es? Esh? esh? eult? <-GABuck?> <-DMBuck?> <-CMFatigue?> <-MPCurveParams?> <-IsoHard?>" << endln;
	return TCL_ERROR;
  }
  
  int tag;
  double fy, fu, Es, Esh, esh, eult;
  double slen = 0.0;
  double Cf = 0.0;
  double alpha = -4.46;
  double Cd = 0.0;
  double beta = 1.0;
  double r = 1.0;
  double gama = 0.5;
  int buckModel = 0;
  double RC1 = 1.0/3.0;
  double RC2 = 18.0;
  double RC3 = 4.0;
  double a1 = 0.0;
  double hardLim = 0.01;
  
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial ReinforcingSteel tag" << endln;
    return TCL_ERROR;		
  }
  
  if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
    opserr << "WARNING invalid fy\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;	
  }
  
  if (Tcl_GetDouble(interp, argv[4], &fu) != TCL_OK) {
    opserr << "WARNING invalid fu\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;
  }
  
  if (Tcl_GetDouble(interp, argv[5], &Es) != TCL_OK) {
    opserr << "WARNING invalid Es\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;	
  }
  
  if (Tcl_GetDouble(interp, argv[6], &Esh) != TCL_OK) {
    opserr << "WARNING invalid Esh\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;	
  }
  
  if (Tcl_GetDouble(interp, argv[7], &esh) != TCL_OK) {
    opserr << "WARNING invalid esh\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;	
  }
  
  if (Tcl_GetDouble(interp, argv[8], &eult) != TCL_OK) {
    opserr << "WARNING invalid eult\n";
    opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
    return TCL_ERROR;	
  }
  int argLoc = 9;
  while (argc > argLoc) {
	  if (strcmp(argv[argLoc],"-GABuck") == 0) {
	    if (argc < ++argLoc+4)  {
	      opserr << "WARNING insufficient optional arguments for -GABuck\n";
		    opserr << "Want: <-GABuck lsr? beta? r? gama?>" << endln;
		    return TCL_ERROR;
	    }
      
	    buckModel = 1;
	    if (Tcl_GetDouble(interp, argv[argLoc++], &slen) != TCL_OK) {
	      opserr << "WARNING invalid lsr\n";
	      opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
	      return TCL_ERROR;	
	    }
	    if (Tcl_GetDouble(interp, argv[argLoc++], &beta) != TCL_OK) {
        opserr << "WARNING invalid beta\n";
        opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
        return TCL_ERROR;	
      }
    
      if (Tcl_GetDouble(interp, argv[argLoc++], &r) != TCL_OK) {
        opserr << "WARNING invalid r\n";
        opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
        return TCL_ERROR;	
      }
    
      if (Tcl_GetDouble(interp, argv[argLoc++], &gama) != TCL_OK) {
        opserr << "WARNING invalid gama\n";
        opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
        return TCL_ERROR;	
      }
    }
	  
    else if (strcmp(argv[argLoc],"-DMBuck") == 0) {
	    if (argc < ++argLoc+1)  {
	      opserr << "WARNING insufficient optional arguments for -DMBuck\n";
		    opserr << "Want: <-DMBuck lsr? <alpha?>>" << endln;
		    return TCL_ERROR;
	    }

	    buckModel = 2;
	    if (Tcl_GetDouble(interp, argv[argLoc++], &slen) != TCL_OK) {
	      opserr << "WARNING invalid lsr\n";
	      opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
	      return TCL_ERROR;	
	    }
      if (argc <= argLoc)  {
        beta = 1.0;
      } else if (argv[argLoc][0]!=45) {
        if (Tcl_GetDouble(interp, argv[argLoc++], &beta) != TCL_OK) {
          opserr << "WARNING invalid alpha\n";
          opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
          return TCL_ERROR;	
        }
      }
      if (beta<0.75 || beta>1.0)
        opserr << "WARNING alpha usually is between 0.75 and 1.0\n";
	  }
	
	  else if (strcmp(argv[argLoc],"-CMFatigue") == 0) {
	    if (argc < ++argLoc+3)  {
	      opserr << "WARNING insufficient optional arguments for -CMFatigue\n";
		    opserr << "Want: <-CMFatigue Cf? alpha? Cd?>" << endln;
		    return TCL_ERROR;
	    }
	    if (Tcl_GetDouble(interp, argv[argLoc++], &Cf) != TCL_OK) {
        opserr << "WARNING invalid Cf\n";
        opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
        return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[argLoc++], &alpha) != TCL_OK) {
        opserr << "WARNING invalid alpha\n";
        opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
        return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[argLoc++], &Cd) != TCL_OK) {
        opserr << "WARNING invalid Cd\n";
        opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
        return TCL_ERROR;	
	    }
	  }
    else if (strcmp(argv[argLoc],"-MPCurveParams") == 0) {
      if (argc < ++argLoc+3)  {
	      opserr << "WARNING insufficient optional arguments for -MPCurveParams\n";
		    opserr << "Want: <-CMFatigue R1? R2? R3?>" << endln;
		    return TCL_ERROR;
	    }
      if (Tcl_GetDouble(interp, argv[argLoc++], &RC1) != TCL_OK) {
        opserr << "WARNING invalid RC1\n";
        opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
        return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[argLoc++], &RC2) != TCL_OK) {
        opserr << "WARNING invalid RC2\n";
        opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
        return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[argLoc++], &RC3) != TCL_OK) {
        opserr << "WARNING invalid RC3\n";
        opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
        return TCL_ERROR;	
	    }
    }
    else if (strcmp(argv[argLoc],"-IsoHard") == 0) {
      if (argc < ++argLoc+1) {
        a1 = 4.3;
        opserr << "uniaxialMaterial ReinforcingSteel -IsoHard: defaut values used\n";
      } else {
        if (argv[argLoc][0]==45) {
          a1 = 4.3;
        } else {
          if (Tcl_GetDouble(interp, argv[argLoc++], &a1) != TCL_OK) {
            opserr << "WARNING invalid a1\n";
            opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
            return TCL_ERROR;	
          }
        }
        if (argc > argLoc) {
          if (argv[argLoc][0]==45) {
            a1 = 4.3;
          } else {
            if (Tcl_GetDouble(interp, argv[argLoc++], &hardLim) != TCL_OK) {
              opserr << "WARNING invalid hardening limit\n";
              opserr << "uniaxialMaterial ReinforcingSteel: " << tag << endln;
              return TCL_ERROR;	
            }
          }
        }
      }
    }
    
    else {
	    opserr << "WARNING did not recognize optional flag\n";
	    opserr << "Possible Optional Flags: <-GABuck?> <-DMBuck?> <-CMFatigue?> <-MPCurveParams?> <-IsoHard?>" << endln;
	    return TCL_ERROR;
    }
  }

  // Parsing was successful, allocate the material
  theMaterial = new ReinforcingSteel(tag, fy, fu, Es, Esh, esh, eult, buckModel, slen, beta, r, gama, Cf, alpha, Cd, RC1, RC2, RC3, a1, hardLim);

  if (theMaterial != 0) 
    return theTclBuilder->addUniaxialMaterial(*theMaterial);
  else
    return -1;
}

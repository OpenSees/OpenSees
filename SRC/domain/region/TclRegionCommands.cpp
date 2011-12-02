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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-25 23:32:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/region/TclRegionCommands.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the function that is invoked
// by the interpreter when the comand 'region' is invoked by the 
// user.
//
// What: "@(#) Region.h, revA"

#include <tcl.h>
#include <tk.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Domain.h>

#include <MeshRegion.h>
#include <ID.h>

int
TclAddMeshRegion(ClientData clientData, Tcl_Interp *interp, int argc, 
		 TCL_Char **argv, Domain &theDomain)
{
  int loc = 1;
  int tag;
  double alphaM = 0.0;
  double betaK  = 0.0;
  double betaK0 = 0.0;
  double betaKc = 0.0;

  ID *theNodes = 0;
  ID *theElements = 0;
  int numNodes = 0;
  int numElements = 0;

  // first get tag for region
  if (argc < 2) {
    opserr << "WARNING region tag? - no tag specified\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
    opserr << "WARNING region tag? .. - invalid tag " << argv[loc] << endln;
    return TCL_ERROR;
  }

  loc++;

  // now contine until end of command
  while (loc < argc) {

    if (strcmp(argv[loc],"-ele") == 0) {
      
      // ensure no segmentation fault if user messes up
      if (argc < loc+2) {
	opserr << "WARNING region tag? .. -ele tag1? .. - no ele tags specified\n";
	return TCL_ERROR;
      }

      //
      // read in a list of ele until end of command or other flag
      //
      loc++;

      if (theElements == 0)
	theElements = new ID(0, 64);
      int eleTag;

      while (loc < argc && Tcl_GetInt(interp, argv[loc++], &eleTag) == TCL_OK) {
	(*theElements)[numElements++] = eleTag;
      }
      if (loc < argc) loc--;

    } else if (strcmp(argv[loc],"-eleRange") == 0) {

      // ensure no segmentation fault if user messes up
      if (argc < loc+3) {
	opserr << "WARNING region tag? .. -eleRange start? end?  .. - no ele tags specified\n";
	return TCL_ERROR;
      }

      //
      // read in start and end tags of two elements & add set [start,end]
      //

      int start, end;
      if (Tcl_GetInt(interp, argv[loc+1], &start) != TCL_OK) {
	opserr << "WARNING region tag? -eleRange start? end? - invalid start " << argv[loc+1] << endln;
	return TCL_ERROR;
      }      
      if (Tcl_GetInt(interp, argv[loc+2], &end) != TCL_OK) {
	opserr << "WARNING region tag? -eleRange start? end? - invalid end " << argv[loc+2] << endln;
	return TCL_ERROR;
      }      
      if (start > end) {
	int swap = end;
	end = start;
	start = swap;
      }
      int numEle = end-start+1;
      
      if (theElements == 0)
	theElements = new ID(0, numEle);
      for (int i=start; i<=end; i++)
	(*theElements)[numElements++] = i;

      loc += 3;
      
    } else if (strcmp(argv[loc],"-node") == 0) {

      // ensure no segmentation fault if user messes up
      if (argc < loc+2) {
	opserr << "WARNING region tag? .. -node tag1? .. - no node tags specified\n";
	return TCL_ERROR;
      }

      loc++;

      // read in list of nodes

      if (theNodes == 0)
	theNodes = new ID(0, 64);
      int nodTag;
      while (loc < argc && Tcl_GetInt(interp, argv[loc++], &nodTag) == TCL_OK) {
	(*theNodes)[numNodes++] = nodTag;
      }      
      
      if (loc < argc) loc--;

    } else if (strcmp(argv[loc],"-nodeRange") == 0) {

      // ensure no segmentation fault if user messes up
      if (argc < loc+3) {
	opserr << "WARNING region tag? .. -nodeRange start? end?  .. - no node tags specified\n";
	return TCL_ERROR;
      }

      // read in start and end ele tags
      int start, end;
      if (Tcl_GetInt(interp, argv[loc+1], &start) != TCL_OK) {
	opserr << "WARNING region tag? -eleRange start? end? - invalid start " << argv[loc+1] << endln;
	return TCL_ERROR;
      }      
      if (Tcl_GetInt(interp, argv[loc+2], &end) != TCL_OK) {
	opserr << "WARNING region tag? -eleRange start? end? - invalid end " << argv[loc+1] << endln;
	return TCL_ERROR;
      }      
      if (start > end) {
	int swap = end;
	end = start;
	start = swap;
      }
      int numNode = end-start+1;
      
      if (theNodes == 0)
	theNodes = new ID(0, numNode);
      for (int i=start; i<=end; i++)
	(*theNodes)[numNodes++] = i;

      loc += 3;
      
    } else if (strcmp(argv[loc],"-rayleigh") == 0) {

      // ensure no segmentation fault if user messes up
      if (argc < loc+5) {
	opserr << "WARNING region tag? .. -rayleigh aM? bK? bK0?  .. - not enough factors\n";
	return TCL_ERROR;
      }

      // read in rayleigh damping factors
      if (Tcl_GetDouble(interp, argv[loc+1], &alphaM) != TCL_OK) {
	opserr << "WARNING region tag? .. -rayleigh aM bK bK0 - invalid aM " << argv[loc+1] << endln;
	return TCL_ERROR;
      }      
      if (Tcl_GetDouble(interp, argv[loc+2], &betaK) != TCL_OK) {
	opserr << "WARNING region tag? .. -rayleigh aM bK bK0 - invalid bK " << argv[loc+2] << endln;
	return TCL_ERROR;
      }      
      if (Tcl_GetDouble(interp, argv[loc+3], &betaK0) != TCL_OK) {
	opserr << "WARNING region tag? .. -rayleigh aM bK bK0 - invalid bK0 " << argv[loc+3] << endln;
	return TCL_ERROR;
      }      
      if (Tcl_GetDouble(interp, argv[loc+4], &betaKc) != TCL_OK) {
	opserr << "WARNING region tag? .. -rayleigh aM bK bK0 - invalid bKc " << argv[loc+4] << endln;
	return TCL_ERROR;
      }      
      loc += 5;

    } else
      loc++;

  }

  MeshRegion *theRegion = new MeshRegion(tag);

  if ((theRegion == 0) || (theDomain.addRegion(*theRegion)) < 0) {
    opserr << "WARNING could not add to domain - region " << tag << endln;
    if (theRegion == 0) 
      opserr << "could not create region\n";
    else
      delete theRegion;
    return TCL_ERROR;
  } 

  // if elements or nodes have been set, set them in the Region
  if (theElements != 0)
    theRegion->setElements(*theElements);

  if (theNodes != 0) {
    if (theElements == 0)
      theRegion->setNodes(*theNodes);  
    else
      opserr << "WARNING region - both elements & nodes set, ONLY set using elements\n";
  }

  // if damping has been specified set the damping factors
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    theRegion->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);

  return TCL_OK;
}


/* ****************************************************************** **
**    Opensees - Open System for Earthquake Engineering Simulation    **
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
                                                                        
// $Revision: 1.0 $
// $Date: 2019-01-28 17:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/damping/TclDampingCommand.cpp,v $

// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 02/2020
// Revision: A

#include <string.h>
#include <TclModelBuilder.h>

#include <UniformDamping.h>
#include <SecStifDamping.h>

static Domain *theTclModelBuilderDomain = 0;
static TclModelBuilder *theTclModelBuilder = 0;

// 
// to create a coordinate transformation 
//
int
TclCommand_addDamping(ClientData clientData, Tcl_Interp *interp,
			 int argc, TCL_Char **argv,
			 Domain *theDomain,
			 TclModelBuilder *theBuilder)
  
{
  // Make sure there is a minimum number of arguments
  if (argc < 2) {
    opserr << "WARNING insufficient number of damping arguments\n";
    opserr << "Want: damping type? tag? <specific transf args>" << endln;
    return TCL_ERROR;
  }   
  
  theTclModelBuilderDomain = theDomain;
  theTclModelBuilder = theBuilder;
  
    int dampingTag;

    if (argc < 3) {
      opserr << "WARNING insufficient arguments - want: damping type? tag? <specific transf args>\n"; 
      return TCL_ERROR;
    }
	    
    if (Tcl_GetInt(interp, argv[2], &dampingTag) != TCL_OK) {	
      opserr << "WARNING invalid tag - want: damping type? tag? <specific transf args>\n";
      return  TCL_ERROR;
    }

    // construct the transformation object
    
    Damping *Damping =0;
    
    if (strcmp(argv[1],"Uniform") == 0)
    {
      double zeta, freq1, freq2;
      if (Tcl_GetDouble(interp, argv[3], &zeta) != TCL_OK)
      {	
        opserr << "WARNING invalid tag - want: damping type? tag? <specific damping args>\n";
        return  TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[4], &freq1) != TCL_OK)
      {	
        opserr << "WARNING invalid tag - want: damping type? tag? <specific damping args>\n";
        return  TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[5], &freq2) != TCL_OK)
      {	
        opserr << "WARNING invalid tag - want: damping type? tag? <specific damping args>\n";
        return  TCL_ERROR;
      }
      Damping = new UniformDamping(dampingTag, zeta*2.0, freq1, freq2);
      
    }
    else if (strcmp(argv[1],"SecStif") == 0 || strcmp(argv[1],"SecStiff") == 0)
    {
      double beta;
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK)
      {	
        opserr << "WARNING invalid tag - want: damping type? tag? <specific damping args>\n";
        return  TCL_ERROR;
      }
      Damping = new SecStifDamping(dampingTag, beta);
      
    }
    else {
      opserr << "WARNING TclModelBuilder - damping - invalid Type\n";
      opserr << argv[1] << endln;
      return TCL_ERROR;
    }
    
    if (Damping == 0) {
      opserr << "WARNING TclModelBuilder - damping - ran out of memory to create damping object\n";
      return TCL_ERROR;
    }
    
    // add the transformation to the modelBuilder
    if (OPS_addDamping(Damping) != true) {
      opserr << "WARNING TclElmtBuilder - damping  - could not add damping to model Builder\n";
      return TCL_ERROR;
    }
  
  return TCL_OK;
}

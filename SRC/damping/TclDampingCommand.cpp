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
// Created: 01/2020
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
      opserr << "WARNING insufficient arguments - want: damping type? tag? <specific damping args>\n"; 
      return TCL_ERROR;
    }
	    
    if (Tcl_GetInt(interp, argv[2], &dampingTag) != TCL_OK) {	
      opserr << "WARNING invalid tag - want: damping type? tag? <specific damping args>\n";
      return  TCL_ERROR;
    }

    // construct the transformation object
    
    Damping *Damping = 0;
    
    if (strcmp(argv[1],"Uniform") == 0)
    {
      double zeta, freq1, freq2;
      double ta = 0.0, td = 1e20;
      TimeSeries *facSeries = 0;
      if (Tcl_GetDouble(interp, argv[3], &zeta) != TCL_OK)
      {	
        opserr << "WARNING invalid damping ratio - want: damping Uniform tag zeta freq1 freq2\n";
        return  TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[4], &freq1) != TCL_OK)
      {	
        opserr << "WARNING invalid frequency range - want: damping Uniform tag zeta freq1 freq2\n";
        return  TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[5], &freq2) != TCL_OK)
      {	
        opserr << "WARNING invalid frequency range - want: damping Uniform tag zeta freq1 freq2\n";
        return  TCL_ERROR;
      }
      int count = 6;
      while (argc > count)
      {
        if ((strcmp(argv[count],"-activateTime") == 0) || (strcmp(argv[count],"-ActivateTime") == 0))
        {
          if (Tcl_GetDouble(interp, argv[count+1], &ta) != TCL_OK)
          {	
            opserr << "WARNING invalid activation time - want: damping Uniform tag zeta freq1 freq2 <-activateTime ta> <-deactivateTime td> <-fact tsTag>\n";
            return  TCL_ERROR;
          }
          count++;
        }
        else if ((strcmp(argv[count],"-deactivateTime") == 0) || (strcmp(argv[count],"-DeactivateTime") ==0 ))
        {
          if (Tcl_GetDouble(interp, argv[count+1], &td) != TCL_OK)
          {	
            opserr << "WARNING invalid deactivation time - want: damping Uniform tag zeta freq1 freq2 <-activateTime ta> <-deactivateTime td> <-fact tsTag>\n";
            return  TCL_ERROR;
          }
          count++;
        }
        else if ((strcmp(argv[count],"-fact") == 0) || (strcmp(argv[count],"-factor") ==0 ))
        {
          int tsTag;
          if (Tcl_GetInt(interp, argv[count+1], &tsTag) != TCL_OK)
          {	
            opserr << "WARNING invalid factor series - want: damping Uniform tag zeta freq1 freq2 <-activateTime ta> <-deactivateTime td> <-fact tsTag>\n";
            return  TCL_ERROR;
          }
          facSeries = OPS_getTimeSeries(tsTag);
          count++;
        }
        count++;
      }
      Damping = new UniformDamping(dampingTag, zeta*2.0, freq1, freq2, ta, td, facSeries);
    }
    else if (strcmp(argv[1],"SecStif") == 0 || strcmp(argv[1],"SecStiff") == 0)
    {
      double beta;
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK)
      {	
        opserr << "WARNING invalid damping factor - want: damping SecStiff tag beta\n";
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

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

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <string.h>
#include <TclModelBuilder.h>

#include <UniformDamping.h>
#include <SecStifDamping.h>
#include <URDDamping.h>
#include <URDDampingbeta.h>

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
      double ta = 0.0, td = 1e20;
      TimeSeries *facSeries = 0;
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK)
      {	
        opserr << "WARNING invalid damping factor - want: damping SecStiff tag beta\n";
        return  TCL_ERROR;
      }
      int count = 4;
      while (argc > count)
      {
        if ((strcmp(argv[count],"-activateTime") == 0) || (strcmp(argv[count],"-ActivateTime") == 0))
        {
          if (Tcl_GetDouble(interp, argv[count+1], &ta) != TCL_OK)
          {	
            opserr << "WARNING invalid activation time - want: damping SecStiff tag beta <-activateTime ta> <-deactivateTime td> <-fact tsTag>\n";
            return  TCL_ERROR;
          }
          count++;
        }
        else if ((strcmp(argv[count],"-deactivateTime") == 0) || (strcmp(argv[count],"-DeactivateTime") ==0 ))
        {
          if (Tcl_GetDouble(interp, argv[count+1], &td) != TCL_OK)
          {	
            opserr << "WARNING invalid deactivation time - want: damping SecStiff tag beta <-activateTime ta> <-deactivateTime td> <-fact tsTag>\n";
            return  TCL_ERROR;
          }
          count++;
        }
        else if ((strcmp(argv[count],"-fact") == 0) || (strcmp(argv[count],"-factor") ==0 ))
        {
          int tsTag;
          if (Tcl_GetInt(interp, argv[count+1], &tsTag) != TCL_OK)
          {	
            opserr << "WARNING invalid factor series - want: damping SecStiff tag beta <-activateTime ta> <-deactivateTime td> <-fact tsTag>\n";
            return  TCL_ERROR;
          }
          facSeries = OPS_getTimeSeries(tsTag);
          count++;
        }
        count++;
      }
      Damping = new SecStifDamping(dampingTag, beta, ta, td, facSeries);
      
    }
    else if (strcmp(argv[1], "URD") == 0) {
        int numfreq;
        int prttag = 0;
        int maxiter = 100;
        double tmpetafeq;
        double ta = 0.0, td = 1e20;
        double dptol = 0.05;
        TimeSeries* facSeries = 0;
        if (Tcl_GetInt(interp, argv[3], &numfreq) != TCL_OK)
        {
            opserr << "WARNING invalid number - want: damping URD tag n freq1 zeta1 ... freqn zetan\n";
            return  TCL_ERROR;
        }
        if (numfreq < 2) {
            opserr << "WARNING - n needs to be larger than 1\n ";
            return  TCL_ERROR;
        }
        Matrix* etaFreq = new Matrix(numfreq, 2);
        for (int i = 0; i < numfreq; i++) {
            for (int j = 0; j < 2; j++) {
                if (Tcl_GetDouble(interp, argv[4 + i*2 +j], &tmpetafeq) != TCL_OK) {
                    opserr << "WARNING invalid factor series - want: damping URD tag n freq1 zeta1 ... freqn zetan\n ";
                    return TCL_ERROR;
                }
                if (tmpetafeq <= 0.0) opserr << "URDDamping::URDDamping:  Invalid frequency or damping ratio\n";
                (*etaFreq)(i, j) = tmpetafeq;
            }  
            //(*etaFreq)(i, 1) *= (2.0);
        }

        int count = 4 + numfreq * 2;
        while (argc > count)
        {
            if ((strcmp(argv[count], "-activateTime") == 0) || (strcmp(argv[count], "-ActivateTime") == 0))
            {
                if (Tcl_GetDouble(interp, argv[count + 1], &ta) != TCL_OK)
                {
                    opserr << "WARNING invalid activation time - want: damping URD tag n freq1 zeta1 ... freqn zetan <-activateTime ta> <-deactivateTime td> <-fact tsTag> <-tol dptol> <-iter maxiter>\n";
                    return  TCL_ERROR;
                }
                count++;
            }
            else if ((strcmp(argv[count], "-deactivateTime") == 0) || (strcmp(argv[count], "-DeactivateTime") == 0))
            {
                if (Tcl_GetDouble(interp, argv[count + 1], &td) != TCL_OK)
                {
                    opserr << "WARNING invalid deactivation time - want: damping URD tag n freq1 zeta1 ... freqn zetan <-activateTime ta> <-deactivateTime td> <-fact tsTag> <-tol dptol> <-iter maxiter>\n";
                    return  TCL_ERROR;
                }
                count++;
            }
            else if ((strcmp(argv[count], "-fact") == 0) || (strcmp(argv[count], "-factor") == 0))
            {
                int tsTag;
                if (Tcl_GetInt(interp, argv[count + 1], &tsTag) != TCL_OK)
                {
                    opserr << "WARNING invalid factor series - want: damping URD tag n freq1 zeta1 ... freqn zetan <-activateTime ta> <-deactivateTime td> <-fact tsTag> <-tol dptol> <-iter maxiter>\n";
                    return  TCL_ERROR;
                }
                facSeries = OPS_getTimeSeries(tsTag);
                count++;
            }
            else if ((strcmp(argv[count], "-tol") == 0) || (strcmp(argv[count], "-tolerence") == 0))
            {
                if (Tcl_GetDouble(interp, argv[count + 1], &dptol) != TCL_OK)
                {
                    opserr << "WARNING invalid tolerance - want: damping URD tag n freq1 zeta1 ... freqn zetan <-activateTime ta> <-deactivateTime td> <-fact tsTag> <-tol dptol> <-iter maxiter>\n";
                    return  TCL_ERROR;
                }
                count++;
            }
            else if ((strcmp(argv[count], "-maxiter") == 0) || (strcmp(argv[count], "-iter") == 0))
            {
                if (Tcl_GetInt(interp, argv[count + 1], &maxiter) != TCL_OK)
                {
                    opserr << "WARNING invalid interation number - want: damping URD tag n freq1 zeta1 ... freqn zetan <-activateTime ta> <-deactivateTime td> <-fact tsTag> <-tol dptol> <-iter maxiter>\n";
                    return  TCL_ERROR;
                }
                count++;
            }
            else if ((strcmp(argv[count], "-prttag") == 0) || (strcmp(argv[count], "-print") == 0))
            {
                prttag = 1;
            }
            count++;
        }
        Damping = new URDDamping(dampingTag, numfreq, etaFreq, dptol, ta, td, facSeries, prttag, maxiter);
        delete etaFreq;
    }
    else if (strcmp(argv[1], "URDbeta") == 0) {
    int numfreq;
    double ta = 0.0, td = 1e20;
    //double tmpdata;
    TimeSeries* facSeries = 0;
    if (Tcl_GetInt(interp, argv[3], &numfreq) != TCL_OK)
    {
        opserr << "WARNING invalid number - want: damping URDbeta tag n freq1 ... freqn zeta1 ... zetan alpha1 ... alphan \n";
        return  TCL_ERROR;
    }
    if (numfreq < 2) {
        opserr << "WARNING - n needs to be larger than 1\n ";
        return  TCL_ERROR;
    }
    Vector* tmpbeta = new Vector(numfreq);
    Vector* tmpomegac = new Vector(numfreq);
    for (int i = 0; i < numfreq; i++) {
        if (Tcl_GetDouble(interp, argv[4 + 2 * i], &(*tmpomegac)(i)) != TCL_OK) {
            opserr << "WARNING invalid factor series - want: damping URDbeta tag n freqc1 beta1 ... freqcn betan\n ";
            return TCL_ERROR;
        }
        if ((*tmpomegac)(i) <= 0.0) opserr << "URDDamping::URDDamping:  Invalid frequency\n";
        (*tmpomegac)(i) *= (6.28318530718);
        if (Tcl_GetDouble(interp, argv[5 + 2 * i], &(*tmpbeta)(i)) != TCL_OK) {
            opserr << "WARNING invalid factor series - want: damping URDbeta tag n freqc1 beta1 ... freqcn betan\n ";
            return TCL_ERROR;
        }
    }
   
    int count = 4 + numfreq * 2;
    while (argc > count)
    {
        if ((strcmp(argv[count], "-activateTime") == 0) || (strcmp(argv[count], "-ActivateTime") == 0))
        {
            if (Tcl_GetDouble(interp, argv[count + 1], &ta) != TCL_OK)
            {
                opserr << "WARNING invalid activation time - want: damping URDbeta tag n freqc1 beta1 ... freqcn betan <-activateTime ta> <-deactivateTime td> <-fact tsTag> <-tol dptol>\n";
                return  TCL_ERROR;
            }
            count++;
        }
        else if ((strcmp(argv[count], "-deactivateTime") == 0) || (strcmp(argv[count], "-DeactivateTime") == 0))
        {
            if (Tcl_GetDouble(interp, argv[count + 1], &td) != TCL_OK)
            {
                opserr << "WARNING invalid deactivation time - want: damping URDbeta tag n freqc1 beta1 ... freqcn betan <-activateTime ta> <-deactivateTime td> <-fact tsTag> <-tol dptol>\n";
                return  TCL_ERROR;
            }
            count++;
        }
        else if ((strcmp(argv[count], "-fact") == 0) || (strcmp(argv[count], "-factor") == 0))
        {
            int tsTag;
            if (Tcl_GetInt(interp, argv[count + 1], &tsTag) != TCL_OK)
            {
                opserr << "WARNING invalid factor series - want: damping URDbeta tag n freqc1 beta1 ... freqcn betan <-activateTime ta> <-deactivateTime td> <-fact tsTag> <-tol dptol>\n";
                return  TCL_ERROR;
            }
            facSeries = OPS_getTimeSeries(tsTag);
            count++;
        }
        
        count++;
    }
    Damping = new URDDampingbeta(dampingTag, numfreq, tmpomegac, tmpbeta, ta, td, facSeries);
    delete tmpbeta;
    delete tmpomegac;
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
      opserr << "WARNING TclElmtBuilder - damping - could not add damping to model Builder\n";
      return TCL_ERROR;
    }
  
  return TCL_OK;
}

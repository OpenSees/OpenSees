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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-04-18 16:31:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/element/TclNLBeamColumnCommand.cpp,v $
                                                                                                                                 
// Written: Remo M. de Souza (rmsouza@ce.berkeley.edu)
// Created: 08/99
//
// Description: This file contains the implementation of the commands used 
// to add coordinate transformation objects  and nonlinear frame elements to the model.

#include <stdlib.h>
#include <string.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <SectionForceDeformation.h>

#include <NLBeamColumn2d.h>
#include <NLBeamColumn3d.h>

#include <LobattoBeamIntegration.h>
#include <ForceBeamColumn2d.h>
#include <ForceBeamColumn3d.h>

// #include <LargeDispBeamColumn3d.h>

#include <GaussLobattoQuadRule1d01.h>
#include <TclModelBuilder.h>

//
// some static variables used in the functions
//

static Domain *theTclModelBuilderDomain = 0;
static TclModelBuilder *theTclModelBuilder =0;

// 
// to create a NL frame element and add to the domain
//
int
TclModelBuilder_addNLBeamColumn(ClientData clientData, Tcl_Interp *interp,
				int inArgc, TCL_Char **inArgv,
				Domain *theDomain,
				TclModelBuilder *theBuilder)
				
{
  theTclModelBuilderDomain = theDomain;
  theTclModelBuilder = theBuilder;
    
  int NDM, NDF;
     
  NDM = theTclModelBuilder->getNDM();   // dimension of the structure (1d, 2d, or 3d)
  NDF = theTclModelBuilder->getNDF();   // number of degrees of freedom per node

  // split possible lists present in argv
  char *List;

  List = Tcl_Merge (inArgc, inArgv);
  if (List == 0)
  {
    opserr << "WARNING - TclModelBuilder_addFrameElement - problem merging list\n";
    return TCL_ERROR;
  }

//  opserr << "List :" << List << endln;

  // remove braces from list
  for (int i = 0; List[i] != '\0'; i++)
  {
    if ((List[i] == '{')  ||  (List[i] == '}'))
      List[i] = ' ';
  }
  
  int argc;
  TCL_Char **argv;
       
  if (Tcl_SplitList(interp, List, &argc, &argv) != TCL_OK)
  {
    opserr <<  "WARNING - TclModelBuilder_addFrameElement - problem spliting list\n";
    return TCL_ERROR;
  }
      
  Tcl_Free (List);
  
//  opserr << "argc : " << argc; 
//  for (int i=0; i<argc; i++)
//  {
//    opserr <<"string " << i << " : " << argv[i] << endln;
//  }


  // create plane frame elements
  if ((NDM == 2 && NDF == 3) || (NDM == 3 && NDF == 6)) {
    
    int eleTag, iNode, jNode, numIntgrPts, transfTag;
    int secTag[10]; // Max size of integration rule ... can change if needed
    
    if (argc < 8) {
      opserr << "WARNING bad command - want: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
      return TCL_ERROR;
    }
    int argi = 2;  
    if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK) {
      opserr << "WARNING invalid eleTag: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n"; 
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK) {
      opserr << "WARNING invalid iNode:  element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numIntgrPts) != TCL_OK) {
      opserr << "WARNING invalid numIntgrPts: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
      return TCL_ERROR;
    }

    if (strcmp(argv[argi], "-sections") == 0) {
      argi++;
      if (argi+numIntgrPts > argc) {
	opserr << "WARNING insufficient number of section tags - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
	return TCL_ERROR;
      }
      int section;
      for (int i = 0; i < numIntgrPts; i++) {
	if (Tcl_GetInt(interp, argv[argi+i], &section) != TCL_OK) {
	  opserr << "WARNING invalid secTag - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
	  return TCL_ERROR;
	}
	secTag[i] = section;
      }
      argi += numIntgrPts;
    }

    else {
      int section;
      if (Tcl_GetInt(interp, argv[argi++], &section) != TCL_OK) {
	opserr << "WARNING invalid secTag - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
	return TCL_ERROR;
      }
      for (int i = 0; i < numIntgrPts; i++)
	secTag[i] = section;
    }

    if (argi >= argc || Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK) {
      opserr << "WARNING invalid transfTag? - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
      return TCL_ERROR;
    }

    // some  additional options at end of command .. setting defaults of 10 and 1.0e-10
    double massDens = 0.0;
    int    nMaxLocIters = 10;
    double locToler = 1e-08;
    
    while (argi != argc) {
      if (strcmp(argv[argi],"-mass") == 0) {
	// allow user to specify mass (per unit length)
	argi++;
	if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK) {
	  opserr << "WARNING invalid massDens - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
	  return TCL_ERROR;
	} 
      }

      else if (strcmp(argv[argi],"-iter") == 0) {
	// allow user to specify maximum number of local iterations
	argi++;
	if (argi == argc || Tcl_GetInt(interp, argv[argi++], &nMaxLocIters) != TCL_OK) {
	  opserr << "WARNING invalid nMaxLocIters - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
	  return TCL_ERROR;
	} 

	// specify local tolerance 
	if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &locToler) != TCL_OK) {
	  opserr << "WARNING invalid locToler - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
	  return TCL_ERROR;
	} 
      }
      else {
	opserr << "WARNING bad command  - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>\n";
	opserr << "invalid: " << argv[argi] << endln;
	return TCL_ERROR;
      }
    }
    
    // create the element

    // get pointer to the sections for the whole beam

    SectionForceDeformation **sections = new SectionForceDeformation* [numIntgrPts];
    
    if (!sections) {
      opserr << "WARNING TclElmtBuilder - addFrameElement - Insufficient memory to create sections\n";
      return TCL_ERROR;
    }

    for (int j=0; j<numIntgrPts; j++) {
      SectionForceDeformation *theSection = theTclModelBuilder->getSection(secTag[j]);

      if (theSection == 0) {
	opserr << "WARNING TclElmtBuilder - frameElement - no Section found with tag ";
	opserr << secTag[j] << endln;
	delete [] sections;
	return TCL_ERROR;
      }

      sections[j] = theSection;
    }

    // opserr << "massDens " << massDens << endln;
     
    // construct the element

    Element *element = 0;
    if (NDM == 2) {
      CrdTransf2d *theCrdTransf = theTclModelBuilder->getCrdTransf2d(transfTag);
      
      if (theCrdTransf == 0) {
	opserr << "WARNING TclElmtBuilder - frameElement - no geometric transformation found with tag ";
	opserr << transfTag << endln;
	return TCL_ERROR;
      }
      
      //element = new NLBeamColumn2d(eleTag, iNode, jNode, numIntgrPts, sections,
      //			   *theCrdTransf, massDens, nMaxLocIters, locToler, 10);
      LobattoBeamIntegration beamIntegr;
      element = new ForceBeamColumn2d(eleTag, iNode, jNode, numIntgrPts,
				      sections, beamIntegr, *theCrdTransf,
				      massDens, nMaxLocIters, locToler);
      
      delete [] sections;
    }
    else {
      CrdTransf3d *theCrdTransf = theTclModelBuilder->getCrdTransf3d(transfTag);
      
      if (theCrdTransf == 0) {
	opserr << "WARNING TclElmtBuilder - frameElement - no geometric transformation found with tag ";
	opserr << transfTag << endln;
	return TCL_ERROR;
      }
      
      //element = new NLBeamColumn3d(eleTag, iNode, jNode, numIntgrPts, sections,
      //			   *theCrdTransf, massDens, nMaxLocIters, locToler);
        
      LobattoBeamIntegration beamIntegr;
      element = new ForceBeamColumn3d(eleTag, iNode, jNode, numIntgrPts,
				      sections, beamIntegr, *theCrdTransf,
				      massDens, nMaxLocIters, locToler);

      delete [] sections;
    }

    if (element == 0) {
      opserr << "WARNING  TclElmtBuilder - addFrameElement - ran out of memory to create element\n";
      return TCL_ERROR;
    }
   
    if (theTclModelBuilderDomain->addElement(element) == false) {
      opserr << "WARNING TclElmtBuilder - addFrameElement - could not add element to domain ";
      opserr << eleTag << endln;
      return TCL_ERROR;
    } 
    
  }
  else {
    opserr << "WARNING NDM = " << NDM << " and NDF = " << NDF << "is imcompatible with available frame elements\n";
    return TCL_ERROR;
  }      

  Tcl_Free ((char *)argv);

  // if get here we have sucessfully created the element and added it to the domain
  
  return TCL_OK;
}

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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/tcl/TclElmtBuilder.cpp,v $
                                                                                                                                 
// File: ~/tcl/TclElmtBuilder.C
// 
// Written: Remo M. de Souza (rmsouza@ce.berkeley.edu)
// Created: 08/99
//
// Description: This file contains the implementation of the commands used 
// to add coordinate transformation objects  and nonlinear frame elements to the model.

#include <stdlib.h>
#include <string.h>
#include <iostream.h>

#include <Domain.h>
#include <Node.h>
#include <ArrayOfTaggedObjects.h>
#include <Matrix.h>

#include <SectionForceDeformation.h>

#include <NLBeamColumn2d.h>
#include <NLBeamColumn3d.h>
//#include <LargeDispBeamColumn3d.h>

#include <LinearCrdTransf2d.h>
//#include <CorotCrdTransf2d.h>
#include <LinearCrdTransf3d.h>
//#include <CorotCrdTransf3d.h>

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
TclModelBuilder_addFrameElement(ClientData clientData, Tcl_Interp *interp,
				int inArgc, char **inArgv,
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
      interp->result = "WARNING - TclModelBuilder_addFrameElement - problem merging list";
      return TCL_ERROR;
  }

//  cerr << "List :" << List << endl;

  // remove braces from list
  for (int i = 0; List[i] != '\0'; i++)
  {
    if ((List[i] == '{')  ||  (List[i] == '}'))
      List[i] = ' ';
  }
  
  int argc;
  char **argv;
       
  if (Tcl_SplitList(interp, List, &argc, &argv) != TCL_OK)
  {
     interp->result = "WARNING - TclModelBuilder_addFrameElement - problem spliting list";
     return TCL_ERROR;
  }
      
  Tcl_Free (List);
  
//  cerr << "argc : " << argc; 
//  for (int i=0; i<argc; i++)
//  {
//    cerr <<"string " << i << " : " << argv[i] << endl;
//  }


  // create plane frame elements
  if (NDM == 2 && NDF == 3)     
  {
       
    int eleTag, iNode, jNode, numIntgrPts, secTag, transfTag;

    if (argc < 8)
    {
      interp->result = "WARNING bad command - want: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }
    int argi = 2;  
    if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK)
    {
      interp->result = "WARNING invalid eleTag: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>"; 
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK)
    {
      interp->result = "WARNING invalid iNode:  element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK)
    {
      interp->result = "WARNING invalid jNode: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numIntgrPts) != TCL_OK)
    {
      interp->result = "WARNING invalid numIntgrPts: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &secTag) != TCL_OK)
    {
      interp->result = "WARNING invalid secTag - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK)
    {
      interp->result = "WARNING invalid transfTag? - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }

    
    // allow some additional options at end of command
    double massDens = 0;
    int    nMaxLocIters = 1;
    double locToler = 1e-16;
 
    while (argi != argc) 
    {
      if (strcmp(argv[argi],"-mass") == 0) 
      {
         // allow user to specify mass (per unit length)
         argi++;
         if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK)
         {
            interp->result = "WARNING invalid massDens - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
            return TCL_ERROR;
         } 
      }

      else if (strcmp(argv[argi],"-iter") == 0) 
      {
         // allow user to specify maximum number of local iterations
         argi++;
         if (argi == argc || Tcl_GetInt(interp, argv[argi++], &nMaxLocIters) != TCL_OK)
         {
	    interp->result = "WARNING invalid nMaxLocIters - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
            return TCL_ERROR;
         } 

         // specify local tolerance 
         if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &locToler) != TCL_OK)
         {
            interp->result = "WARNING invalid locToler - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
            return TCL_ERROR;
         } 
      }
      else
      {
         interp->result = "WARNING bad command  - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
         cerr << "invalid: " << argv[argi] << endl;
         return TCL_ERROR;
      }
    }

    
    CrdTransf2d *theCrdTransf = theTclModelBuilder->getCrdTransf2d(transfTag);

    if (theCrdTransf == 0) 
    {
       cerr << "WARNING TclElmtBuilder - frameElement - no geometric transformation found with tag ";
       cerr << transfTag << endl;
       return TCL_ERROR;
    }
    
    
    SectionForceDeformation *theSection = theTclModelBuilder->getSection(secTag);

    if (theSection == 0) 
    {
       cerr << "WARNING TclElmtBuilder - frameElement - no Section found with tag ";
       cerr << secTag << endl;
       return TCL_ERROR;
    }

    // create the element

    // get pointer to the sections for the whole beam

    SectionForceDeformation **sections = new SectionForceDeformation* [numIntgrPts];
    
    if (!sections)
    {
      interp->result = "WARNING TclElmtBuilder - addFrameElement - Insufficient memory to create sections";
      return TCL_ERROR;
    }

    for (int j=0; j<numIntgrPts; j++)
       sections[j] = theSection;

    // cerr << "massDens " << massDens << endl;
     
    // construct the element

    // TEMPORARY - CHANGE LATER
    // GaussLobattoQuadRule1d01 *theQuadRule = new GaussLobattoQuadRule1d01 (5);
    
    Element *element;
    if (strcmp(argv[1],"nonlinearBeamColumn") == 0)
    {
       element = new NLBeamColumn2d (eleTag, iNode, jNode, numIntgrPts, sections, *theCrdTransf, massDens, nMaxLocIters, locToler);
       
     // REMO !!!!!!!!!!!!!!!!
     //element = new NLBeamColumn2d (eleTag, iNode, jNode, theQuadRule, sections, *theCrdTransf, massDens, nMaxLocIters, locToler);

       delete [] sections;
    }   
    else
    {
      interp->result = "WARNING TclElmtBuilder - addFrameElement - invalid elemType";
      return TCL_ERROR;
    }
     
    if (element == 0)
    {
      interp->result = "WARNING  TclElmtBuilder - addFrameElement - ran out of memory to create element";
      return TCL_ERROR;
    }
   
   if (theTclModelBuilderDomain->addElement(element) == false) 
   {
      cerr << "WARNING TclElmtBuilder - addFrameElement - could not add element to domain ";
      cerr << eleTag << endl;
      return TCL_ERROR;
    } 
  }
  
  // create 3d frame element 

  else if (NDM == 3 && NDF == 6)
  {
    int   eleTag, iNode, jNode, numIntgrPts, secTag, transfTag, PdeltaFlag = 0;
        
    if (argc < 8)
    {
      interp->result = "WARNING bad command - want: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }
    int argi = 2;
    if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK)
    {
      interp->result = "WARNING invalid eleTag: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK)
    {
      interp->result = "WARNING invalid iNode: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK)
    {
      interp->result = "WARNING invalid jNode: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numIntgrPts) != TCL_OK)
    {
      interp->result = "WARNING invalid numIntgrPts: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &secTag) != TCL_OK)
    {
      interp->result = "WARNING invalid secTag: element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }
   
    if (Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK)
    {
      interp->result = "WARNING invalid transfTag - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
      return TCL_ERROR;
    }

    // allow some additional options at end of command
    double massDens = 0;
    int    nMaxLocIters = 1;
    double locToler = 1e-16;
 
    while (argi != argc) 
    {
      if (strcmp(argv[argi],"-mass") == 0) 
      {
         // allow user to specify mass (per unit length)
	 argi++;
         if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK)
         {
            interp->result = "WARNING invalid massDens - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
            return TCL_ERROR;
         } 
      }
      else if (strcmp(argv[argi],"-Pdelta") == 0) 
      {
         // allow user to include second order effects along the element
	 argi++;
	 PdeltaFlag = 1;  
	 cerr << "PdeltaFlag: " << PdeltaFlag;
      }
      else if (strcmp(argv[argi],"-iter") == 0) 
      {
         // allow user to specify maximum number of local iterations
         argi++;
         if (argi == argc || Tcl_GetInt(interp, argv[argi++], &nMaxLocIters) != TCL_OK)
         {
	    interp->result = "WARNING invalid nMaxLocIters - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
            return TCL_ERROR;
         } 

         // specify local tolerance 
         if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &locToler) != TCL_OK)
         {
            interp->result = "WARNING invalid locToler - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>";
            return TCL_ERROR;
         } 
      }
      else
      {
         interp->result = "WARNING bad command  - element nonlinearBeamColumn eleTag? iNode? jNode? numIntgrPts? secTag? transfTag? <-mass massDens?> <-iter nMaxLocIters? locToler?>"; 
         cerr << "invalid: " << argv[argi] << endl;
         return TCL_ERROR;
      }
    }
   
    CrdTransf3d *theCrdTransf = theTclModelBuilder->getCrdTransf3d(transfTag);

    if (theCrdTransf == 0) 
    {
       cerr << "WARNING TclElmtBuilder - frameElement - no geometric transformation found with tag ";
       cerr << transfTag << endl;
       return TCL_ERROR;
    }
    
    SectionForceDeformation *theSection = theTclModelBuilder->getSection(secTag);

    if (theSection == 0) 
    {
       cerr << "WARNING Tcl3dFrame - frameElement - no Section found with tag";
       cerr << secTag << endl;
       return TCL_ERROR;
    }

    // create the element

    // get pointer to the sections for the whole beam

    SectionForceDeformation **sections = new SectionForceDeformation* [numIntgrPts];
    
    if (!sections)
    {
      interp->result = "WARNING TclElmtBuilder - addFrameElement - Insufficient memory to create sections";
      return TCL_ERROR; 
    }

    for (int j=0; j<numIntgrPts; j++)
       sections[j] = theSection;
    
    // construct the element
    
    Element *element;

    if (strcmp(argv[1],"nonlinearBeamColumn") == 0)
    {
//	if (PdeltaFlag)
//	    element = new LargeDispBeamColumn3d (eleTag, iNode, jNode, numIntgrPts, sections, *theCrdTransf, massDens, nMaxLocIters, locToler);
//	else
	    element = new NLBeamColumn3d (eleTag, iNode, jNode, numIntgrPts, sections, *theCrdTransf, massDens, nMaxLocIters, locToler);
	 
	delete [] sections;
    }
    else
    {
       interp->result = "WARNING TclElmtBuilder - addFrameElement - invalid elemType";
       return TCL_ERROR;
    }
     
    if (element == 0)
    {
      interp->result = "WARNING  TclElmtBuilder - addFrameElement - ran out of memory to create element";
      return TCL_ERROR;
    }
   
    if (theTclModelBuilderDomain->addElement(element) == false) 
    {
      cerr << "WARNING TclElmtBuilder - addFrameElement - could not add element to domain ";
      cerr << eleTag << endl;
      return TCL_ERROR;
    }
    
  }
  
  else
  {
     cerr << "WARNING NDM = " << NDM << " and NDF = " << NDF << "is imcompatible with available frame elements";
     return TCL_ERROR;
  }      

  Tcl_Free ((char *)argv);
        
  // if get here we have sucessfully created the element and added it to the domain

  return TCL_OK;
}






// 
// to create a coordinate transformation 
//
int
TclModelBuilder_addGeomTransf(ClientData clientData, Tcl_Interp *interp,
				int argc, char **argv,
				Domain *theDomain,
				TclModelBuilder *theBuilder)
				
{
	// Make sure there is a minimum number of arguments
    if (argc < 2) {
		cerr << "WARNING insufficient number of geomTransf arguments\n";
		cerr << "Want: geomTransf type? tag? <specific transf args>" << endl;
		return TCL_ERROR;
    }   
	
	theTclModelBuilderDomain = theDomain;
   theTclModelBuilder = theBuilder;
    
   int NDM, NDF;
     
   NDM = theTclModelBuilder->getNDM();   // dimension of the structure (1d, 2d, or 3d)
   NDF = theTclModelBuilder->getNDF();   // number of degrees of freedom per node

   // create 2d coordinate transformation
   if (NDM == 2 && NDF == 3)     
   {
      if ((strcmp(argv[1],"Linear") == 0) || (strcmp(argv[1],"LinearWithPDelta") == 0) || (strcmp(argv[1],"Corotational") == 0))
      {
	 int crdTransfTag;
         Vector jntOffsetI(2), jntOffsetJ(2);
	 
	 if (argc < 3) 
	 {
       	    interp->result = "WARNING insufficient arguments - want: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>"; 
	    return TCL_ERROR;
	 }
	    
         int argi = 2;  
         if (Tcl_GetInt(interp, argv[argi++], &crdTransfTag) != TCL_OK)
	 {	
	    interp->result = "WARNING invalid tag - want: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>";
	    return  TCL_ERROR;
	 }

	 // allow additional options at end of command
	 int i;

	 while (argi != argc) 
         {
	    if (strcmp(argv[argi],"-jntOffset") == 0) 
            {
	       argi++;
               for (i = 0; i < 2; i++)
	       {
                  if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &jntOffsetI(i)) != TCL_OK) 
                  {
                     interp->result = "WARNING invalid jntOffset value - want: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>";
		     return TCL_ERROR;
                  }
	       }
 
	       for (i = 0; i < 2; i++)
               {
	          if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &jntOffsetJ(i)) != TCL_OK) 
		  {
                     interp->result = "WARNING invalid jntOffset value - want: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>";
		     return TCL_ERROR;
		  }
               }
	    }
	 
	    else
	    {
               interp->result = "WARNING bad command - want: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>";
               cerr << "invalid: " << argv[argi] << endl;
               return TCL_ERROR;
            }
         }

	 // construct the transformation object
    
         CrdTransf2d *crdTransf2d;

	 if (strcmp(argv[1],"Linear") == 0)
     	    crdTransf2d = new LinearCrdTransf2d(crdTransfTag, jntOffsetI, jntOffsetJ, 0);

	 else if (strcmp(argv[1],"LinearWithPDelta") == 0)
     	    crdTransf2d = new LinearCrdTransf2d(crdTransfTag, jntOffsetI, jntOffsetJ, 1);

	 //else if (strcmp(argv[1],"Corotational") == 0)
     //	    crdTransf2d = new CorotCrdTransf2d(crdTransfTag, jntOffsetI, jntOffsetJ);
	 else
         {
            interp->result = "WARNING TclElmtBuilder - addGeomTransf - invalid Type";
	    return TCL_ERROR;
	 }
     
	 if (crdTransf2d == 0)
	 {
            interp->result = "WARNING TclElmtBuilder - addGeomTransf - ran out of memory to create geometric transformation object";
	    return TCL_ERROR;
	 }

	 // add the transformation to the modelBuilder
	 if (theTclModelBuilder->addCrdTransf2d(*crdTransf2d)) 
	 {
             interp->result = "WARNING TclElmtBuilder - addGeomTransf  - could not add geometric transformation to model Builder";
             return TCL_ERROR;
         }
      }
      else
      {
         interp->result = "WARNING TclElmtBuilder - addGeomTransf - invalid geomTransf type";
         return TCL_ERROR;
      }
   }
   else if  (NDM == 3 && NDF == 6)
   {
      if ((strcmp(argv[1],"Linear") == 0) || (strcmp(argv[1],"LinearWithPDelta") == 0) || (strcmp(argv[1],"Corotational") == 0))
      { 
	 int crdTransfTag;
         Vector vecxzPlane(3);                  // vector that defines local xz plane
         Vector jntOffsetI(3), jntOffsetJ(3);   // joint offsets in global coordinates

	 if (argc < 6) 
	 {
       	    interp->result = "WARNING insufficient arguments - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";
	    return TCL_ERROR;
	 }
	    
         int argi = 2;  
         if (Tcl_GetInt(interp, argv[argi++], &crdTransfTag) != TCL_OK)
	 {	
	    interp->result = "WARNING invalid tag - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";
	    return  TCL_ERROR;
	 }

	 if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(0)) != TCL_OK)
	 {
             interp->result = "WARNING invalid vecxzPlaneX - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";
	     return TCL_ERROR;
	 }
   
         if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(1)) != TCL_OK)
	 {
             interp->result = "WARNING invalid vecxzPlaneY - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";
	     return TCL_ERROR;
	 }
  
         if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(2)) != TCL_OK)
	 {
             interp->result = "WARNING invalid vecxzPlaneZ - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";
	     return TCL_ERROR;
	 }
  
	 // allow additional options at end of command
	 int i;

	 while (argi != argc) 
         {
	    if (strcmp(argv[argi],"-jntOffset") == 0) 
            {
               argi++;
               for (i = 0; i < 3; i++)
	       {
                  if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &jntOffsetI(i)) != TCL_OK) 
                  {
                     interp->result = "WARNING invalid jntOffset value - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >"; 
		     return TCL_ERROR;
                  }
	       }
 
	       for (i = 0; i < 3; i++)
               {
	          if (argi == argc || Tcl_GetDouble(interp, argv[argi++], &jntOffsetJ(i)) != TCL_OK) 
		  {
                     interp->result = "WARNING invalid jntOffset value - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >";  
		     return TCL_ERROR;
		  }
               }
	    }
	 
	    else
	    {
               interp->result = "WARNING bad command - want: geomTransf type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >"; 
               cerr << "invalid: " << argv[argi] << endl;
               return TCL_ERROR;
            }
         }

	 // construct the transformation object
    
         CrdTransf3d *crdTransf3d;

	 if (strcmp(argv[1],"Linear") == 0)
     	    crdTransf3d = new LinearCrdTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ, 0);

	 else if (strcmp(argv[1],"LinearWithPDelta") == 0)
     	    crdTransf3d = new LinearCrdTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ, 1);

	 //else if (strcmp(argv[1],"Corotational") == 0)
     //	    crdTransf3d = new CorotCrdTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ);
	 else
         {
            interp->result = "WARNING TclElmtBuilder - addGeomTransf - invalid Type";
	    return TCL_ERROR;
	 }
     
	 if (crdTransf3d == 0)
	 {
            interp->result = "WARNING TclElmtBuilder - addGeomTransf - ran out of memory to create geometric transformation object";
	    return TCL_ERROR;
	 }

	 // add the transformation to the modelBuilder
	 if (theTclModelBuilder->addCrdTransf3d(*crdTransf3d)) 
	 {
             interp->result = "WARNING TclElmtBuilder - addGeomTransf  - could not add geometric transformation to model Builder";

             return TCL_ERROR;
         }
      }
      else 
      {
         interp->result = "WARNING TclElmtBuilder - addGeomTransf - invalid geomTransf type ";
         return TCL_ERROR;
      }
  }
  else
  {
     cerr << "WARNING NDM = " << NDM << " and NDF = " << NDF << "is imcompatible with available frame elements";
     return TCL_ERROR;
  }      

   //  Tcl_Free ((char *)argv);
        
  // if get here we have sucessfully created the element and added it to the domain

  return TCL_OK;
}




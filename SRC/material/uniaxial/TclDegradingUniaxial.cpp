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
#include <stdlib.h>
#include <Steel02.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <elementAPI.h>
#include <OPS_Globals.h>

#if defined(OPSDEF_LIBDMG) || 1
#include <fedeas/LibDmg.h>
#else
// template<typename T> class Damage1D: public T {using T::T;};
#endif


static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
  opserr << argv[i] << " ";
    opserr << endln;
} 

// UniaxialMaterial *
Damage1D<Steel02> *
TclModelBuilder_addDegradingMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv) {
  
  UniaxialMaterial *theMaterial = nullptr;
  const char *dmg_tag = nullptr;
  Damage1D<Steel02> *dmgMat = nullptr;
  int argi = 2;

  if (strcmp(argv[argc-2], "-dmg")==0){
    dmg_tag = argv[argc-1];
    opserr << "Wrapping damage from tag" << dmg_tag << "\n";
  }

  // for (int argi=0; argi<argc; argi++){
  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial tag\n";
      printCommand(argc, argv);
      return 0;
  }

  if (strcmp(argv[1],"DmgDegradingSteel02") == 0) {

    double fy, E, b;
    double R0, cR1=0.925, cR2=0.15;
    double a1, a2=1.0, a3=0.0, a4=1.0;
    

    if (Tcl_GetDouble(interp, argv[++argi], &fy) != TCL_OK) {
      opserr << "WARNING invalid fy\n";
      printCommand(argc, argv);
      return 0;  
    }
    if (Tcl_GetDouble(interp, argv[++argi], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      printCommand(argc, argv);
      return 0;  
    }
    if (Tcl_GetDouble(interp, argv[++argi], &b) != TCL_OK) {
      opserr << "WARNING invalid b\n";
      printCommand(argc, argv);
      return 0;  
    }
    if (argc > 8) {
      if (Tcl_GetDouble(interp, argv[++argi], &R0) != TCL_OK) {
        opserr << "WARNING invalid R0\n";
        printCommand(argc, argv);
        return 0;  
      }
      if (Tcl_GetDouble(interp, argv[++argi], &cR1) != TCL_OK) {
        opserr << "WARNING invalid cR1\n";
        printCommand(argc, argv);
        return 0;  
      }
      if (Tcl_GetDouble(interp, argv[++argi], &cR2) != TCL_OK) {
        opserr << "WARNING invalid cR2\n";
        printCommand(argc, argv);
        return 0;  
      }
      if (argc > 12) {
        if (Tcl_GetDouble(interp, argv[++argi], &a1) != TCL_OK) {
          opserr << "WARNING invalid a1\n";
          printCommand(argc, argv);
          return 0;  
        }
        if (Tcl_GetDouble(interp, argv[++argi], &a2) != TCL_OK) {
          opserr << "WARNING invalid a2\n";
          printCommand(argc, argv);
          return 0;  
        }
        if (Tcl_GetDouble(interp, argv[++argi], &a3) != TCL_OK) {
          opserr << "WARNING invalid a3\n";
          printCommand(argc, argv);
          return 0;  
        }
        if (Tcl_GetDouble(interp, argv[++argi], &a4) != TCL_OK) {
          opserr << "WARNING invalid a4\n";
          printCommand(argc, argv);
          return 0;  
        }
        dmgMat = new Damage1D<Steel02>(tag, fy, E, b, R0, cR1, cR2, a1, a2, a3, a4);
      } else {
        dmgMat = new Damage1D<Steel02>(tag, fy, E, b, R0, cR1, cR2);
      }
    } else {
      dmgMat = new Damage1D<Steel02>(tag, fy, E, b);
    }
   //  printf("fy=%lf, E=%lf, b=%lf, b=%lf\n", fy, E, b);
  }

  if (dmg_tag)
    dmgMat->setDamageWrapper(interp, dmg_tag);

  return dmgMat;
}

/*
void *
OPS_Steel02D(Tcl_Interp *interp)
{
  // Pointer to a uniaxial material that will be returned
  Damage1D<Steel02> *theMaterial = nullptr;

  int    iData[1];
  double dData[12];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel02 tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 6 && numData != 10 && numData != 11) {
    opserr << "Invalid #args, want: uniaxialMaterial Steel02 " << iData[0] << 
      " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 3) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid double: uniaxialMaterial Steel02 " << iData[0] << 
  " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Damage1D<Steel02>(iData[0], dData[0], dData[1], dData[2]);    

  } else if (numData == 6) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid int: uniaxialMaterial Steel02 " << iData[0] << 
  " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Damage1D<Steel02>(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);    

  } else if (numData == 10) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02 " << iData[0] << 
  " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Damage1D<Steel02>(iData[0], dData[0], dData[1], dData[2], 
            dData[3], dData[4], dData[5], dData[6], 
            dData[7], dData[8], dData[9]);    

  } else if (numData == 11) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial Steel02 " << iData[0] << 
  " fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new Damage1D<Steel02>(iData[0], dData[0], dData[1], dData[2], 
            dData[3], dData[4], dData[5], dData[6], 
            dData[7], dData[8], dData[9], dData[10]);    

  }   

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel02 Material\n";
    return 0;
  }

  return theMaterial;
}
*/

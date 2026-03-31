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
// $Date: 2010-02-04 20:12:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SAWSMaterial.cpp,v $

// Written: Patxi (Converted from FORTRAN code originally written by Bryan Folz)
// Created: June 2006
//
// Description: This file contains the class definition for 
// SAWSMaterial.  SAWSMaterial provides the implementation
// of a one-dimensional hysteretic model develeped as part of 
// the CUREe Caltech wood frame project.
// Reference: Folz, B. and Filiatrault, A. (2001). "SAWS - Version 1.0, 
// A Computer Program for the Seismic Analysis of Woodframe Structures", 
// Structural Systems Research Project Report No. SSRP-2001/09, 
// Dept. of Structural Engineering, UCSD, La Jolla, CA .
#include <stdlib.h>
#include <SAWSMaterial.h>
#include <OPS_Globals.h>
#include <math.h>
#include <float.h>
#include <Channel.h>
#include <classTags.h>

#include <elementAPI.h>
#define OPS_Export 

static int numSAWSMaterials = 0;

OPS_Export void *
OPS_SAWSMaterial(void)
{
  if (numSAWSMaterials == 0) {
    numSAWSMaterials++;
    //OPS_Error("SAWSMaterial unaxial material - Written by Paxti Uriz, Exponent 2009\n", 1);
    opserr << "SAWSMaterial unaxial material - Written by Paxti Uriz, Exponent 2009\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[10];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial SAWSMaterial tag" << endln;
    return 0;
  }

  numData = 10;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial SAWS tag? F0? FI? dU? S0?" << endln;
    opserr << "    R1? R2? R3? R4? alpha? beta?" << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  theMaterial = new SAWSMaterial(iData[0], 
				 dData[0], dData[1], dData[2],
				 dData[3], dData[4], dData[5],
				 dData[6], dData[7], dData[8],
				 dData[9]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type SAWSMaterial\n";
    return 0;
  }

  return theMaterial;
}

SAWSMaterial::SAWSMaterial(int tag,
			   double f0, double fI, double dU, double s0,
			   double r1, double r2, double r3, double r4,
			   double A, double B):
UniaxialMaterial(tag, MAT_TAG_SAWSMaterial), 
F0(f0), FI(fI), DU(dU), S0(s0), R1(r1), R2(r2), R3(r3), R4(r4),
ALPHA(A), BETA(B)
{
  
  //TOL    = 1.0E-6 ;  //  Tolerance for bi-section search.

  // Initialize history variables
  this->revertToStart();
  this->revertToLastCommit();

}


SAWSMaterial::SAWSMaterial()
  :UniaxialMaterial(0, MAT_TAG_SAWSMaterial), 
   F0(0.0), FI(0.0), DU(0.0), S0(0.0), R1(0.0), R2(0.0), R3(0.0), R4(0.0),
   ALPHA(0.0), BETA(0.0)
{
  //TOL    = 1.0E-6 ;  //  Tolerance for bi-section search.

  // Initialize history variables
  this->revertToStart();
  this->revertToLastCommit();
}

SAWSMaterial::~SAWSMaterial()
{
  // Nothing to do
}

int
SAWSMaterial::setTrialStrain(double strain, double strainRate)
{

  // Set the initial parameters from the committed stage
  DISPL  = strain;
  LPATH  = cLPATH;
  LPPREV = cLPPREV;
  IYPLUS = cIYPLUS;
  IYMINS = cIYMINS;
  DOLD   = cDOLD;
  DUNP   = cDUNP;
  FUNP   = cFUNP;
  DUNM   = cDUNM;
  FUNM   = cFUNM;
  DMAXP  = cDMAXP;
  FMAXP  = cFMAXP;
  DMAXM  = cDMAXM;
  FMAXM  = cFMAXM;
  SP     = cSP;

  
  int foundStateFlag = 0 ;
  int checkThisStateFlag = 0;

  // The following has been converted from SAWS (originally written in 
  //     FORTRAN by Bryan Folz), converted by Patxi Uriz 6/12/2006
  /* -----------------------------------------------------------------------
   *      SUBROUTINE HYSTR CALCULATES THE FORCE-DISPLACEMENT RESPONSE
   *      BASED ON A MODIFIED STEWART HYSTERETIC MODEL. IN THIS MODEL THE
   *      ENVELOPE CURVE IS THE FOSCHI EXPONENTIAL CURVE WITH A LINEAR
   *      SOFTENING BRANCH. FORCE AND STIFFNESS ARE RETURNED FOR AN
   *      INPUTTED DISPLACEMENT.
   *-----------------------------------------------------------------------
   */
  
  /*  SUBROUTINE HYSTR (DISPL, FORCE, STIFF, JFLAG)
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      COMMON /H01/ LPATH, LPPREV, IMODE, IYPLUS, IYMINS
      COMMON /H02/ F0, FI, DU, S0, R1, R2, R3, R4, ALPHA, BETA
      COMMON /H03/ DOLD, DUNP, FUNP, DUNM, FUNM, DMAXP, FMAXP, DMAXM,
      FMAXM, SP
  */
  
  int IFLAG = 0    ;   //  Looping flag to capture infinite loop.
  int JFLAG = 0    ;
  int IMODE = 3    ;   //  Full hysteretic response permitted.

  // Determine ultimate force FU, corresponding to DU.

  FAC1 = F0 + (R1 * S0 * DU);
  FAC2 = 1.0 - exp(-S0*DU/F0);
  FU = FAC1 * FAC2;

  /* Check if DISPL > DF.  If so, there is no capacity left in the
     connector: set FORCE = 0, STIFF = 0 and LPATH = 0.
     
     Corrected DF1 failure displacement criterion
     McDonald 9/20/05
     
     DF1 = (FU - FI - (R2*S0*DU)) / (S0*(R4 - R2))
  */
  
  DF1 = (FU + FI - (R2*S0*DU)) / (S0*(R4 - R2));
  DF = DF1;
  DF2 = DU - (FU/(R2*S0));
  if (DF2 < DF1 ) 
    DF = DF2 ;
  if ( fabs(DISPL) >= DF || LPATH == 0) {
    // Reduce the force and stiffness to a small number
    FORCE = 1.0e-8*DISPL;
    STIFF = 1.0e-8;
    LPATH = 0 ;
    opserr << "Strain too large" << endln;
    return 0;
  } 
  
  /*  For the static analysis option (IMODE = 1 or 2) the load
      displacement response remains on the envelope curve.
      The option of IMODE set to 1 or 2 is used in the CASHEW program.
  */

  if (IMODE ==  2) {
    if (DISPL >= 0.0) {
      if  (DISPL <= DU) {
	FAC1 = F0 + (R1 * S0 * fabs(DISPL));
	FAC2 = 1.0 - exp(-S0*fabs(DISPL)/F0);
	FAC3 = 1 - FAC2;
	FORCE = FAC1 * FAC2;
	STIFF = (FAC1 * (S0/F0) * FAC3) + (R1 * S0 * FAC2);
	return 0;
      } else if  (DISPL > DU) {
	D0 = DU - (FU/(R2*S0));
	if (DISPL >= D0) {
	  // Reduce the force and stiffness to a small number
	  FORCE = 1.0e-8*DISPL;
	  STIFF = 1.0e-8;
	  return 0;
	}
	FORCE = FU + (R2*S0*(DISPL - DU));
	STIFF = R2 * S0;
	return 0;
      }

    } else if (DISPL <  0.0) {

      if (DISPL >= -DU) {
	FAC1 = F0 + (R1 * S0 * fabs(DISPL));
	FAC2 = 1.0 - exp(-S0*fabs(DISPL)/F0);
	FAC3 = 1 - FAC2;
	FORCE = -FAC1 * FAC2;
	STIFF = (FAC1 * (S0/F0) * FAC3) + (R1 * S0 * FAC2);
	return 0;
      } else if (DISPL < -DU) {
	D0 = (FU/(R2*S0)) - DU;
	if (DISPL <= D0) {
	  // Reduce the force and stiffness to a small number
	  FORCE = 1.0e-8*DISPL;
	  STIFF = 1.0e-8;
	  return 0;
	}
	FORCE = -FU + (R2*S0*(DISPL + DU));
	STIFF = R2 * S0;
	return 0;
      }
    }
  }

  // Find DINT2 (using bi-section search).

  DA = 0.0;
  DB = 2.0 * DU;

  /* 
     NEED TO DEAL WITH THIS PESKY 'GO TO' BUSINESS,
     TRY WHILE LOOP, INSTEAD
     Edited by Patxi 6/12/2006
  */

  const double TOL = 1e-6;
  DIFF = TOL+1.0;

  while ( fabs(DIFF) > TOL ) {
    //10 DINT2 = (DA + DB) / 2.0D0;
    DINT2 = (DA + DB) / 2.0;
    FAC1 = F0 + (R1 * S0 * fabs(DINT2));
    FAC2 = 1.0 - exp(-S0*fabs(DINT2)/F0);
    FAC3 = FI + (R4 * S0 * DINT2);
    DIFF = (FAC1 * FAC2) - FAC3;
    
    if (DIFF >= 0.0) {
      DB = DINT2;
    } else { 
      DA = DINT2;
    }
    
    // if (DABS(DIFF) .GT.  TOL)  GO TO 10;
  }

  
  DINT2 = fabs(DINT2);
  
  if (LPATH == 5) 
    DINT2 = - DINT2;
  
  DY = F0 / S0;
  
  // -------------------------------------------
  // NEED TO DEAL WITH THIS PESKY 'GO TO' BUSINESS,
  // TRY WHILE LOOP, INSTEAD
  // Edited by Patxi 6/12/2006
  // 100 CONTINUE;
  while (foundStateFlag == 0 ) {
    
    // if we are looping through, and looking for the right state, 
    // we need to stop at all LPATH checks
    checkThisStateFlag = 0;
    
    IFLAG = IFLAG + 1;
    
    if (IFLAG > 10) {
      JFLAG = 1;
      foundStateFlag = 1;
      return 0;
    }
    
    //================================================================
    /*
     *     LPATH = 1: Load-displacement response is on the non-linear
     *                (exponential) segment of the backbone curve.
     *                Unloading occurs along LPATH 1.
     */
    if (LPATH == 1  && foundStateFlag == 0 && checkThisStateFlag == 0) {
      
      DLIM = 1.05 * fabs(DINT2);
      
      if ( (DISPL >= 0.0) && (DISPL <= DLIM)) {
	
	FAC1 = F0 + (R1 * S0 * fabs(DISPL));
	FAC2 = 1.0 - exp(-S0*fabs(DISPL)/F0);
	FAC3 = 1 - FAC2;
	FORCE = FAC1 * FAC2;
	STIFF = (FAC1 * (S0/F0) * FAC3) + (R1 * S0 * FAC2);
	DOLD = DISPL;
	foundStateFlag = 1;
	return 0;
	
      } else if  ( (DISPL < 0.0) && (DISPL >= (-DLIM)) ) {
	
	FAC1 = F0 + (R1 * S0 * fabs(DISPL));
	FAC2 = 1.0 - exp(-S0*fabs(DISPL)/F0);
	FAC3 = 1 - FAC2;
	FORCE = -FAC1 * FAC2;
	STIFF = (FAC1 * (S0/F0) * FAC3) + (R1 * S0 * FAC2);
	DOLD = DISPL;
	foundStateFlag = 1;
	return 0;
	
      } else {
	
	LPATH = 2;
	LPPREV = 1;
	
      }
      
    }
    
    //====================================================================
    /*
     *     LPATH = 2: Load-displacement response is on the non-linear
     *                (exponential) segment of the backbone curve.
     */
    if ( LPATH == 2 && foundStateFlag == 0 && checkThisStateFlag == 0) {
      
      if (fabs(DISPL) <= DU) {
	
	if (fabs(DISPL) >= fabs(DOLD)) {
	  
	  FAC1 = F0 + (R1 * S0 * fabs(DISPL));
	  FAC2 = 1.0 - exp(-S0*fabs(DISPL)/F0);
	  FAC3 = 1 - FAC2;
	  
	  if (DISPL >= 0.0) {
	    
	    FORCE = FAC1 * FAC2;
	    IYPLUS = 1;
	    DUNP = DISPL;
	    FUNP = FORCE;
	    DMAXP = BETA * DUNP;
	    FAC1 = F0 + (R1 * S0 * DMAXP);
	    FAC2 = 1.0 - exp(-S0*DMAXP/F0);
	    FMAXP = FAC1 * FAC2;
	    
	    if (FMAXP > FU) {
	      
	      FMAXP = FU;
	      
	    }
	    
	  } else {
	    
	    FORCE = -FAC1 * FAC2;
	    IYMINS = 1;
	    DUNM = DISPL;
	    FUNM = FORCE;
	    DMAXM = BETA * DUNM;
	    FAC1 = F0 + (R1 * S0 * fabs(DMAXM));
	    FAC2 = 1.0 - exp(-S0*fabs(DMAXM)/F0);
	    FMAXM = -FAC1 * FAC2;
	    
	    if (FMAXM < (-FU)) {
	      
	      FMAXM = -FU;
	      
	    }

	  }
	  
	  STIFF = (FAC1 * (S0/F0) * FAC3) + (R1 * S0 * FAC2);
	  DOLD = DISPL;
	  LPPREV = 2;
	  foundStateFlag = 1;
	  return 0;

	} else {
	  
	  LPATH = 4 ;  //Unloading off LPATH 2.
	  
	}
	
      } else {

	LPATH = 3;     //Loading continues on LPATH 3.

      }

    }
    
    
    // ==================================================================
    /*
     *     LPATH = 3: Load-displacement response is on the second segment
     *                of the backbone curve
     *                (hardening or softening response).
     */
    if (LPATH == 3 && foundStateFlag == 0 && checkThisStateFlag == 0 ) {
      
      /* Corrected D6 - displacement where hysteresis slope intersects R2*S0
	 McDonald 9/20/05
	 
	 D6 = ((FU - FI) - (R2*S0*DU)) / (S0*(R4-R2))
      */
      
      D6 = ((FU + FI) - (R2*S0*DU)) / (S0*(R4-R2));
      
      if (DISPL > D6) {
	
	if (fabs(DISPL) >= fabs(DOLD)) {
	  
	  FORCE = FU + (R2*S0*(DISPL - DU));
	  STIFF = R2 * S0;
	  DOLD = DISPL;
	  DUNP = DISPL;
	  FUNP = FORCE;
	  DMAXP = BETA * DUNP;
	  FMAXP = FU + (R2*S0*(DMAXP - DU));
	  LPATH  = 3;
	  LPPREV = 3;
	  foundStateFlag = 1;
	  return 0;
	  
	} else {
	  
	  LPATH = 4;
	  LPPREV = 3;
	  checkThisStateFlag = 1;
	  // GO TO 100;
	  
	}
	
      } else if  (DISPL < -D6) {
	
	if (fabs(DISPL) >= fabs(DOLD)) {
	  
	  FORCE = -FU + (R2*S0*(DISPL + DU));
	  STIFF = R2 * S0;
	  DOLD = DISPL;
	  DUNM = DISPL;
	  FUNM = FORCE;
	  DMAXM = BETA * DUNM;
	  FMAXM = -FU + (R2*S0*(DMAXM + DU));
	  LPATH  = 3;
	  LPPREV = 3;
	  foundStateFlag = 1;
	  return 0;
	  
	} else {
	  
	  LPATH = 4;
	  LPPREV = 3;
	  checkThisStateFlag = 1;
	  //GO TO 100;
	  
	}
	
      }
      
      if (fabs(DISPL) >= fabs(DOLD) && checkThisStateFlag == 0) {
	
	if (DISPL > 0.0) {
	  
	  D0 = DU - (FU/(R2*S0));
	  
	  if (DISPL >= D0) {
	    
	    LPATH = 3;
	    LPPREV = 3;
	    // Reduce the force and stiffness to a small number
	    FORCE = 1.0e-8*DISPL;
	    STIFF = 1.0e-8;
	    foundStateFlag = 1;
	    return 0;
	    
	  }
	  
	  IYPLUS = 1;
	  FORCE = FU + (R2*S0*(DISPL - DU));
	  STIFF = R2 * S0;
	  DOLD = DISPL;
	  DUNP = DISPL;
	  FUNP = FORCE;
	  DMAXP = BETA * DUNP;
	  FMAXP = FU + (R2*S0*(DMAXP - DU));
	  LPPREV = 3;
	  foundStateFlag = 1;
	  return 0;
	  
	} else if  (DISPL < 0.0 && checkThisStateFlag == 0) {
	  
	  D0 = (FU/(R2*S0)) - DU;
	  
	  if (DISPL <= D0) {
	    
	    LPATH = 3;
	    // Reduce the force and stiffness to a small number
	    FORCE = 1.0e-8*DISPL;
	    STIFF = 1.0e-8;
	    foundStateFlag = 1;
	    return 0;
	    
	  }
	  
	  IYMINS = 1;
	  FORCE = -FU + (R2*S0*(DISPL + DU));
	  STIFF = R2 * S0;
	  DOLD = DISPL;
	  DUNM = DISPL;
	  FUNM = FORCE;
	  DMAXM = BETA * DUNM;
	  FMAXM = -FU + (R2*S0*(DMAXM + DU));
	  LPPREV = 3;
	  foundStateFlag = 1;
	  return 0;
	  
	}
	
      } else if ( checkThisStateFlag == 0) {
	
	LPATH = 4;  // Unloading off LPATH 4
	
      }
      
    }

  
    // ======================================================
    /*
     * LPATH = 4: Load-displacement response is unloading off
     *            the backbone curve (LPATH 2 or 3).
    */
    if (LPATH == 4 && foundStateFlag == 0 && checkThisStateFlag == 0) {

      if ((DOLD < 0.0) && (DISPL >= 0.0)) {

	LPPREV = 4;
	LPATH = 13;
	checkThisStateFlag = 1;
	//GO TO 100;

      }
      
      if ((DOLD > 0.0) && (DISPL <= 0.0) && checkThisStateFlag == 0) {

	LPPREV = 4;
	LPATH = 14;
	checkThisStateFlag = 1;
	//GO TO 100

      }

      if (DISPL >= 0.0 && checkThisStateFlag == 0) {

	dZero = DUNP - (FUNP/(R3*S0));
	DINT1P = (FI - (R3*S0*dZero))/(S0*(R4 - R3));

	if (DISPL >= DINT1P) {

	  FORCE = R3*S0*(DISPL - dZero);

	  if (FORCE > FUNP) {

	    LPPREV = 4;
	    LPATH = 2;
	    checkThisStateFlag = 1;
	    //GO TO 100

	  }

	  STIFF = R3 * S0;
	  DOLD = DISPL;
	  foundStateFlag = 1;
	  return 0;
          
	} else {

	  LPATH = 5;

	}
	
      } else if ( checkThisStateFlag == 0) {

	dZero = DUNM - (FUNM/(R3*S0));
	DINT1M = (-FI - (R3*S0*dZero))/(S0*(R4 - R3));

	if (DISPL <= DINT1M) {

	  FORCE = R3*S0*(DISPL - dZero);

	  if (FORCE < FUNM) {

	    LPPREV = 4;
	    LPATH = 2;
	    checkThisStateFlag = 1;
	    //GO TO 100

	  }

	  STIFF = R3 * S0;
	  DOLD = DISPL;
	  foundStateFlag = 1;
	  return 0;

	} else {

	  LPATH = 7;

	}
      }
    }

    //====================================================
    /*
     * LPATH = 5: Load-displacement response is unloading or
     *            reloading on the pinched curve
     *           (for negative forces).
    */
    if (LPATH == 5 && foundStateFlag == 0 && checkThisStateFlag == 0) {

      DINT2 = -fabs(DINT2);
      //c Corrected McDonald 9/20/05
      //c        DINT4 = (-FU + FI + (R2*S0*DU))/(S0*(R4 - R2))
      DINT4 = (-FU - FI + (R2*S0*DU))/(S0*(R4 - R2));

      if (DISPL <= DINT4) {

	LPATH = 3;
	LPPREV = 5;
	checkThisStateFlag = 1;
	//GOTO 100

      }

      if ((LPPREV == 5) && (DISPL > DOLD)  && checkThisStateFlag == 0 ) {

	LPATH = 9;

      } else if ( checkThisStateFlag == 0) {

	if (IYMINS == 1) {

	  if (DMAXM != 0.0) 
	    SP = S0 * pow((DY/fabs(DMAXM)),ALPHA);

	} else {

	  SP = S0;

	}

	DINT3 = (-FI - FMAXM + (SP*DMAXM))/(SP - (R4*S0));
	
	if (DISPL >= DINT2) {

	  FORCE = -FI + (R4*S0*DISPL);
	  STIFF = R4 * S0;
	  DOLD = DISPL;
	  LPPREV = 5;
	  foundStateFlag = 1;
	  return 0;

	} else if  ((DISPL < DINT2) && (IYMINS == 0)) {

	  LPPREV = 5;
	  LPATH = 1;
	  checkThisStateFlag = 1;
	  //GOTO 100

	} else if  ((DISPL < DINT2) &&  (DISPL >= DINT3)) {

	  FORCE = -FI + (R4*S0*DISPL);
	  STIFF = R4 * S0;
	  DOLD = DISPL;
	  LPPREV = 5;
	  foundStateFlag = 1;
	  return 0;

	} else {

	  LPATH = 6;

	}

      }

    }

    //==================================================
    /*
     *     LPATH = 6: Load-displacement response is reloading from the
     *                pinched curve (LPATH 5) to the backbone curve
     *                (LPATH 2 or 3).
    */
    
    if (LPATH == 6 && foundStateFlag == 0 && checkThisStateFlag == 0) {
      
      if ((LPPREV == 6) && (DISPL > DOLD)) {

	LPATH = 11;

      } else {
	
	if (DISPL >= DMAXM) {

	  FORCE = FMAXM + (SP*(DISPL - DMAXM));
	  STIFF = SP;
	  DOLD = DISPL;
	  LPPREV = 6;
	  foundStateFlag = 0;
	  return 0;

	} else {

	  LPATH = 2;
	  checkThisStateFlag = 1;
	  //GOTO 100

	}

      }

    }

    //========================================================
    /*
     *     LPATH = 7: Load-displacement response is unloading or
     *                reloading on the pinched curve
     *                (for positive forces).
    */

    if (LPATH == 7 && foundStateFlag == 0 && checkThisStateFlag == 0) {

      //c Corrected McDonald 9/20/05
      //c        DINT4 = (FU - FI - (R2*S0*DU))/(S0*(R4 - R2))
      DINT4 = (FU + FI - (R2*S0*DU))/(S0*(R4 - R2));

      if (DISPL >= DINT4) {
	
	LPPREV = 7;
	LPATH = 3;
	checkThisStateFlag = 1;
	//GOTO 100
        
      }

      if ((LPPREV == 7) && (DISPL < DOLD) && checkThisStateFlag == 0 ) {

	LPATH = 10;

      } else if ( checkThisStateFlag == 0 ) {

	if (IYPLUS == 1) {

	  if (DMAXP != 0.0) 
	    SP = S0 * pow((DY/DMAXP),ALPHA);

	} else {

	  SP = S0 ;

	}

	DINT3 = (FI - FMAXP + (SP*DMAXP))/(SP - (R4*S0));

	if (DISPL <= DINT2) {

	  FORCE = FI + (R4*S0*DISPL);
	  STIFF = R4 * S0;
	  DOLD = DISPL;
	  LPPREV = 7;
	  foundStateFlag = 1;
	  return 0;

	} else if  ((DISPL > DINT2) &&  (IYPLUS == 0)) {

	  LPPREV = 7;
	  LPATH = 1;
	  checkThisStateFlag = 1;
	  //GOTO 100

	} else if ((DISPL > DINT2) && (DISPL <= DINT3)) {

	  FORCE = FI + (R4*S0*DISPL);
	  STIFF = R4 * S0;
	  DOLD = DISPL;
	  LPPREV = 7;
	  foundStateFlag = 1;
	  return 0;
	  
	} else {

	  LPPREV = 7;
	  LPATH = 8;

	}

      }

    }

    //===============================================================
    /*
     *     LPATH = 8: Load-displacement response is reloading from the
     *                pinched curve (LPATH 7) to the backbone curve
     *                (LPATH 2 or 3).
     */

    if (LPATH == 8 && foundStateFlag == 0 && checkThisStateFlag == 0) {

      if ((LPPREV == 8) && (DISPL < DOLD)) {

	LPATH = 12;

        } else {

          if  (DISPL <= DMAXP) {

            FORCE = FMAXP + (SP*(DISPL - DMAXP));
            STIFF = SP;
            DOLD = DISPL;
            LPPREV = 8;
	    foundStateFlag = 1;
            return 0;
	    
	  } else {

            LPATH = 2;
	    checkThisStateFlag = 1;
            // GOTO 100

	    }
	  }
    }
     
    //====================================================
    /*
     *     LPATH = 9: Load-dispalcement response is between the two
     *                pinched curves (LPATH 5 and 7).
    */
    if (LPATH == 9 && foundStateFlag == 0 && checkThisStateFlag == 0) {

      FOLD = -FI + (R4*S0*DOLD);
      DUPPER = (FOLD - FI - (R3*S0*DOLD))/(S0*(R4 - R3));
      DLOWER = DOLD;
      
      if (DISPL <= DLOWER) {

	LPATH = 5;
	LPPREV = 9;
	checkThisStateFlag = 1;
	//GO TO 100

      } else if  ((DISPL > DLOWER) && (DISPL < DUPPER)) {

	FORCE = FOLD + (R3*S0*(DISPL - DLOWER));
	STIFF = R3 * S0;
	foundStateFlag = 1;
	return 0;

        } else {

          LPATH = 7;
	  LPPREV = 9;
	  checkThisStateFlag = 1;
          //GOTO 100;

        }

    }

    //==============================================================
    /*
     *     LPATH = 10: Load-displacement response is between the
     *                two pinched curves (LPATH 7 and 5).
     */

    if (LPATH == 10 && foundStateFlag == 0 && checkThisStateFlag == 0) {

      FOLD = FI + (R4*S0*DOLD);
      DUPPER = (FOLD - FI - (R3*S0*DOLD))/(S0*(R4 - R3));
      DLOWER = (FOLD + FI - (R3*S0*DOLD))/(S0*(R4 - R3));

      if (DISPL <= DLOWER) {

	LPATH = 5;
	LPPREV = 10;
	checkThisStateFlag = 1;
	//GO TO 100

      } else if  ((DISPL > DLOWER) && (DISPL < DUPPER)) {

	FORCE = FOLD + (R3*S0*(DISPL - DUPPER));
	STIFF = R3 * S0;
	foundStateFlag = 1;
	return 0;

      } else {

	LPATH = 7;
	LPPREV = 10;
	checkThisStateFlag = 1;
	//GOTO 100
        
      }

    }

    //==============================================================
    /*
     *     LPATH = 11: Load-displacement response is unloading off the
     *                 reloading curve (LPATH 6).
     */

    if (LPATH == 11 && foundStateFlag == 0 && checkThisStateFlag == 0) {

      SP = S0 * pow((DY/fabs(DMAXM)),ALPHA);
      DINT3 = (-FI - FMAXM + (SP*DMAXM))/(SP - (R4*S0));

      if ((LPPREV == 4) && (DISPL <= DINT3)) {

	LPATH = 6;
	LPPREV = 11;
	checkThisStateFlag = 1;
	//GO TO 100

      } else {

	FOLD = FMAXM + (SP*(DOLD - DMAXM));
	DLOWER = DOLD;
	DUPPER = (FOLD - FI - (R3*S0*DLOWER))/(S0*(R4 - R3));

	if (DISPL <= DLOWER) {

	  if (DISPL <= DMAXM) {

	    LPATH = 2;

	  } else {

	    LPATH = 6;
            
	  }
           
	  LPPREV = 11;
	  checkThisStateFlag = 1;
	  //GO TO 100

	} else if  ((DISPL > DLOWER) && (DISPL < DUPPER)) {

	  FORCE = FOLD + (R3*S0*(DISPL - DLOWER));
	  STIFF = R3 * S0;
	  foundStateFlag = 1;
	  return 0;
	  
          } else {

            LPATH = 7;
            LPPREV = 11;
	    checkThisStateFlag = 1;
            //GOTO 100

	  }
      }
    }

    //============================================================
    /*
     *     LPATH = 12: Load-displacement response is unloading off the
     *                 reloading curve (LPATH 8).
     */
    if (LPATH == 12 && foundStateFlag == 0 && checkThisStateFlag == 0) {

      SP = S0 * pow((DY/DMAXP),ALPHA);
      DINT3 = (FI - FMAXP + (SP*DMAXP))/(SP - (R4*S0));

      if ((LPPREV == 7) && (DISPL >= DINT3)) {

	LPATH = 8;
	LPPREV = 12;
	checkThisStateFlag = 1;
	//GO TO 100

      } else {

	FOLD = FMAXP + (SP*(DOLD - DMAXP));
	DUPPER = DOLD;
	DLOWER = (FOLD + FI - (R3*S0*DUPPER))/(S0*(R4 - R3));

	if (DISPL <= DLOWER) {
	  LPATH = 5;
	  LPPREV = 12;
	  checkThisStateFlag = 1;
	  //GO TO 100

	} else if ((DISPL > DLOWER) && (DISPL < DUPPER)) {

	  FORCE = FOLD + (R3*S0*(DISPL - DUPPER));
	  STIFF = R3 * S0;
	  foundStateFlag = 1;
	  return 0;

          } else { 

            LPATH = 8;
            LPPREV = 12;
	    checkThisStateFlag = 1;
            //GOTO 100

	  }

      }

    }

    //============================================================
    /*
     *     LPATH = 13: Load-displacement response following LPATH 3 when
     *                 displacement goes from negative to positive.
    */
    if (LPATH == 13 && foundStateFlag == 0 && checkThisStateFlag == 0) {

      dZero = DUNM - (FUNM/(R3*S0));
      DINT5= (-FI - (R3*S0*dZero))/(S0*(R4 - R3));

      if (DISPL < DINT5) {

	LPPREV = 13;
	FORCE = R3*S0*(DISPL - dZero);

	if (FORCE < FMAXM) {

	  LPATH = 1;
	  checkThisStateFlag = 1;
	  //GO TO 100

	}

	STIFF = R3 * S0;
	DOLD = DISPL;
	foundStateFlag = 1;
	return 0;

      } else {
	
	LPPREV = 13;
	LPATH = 7;
	checkThisStateFlag = 1;
	//GO TO 100

      }

    }

    //  ================================================================
    /*
     *     LPATH = 14: Load-displacement response following LPATH 4 when
     *                 displacement goes from positive to negative.
     */
    if (LPATH == 14 && foundStateFlag == 0 && checkThisStateFlag == 0) {

      dZero = DUNP - (FUNP/(R3*S0));
      DINT5= (FI - (R3*S0*dZero))/(S0*(R4 - R3));

      if (DISPL > DINT5) {

	LPPREV = 14;
	FORCE = R3*S0*(DISPL - dZero);

	if (FORCE > FMAXP) {

	  LPATH = 1;
	  checkThisStateFlag = 1;
	  //GO TO 100
          
	}
         
	STIFF = R3 * S0;
	DOLD = DISPL;
	foundStateFlag = 1;
	return 0;

      } else {

	LPPREV = 14;
	LPATH = 5;
	checkThisStateFlag = 1;
	//GO TO 100
        
      }

    }

  }
return 0;

}


double
SAWSMaterial::getStrain(void)
{
  return DISPL;
}

double
SAWSMaterial::getStress(void)
{
  return FORCE;
}

double
SAWSMaterial::getTangent(void)
{
  return STIFF;
}


double
SAWSMaterial::getInitialTangent(void)
{
  return S0;
}

int
SAWSMaterial::commitState(void)
{
  cDISPL  = DISPL;
  cFORCE  = FORCE;
  cSTIFF  = STIFF;
  cLPATH  = LPATH;
  cLPPREV = LPPREV;
  cIYPLUS = IYPLUS;
  cIYMINS = IYMINS;
  cDOLD   = DOLD;
  cDUNP   = DUNP;
  cFUNP   = FUNP;
  cDUNM   = DUNM;
  cFUNM   = FUNM;
  cDMAXP  = DMAXP;
  cFMAXP  = FMAXP;
  cDMAXM  = DMAXM;
  cFMAXM  = FMAXM;
  cSP     = SP;
  return 0;
}

int
SAWSMaterial::revertToLastCommit(void)
{
  DISPL = cDISPL;
  FORCE = cFORCE;
  STIFF = cSTIFF;
  LPATH = cLPATH;
  LPPREV= cLPPREV;
  IYPLUS = cIYPLUS;
  IYMINS = cIYMINS;
  DOLD   = cDOLD;
  DUNP   = cDUNP;
  FUNP   = cFUNP;
  DUNM   = cDUNM;
  FUNM   = cFUNM;
  DMAXP  = cDMAXP;
  FMAXP  = cFMAXP;
  DMAXM  = cDMAXM;
  FMAXM  = cFMAXM;
  SP     = cSP;

  return 0;
}

int
SAWSMaterial::revertToStart(void)
{
  cDISPL  = 0.0;
  cFORCE  = 0.0;
  cSTIFF  = S0 ;
  cLPATH  = 1  ;
  cLPPREV = 1  ;
  cIYPLUS = 0  ;
  cIYMINS = 0  ;
  cDOLD   = 0.0;
  cDUNP   = 0.0;
  cFUNP   = 0.0;
  cDUNM   = 0.0;
  cFUNM   = 0.0;
  cDMAXP  = 0.0;
  cFMAXP  = 0.0;
  cDMAXM  = 0.0;
  cFMAXM  = 0.0;
  cSP     = 0.0;
  return 0;
}

UniaxialMaterial*
SAWSMaterial::getCopy(void)
{
  SAWSMaterial *theCopy = new SAWSMaterial (this->getTag(),
					    F0, FI, DU, S0, 
					    R1, R2, R3, R4,
					    ALPHA, BETA );
  
  theCopy->cDISPL  = cDISPL;
  theCopy->cFORCE  = cFORCE;
  theCopy->cSTIFF  = cSTIFF;
  theCopy->cLPATH  = cLPATH;
  theCopy->cIYPLUS = cIYPLUS;
  theCopy->cIYMINS = cIYMINS;
  theCopy->cDOLD   = cDOLD;
  theCopy->cDUNP   = cDUNP;
  theCopy->cFUNP   = cFUNP;
  theCopy->cDUNM   = cDUNM;
  theCopy->cFUNM   = cFUNM;
  theCopy->cDMAXP  = cDMAXP;
  theCopy->cFMAXP  = cFMAXP;
  theCopy->cDMAXM  = cDMAXM;
  theCopy->cFMAXM  = cFMAXM;
  theCopy->cSP     = cSP;

  return theCopy;
}

int
SAWSMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector dataVec(28);
  
  dataVec(0) = this->getTag();
  dataVec(1) = F0;
  dataVec(2) = FI;
  dataVec(3) = DU;
  dataVec(4) = S0;
  dataVec(5) = R1;
  dataVec(6) = R2;
  dataVec(7) = R3;
  dataVec(8) = R4;
  dataVec(9) = ALPHA;
  dataVec(10) = BETA;
  dataVec(11) = cDISPL;
  dataVec(12) = cSTIFF;
  dataVec(13) = cFORCE;
  dataVec(14) = cLPATH;
  dataVec(15) = cLPPREV;
  dataVec(16) = cIYPLUS;
  dataVec(17) = cIYMINS;
  dataVec(18) = cDOLD;
  dataVec(19) = cDUNP;
  dataVec(20) = cFUNP;
  dataVec(21) = cDUNM;
  dataVec(22) = cFUNM;
  dataVec(23) = cDMAXP;
  dataVec(24) = cFMAXP;
  dataVec(25) = cDMAXM;
  dataVec(26) = cFMAXM;
  dataVec(27) = cSP;



  res = theChannel.sendVector(this->getDbTag(), commitTag, dataVec);
  if (res < 0) 
    opserr << "SAWSMaterial::sendSelf() - failed to send data\n";


  return res;
}

int
SAWSMaterial::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector dataVec(28);
  res = theChannel.recvVector(this->getDbTag(), commitTag, dataVec);
  
  if (res < 0) {
      opserr << "SAWSMaterial::recvSelf() - failed to receive data\n";
      return res;
  }
  else {
    this->setTag((int)dataVec(0));
    F0      = dataVec(1);
    FI      = dataVec(2);
    DU      = dataVec(3);
    S0      = dataVec(4);
    R1      = dataVec(5);
    R2      = dataVec(6);
    R3      = dataVec(7);
    R4      = dataVec(8);
    ALPHA   = dataVec(9);
    BETA    = dataVec(10);
    cDISPL  = dataVec(11);
    cSTIFF  = dataVec(12);
    cFORCE  = dataVec(13);
    cLPATH  = int(dataVec(14));
    cLPPREV = int(dataVec(15));
    cIYPLUS = int(dataVec(16));
    cIYMINS = int(dataVec(17));
    cDOLD   = dataVec(18);
    cDUNP   = dataVec(19);
    cFUNP   = dataVec(20);
    cDUNM   = dataVec(21);
    cFUNM   = dataVec(22);
    cDMAXP  = dataVec(23);
    cFMAXP  = dataVec(24);
    cDMAXM  = dataVec(25);
    cFMAXM  = dataVec(26);
    cSP     = dataVec(27);

    // set the trial values
    FORCE  = cFORCE;
    DISPL  = cDISPL;
    LPATH  = cLPATH;
    LPPREV = cLPPREV;

  }

  return 0;
}
    
void
SAWSMaterial::Print(OPS_Stream &s, int flag)
{
  
	s << "SAWSMaterial, tag: " << this->getTag() << endln;
	s << "F0: " << F0 << endln;
	s << "FI: " << FI << endln;
	s << "DU: " << DU << endln;
	s << "S0: " << S0 << endln;
	s << "R1: " << R1 << endln;
	s << "R2: " << R2 << endln;
	s << "R3: " << R3 << endln;
	s << "R4: " << R4 << endln;
	s << "ALPHA: " << ALPHA << endln;
	s << "BETA: " << BETA << endln;

}


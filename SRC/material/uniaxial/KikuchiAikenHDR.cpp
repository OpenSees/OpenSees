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

// $Revision: 1.0 $
// $Date: 2013-05-29 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/KikuchiAikenHDR.cpp,v $

// Written: Ken Ishii
// Created: June 2012
//
// Kikuchi&Aiken model for high-damping rubber bearing
//
// Description: This file contains the class definition for KikuchiAikenHDR.
//              This file contains the function to parse the TCL input
// uniaxialMaterial KikuchiAikenHDR matTag? tp? ar? hr? <-coGHU cg? ch? cu?> <-coMSS rs? rf?>



#include <KikuchiAikenHDR.h>
#include <Vector.h>
#include <Channel.h>

#include <TclModelBuilder.h>
#include <string.h>
#include <tcl.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>



int
TclCommand_KikuchiAikenHDR(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  //arguments (necessary)
  int tag;
  int tp;
  double ar;
  double hr;

  //arguments (optional)
  double cg = 1.0;
  double ch = 1.0;
  double cu = 1.0;
  double rs = 1.0;
  double rf = 1.0;
      
  //
  UniaxialMaterial *theMaterial = 0;


  //error flag
  bool ifNoError = true;
  
  if (argc < 6) { // uniaxialMaterial KikuchiAikenHDR matTag? tp? ar? hr?

    opserr << "WARNING invalid number of arguments\n";
    ifNoError = false;

  } else {

    //argv[2~5]
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid KikuchiAikenHDR tag" << endln;
      ifNoError = false;
    }
  
    if ((strcmp(argv[3],"X0.6") == 0) || (strcmp(argv[3],"1") == 0)) {
      tp = 1;
    } else if ((strcmp(argv[3],"X0.6-0MPa") == 0) || (strcmp(argv[3],"2") == 0)) {
      tp = 2;
    } else {
      opserr << "WARNING invalid KikuchiAikenHDR tp" << endln;
      ifNoError = false;
    }
  
    if (Tcl_GetDouble(interp, argv[4], &ar) != TCL_OK || ar <= 0.0) {
      opserr << "WARNING invalid ar\n";
      ifNoError = false;
    }
    
    if (Tcl_GetDouble(interp, argv[5], &hr) != TCL_OK || hr <= 0.0) {
      opserr << "WARNING invalid hr\n";
      ifNoError = false;
    }
    
    
    //argv[6~]
    for (int i=6; i<=(argc-1); i++) {
      
      if (strcmp(argv[i],"-coGHU")==0 && (i+3)<=(argc-1)) { // <-coGHU cg? ch? cu?>

	if (Tcl_GetDouble(interp,argv[i+1], &cg) != TCL_OK || cg < 0.0) {
	  opserr << "WARNING invalid cg\n";
	  ifNoError = false;
	}
	
	if (Tcl_GetDouble(interp,argv[i+2], &ch) != TCL_OK || ch < 0.0) {
	  opserr << "WARNING invalid ch\n";
	  ifNoError = false;
	}
	
	if (Tcl_GetDouble(interp,argv[i+3], &cu) != TCL_OK || cu < 0.0) {
	  opserr << "WARNING invalid cu\n";
	  ifNoError = false;
	}
	
	i += 3;
	
      } else if (strcmp(argv[i],"-coMSS")==0 && (i+2)<=(argc-1)) { // <-coMSS rs? rf?>
	
	if (Tcl_GetDouble(interp,argv[i+1], &rs) != TCL_OK || rs < 0.0) {
	  opserr << "WARNING invalid rs\n";
	  ifNoError = false;
	}
	
	if (Tcl_GetDouble(interp,argv[i+2], &rf) != TCL_OK || rf < 0.0) {
	  opserr << "WARNING invalid rf\n";
	  ifNoError = false;
	} 

	i += 2;
	
      } else { // invalid option
	opserr << "WARNING invalid optional arguments \n";
	ifNoError = false;
	break;
      }

    }

  }

  //if error detected
  if (!ifNoError) {
    //input:
    opserr << "Input command: ";
    for (int i=0; i<argc; i++){
      opserr << argv[i] << " ";
    }
    opserr << endln;
    
    //want:
    opserr << "Want: uniaxialMaterial KikuchiAikenHDR matTag? tp? ar? hr? <-coGHU cg? ch? cu?> <-coMSS rs? rf?>" << endln;
    return TCL_ERROR;
  }

  //regard 0.0 input as mistake (substitute 1.0 for 0.0)
  if (cg == 0.0) cg = 1.0;
  if (ch == 0.0) ch = 1.0;
  if (cu == 0.0) cu = 1.0;
  if (rs == 0.0) rs = 1.0;
  if (rf == 0.0) rf = 1.0;


  // Parsing was successful, allocate the material
  theMaterial = new KikuchiAikenHDR(tag, tp, ar, hr, cg, ch, cu, rs, rf);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial " << argv[1] << endln;
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (OPS_addUniaxialMaterial(theMaterial) == false) {
    opserr << "WARNING could not add uniaxialMaterial to the modelbuilder\n";
    opserr << *theMaterial << endln;
    delete theMaterial; // invoke the material objects destructor, otherwise mem leak
    return TCL_ERROR;
  } 

  // succeeded
  return TCL_OK;

}



KikuchiAikenHDR::KikuchiAikenHDR(int tag, int tp, double ar, double hr, 
			 double cg, double ch, double cu, double rs, double rf)
  :UniaxialMaterial(tag,MAT_TAG_KikuchiAikenHDR),Tp(tp),Ar(ar),Hr(hr),
   Cg(cg),Ch(ch),Cu(cu),Rs(rs),Rf(rf)
{
  
  //parameter function for each rubber
  switch (Tp) {
    
  case 1: // Bridgestone HDR X0.6 (standard compressive stress)
    trgStrain = 0.05;
    lmtStrain = 4.10;
    calcGeq = KikuchiAikenHDR::calcGeqTp1;
    calcHeq = KikuchiAikenHDR::calcHeqTp1;
    calcU   = KikuchiAikenHDR::calcUTp1;
    calcN   = KikuchiAikenHDR::calcNTp1;
    calcA   = KikuchiAikenHDR::calcATp1;
    calcB   = KikuchiAikenHDR::calcBTp1;
    calcC   = KikuchiAikenHDR::calcCTp1;
    break;
    
  case 2: // Bridgestone HDR X0.6 (zero compressive stress)
    trgStrain = 0.05;
    lmtStrain = 4.10;
    calcGeq = KikuchiAikenHDR::calcGeqTp2;
    calcHeq = KikuchiAikenHDR::calcHeqTp2;
    calcU   = KikuchiAikenHDR::calcUTp2;
    calcN   = KikuchiAikenHDR::calcNTp2;
    calcA   = KikuchiAikenHDR::calcATp2;
    calcB   = KikuchiAikenHDR::calcBTp2;
    calcC   = KikuchiAikenHDR::calcCTp2;
    break;
  }
  

  // intitialize
  initialStiff = (this->calcGeq)(trgStrain)*Cg*Ar/Hr;


  //
  numIdx = 500;
  revXBgn   = new double [numIdx];
  revQ2Bgn  = new double [numIdx];
  revXEnd   = new double [numIdx];
  revQ2End  = new double [numIdx];
  revB      = new double [numIdx];
  revAlpha  = new double [numIdx];

  trialDeform  = 0.0;
  trialForce        = 0.0;
  trialStiff    = initialStiff;
  trialStrain = 0.0;
  trialStress = 0.0;
  trialTangent = initialStiff*Hr/Ar;
  trialIfElastic = true;
  trialQ1 = 0.0;
  trialQ2 = 0.0;
  trialMaxStrain = 0.0;
  trialDStrain = 0.0;
  trialDStrainLastSign = 0;
  trialIdxRev=0;

  commitDeform  = 0.0;
  commitForce        = 0.0;
  commitStiff    = initialStiff;
  commitStrain = 0.0;
  commitStress = 0.0;
  commitTangent = initialStiff*Hr/Ar;
  commitIfElastic = true;
  commitQ1 = 0.0;
  commitQ2 = 0.0;
  commitMaxStrain = 0.0;
  commitDStrain = 0.0;
  commitDStrainLastSign = 0;
  commitIdxRev=0;

  revB[0] = 0.0;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenHDR::KikuchiAikenHDR\n";
  // -------------------------------------------------------

}

KikuchiAikenHDR::KikuchiAikenHDR()
  :UniaxialMaterial(0,MAT_TAG_KikuchiAikenHDR)
{

  //

  trialDeform  = 0.0;
  trialForce   = 0.0;
  trialStiff   = 0.0;
  commitDeform = 0.0;
  commitForce  = 0.0;
  commitStiff  = 0.0;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenHDR::KikuchiAikenHDR()\n";
  // -------------------------------------------------------

}

KikuchiAikenHDR::~KikuchiAikenHDR()
{

  if (revXBgn != 0)
    delete [] revXBgn;

  if (revQ2Bgn != 0)
    delete [] revQ2Bgn;

  if (revXEnd != 0)
    delete [] revXEnd;

  if (revQ2End != 0)
    delete [] revQ2End;

  if (revB != 0)
    delete [] revB;

  if (revAlpha != 0)
    delete [] revAlpha;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenHDR::~KikuchiAikenHDR\n";
  // -------------------------------------------------------

}

int 
KikuchiAikenHDR::setTrialStrain(double strain, double strainRate)
{

  // (input)                             (output)
  // deformation -> strain -> stress  -> force
  //                strain -> modulus -> spring constant



  //(input: deformation->strain)
  trialDeform = strain;
  trialDeform = trialDeform * (Rs/Rf); //for MSS model
  trialStrain  = trialDeform/Hr;


  //incremental strain
  trialDStrain = trialStrain - commitStrain;

  //sign of incremental strain (used in next caluculation step)
  if (trialDStrain > 0) {
    trialDStrainLastSign = +1;
  } else if (trialDStrain < 0) {
    trialDStrainLastSign = -1;
  } else {
    trialDStrainLastSign = commitDStrainLastSign;
  }



  // if incremental strain is zero
  if ( fabs(trialDStrain) < DBL_EPSILON ) { // if ( trialDStrain == 0.0 )

    //(strain->stress, strain->modulus)
    trialStress  = commitStress;
    trialTangent = commitTangent;

    //(stress->force, modulus->spring constant)
    trialForce     = commitForce;
    trialStiff = commitStiff;

    return 0;
  }


  //application limit strain
  if ( fabs(trialStrain) > lmtStrain ) {
    opserr << "uniaxialMaterial KikuchiAikenHDR: \n";
    opserr << "   Response value exceeded limited strain.\n";
    // return -1;
    // warning, extrapolation
  }

  //elastic limit strain
  if ( fabs(trialStrain) > trgStrain ) {
      trialIfElastic = false;
  }

  //maximum strain
  if ( fabs(trialStrain) > commitMaxStrain ) {
    trialMaxStrain = fabs(trialStrain);
  }


  //parameters for Kikuchi&Aiken model

  if ( trialIfElastic || fabs(trialStrain) == trialMaxStrain ) {

    //if elastic, tmpstrain is max(fabs(trialStrain),trgStrain)
    tmpStrain = (fabs(trialStrain)>trgStrain) ? fabs(trialStrain) : trgStrain ; //max(fabs(trialStrain),trgStrain)
    
    geq = (this->calcGeq)(tmpStrain)*Cg;
    heq = (this->calcHeq)(tmpStrain)*Ch;
    u = (this->calcU)(tmpStrain)*Cu;

    n = (this->calcN)(fabs(trialStrain));
    c = (this->calcC)(fabs(trialStrain));
    a = (this->calcA)(fabs(trialStrain),heq,u);

    xm = fabs(trialStrain);
    fm = geq*xm;

  }
  
  // normalized strain
  x = trialStrain/xm;


  // reversal points
  if (!trialIfElastic) {

    // unload or reverse
    if (trialDStrain*commitDStrainLastSign < 0) {

            
      if ( trialIdxRev == 0 ) { // unload

	trialIdxRev = 1;

	revXBgn[1]  = commitStrain/xm;
	revQ2Bgn[1] = commitQ2;
	
	//b
	b = (this->calcB)(fabs(commitStrain),a,c,(this->calcHeq)(fabs(commitStrain))*Ch,(this->calcU)(fabs(commitStrain))*Cu);

	revB[1] = b;
	revAlpha[1] = 1.0;
                
            
      } else { //reverse

	trialIdxRev++;


	if (trialIdxRev >= numIdx ) { //memorize reversal points of hysteretic curve
	  int newIdx;
	  double *newXBgn;
	  double *newQ2Bgn;
	  double *newXEnd;
	  double *newQ2End;
	  double *newB;
	  double *newAlpha;
	  
	  newIdx = numIdx + 500;
	  newXBgn  = new double [newIdx];
	  newQ2Bgn = new double [newIdx];
	  newXEnd  = new double [newIdx];
	  newQ2End = new double [newIdx];
	  newB     = new double [newIdx];
	  newAlpha = new double [newIdx];
	  
	  for(int i = 0; i<numIdx; i++){
	    newXBgn[i]  = revXBgn[i];
	    newQ2Bgn[i] = revQ2Bgn[i];
	    newXEnd[i]  = revXEnd[i];
	    newQ2End[i] = revQ2End[i];
	    newB[i]     = revB[i];
	    newAlpha[i] = revAlpha[i];
	  }
	  
	  numIdx = newIdx;

	  delete [] revXBgn;
	  delete [] revQ2Bgn;
	  delete [] revXEnd;
	  delete [] revQ2End;
	  delete [] revB;
	  delete [] revAlpha;

	  revXBgn  = newXBgn;
	  revQ2Bgn = newQ2Bgn;
	  revXEnd  = newXEnd;
	  revQ2End = newQ2End;
	  revB     = newB;
	  revAlpha = newAlpha;
	}


	revXEnd[trialIdxRev]  = revXBgn[trialIdxRev-1];
	revQ2End[trialIdxRev] = revQ2Bgn[trialIdxRev-1];
	revXBgn[trialIdxRev]  = commitStrain/xm;
	revQ2Bgn[trialIdxRev] = commitQ2;
	
	//b
	if ( revB[trialIdxRev-1] == 0 ) { // after reversal point with b=0
	  b = 0;
	} else if (revXEnd[trialIdxRev]*revXBgn[trialIdxRev] > 0) { //consecutive reverse at same sign of "x"
	  b = 0;
	} else { // reverse at different sign of "x"
	  b = (this->calcB)(fabs(commitStrain),a,c,(this->calcHeq)(fabs(commitStrain))*Ch,(this->calcU)(fabs(commitStrain))*Cu);
	}

	
	//alpha
	if (trialDStrain > 0) {
	  alpha = (this->compAlpha)(a,revB[trialIdxRev-1],b,c,revXEnd[trialIdxRev],revXBgn[trialIdxRev],revAlpha[trialIdxRev-1]);
	} else {
	  alpha = (this->compAlpha)(a,revB[trialIdxRev-1],b,c,-revXEnd[trialIdxRev],-revXBgn[trialIdxRev],revAlpha[trialIdxRev-1]);
	}
		  
	revB[trialIdxRev] = b;
	revAlpha[trialIdxRev] = alpha;
                
      }

    }

    
    //forget reversal points
    if  (fabs(trialStrain) == trialMaxStrain) { //if reach skeleton curve
      trialIdxRev = 0;
    } else if (trialIdxRev >= 2 ) { //if pass through reversal point
      while (trialIdxRev >= 2 && (x-revXBgn[trialIdxRev])*(x-revXEnd[trialIdxRev]) > 0) {
	trialIdxRev--;
      }
    }
    
  }   


  // calculate stress

  // Q1 component (nonlinear elastic)
  if (trialStrain > 0) {
    trialQ1 = (this->compQ1)(u,n,fm,x);
    q1Tan   = (this->compQ1Derivertive)(u,n,geq,x);
  } else {
    trialQ1 = (this->compQ1)(u,n,-fm,-x);    
    q1Tan   = (this->compQ1Derivertive)(u,n,geq,-x);
  }
        
  // Q2 component (hysteretic)
  if (trialIdxRev == 0) {// skeleton curve
    
    if (trialDStrain > 0) {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,-fm,-x);
      q2Tan   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,geq,x);
    } else {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,fm,x);
      q2Tan   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,geq,-x);
    }

  } else if (trialIdxRev == 1) { //unload

    if (trialDStrain > 0) {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,fm,x);
      q2Tan   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,geq,x);
    } else {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,-fm,-x);
      q2Tan   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,geq,-x);
    }

  } else { //reverse

    if (trialDStrain > 0) {
      trialQ2 = (this->compQ2Masing)(u,a,revB[trialIdxRev],c,fm,x,revXBgn[trialIdxRev],revQ2Bgn[trialIdxRev],revAlpha[trialIdxRev]);
      q2Tan   = (this->compQ2MasingDerivertive)(u,a,revB[trialIdxRev],c,geq,x,revXBgn[trialIdxRev],revAlpha[trialIdxRev]);
    } else {
      trialQ2 = (this->compQ2Masing)(u,a,revB[trialIdxRev],c,-fm,-x,-revXBgn[trialIdxRev],revQ2Bgn[trialIdxRev],revAlpha[trialIdxRev]);
      q2Tan   = (this->compQ2MasingDerivertive)(u,a,revB[trialIdxRev],c,geq,-x,-revXBgn[trialIdxRev],revAlpha[trialIdxRev]);
    }

  }


  //(strain->stress, strain->modulus)
  trialStress  = trialQ1 + trialQ2;

  //trialTangent = (trialStress-commitStress)/trialDStrain;
  if (trialIfElastic) {
    trialTangent = initialStiff*Hr/Ar;
  } else {
    trialTangent = q1Tan + q2Tan;
  }

  
  //(output: stress->force, modulus->spring constant)
  trialForce = trialStress * Ar;
  trialStiff = trialTangent * Ar / Hr;

  trialForce = trialForce * Rf; //for MSS model
  trialStiff = trialStiff * Rs; //for MSS model

  // check print--------------------------------------------
  // opserr << "KikuchiAikenHDR::setTrialStrain\n";
  // -------------------------------------------------------

  return 0;

}



double 
KikuchiAikenHDR::getStress(void)
{
  return trialForce;
}

double 
KikuchiAikenHDR::getTangent(void)
{
  return trialStiff;
}

double 
KikuchiAikenHDR::getInitialTangent(void)
{
  return initialStiff;
}

double 
KikuchiAikenHDR::getStrain(void)
{
  return trialDeform;
}

int 
KikuchiAikenHDR::commitState(void)
{
  commitDeform  = trialDeform;
  commitForce   = trialForce;
  commitStiff   = trialStiff;
  commitStrain  = trialStrain;
  commitStress  = trialStress;
  commitTangent = trialTangent;
  commitIfElastic = trialIfElastic;
  commitQ1 = trialQ1;
  commitQ2 = trialQ2;
  commitMaxStrain = trialMaxStrain;
  commitDStrain   = trialDStrain;
  commitDStrainLastSign = trialDStrainLastSign;
  commitIdxRev = trialIdxRev;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenHDR::commitState\n";
  // -------------------------------------------------------
  return 0;
}

int 
KikuchiAikenHDR::revertToLastCommit(void)
{

  trialDeform  = commitDeform;
  trialForce   = commitForce;
  trialStiff   = commitStiff;
  trialStrain  = commitStrain;
  trialStress  = commitStress;
  trialTangent = commitTangent;
  trialIfElastic = commitIfElastic;
  trialQ1 = commitQ1;
  trialQ2 = commitQ2;
  trialMaxStrain = commitMaxStrain;
  trialDStrain = commitDStrain;
  trialDStrainLastSign = commitDStrainLastSign;
  trialIdxRev = commitIdxRev;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenHDR::revertToLastCommit\n";
  // -------------------------------------------------------
  return 0;
}

int 
KikuchiAikenHDR::revertToStart(void)
{

  trialDeform  = 0.0;
  trialForce   = 0.0;
  trialStiff   = initialStiff;
  trialStrain  = 0.0;
  trialStress  = 0.0;
  trialTangent = initialStiff*Hr/Ar;
  trialIfElastic = true;
  trialQ1 = 0.0;
  trialQ2 = 0.0;
  trialMaxStrain = 0.0;
  trialDStrain = 0.0;
  trialDStrainLastSign = 0;
  trialIdxRev=0;

  commitDeform  = 0.0;
  commitForce   = 0.0;
  commitStiff   = initialStiff;
  commitStrain  = 0.0;
  commitStress  = 0.0;
  commitTangent = initialStiff*Hr/Ar;
  commitIfElastic = true;
  commitQ1 = 0.0;
  commitQ2 = 0.0;
  commitMaxStrain = 0.0;
  commitDStrain = 0.0;
  commitDStrainLastSign = 0;
  commitIdxRev=0;

  revB[0] = 0.0;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenHDR::revertToStart\n";
  // -------------------------------------------------------
  return 0;
}

UniaxialMaterial *
KikuchiAikenHDR::getCopy(void)
{

  KikuchiAikenHDR *theCopy = new KikuchiAikenHDR(this->getTag(), Tp, Ar, Hr, 
					 Cg, Ch, Cu, Rs, Rf );

  // Copy temporary variables
  theCopy->tmpStrain = tmpStrain;
  theCopy->geq = geq;
  theCopy->heq = heq;
  theCopy->u = u;
  theCopy->n = n;
  theCopy->a = a;
  theCopy->b = b;
  theCopy->c = c;
  theCopy->xm = xm;
  theCopy->fm = fm;
  theCopy->x = x;
  theCopy->alpha = alpha;

  // Copy trial variables
  theCopy->trialDeform = trialDeform;
  theCopy->trialForce = trialForce;
  theCopy->trialStiff = trialStiff;
  theCopy->trialStrain = trialStrain;
  theCopy->trialStress = trialStress;
  theCopy->trialTangent = trialTangent;
  theCopy->trialIfElastic = trialIfElastic;
  theCopy->trialQ1 = trialQ1;
  theCopy->trialQ2 = trialQ2;
  theCopy->trialMaxStrain = trialMaxStrain;
  theCopy->trialDStrain = trialDStrain;
  theCopy->trialDStrainLastSign = trialDStrainLastSign;
  theCopy->trialIdxRev = trialIdxRev;

  // Copy commit variables
  theCopy->commitDeform = commitDeform;
  theCopy->commitForce = commitForce;
  theCopy->commitStiff = commitStiff;
  theCopy->commitStrain = commitStrain;
  theCopy->commitStress = commitStress;
  theCopy->commitTangent = commitTangent;
  theCopy->commitIfElastic = commitIfElastic;
  theCopy->commitQ1 = commitQ1;
  theCopy->commitQ2 = commitQ2;
  theCopy->commitMaxStrain = commitMaxStrain;
  theCopy->commitDStrain = commitDStrain;
  theCopy->commitDStrainLastSign = commitDStrainLastSign;
  theCopy->commitIdxRev = commitIdxRev;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenHDR::getCopy\n";
  // -------------------------------------------------------

  return theCopy;
}

int 
KikuchiAikenHDR::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
KikuchiAikenHDR::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
KikuchiAikenHDR::Print(OPS_Stream &s, int flag)
{
  s << "KikuchiAikenHDR : " << this->getTag() << endln;
  s << "  Tp: " << Tp << endln;
  s << "  Ar: " << Ar << endln;
  s << "  Hr: " << Hr << endln;
  s << "  Cg: " << Cg << endln;
  s << "  Ch: " << Ch << endln;
  s << "  Cu: " << Cu << endln;
  s << "  Rs: " << Rs << endln;
  s << "  Rf: " << Rf << endln;

  return;
}



//------------------------------------------------------------
//formulae of Kikuchi&Aiken model
//------------------------------------------------------------

double
KikuchiAikenHDR::compQ1(double u, double n, double fm, double x)
{
  return 0.5*(1-u)*fm*(x+pow(x,n));
}


double
KikuchiAikenHDR::compQ2Unload(double u, double a, double b, double c, double fm, double x)
{
  return u*fm*(1-2*exp(-a*(1+x))+b*(1+x)*exp(-c*(1+x)));
}

double
KikuchiAikenHDR::compQ2Masing(double u, double a, double b, double c, double fm, double x1, double x2, double q2i, double alpha)
{
  return q2i + alpha*u*fm*(2-2*exp(-a*(x1-x2))+b*(x1-x2)*exp(-c*(x1-x2)));
}


double
KikuchiAikenHDR::compAlpha(double a, double b1, double b2, double c, double x1, double x2, double alpha0)
{
  return alpha0*(2-2*exp(-a*(x1-x2))+b1*(x1-x2)*exp(-c*(x1-x2)))/(2-2*exp(-a*(x1-x2))+b2*(x1-x2)*exp(-c*(x1-x2)));
}


double
KikuchiAikenHDR::compQ1Derivertive(double u, double n, double geq, double x)
{
  return 0.5*(1-u)*geq*(1+n*pow(x,n-1));
}

double
KikuchiAikenHDR::compQ2UnloadDerivertive(double u, double a, double b, double c, double geq, double x)
{
  return u*geq*(2*a*exp(-a*(1+x))+b*exp(-c*(1+x))-b*c*(1+x)*exp(-c*(1+x)));
}


double
KikuchiAikenHDR::compQ2MasingDerivertive(double u, double a, double b, double c, double geq, double x1, double x2,double alpha)
{
  return alpha*u*geq*(2*a*exp(-a*(x1-x2))+b*exp(-c*(x1-x2))-b*c*(x1-x2)*exp(-c*(x1-x2)));
}

//bisection method to calculate parameter "a"
double
KikuchiAikenHDR::compABisection(double heq, double u, double min, double max, double tol, double lim)
{
  double aMin, aMax, aTmp ;
  double RHS,LHS ;

  RHS = (2*u-M_PI*heq)/(2*u);

  aMin = min;
  aMax = max;

  while (true) {
    aTmp = (aMin+aMax)/2;
    LHS = (1-exp(-2*aTmp))/(aTmp);
    if (fabs((LHS-RHS)/RHS) < tol) {
      break;
    } else if (LHS < RHS) {
      aMax = aTmp;
    } else {
      aMin = aTmp;
    }
  }

  return (aTmp < lim) ? aTmp : lim ; //min(aTmp,lim)

}



//------------------------------------------------------------
//parameters for each rubber
//------------------------------------------------------------

// tp=1 : Bridgestone HDR X0.6 (standard compressive stress)
double KikuchiAikenHDR::calcGeqTp1(double gm)
{if (gm<1.15) {return (0.68385*pow(gm,-0.39964))*1e6;}
         else {return (1.1015-0.61971*gm+0.22056*gm*gm-0.022200*gm*gm*gm)*1e6;}}

double KikuchiAikenHDR::calcHeqTp1(double gm)
{return 0.23834-0.028108*gm-0.000087664*gm*gm;}

double KikuchiAikenHDR::calcUTp1(double gm)
{return 0.44236-0.062943*gm-0.0037571*gm*gm;}

double KikuchiAikenHDR::calcNTp1(double gm)
{if (gm<1.5) {return 1.0;} 
        else {return 0.91173-0.59184*gm+0.43379*gm*gm;}}

double KikuchiAikenHDR::calcATp1(double gm, double heq, double u)
{if (gm<1.4) {return compABisection(heq,u,0.0,20.0,1e-6,10.19171);}
        else {return 10.19171;}}

double KikuchiAikenHDR::calcBTp1(double gm, double a, double c, double heq, double u)
{if (gm<1.4) {return 0.0;} 
        else {return c*c*(M_PI*heq/u-(2+2/a*(exp(-2*a)-1)));}}

double
KikuchiAikenHDR::calcCTp1(double gm)
{return 5.0;}



// tp=2 : Bridgestone HDR X0.6  (zero compressive stress)
double KikuchiAikenHDR::calcGeqTp2(double gm)
{if (gm<1.1) {return (0.67148*pow(gm,-0.41264))*1e6;}
        else {return (0.99225-0.45820*gm+0.14236*gm*gm-0.011238*gm*gm*gm)*1e6;}}

double KikuchiAikenHDR::calcHeqTp2(double gm)
{return 0.21907-0.021302*gm-0.0013313*gm*gm;}

double KikuchiAikenHDR::calcUTp2(double gm)
{return 0.44226-0.091218*gm+0.0030028*gm*gm;}

double KikuchiAikenHDR::calcNTp2(double gm)
{if (gm<1.5) {return 1.0;} 
        else {return 2.10973-1.6166*gm+0.58452*gm*gm;}}

double KikuchiAikenHDR::calcATp2(double gm, double heq, double u)
{if (gm<1.3) {return compABisection(heq,u,0.0,20.0,1e-6,10.38035);}
        else {return 10.38035;}}

double KikuchiAikenHDR::calcBTp2(double gm, double a, double c, double heq, double u)
{if (gm<1.3) {return 0.0;} 
        else {return c*c*(M_PI*heq/u-(2+2/a*(exp(-2*a)-1)));}}

double
KikuchiAikenHDR::calcCTp2(double gm)
{return 6.0;}



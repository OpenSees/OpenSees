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

// Written: Akiri Seki
// Created: Oct 2026
//
// Description: This file contains the class implementation for
// FatigueMemberMaterial. See FatigueMemberMaterial.h. The modified
// rainflow counting and damage accumulation follow FatigueMaterial
// (Uriz and Mahin), applied to the normalized node-level deformation
// instead of the strain of the wrapped material.
//
//   Seki A. and Simpson B. (2025) "Consistent Estimate of Brace 
//   Fatigue Life  Using Axial Deformations Evaluated from Nodal Displacements"
//   Journal of Structural Engineering, 151(8).


#include <math.h>

#include <stdlib.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <string.h>

#include <FatigueMemberMaterial.h>
#include <ID.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <Node.h>

#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <elementAPI.h>

void* OPS_FatigueMemberMaterial()
{
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial FatigueMember tag? matTag? node1? node2? <-k k?>" << endln;
    opserr << " <-Dmax dmax?> <-dbar0 dbar0?> <-m m?>";
    opserr << " <-min min?> <-max max?>" << endln;
    return 0;
  }

  // tag, matTag, node1, node2
  int idata[4];
  int numdata = 4;
  if (OPS_GetIntInput(&numdata,idata)<0) {
    opserr<<"WARNING invalid int inputs\n";
    return 0;
  }

  int nodeI = idata[2];
  int nodeJ = idata[3];
  double k       =  1.0;
  double Dmax    =  1.0;
  double dbar0   =  0.0307; // calibrated value for HSS braces
  double m       = -0.481;
  double dbarmin = NEG_INF_STRAIN;
  double dbarmax = POS_INF_STRAIN;
  numdata = 1;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* type = OPS_GetString();
    if (OPS_GetNumRemainingInputArgs() < 1) {
      opserr << "WARNING missing value for " << type
	     << " in uniaxialMaterial FatigueMember " << idata[0] << endln;
      return 0;
    }
    if (strcmp(type,"-k")==0) {
      if (OPS_GetDoubleInput(&numdata,&k)<0) {
	opserr<<"WARNING invalid double inputs\n";
	return 0;
      }
    } else if (strcmp(type,"-Dmax")==0) {
      if (OPS_GetDoubleInput(&numdata,&Dmax)<0) {
	opserr<<"WARNING invalid double inputs\n";
	return 0;
      }
    } else if (strcmp(type,"-dbar0")==0) {
      if (OPS_GetDoubleInput(&numdata,&dbar0)<0) {
	opserr<<"WARNING invalid double inputs\n";
	return 0;
      }
    } else if (strcmp(type,"-m")==0) {
      if (OPS_GetDoubleInput(&numdata,&m)<0) {
	opserr<<"WARNING invalid double inputs\n";
	return 0;
      }
    } else if (strcmp(type,"-min")==0) {
      if (OPS_GetDoubleInput(&numdata,&dbarmin)<0) {
	opserr<<"WARNING invalid double inputs\n";
	return 0;
      }
    } else if (strcmp(type,"-max")==0) {
      if (OPS_GetDoubleInput(&numdata,&dbarmax)<0) {
	opserr<<"WARNING invalid double inputs\n";
	return 0;
      }
    } else {
      opserr << "WARNING unknown option " << type
	     << " in uniaxialMaterial FatigueMember " << idata[0] << endln;
      return 0;
    }
  }

  if (nodeI == nodeJ) {
    opserr << "WARNING node1 and node2 must be different for uniaxialMaterial FatigueMember "
	   << idata[0] << endln;
    return 0;
  }
  if (k < 0.5 || k > 1.0) {
    opserr << "WARNING k must be between 0.5 and 1.0 for uniaxialMaterial FatigueMember "
	   << idata[0] << ", k = " << k << endln;
    return 0;
  }

  UniaxialMaterial* mat = OPS_getUniaxialMaterial(idata[1]);
  if (mat == 0) {
    opserr << "WARNING component material does not exist\n";
    opserr << "Component material: " << idata[1];
    opserr << "\nuniaxialMaterial FatigueMember: " << idata[0] << endln;
    return 0;
  }

  UniaxialMaterial* theMat = new FatigueMemberMaterial(idata[0], *mat, nodeI, nodeJ, k,
						       Dmax, dbar0, m, dbarmin, dbarmax);
  if (theMat == 0) {
    opserr << "WARNING: failed to create FatigueMemberMaterial material\n";
    return 0;
  }

  return theMat;
}

FatigueMemberMaterial::FatigueMemberMaterial(int tag, UniaxialMaterial &material,
					     int nodeI, int nodeJ, double k_end,
					     double dmax, double dbar_0, double slope_m,
					     double dbarmin, double dbarmax)
  :UniaxialMaterial(tag,MAT_TAG_FatigueMember), theMaterial(0),
   k(k_end), L0(0.0), Cdbar(0.0),
   Cfailed(false), trialStrain(0)
{
  nodeTags[0] = nodeI;
  nodeTags[1] = nodeJ;
  theNodes[0] = 0;
  theNodes[1] = 0;

  DI  = 0; //Damage index
  X   = 0; //Range in consideration
  Y   = 0; //Previous Adjacent Range
  A   = 0; //Peak or valley 1
  B   = 0; //Peak or valley 2
  C   = 0; //Peak or valley 3
  D   = 0; //Peak or valley 4
  PCC = 0; /*Previous Cycle counter flag if >1 then previous 'n'
	     cycles did not flag a complete cycle */
  R1F = 0; //Flag for first  peak count
  R2F = 0; //Flag for second peak count
  cSlope  = 0; //Current Slope
  PS  = 0; //Previous slope
  EP  = 0; //Previous normalized deformation
  SF  = 0; /*Start Flag = 0 if very first deformation, (i.e. when initializing)
	     = 1 otherwise */
  DL  = 0; //Damage if current deformation was last peak.

  SR1 = 0;
  NC1 = 0;
  SR2 = 0;
  NC2 = 0;
  SR3 = 0;
  NC3 = 0;

  if ( dmax > 1.0 || dmax < 0.0 ) {
    opserr << "FatigueMemberMaterial::FatigueMemberMaterial " <<
      "- Dmax must be between 0 and 1, assuming Dmax = 1\n" ;
    Dmax    = 1;
  } else
    Dmax    = dmax;

  dbar0   = dbar_0;
  m       = slope_m;
  minDbar = dbarmin;
  maxDbar = dbarmax;

  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "FatigueMemberMaterial::FatigueMemberMaterial " <<
      " -- failed to get copy of material\n" ;
    exit(-1);
  }

  energy = 0;
  CStress = 0;
}

FatigueMemberMaterial::FatigueMemberMaterial()
  :UniaxialMaterial(0,MAT_TAG_FatigueMember), theMaterial(0),
   k(1.0), L0(0.0), Cdbar(0.0),
   Cfailed(false), trialStrain(0)
{
  nodeTags[0] = 0;
  nodeTags[1] = 0;
  theNodes[0] = 0;
  theNodes[1] = 0;

  DI  = 0;
  X   = 0;
  Y   = 0;
  A   = 0;
  B   = 0;
  C   = 0;
  D   = 0;
  PCC = 0;
  R1F = 0;
  R2F = 0;
  cSlope  = 0;
  PS  = 0;
  EP  = 0;
  SF  = 0;
  DL  = 0;

  Dmax    = 0;
  dbar0   = 0;
  m       = 0;
  minDbar = 0;
  maxDbar = 0;

  SR1 = 0;
  NC1 = 0;
  SR2 = 0;
  NC2 = 0;
  SR3 = 0;
  NC3 = 0;

  energy = 0;
  CStress = 0;
}

FatigueMemberMaterial::~FatigueMemberMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

// Look up the two nodes in the Domain and compute the initial length.
// Called the first time the nodes are needed, and again after getCopy()
// or recvSelf(), since only the node tags are copied or sent.
int
FatigueMemberMaterial::setNodes(void)
{
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0) {
    opserr << "FatigueMemberMaterial::setNodes() - material tag " << this->getTag()
	   << " - no Domain available\n";
    return -1;
  }

  for (int i = 0; i < 2; i++) {
    theNodes[i] = theDomain->getNode(nodeTags[i]);
    if (theNodes[i] == 0) {
      opserr << "FatigueMemberMaterial::setNodes() - material tag " << this->getTag()
	     << " - node " << nodeTags[i] << " does not exist in the Domain"
	     << " (in parallel runs both nodes must be on the same process)\n";
      theNodes[0] = 0;
      theNodes[1] = 0;
      return -1;
    }
  }

  const Vector &crdI = theNodes[0]->getCrds();
  const Vector &crdJ = theNodes[1]->getCrds();
  int ndm = crdI.Size();
  if (crdJ.Size() != ndm) {
    opserr << "FatigueMemberMaterial::setNodes() - material tag " << this->getTag()
	   << " - nodes " << nodeTags[0] << " and " << nodeTags[1]
	   << " have different numbers of coordinates\n";
    theNodes[0] = 0;
    theNodes[1] = 0;
    return -1;
  }

  double L2 = 0.0;
  for (int i = 0; i < ndm; i++) {
    double d = crdJ(i) - crdI(i);
    L2 += d*d;
  }
  L0 = sqrt(L2);

  if (L0 <= 0.0) {
    opserr << "FatigueMemberMaterial::setNodes() - material tag " << this->getTag()
	   << " - nodes " << nodeTags[0] << " and " << nodeTags[1]
	   << " have the same coordinates\n";
    theNodes[0] = 0;
    theNodes[1] = 0;
    return -1;
  }

  return 0;
}

// Normalized axial deformation between the two nodes, dbar = (L-L0)/(k*L0).
// Uses the translational displacements (the first ndm components).
int
FatigueMemberMaterial::computeDbar(double &dbar)
{
  if (theNodes[0] == 0 || theNodes[1] == 0)
    if (this->setNodes() < 0)
      return -1;

  const Vector &crdI  = theNodes[0]->getCrds();
  const Vector &crdJ  = theNodes[1]->getCrds();
  const Vector &dispI = theNodes[0]->getTrialDisp();
  const Vector &dispJ = theNodes[1]->getTrialDisp();

  int ndm = crdI.Size();
  if (dispI.Size() < ndm || dispJ.Size() < ndm) {
    opserr << "FatigueMemberMaterial::computeDbar() - material tag " << this->getTag()
	   << " - nodes need at least " << ndm << " translational dofs\n";
    return -1;
  }

  double L2 = 0.0;
  for (int i = 0; i < ndm; i++) {
    double d = (crdJ(i) + dispJ(i)) - (crdI(i) + dispI(i));
    L2 += d*d;
  }

  dbar = (sqrt(L2) - L0)/(k*L0);
  return 0;
}

static int sign(double a) {
  if (a < 0)
    return -1;
  else if (a == 0)
    return 0;
  else
    return 1;
}

int
FatigueMemberMaterial::setTrialStrain(double strain, double strainRate)
{
  trialStrain = strain;
  return theMaterial->setTrialStrain(strain, strainRate);
}

double
FatigueMemberMaterial::getStress(void)
{
  if (Cfailed)
    // Reduce stress to nearly zero
    return theMaterial->getStress()*1.0e-8;
  else
    return theMaterial->getStress();
}

double
FatigueMemberMaterial::getTangent(void)
{
  if (Cfailed)
    // Reduce tangent to nearly zero
    return 1.0e-8*theMaterial->getInitialTangent();
  else
    return theMaterial->getTangent();
}

double
FatigueMemberMaterial::getDampTangent(void)
{
  if (Cfailed)
    return 0.0;
  else
    return theMaterial->getDampTangent();
}

double
FatigueMemberMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

double
FatigueMemberMaterial::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int
FatigueMemberMaterial::commitState(void)
{
  // The rainflow counting below is the same as FatigueMaterial::commitState(),
  // applied to the normalized node-level deformation dbar instead of the
  // strain of the wrapped material. Domain::commit() commits the nodes before
  // the elements, so the node displacements are those of the current step.

  // for recorder
  SR1 = 0;
  NC1 = 0;

  // A missing node is a modeling error; stop rather than continue with an
  // uncommitted base material (Domain::commit() ignores the return value).
  double dbar;
  if (this->computeDbar(dbar) < 0) {
    opserr << "FatigueMemberMaterial::commitState() - material tag " << this->getTag()
	   << " - cannot compute the node-level deformation, exiting\n";
    exit(-1);
  }
  Cdbar = dbar;

  // No need to continue if the uniaxial material copy
  // has already failed.
  if (Cfailed) {
    return 0;
  }

  //Simple check to see if we reached max deformation capacities
  if (dbar >= maxDbar || dbar <= minDbar) {
      Cfailed = true;
      opserr << "FatigueMemberMaterial: material tag " << this->getTag() << " failed from excessive deformation\n";
      DI = Dmax;
      DL = Dmax;
      return 0;
  }

  //Initialize the fatigue parameters if they have
  // not been initialized yet
  if (SF == 0) {

    A   = dbar;
    SF  = 1  ;
    EP  = dbar;
    // Initialize other params if not done so already
    PCC = 0;
    B   = 0;
    C   = 0;
    D   = 0;

  }

  /* Now we need to determine if we are at a peak or not
     If we are, then we need to do some calcs to determine
     the amount of damage suffered. If we are not at a peak, we need to
     pretend like we are at a peak, so that we can calculate the damage
     as if it WERE a peak.
  */

  // Determine the slope of the deformation history
  if ( EP == dbar ) {
    cSlope = PS;         // No real slope here....
  } else {
    cSlope = dbar - EP;   // Determine Current Slope
  }


  // If we are at a peak or a valley, then check for damage
  if (sign(PS) != sign(cSlope) && sign(PS) != 0 ) {

    if ( R1F == 0 )  {    // mark second peak

      B = EP;
      Y = fabs(B-A);
      R1F = 1;

    } else {   // start at least the third peak

      // begin modified Rainflow cycle counting
      if (PCC == 1) {

	D = EP;
	X = fabs(D-C);

      } else {

	C = EP;
	X = fabs(C-B);

      }

      if (X < Y) {

	PCC = PCC + 1;

	if (PCC == 1) {
	  Y = fabs(C-B);
	} else if (PCC == 2 ) {
	  // Count X = |D-C| as a 1.0 cycle
	  DI = DI + 1.0 / fabs(pow( (X/dbar0) , 1/m )) ;
	  SR1 = X;
	  NC1 = 1.0;
	  // Reset parameters
	  D = 0;
	  C = 0;
	  Y = fabs(B-A);
	  PCC = 0;

	}

      } else {

	if (PCC == 1 ) {

	  // Count Y = |C-B| as a 1.0 cycle
	  DI = DI + 1.0 / fabs(pow( (Y/dbar0) , 1/m ));
	  SR1 = Y;
	  NC1 = 1.0;

	  // Reset parameters
	  B = D;
	  C = 0;
	  D = 0;
	  Y = fabs(B-A);
	  PCC = 0;

	} else {

	  // Count Y = |A-B| as a 0.5 cycle
	  DI = DI + 0.5 / fabs(pow( (Y/dbar0), 1/m ));

	  // For recorder
	  SR1 = Y;
	  NC1 = 0.5;

	  // Reset parameters
	  A = B;
	  B = C;
	  C = 0;
	  D = 0;
	  Y = X;
	  PCC = 0;

	}

      }

    }

    // Flag failure if we have reached that point
    if (DI >= Dmax )  {
      // Most likely will not fail at this point, more
      // likely at the pseudo peak. But this step is
      // is important for accumulating damage
      Cfailed = true;
      opserr << "FatigueMemberMaterial: material tag " << this->getTag() << " failed at peak\n";
      DL=DI;
    } else {
      Cfailed = false;
      DL=DI;
    }

  }
  if (Cfailed == false) {

    // Now check for damage, although we may not be at a peak at all.
    // Store temporary damage only as if it were the last peak: DL
    // Commit to DI only if failure occurs.
    if (B == 0 && C == 0 &&  D == 0) {

      // If we have not yet found the second peak
      X = fabs(dbar - A);

      if (fabs(X) < 1e-10) {
	DL = DI ;
	//values for recorder
	SR2 = 0.0;
	NC2 = 0.0;
	SR3 = 0.0;
	NC3 = 0.0;
      } else {
	DL = DI +  0.5 / fabs(pow( (X/dbar0), 1/m ));
	//values for recorder
	SR2 = X;
	NC2 = 0.5;
	SR3 = 0.0;
	NC3 = 0.0;
      }

    } else if (B != 0 && C == 0 &&  D == 0) {

      // On our way to find point C. Range Y defined, no X yet
      X = fabs(dbar - B);

      if (fabs(X) < 1e-10) {
	DL = DI;
	//values for recorder
	SR2 = 0.0;
	NC2 = 0.0;
      } else {
	DL = DI +  0.5 / fabs(pow( (X/dbar0) , 1/m ));
	//values for recorder
	SR2 = X;
	NC2 = 0.5;
      }

      if (fabs(Y) < 1e-10) {
	DL = DL;
	//values for recorder
	SR3 = 0.0;
	NC3 = 0.0;
      } else {
	DL = DL +  0.5 / fabs(pow( (Y/dbar0) , 1/m ));
	//values for recorder
	SR3 = Y;
	NC3 = 0.5;
      }

    } else if (B != 0 && C != 0 &&  D == 0) {

      // Two ranges stored, but no cycles for either stored
      // There are two scenarios that can result:
      // 1.)  |A-D| is the predominate 1/2 cycle,
      //      and |B-C| is a small full cycle.
      //
      // 2.)  One full cycle at |D-C|, 1/2 cycle at |A-B|

      if (fabs(A-dbar)> fabs(A-B)) {

	//   count 1/2 cycle at |D-A|, and one full cycle at |B-C|.
	X = fabs(dbar-A);

	if (fabs(Y) < 1e-10) {
	  DL = DI;
	  //values for recorder
	  SR3 = 0.0;
	  NC3 = 0.0;
	} else {
	  DL = DI +  1.0 / fabs(pow( (Y/dbar0) , 1/m ));
	  //values for recorder
	  SR3 = Y;
	  NC3 = 1.0;
	}

	if (fabs(X) < 1e-10) {
	  DL = DL;
	  //values for recorder
	  SR2 = 0.0;
	  NC2 = 0.0;
	} else {
	  DL = DL +  0.5 / fabs(pow( (X/dbar0) , 1/m ));
	  //values for recorder
	  SR2 = X;
	  NC2 = 0.5;
	}

      } else {

	// One full cycle of |C-D| and 1/2 cycle of |A-B|

	if (fabs(C-dbar) < 1e-10) {
	  DL = DI;
	  //values for recorder
	  SR3 = 0.0;
	  NC3 = 0.0;
	} else {
	  DL = DI +  1.0 / fabs(pow( ( fabs(C-dbar)/dbar0 ) , 1/m ));
	  //values for recorder
	  SR3 = fabs(C-dbar);
	  NC3 = 1.0;
	}

	if (fabs(A-B) < 1e-10) {
	  DL = DL;
	  //values for recorder
	  SR2 = 0.0;
	  NC2 = 0.0;
	} else {
	  DL = DL +  0.5 / fabs(pow( (fabs(A-B) /dbar0) , 1/m ));
	  //values for recorder
	  SR2 = fabs(A-B);
	  NC2 = 0.5;
	}
      }
    }

    // Did we fail before a peak? Same as FatigueMaterial: only if the stress
    // of the wrapped material is positive, so copies in different fibers can
    // fail in different steps.
    double mStress = theMaterial->getStress();
    if (DL > Dmax && mStress > 0.0 ) {
      DI = DL;
      Cfailed = true;
      opserr << "FatigueMemberMaterial: material tag " << this->getTag() << " failed at pseudo peak\n";
    } else {
      Cfailed = false;
    }

  }

  // Same expression as FatigueMaterial (SAJalali). Note that PS is the
  // previous slope of dbar, not the previous strain; kept as-is for now.
  if (!Cfailed)
  {
	  double TStress = getStress();
	  energy += 0.5*(trialStrain - PS)*(TStress + CStress);
	  CStress = TStress;
  }

  PS = cSlope;     // Previous Slope
  EP = dbar;       // Keep track of previous deformation

  // Check if failed at current step
  if (Cfailed) {
    return 0;
  }
  else
    return theMaterial->commitState();
}


int
FatigueMemberMaterial::revertToLastCommit(void)
{
  // Check if failed at last step
  if (Cfailed)
    return 0;
  else
    return theMaterial->revertToLastCommit();
}

int
FatigueMemberMaterial::revertToStart(void)
{
  Cfailed = false;
  DI  = 0;
  X   = 0;
  Y   = 0;
  A   = 0;
  B   = 0;
  C   = 0;
  D   = 0;
  PCC = 0;
  R1F = 0;
  R2F = 0;
  cSlope  = 0;
  PS  = 0;
  EP  = 0;
  SF  = 0;
  DL  = 0;
  Cdbar = 0;

  //values for recorder
  SR1 = 0;
  NC1 = 0;
  SR2 = 0;
  NC2 = 0;
  SR3 = 0;
  NC3 = 0;

  return theMaterial->revertToStart();
}

UniaxialMaterial *
FatigueMemberMaterial::getCopy(void)
{
  FatigueMemberMaterial *theCopy =
    new FatigueMemberMaterial(this->getTag(), *theMaterial, nodeTags[0], nodeTags[1], k,
			      Dmax, dbar0, m, minDbar, maxDbar);

  // copy the full fatigue state; the nodes are looked up again by tag
  theCopy->Cdbar = Cdbar;
  theCopy->DI  = DI;
  theCopy->X   = X;
  theCopy->Y   = Y;
  theCopy->A   = A;
  theCopy->B   = B;
  theCopy->C   = C;
  theCopy->D   = D;
  theCopy->PCC = PCC;
  theCopy->R1F = R1F;
  theCopy->R2F = R2F;
  theCopy->cSlope = cSlope;
  theCopy->PS  = PS;
  theCopy->EP  = EP;
  theCopy->SF  = SF;
  theCopy->DL  = DL;
  theCopy->SR1 = SR1;
  theCopy->NC1 = NC1;
  theCopy->SR2 = SR2;
  theCopy->NC2 = NC2;
  theCopy->SR3 = SR3;
  theCopy->NC3 = NC3;
  theCopy->Cfailed = Cfailed;
  theCopy->trialStrain = trialStrain;
  theCopy->energy  = energy;
  theCopy->CStress = CStress;

  return theCopy;
}

int
FatigueMemberMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static ID dataID(5);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if ( matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  dataID(3) = nodeTags[0];
  dataID(4) = nodeTags[1];
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "FatigueMemberMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(26);
  dataVec(0)  = DI;
  dataVec(1)  = X;
  dataVec(2)  = Y;
  dataVec(3)  = A;
  dataVec(4)  = B;
  dataVec(5)  = C;
  dataVec(6)  = D;
  dataVec(7)  = PCC;
  dataVec(8)  = R1F;
  dataVec(9)  = R2F;
  dataVec(10) = cSlope;
  dataVec(11) = PS;
  dataVec(12) = EP;
  dataVec(13) = SF;
  dataVec(14) = DL;
  dataVec(15) = Dmax;
  dataVec(16) = dbar0;
  dataVec(17) = m;
  dataVec(18) = minDbar;
  dataVec(19) = maxDbar;

  if (Cfailed == true)
    dataVec(20) = 1.0;
  else
    dataVec(20) = 0.0;

  dataVec(21) = k;
  dataVec(22) = Cdbar;
  dataVec(23) = trialStrain;
  dataVec(24) = energy;
  dataVec(25) = CStress;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "FatigueMemberMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "FatigueMemberMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int
FatigueMemberMaterial::recvSelf(int cTag, Channel &theChannel,
				FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(5);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "FatigueMemberMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "FatigueMemberMaterial::recvSelf() - failed to create Material with classTag "
	   << dataID(1) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  // nodes are looked up again by tag the next time they are needed
  nodeTags[0] = dataID(3);
  nodeTags[1] = dataID(4);
  theNodes[0] = 0;
  theNodes[1] = 0;

  static Vector dataVec(26);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "FatigueMemberMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  DI   = dataVec(0);
  X    = dataVec(1);
  Y    = dataVec(2);
  A    = dataVec(3);
  B    = dataVec(4);
  C    = dataVec(5);
  D    = dataVec(6);
  PCC  = int(dataVec(7));
  R1F  = int(dataVec(8));
  R2F  = int(dataVec(9));
  cSlope   = dataVec(10);
  PS   = dataVec(11);
  EP   = dataVec(12);
  SF   = int(dataVec(13));
  DL   = dataVec(14);
  Dmax = dataVec(15);
  dbar0 = dataVec(16);
  m    = dataVec(17);
  minDbar = dataVec(18);
  maxDbar = dataVec(19);

  if (dataVec(20) == 1.0)
    Cfailed = true;
  else
    Cfailed = false;

  k       = dataVec(21);
  Cdbar   = dataVec(22);
  trialStrain = dataVec(23);
  energy  = dataVec(24);
  CStress = dataVec(25);

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "FatigueMemberMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

//Printing of current damage.  NOTE:  The damage that is returned
//  is damage at the pseudo peak, DL, ( if not at a peak when the
//  print function  is called). The generic print will print both
//  the damage recorded at the last peak, along with current damage

void
FatigueMemberMaterial::Print(OPS_Stream &s, int flag)
{
	if (flag == 100) {
		s << DL << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "FatigueMemberMaterial tag: " << this->getTag() << endln;
		s << "\tMaterial: " << theMaterial->getTag() << endln;
		s << "\tNodes: " << nodeTags[0] << " " << nodeTags[1] << " k: " << k << endln;
		s << "\tDI: " << DI << " Dmax: " << Dmax << endln;
		s << "\tdbar0: " << dbar0 << " m: " << m << endln;
		s << "\tDL: " << DL << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"FatigueMemberMaterial\", ";
		s << "\"material\": \"" << theMaterial->getTag() << "\", ";
		s << "\"nodes\": [" << nodeTags[0] << ", " << nodeTags[1] << "], ";
		s << "\"k\": " << k << ", ";
		s << "\"tDI\": " << DI << ", ";
		s << "\"Dmax\": " << Dmax << ", ";
		s << "\"dbar0\": " << dbar0 << ", ";
		s << "\"m\": " << m << ", ";
		s << "\"tDL\": " << DL << "}";
	}
}

Response*
FatigueMemberMaterial::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{
  if (argc == 0)
    return 0;

  Response *theResponse = 0;

  theOutput.tag("UniaxialMaterialOutput");
  theOutput.attr("matType", this->getClassType());
  theOutput.attr("matTag", this->getTag());


  // stress
  if (strcmp(argv[0],"stress") == 0) {
    theOutput.tag("ResponseType", "sigma11");
    theResponse =  new MaterialResponse(this, 1, this->getStress());
  }
  // tangent
  else if (strcmp(argv[0],"tangent") == 0) {
    theOutput.tag("ResponseType", "C11");
    theResponse =  new MaterialResponse(this, 2, this->getTangent());
  }

  // strain
  else if (strcmp(argv[0],"strain") == 0) {
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 3, this->getStrain());
  }

  // stress and strain
  else if ((strcmp(argv[0],"stressStrain") == 0) ||
	   (strcmp(argv[0],"stressANDstrain") == 0)) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 4, Vector(2));
  }

  // damage
  else if (strcmp(argv[0],"damage") == 0) {
    theResponse =  new MaterialResponse(this, 5, DI);
    theOutput.tag("ResponseType", "DI");
  }

  else if (strcmp(argv[0],"cyclesAndRange") == 0) {
    theOutput.tag("ResponseType", "UnknownResponse");
    theOutput.tag("ResponseType", "UnknownResponse");
    theOutput.tag("ResponseType", "UnknownResponse");
    theOutput.tag("ResponseType", "UnknownResponse");
    theOutput.tag("ResponseType", "UnknownResponse");
    theOutput.tag("ResponseType", "UnknownResponse");
    theResponse =  new MaterialResponse(this, 6, Vector(6));
  }

  else if (strcmp(argv[0],"failure") == 0) {
    int res = 0;
    theResponse =  new MaterialResponse(this, 7, res);
    theOutput.tag("ResponseType", "Failure");
  }

  // normalized node-level deformation
  else if (strcmp(argv[0],"dbar") == 0) {
    theOutput.tag("ResponseType", "dbar");
    theResponse =  new MaterialResponse(this, 8, Cdbar);
  }

  theOutput.endTag();
  return theResponse;
}

int
FatigueMemberMaterial::getResponse(int responseID, Information &matInfo)
{
  static Vector stressStrain(2);
  static Vector cyclesAndRange(6);

  switch (responseID) {
  case 1:
    matInfo.setDouble(this->getStress());
    return 0;

  case 2:
    matInfo.setDouble(this->getTangent());
    return 0;

  case 3:
    matInfo.setDouble(this->getStrain());
    return 0;

  case 4:
    stressStrain(0) = this->getStress();
    stressStrain(1) = this->getStrain();
    matInfo.setVector(stressStrain);
    return 0;

  case 5:
    matInfo.setDouble(DI);
    return 0;

  case 6:
    cyclesAndRange(0) = NC1;
    cyclesAndRange(1) = SR1;
    cyclesAndRange(2) = NC2;
    cyclesAndRange(3) = SR2;
    cyclesAndRange(4) = NC3;
    cyclesAndRange(5) = SR3;
    matInfo.setVector(cyclesAndRange);
    return 0;

  case 7:
    if (Cfailed == true)
      matInfo.setInt(1);
    else
      matInfo.setInt(0);
    return 0;

  case 8:
    matInfo.setDouble(Cdbar);
    return 0;

  default:
    return -1;

  }
}

bool
FatigueMemberMaterial::hasFailed(void) {
  return Cfailed;
}

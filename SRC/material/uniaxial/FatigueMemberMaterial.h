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
// Description: This file contains the class definition for
// FatigueMemberMaterial. FatigueMemberMaterial wraps a UniaxialMaterial
// and imposes fatigue limits based on the normalized axial deformation
// between two nodes ("fatigue-member"), rather than the strain of the
// wrapped material ("fatigue-material"). Damage is accumulated with 
// the Coffin-Manson relationship, Miner's rule and the modified rainflow 
// counting scheme of FatigueMaterial (Uriz and Mahin).
//
//   dbar = (L - L0) / (k * L0)
//
// where L and L0 are the current and initial distances between the two
// nodes and k is the effective length factor. Reference:
//
//   Seki, A., and B. Simpson (2025). "Consistent Estimate of Brace
//      Fatigue Life Using Axial Deformations Evaluated from Nodal
//      Displacements." J. Struct. Eng., 151(8): 04025108.
//      https://doi.org/10.1061/JSENDH.STENG-13866
//
// The two nodes are looked up in the Domain by tag the first time they
// are needed (and again after getCopy() or recvSelf()), so the nodes may
// be defined before or after this material.

#ifndef FatigueMemberMaterial_h
#define FatigueMemberMaterial_h

#include <UniaxialMaterial.h>

class Node;

class FatigueMemberMaterial : public UniaxialMaterial
{
 public:
  FatigueMemberMaterial(int tag, UniaxialMaterial &material,
			int nodeI, int nodeJ,
			double k        =  1.0,
			double Dmax     =  1.0,
			double dbar0    =  0.0307,
			double m        = -0.481,
			double minDbar  = NEG_INF_STRAIN,
			double maxDbar  = POS_INF_STRAIN);
  FatigueMemberMaterial();
  ~FatigueMemberMaterial();

  const char *getClassType(void) const {return "FatigueMemberMaterial";};

  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void);
  double getStrainRate(void);
  double getStress(void);
  double getTangent(void);
  double getDampTangent(void);
  double getInitialTangent(void) {return theMaterial->getInitialTangent();}

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  UniaxialMaterial *getCopy(void);

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel,
	       FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag =0);

  Response *setResponse (const char **argv, int argc, OPS_Stream &s);
  int getResponse (int responseID, Information &matInformation);
  bool hasFailed(void);

  virtual double getEnergy(void) { return energy; }

 protected:

 private:
  int setNodes(void);               // look up the two nodes in the Domain
  int computeDbar(double &dbar);    // normalized node-level deformation

  double energy, CStress;

  UniaxialMaterial *theMaterial;

  // node-level deformation
  int nodeTags[2];      // tags of the two nodes
  Node *theNodes[2];    // resolved from nodeTags
  double k;             // effective length factor
  double L0;            // initial chord length between the two nodes
  double Cdbar;         // committed normalized deformation

  // modified rainflow counting, same as FatigueMaterial
  double DI; //Damage index
  double  X; //Range in consideration
  double  Y; //Previous Adjacent Range
  double  A; //Peak or valley 1
  double  B; //Peak or valley 2
  double  C; //Peak or valley 3
  double  D; //Peak or valley 4
  int   PCC; /*Previous Cycle counter flag if >1 then previous 'n'
	       cycles did not flag a complete cycle */
  int   R1F; //Flag for first  peak count
  int   R2F; //Flag for second peak count
  double cSlope; //Current Slope
  double PS; //Previous slope
  double EP; //Previous normalized deformation
  int    SF; /*Start Flag = 0 if very first deformation,
	       (i.e. when initializing)    = 1 otherwise */
  double DL; //Damage if current deformation was last peak.

  double Dmax;
  double dbar0;   // normalized deformation at which failure occurs in one cycle
  double m;

  double minDbar;
  double maxDbar;

  bool Cfailed;
  double trialStrain;   // strain of the wrapped material

  // For recording ranges (SRXX) and Number of Cycles (NCXX)
  double SR1;  // Committed range at peak
  double NC1;  // Committed number of cycles at SR1 (i.e. 1.0 or 0.5)
  double SR2;  // Committed range 2 at PSEUDO peak
  double NC2;  // Committed number of cycles at SR2 (at PSEUDO peak)
  double SR3;  // Committed range 3 at PSEUDO peak
  double NC3;  // Committed number of cycles at SR3 (at PSEUDO peak)
};


#endif

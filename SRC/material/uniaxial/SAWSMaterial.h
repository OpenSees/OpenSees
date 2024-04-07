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
                                                                        
// $Revision: 1.1 $
// $Date: 2009-11-23 23:29:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SAWSMaterial.h,v $

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


#ifndef SAWSMaterial_h
#define SAWSMaterial_h

#include <UniaxialMaterial.h>

class SAWSMaterial : public UniaxialMaterial
{
 public:
  SAWSMaterial(int tag,
	       double F0, double FI, double DU, double S0,
	       double R1, double R2, double R3, double R4,
	       double ALPHA, double BETA );
  SAWSMaterial();
  ~SAWSMaterial();
  
  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void);

  
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
 protected:
  
 private:
  
  double F0;
  double FI;
  double DU;
  double S0;
  double R1;
  double R2;
  double R3;
  double R4;
  double ALPHA;
  double BETA;

  // History Variables
  //double TOL     ;
  int    IFLAG   ;
  int    JFLAG   ;
  int    IMODE   ;
  double FU      ;
  int    LPATH   ;
  int    LPPREV  ;
  double FAC1    ;
  double FAC2    ;
  double FAC3    ;
  double DF1     ;
  double DF2     ;
  double DF      ;
  double DINT    ;
  double DINT2   ;
  double DA      ;
  double DB      ;
  double DY      ;
  double FORCE   ;
  double STIFF   ;
  double DISPL   ;
  double DOLD    ;
  double DUNP    ;
  double FUNP    ;
  double DUNM    ;
  double FUNM    ;
  double DMAXP   ;
  double FMAXP   ;
  double DMAXM   ;
  double FMAXM   ;
  double SP      ;
  int IYPLUS  ;
  int IYMINS  ;
  double D0      ;
  double DIFF    ;
  double DLIM    ;
  double D6      ;
  double dZero   ;
  double DINT1P  ;
  double DINT1M  ;
  double DINT4   ;
  double DINT3   ;
  double FOLD    ;
  double DUPPER  ;
  double DLOWER  ;
  double DINT5   ;


  // Converged history variables
  double cFORCE;
  double cDISPL;
  double cSTIFF;
  int cLPATH;
  int cLPPREV;
  int cIYPLUS;
  int cIYMINS;
  double cDOLD;
  double cDUNP;
  double cFUNP;
  double cDUNM;
  double cFUNM;
  double cDMAXP;
  double cFMAXP;
  double cDMAXM;
  double cFMAXM;
  double cSP   ;
};

#endif

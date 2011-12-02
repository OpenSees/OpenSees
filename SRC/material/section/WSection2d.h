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
// $Date: 2008-12-03 23:46:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/WSection2d.h,v $
                                                                        
// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// WSection2d.h. WSection2d provides the abstraction of a 
// rectangular section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
// The fiber stresses are the 11, 12, and 13 components of stress, from
// which all six beam stress resultants are obtained.

#ifndef WSection2d_h
#define WSection2d_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

class NDMaterial;

class WSection2d : public SectionForceDeformation
{
 public:
  WSection2d(); 
  WSection2d(int tag, NDMaterial &theMat,
	     double d, double tw, double bf, double tf,
	     int nfdw, int nftf, double shape = 1.0/4.8, double flag = true);
  ~WSection2d();
  
  int   setTrialSectionDeformation(const Vector &deforms); 
  const Vector &getSectionDeformation(void);
  
  const Vector &getStressResultant(void);
  const Matrix &getSectionTangent(void);
  const Matrix &getInitialTangent(void) {return this->getSectionTangent();}
  
  int   commitState(void);
  int   revertToLastCommit(void);    
  int   revertToStart(void);
  
  SectionForceDeformation *getCopy(void);
  const ID &getType(void);
  int getOrder(void) const;
  
  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag = 0);
  
 protected:
  
 private:
  
  NDMaterial **theFibers;
  double *yFibers;
  double *AFibers;

  Vector e;     // section trial deformations     
  
  double d;
  double tw;
  double bf;
  double tf;

  int nfdw;
  int nftf;
  
  double shapeFactor;

  static ID code;
  
  static Vector s;  // section resisting forces
  static Matrix ks; // section stiffness
};

#endif

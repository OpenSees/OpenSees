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
                                                                        
// $Revision: 1.6 $
// $Date: 2007-02-02 01:18:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionGJ.h,v $
                                                                        
// Written: fmk
// Created: 04/01
//
// Description: This file contains the class definition for 
// FiberSectionGJ.h. FiberSectionGJ provides the abstraction of a 
// 3d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

#ifndef FiberSectionGJ_h
#define FiberSectionGJ_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

class UniaxialMaterial;
class Fiber;
class Response;

class FiberSectionGJ : public SectionForceDeformation
{
 public:
  FiberSectionGJ(); 
  FiberSectionGJ(int tag, int numFibers, Fiber **fibers, double GJ = 1.0e10); 
  ~FiberSectionGJ();

    const char *getClassType(void) const {return "FiberSectionGJ";};
  
  int   setTrialSectionDeformation(const Vector &deforms); 
  const Vector &getSectionDeformation(void);
  
  const Vector &getStressResultant(void);
  const Matrix &getSectionTangent(void);
  const Matrix &getInitialTangent(void);
  
  int   commitState(void);
  int   revertToLastCommit(void);    
  int   revertToStart(void);
  
  SectionForceDeformation *getCopy(void);
  const ID &getType (void);
  int getOrder (void) const;
  
  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag = 0);
  
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &info);
  
  int addFiber(Fiber &theFiber);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  
 protected:
  
 private:
  int numFibers;                   // number of fibers in the section
  UniaxialMaterial **theMaterials; // array of pointers to materials
  double *matData;               // data for the materials [yloc and area]
  double kData[6];               // data for ks matrix 
  double sData[3];               // data for s vector 
  
  double yBar;       // Section centroid
  double zBar;
  
  Vector e;          // trial section deformations 
  
  static ID code;
  static Vector s;         // section resisting forces
  static Matrix ks;        // section stiffness
  
  double GJ;
};

#endif

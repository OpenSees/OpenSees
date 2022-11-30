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
// $Date: 2019-01-28 17:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/damping/Damping.h,v $

// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 02/2020
// Revision: A
//
// Description: This file contains the definition for the Damping class.
// Damping provides the abstraction of an elemental damping imposition.
// It is an abstract base class and thus no objects of it's type can be instatiated
// It has pure virtual functions which  must be implemented in it's derived classes.
//
// Reference:
// Yuan Tian, Yuli Huang, Zhe Qu, Yifan Fei, Xinzheng Lu,
// High-performance uniform damping model for response history analysis in OpenSees,
// Journal of Earthquake Engineering,
// 2022,
// https://doi.org/10.1080/13632469.2022.2124557
//
// What: "@(#) Damping.h, revA"

#ifndef Damping_h
#define Damping_h

#include <MovableObject.h>
#include <TaggedObject.h>

class Vector;

// class definition

class Damping: public TaggedObject, public MovableObject
{
public:
  Damping(int tag, int classTag);
  Damping();
  virtual ~Damping();

  virtual Damping *getCopy(void) = 0;
  virtual int setDomain(Domain *domain, int nComp) = 0;
  virtual int update(Vector q) = 0;
  virtual int commitState(void) = 0;
  virtual int revertToLastCommit(void) = 0;
  virtual int revertToStart(void) = 0;
  virtual const Vector &getDampingForce(void) = 0;
  virtual double getStiffnessMultiplier(void) = 0;
  

protected:

private:
};

// some additional methods related to prototypes created for copy constructors
extern bool     OPS_addDamping(Damping *newComponent);
extern Damping *OPS_getDamping(int tag);
extern bool     OPS_removeDamping(int tag);
extern void     OPS_clearAllDamping(void);
extern void     OPS_printDamping(OPS_Stream &s, int flag=0);

#endif

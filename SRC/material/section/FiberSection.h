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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSection.h,v $
                                                                        
// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the class definition for 
// FiberSection.h. FiberSection provides the abstraction of a 
// section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

#ifndef FiberSection_h
#define FiberSection_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

class Fiber;

class FiberSection : public SectionForceDeformation
{
  public:
    FiberSection(); 
    FiberSection(int tag, int estNumFibers = 8); 
    FiberSection(int tag, int numFibers, Fiber **fibers); 
    ~FiberSection();

    int   setTrialSectionDeformation(const Vector &deforms); 
    const Vector &getSectionDeformation(void);

    const Vector &getStressResultant(void);
    const Matrix &getSectionTangent(void);

    int   commitState(void);
    int   revertToLastCommit(void);    
    int   revertToStart(void);
 
    SectionForceDeformation *getCopy(void);
    const ID &getType (void) const;
    int getOrder (void) const;
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    void Print(ostream &s, int flag = 0);
	    
	int addFiber(Fiber &theFiber);

  protected:
    
  private:
    int numFibers;  // number of fibers in the section
    Fiber **theFibers;   // array of pointers to fibers
                                 // that form the section
    int sizeFibers;              // size of the fibers array
    
	int order;
	ID *code;

    Vector *e;       // section trial deformations 
	Vector *eCommit;
    Vector *s;       // section resisting forces  (axial force, bending moment)
    Matrix *ks;       // section stiffness

	int otherDbTag;
};

#endif

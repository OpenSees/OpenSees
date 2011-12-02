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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/GenericSection1d.h,v $
                                                                        
                                                                        
#ifndef GenericSection1d_h
#define GenericSection1d_h

// File: ~/material/GenericSection1d.h
//
// Written: MHS
// Created: Apr 2000
// Revision: A
//
// Description: This file contains the class definition for GenericSection1d.
//
// What: "@(#) GenericSection1d.h, revA"

#include <SectionForceDeformation.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

class Information;
class UniaxialMaterial;

class GenericSection1d : public SectionForceDeformation
{
  public:
    GenericSection1d (int tag, UniaxialMaterial& m, int code);
    GenericSection1d ();

    ~GenericSection1d ();

    int setTrialSectionDeformation (const Vector&);
    const Vector &getSectionDeformation (void);

    const Vector &getStressResultant (void);
    const Matrix &getSectionTangent (void);
    const Matrix &getSectionFlexibility (void);
    
    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);
    
    SectionForceDeformation *getCopy (void);
    const ID &getType (void) const;
    int getOrder (void) const;
    
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel,
			  FEM_ObjectBroker &theBroker);
    
    void Print (ostream &s, int flag = 0);
    
  protected:
    UniaxialMaterial *theModel;
    int code;

  private:
};


#endif

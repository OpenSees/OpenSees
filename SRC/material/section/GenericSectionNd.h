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
// $Date: 2006-08-03 23:49:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/GenericSectionNd.h,v $
                                                                        
                                                                        
#ifndef GenericSectionNd_h
#define GenericSectionNd_h

// Written: MHS
// Created: Apr 2000
// Revision: A
//
// Description: This file contains the class definition for GenericSectionNd.
//
// What: "@(#) GenericSectionNd.h, revA"

#include <SectionForceDeformation.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

class Information;
class NDMaterial;

class GenericSectionNd : public SectionForceDeformation
{
  public:
    GenericSectionNd (int tag, NDMaterial& m, const ID &mCode);
    GenericSectionNd ();

    ~GenericSectionNd ();

    const char *getClassType(void) const {return "GenericSectionNd";};

    int setTrialSectionDeformation (const Vector&);
    const Vector &getSectionDeformation (void);
    const Vector &getStressResultant (void);
    const Matrix &getSectionTangent (void);
    
    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);
    
    SectionForceDeformation *getCopy (void);
    const ID &getType (void);
    int getOrder (void) const;

    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel,
			  FEM_ObjectBroker &theBroker);

    void Print (OPS_Stream &s, int flag = 0);

  protected:
    NDMaterial *theModel;
    ID *code;
	int order;
	int otherDbTag;

  private:
};


#endif

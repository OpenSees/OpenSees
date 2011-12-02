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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/oldElasticSection2d.h,v $
                                                                        
                                                                        
///////////////////////////////////////////////////////
// File:  ~/Src/element/hinge/ElasticSection2d.h
//
// Written by Matthew Peavy
//
// Written:	 Feb 13, 2000
// Debugged: Feb 14, 2000
// Revised:        , 200x
//
//
// Purpose:  This header file contains the prototype
// for the ElasticSection2d class.

#ifndef ElasticSection2d_h
#define ElasticSection2d_h

#define SECTION_TAG_ElasticSection2d  501	//Update this

#include <SectionForceDeformation.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Channel;
class FEM_ObjectBroker;
class Information;

class ElasticSection2d: public SectionForceDeformation
{
  public:	
    ElasticSection2d (void);    
    ElasticSection2d (int tag, double E, double I, double A, 
		      double massDensityPerUnitLength=0);
    ~ElasticSection2d (void);
    
    int commitState (void);
    int revertToLastCommit (void);

    void setTrialDeformation (const Vector&);
    const Vector &getDeformation (void);
    
    const Vector &getResistingForce (void);
    const Vector &getPrevResistingForce (void);
    
    const Matrix &getTangentStiff (void);
    const Matrix &getPrevTangentStiff (void);
    
    const Matrix &getFlexMatrix (void);
    const Matrix &getPrevFlexMatrix (void);
    
    SectionForceDeformation *getCopy (void);
    const ID &getType (void) const;
    int getOrder (void) const;
    
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel,
		  FEM_ObjectBroker &theBroker);
    
    void Print (ostream &s, int flag =0);
    friend ostream &operator<< (ostream &s, ElasticSection2d &H);
    
  protected:
    
  private:
    
    double A, E, I;
    double massDens;
    
    Matrix k;			// section stiffness matrix
    Matrix f;			// section flexibility matrix
    Vector e;			// section trial deformations
    Vector s;			// section resisting forces
		
    static ID code;
};

#endif

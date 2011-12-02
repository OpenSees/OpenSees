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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-25 23:33:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/Bidirectional.h,v $
                                                                        
                                                                        
#ifndef Bidirectional_h
#define Bidirectional_h

// File: ~/material/Bidirectional.h
//
// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: 
//
// What: "@(#) Bidirectional.h, revA"

#include <SectionForceDeformation.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Bidirectional : public SectionForceDeformation
{
  public:
    Bidirectional(int tag, double E, double sigY, double Hiso, double Hkin);
    Bidirectional();
    ~Bidirectional();

    int setTrialSectionDeformation(const Vector &v);
    const Matrix &getSectionTangent(void);
    const Matrix &getInitialTangent(void);
    const Vector &getStressResultant(void);
    const Vector &getSectionDeformation(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    
    SectionForceDeformation *getCopy(void);
    const ID &getType(void);
    int getOrder(void) const;
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

  protected:

  private:
	double E;
	double sigY;
	double Hiso;
	double Hkin;

	double e_n1[2];
	double eP_n[2];
	double eP_n1[2];

	double q_n[2];
	double q_n1[2];

	double alpha_n;
	double alpha_n1;
	
	static Vector s;
	static Matrix ks;
	static ID code;
};


#endif

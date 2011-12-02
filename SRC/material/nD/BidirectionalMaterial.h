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
// $Date: 2000-12-18 09:49:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BidirectionalMaterial.h,v $
                                                                        
                                                                        
#ifndef BidirectionalMaterial_h
#define BidirectionalMaterial_h

// File: ~/material/BidirectionalMaterial.h
//
// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: 
//
// What: "@(#) BidirectionalMaterial.h, revA"

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <Tensor.h>

#define ND_TAG_Bidirectional 1234

class BidirectionalMaterial : public NDMaterial
{
  public:
    BidirectionalMaterial(int tag, double E, double sy, double hi, double hk);
    BidirectionalMaterial();
    ~BidirectionalMaterial();

    int setTrialStrain(const Vector &v);
    int setTrialStrain(const Vector &v, const Vector &r);
    int setTrialStrainIncr(const Vector &v);
    int setTrialStrainIncr(const Vector &v, const Vector &r);
    const Matrix &getTangent(void);
    const Vector &getStress(void);
	const Vector &getStrain(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    
    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *code);
    const char *getType(void) const;
    int getOrder(void) const;
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

	void Print(ostream &s, int flag = 0);

  protected:

  private:
	double E;
	double sigY;
	double Hi;
	double Hk;

	double epsn1[2];
	double epsPn[2];
	double epsPn1[2];

	double backn[2];
	double backn1[2];

	double effn;
	double effn1;
	
	static Vector sigma;
	static Matrix D;
};


#endif

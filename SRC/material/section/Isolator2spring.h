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
// $Date: 2006-01-17 20:44:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/Isolator2spring.h,v $
                                                                        
                                                                        
// Written: K. Ryan
// Created: September 2003
// Updates: November 2005
//
// Description: This file contains the class implementation for a "two-spring isolator" 
// material.  This material is based on the two-spring model originally developed by 
// Koh and Kelly to represent the buckling behavior of an elastomeric bearing.  The 
// material model has been modified to include material nonlinearity and optional 
// strength degradation.

#ifndef Isolator2spring_h
#define Isolator2spring_h
#include <SectionForceDeformation.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Isolator2spring : public SectionForceDeformation
{
  public:
    Isolator2spring(int tag, double tol_in, double k1_in, double Fy_in, double kb_in, double kvo_in, 
	       double hb_in, double Pe_in, double po_in);
    Isolator2spring();
    ~Isolator2spring();

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
        double tol;
	double k1;
	double Fyo;
	double kbo;
	double kvo;
	double h;
	double Pe;
	double po;

	double utpt[2];
	double sP_n;
	double sP_n1;

	double q_n;
	double q_n1;

	double H;
	double pcr;

	Vector x0;
	Matrix ks;
	static Vector f0;
	static Matrix df;
	static Vector s;
	static Vector s3;
	static ID code;
};

#endif

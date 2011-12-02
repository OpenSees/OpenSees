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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:01:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/DrainMaterial.h,v $
                                                                      
// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// DrainMaterial. DrainMaterial wraps a Drain spring element subroutine
// and converts it to the UniaxialMaterial interface for use in
// zero length elements, beam sections, or anywhere else
// UniaxialMaterials may be used.
//
// For more information see the Drain-2DX user guide:
//    Allahabadi, R.; Powell, G. H.
//    UCB/EERC-88/06, Berkeley: Earthquake Engineering Research Center,
//    University of California, Mar. 1988, 1 vol.

#ifndef DrainMaterial_h
#define DrainMaterial_h

#include <UniaxialMaterial.h>

class DrainMaterial : public UniaxialMaterial
{
  public:
    DrainMaterial(int tag, int classTag, int numHV, int numData, double beto = 0.0);
    virtual ~DrainMaterial();

    virtual int setTrialStrain(double strain, double strainRate = 0.0);
    virtual int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0); 
    virtual double getStrain(void);
    virtual double getStrainRate(void);
    virtual double getStress(void);
    virtual double getTangent(void);
    virtual double getDampTangent(void);
    virtual double getInitialTangent(void);

    virtual int commitState(void);
    virtual int revertToLastCommit(void);    
    virtual int revertToStart(void);        

    // WARNING -- if you wish to override any method in this base class, you must
    // also override the getCopy method to return a pointer to the derived class!!!
    virtual UniaxialMaterial *getCopy(void);

    virtual int sendSelf(int commitTag, Channel &theChannel);  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);    
    
    virtual void Print(OPS_Stream &s, int flag = 0);

  protected:
	// Invokes the Drain subroutine
	virtual int invokeSubroutine(void);
	
	double *data;		// Material parameters array
	double *hstv;		// History array: first half is committed, second half is trial

	int numData;		// Number of material parameters
	int numHstv;		// Number of history variables

	double epsilonP;	// Committed strain
	double sigmaP;		// Committed stress
	double tangentP;	// Committed tangent

	double beto;		// Stiffness proportional damping factor
	double initialTangent;  // initial tangent

  private:
	double epsilon;		// Trial strain
	double epsilonDot;	// Trial strain rate
	double sigma;		// Trial stress
	double tangent;		// Trial tangent
};


#endif


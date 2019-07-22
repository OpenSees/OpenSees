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
// $Date: 2006-08-04 18:18:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticOrthotropicPlaneStress2D.h,v $

#ifndef ElasticOrthotropicThreeDimensional_h
#define ElasticOrthotropicThreeDimensional_h

// Written: fmk
// Created: 10/11
//
// Description: 
//
// What: "@(#) ElasticOrthotropicThreeDimesnional.h, revA"

#include <ElasticOrthotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class ElasticOrthotropicThreeDimensional : public ElasticOrthotropicMaterial
{
  public:
    ElasticOrthotropicThreeDimensional(int tag, double Ex, double Ey, 
double Ez, double vxy, double vyz, double vzx, double Gxy, double Gyz, 
double Gzx, double rho = 0.0);
    ElasticOrthotropicThreeDimensional();
    ~ElasticOrthotropicThreeDimensional();

    const char *getClassType(void) const {return "ElasticOrthotropicThreeDimensional";};

    int setTrialStrain (const Vector &v);
    int setTrialStrain (const Vector &v, const Vector &r);
    int setTrialStrainIncr (const Vector &v);
    int setTrialStrainIncr (const Vector &v, const Vector &r);
    const Matrix &getTangent (void);
    const Matrix &getInitialTangent (void);
    
    const Vector &getStress (void);
    const Vector &getStrain (void);
    
    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);
    
    NDMaterial *getCopy (void);
    const char *getType (void) const;
    int getOrder (void) const;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    const Vector& getStressSensitivity(int gradIndex,
				       bool conditional);

 protected:

  private:
    static Vector sigma;	// Stress vector ... class-wide for returns
    static Matrix D;		// Elastic constants
    Vector epsilon;	        // Trial strains
    Vector Cepsilon;	        // Committed strain
};

#endif

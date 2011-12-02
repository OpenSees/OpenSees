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
// $Date: 2002-12-05 22:49:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlaneStress2D.h,v $
                                                                        
                                                                        
#ifndef ElasticIsotropicPlaneStress2D_h
#define ElasticIsotropicPlaneStress2D_h

// File: ~/material/ElasticIsotropicPlaneStress2D.h
//
// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: 
//
// What: "@(#) ElasticIsotropicPlaneStress2D.h, revA"

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <Tensor.h>

class ElasticIsotropicPlaneStress2D : public ElasticIsotropicMaterial
{
  public:
    ElasticIsotropicPlaneStress2D (int tag, double E, double nu, double rho);
    ElasticIsotropicPlaneStress2D ();
    ~ElasticIsotropicPlaneStress2D ();

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
    
  protected:

  private:
    static Vector sigma;	// Stress vector ... class-wide for returns
    static Matrix D;		// Elastic constants
    Vector epsilon;	        // Trial strains
};


#endif

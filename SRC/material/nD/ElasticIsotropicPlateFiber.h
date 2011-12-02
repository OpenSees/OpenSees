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
// $Date: 2002-12-05 22:49:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlateFiber.h,v $
                                                                        
                                                                        
#ifndef ElasticIsotropicPlateFiber_h
#define ElasticIsotropicPlateFiber_h

// File: ~/material/ElasticIsotropicPlateFiber.h
//
// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: 
//
// What: "@(#) ElasticIsotropicPlateFiber.h, revA"

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class ElasticIsotropicPlateFiber : public ElasticIsotropicMaterial
{
  public:
    ElasticIsotropicPlateFiber (int tag, double E, double nu, double rho);
    ElasticIsotropicPlateFiber ();
    ~ElasticIsotropicPlateFiber ();

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
    Vector epsilon;		// Trial strains
};


#endif

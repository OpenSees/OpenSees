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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlaneStress2D.h,v $

#ifndef PlasticDamageConcrete3d_h
#define PlasticDamageConcrete3d_h

// Written: Thanh Do
// Created: 07/16
//
// Description: 
//
// What: "@(#) ElasticIsotropicThreeDimesnional.h, revA"

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class PlasticDamageConcrete3d : public NDMaterial
{
  public:
  PlasticDamageConcrete3d(int tag, 
			  double E, 
			  double nu, 
			  double ft,
			  double fc, 
			  double beta = 0.6, 
			  double Ap = 0.5, 
			  double An = 2.0, 
			  double Bn = 0.75);
    PlasticDamageConcrete3d();
    ~PlasticDamageConcrete3d();

    const char *getClassType(void) const {return "PlasticDamageConcrete3d";};

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

    NDMaterial*getCopy(const char *type);
    NDMaterial *getCopy (void);
    const char *getType (void) const;
    int getOrder (void) const;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    void Print(OPS_Stream &s, int flag =0);       

  protected:

  private:
    // parameters
    double E;     // elastic modulus
    double nu;    // Poisson ratio 
    double ft;    // tensile yield strength
    double fc;    // compressive yield strength
    double beta;  // plastic deformation rate
    double Ap;    // damage parameter
    double An;    // damage parameter
    double Bn;    // damage parameter

    // current state variables
    double rp;    // positive damage threshold
    double rn;    // negative damage threshold
    double dp;    // positive damage variable
    double dn;    // negative damage variable

    Vector eps;   // strain
    Vector sig;   // stress
    Vector sige;  // effective stress
    Vector eps_p; // plastic strain
    Vector sigeP; // effective stress

    // committed state variables
    double rpCommit; 
    double rnCommit; 
    double dpCommit; 
    double dnCommit; 

    Vector epsCommit;
    Vector sigCommit;
    Vector sigeCommit;  
    Vector eps_pCommit; 
    Vector sigePCommit; 

    // tangent matrices
    Matrix Ce; 
    Matrix C; 
    Matrix Ccommit; 
};

#endif

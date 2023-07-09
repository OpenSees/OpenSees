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

#ifndef ContactMaterial2D_h
#define ContactMaterial2D_h

// $Revision: 1.2
// $Date: 2010-11-10
// $Source: /OpenSees/SRC/material/nD/ContactMaterial2D.h,v $
                                                                        
// Written: Kathryn Petek
// Created: February 2004
// Modified: Chris McGann
//           November 2010 -> changes for incorporation into main source code
// Modified: Chris McGann
//           Jan 2011 -> added update for frictional state

// Description: This file contains the class prototype for ContactMaterial2D.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class ContactMaterial2D : public NDMaterial
{
public:
    // Full constructor
    ContactMaterial2D (int tag, double mu, double G, double c, double t);

    // Null constructor
    ContactMaterial2D ();

    // Destructor: clean up memory storage space.
    ~ContactMaterial2D ();

    // Sets the values of the trial strain tensor.
    int setTrialStrain (const Vector &strain_from_element);

    // Unused trialStrain functions
    int setTrialStrain(const Vector &v, const Vector &r);

    // Calculates current tangent stiffness.
    const Matrix &getTangent ();
    const Matrix &getInitialTangent ();

    // Calculates the corresponding stress increment (rate),
    // for a given strain increment. 
    const Vector &getStress ();
    const Vector &getStrain ();

    //Get cohesion function for use in contact element
    double getcohesion ();
    void ScaleCohesion (const double len);

    //Get tensile strength function for use in contact element
    double getTensileStrength ();
    void ScaleTensileStrength (const double len);

    // Accepts the current trial strain values as being on the
    // solution path, and updates all model parameters related
    // to stress/strain states. Return 0 on success.
    int commitState ();

    // Revert the stress/strain states to the last committed states.
    // Return 0 on success.
    int revertToLastCommit ();

    int revertToStart();

    // Return an exact copy of itself.
    NDMaterial *getCopy (void);
    
    // Return a copy of itself if "code"="ContactMaterial2D", 
    // otherwise return null.
    NDMaterial *getCopy (const char *code);
    
    // Return the string "ContactMaterial2D".
    const char *getType () const ;

    // Return ndm.
    int getOrder () const ;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		  FEM_ObjectBroker &theBroker);    

    void Print(OPS_Stream &s, int flag =0);

	// public methods for material stage update
	int setParameter(const char **argv, int argc, Parameter &param);
	int updateParameter(int responseID, Information &eleInformation);

protected:

    //material parameters
    double frictionCoeff;   // friction coefficient tan(phi)
    double stiffness;       // shear stiffness of the interface
    double cohesion;        // interface cohesion (force)
	double tensileStrength;  // tensile strength

    // functions
    void zero();
	int UpdateFrictionalState(void);

	// variables for update of friction coefficient
	static int mFrictFlag;
	int mFlag;
	double mMu;
	double mCo;
	double mTen;

private:

    // state variables
    double s_e_n;           // elastic slip from previous increment
    double s_e_nplus1;      // elastic slip after current increment
    
    double r_nplus1;        // sliding direction

    bool inSlip;            // sliding indicator

    // static vectors and matrices
    Vector strain_vec;      // generalized strain state
        // strain_vec(0) -> gap     ... current gap distance
        // strain_vec(1) -> slip    ... incremental slip
        // strain_vec(2) -> lambda  ... lagrangean multiplier -> t_n

    Vector stress_vec;      // generalized stress state
        // stress_vec(0) -> t_n     ... normal contact force
        // stress_vec(1) -> t_s     ... tangential contact force
        // stress_vec(2) -> gap     ... current gap

    Matrix tangent_matrix;  // material tangent

};

#endif


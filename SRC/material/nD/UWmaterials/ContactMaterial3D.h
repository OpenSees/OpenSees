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

#ifndef ContactMaterial3D_h
#define ContactMaterial3D_h

// $Revision: 1.2
// $Date: 2010-11-10
// $Source: /OpenSees/SRC/material/nD/ContactMaterial3D.h,v $
                                                                        
// Written: Kathryn Petek
// Created: February 2004
// Modified: Chris McGann
//           November 2010 -> changes for incorporation into main source code
// Modified: Chris McGann
//           Jan 2011 -> added update for frictional state

// Description: This file contains the class prototype for ContactMaterial3D.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class ContactMaterial3D : public NDMaterial
{
  public:
  // Full constructor
  ContactMaterial3D (int tag, double mu, double G, double c, double t);
  ContactMaterial3D (const ContactMaterial3D &);
  
  // Null constructor
  ContactMaterial3D ();
  
  // Destructor: clean up memory storage space.
  ~ContactMaterial3D ();
  
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
  
  // get metric tensor for material class
  void setMetricTensor(Matrix &m);
  bool getContactState() {return inSlip;}
  
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
  
  // Return a copy of itself if "code"="ContactMaterial3D", 
  // otherwise return null.
  NDMaterial *getCopy (const char *code);
  
  // Return the string "ContactMaterial3D".
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
  
  double frictionCoeff;
  double stiffness;
  double cohesion;        // interface cohesion (force)
  double tensileStrength;  // interface tensile strength (force)
  
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
  Vector s_e_n;           // elastic slip from previous increment
  Vector s_e_nplus1;      // elastic slip after current increment
  
  Vector r_nplus1;        // sliding direction
  
  double gamma;           // consistency parameter
  double s_e_nplus1_norm;	// norm of trial slip
  
  bool inSlip;            // sliding indicator
  
  Matrix g;				// metric tensor
  Matrix G;				// dual basis metric tensor
  
  // static vectors and matrices
  Vector strain_vec;      // generalized strain state
  // strain_vec(0) -> gap     ... current gap distance
  // strain_vec(1) -> slip(1) ... incremental slip in xi direction
  // strain_vec(2) -> slip(2) ... incremental slip in eta direciton
  // strain_vec(3) -> lambda  ... lagrangean multiplier -> t_n
  
  Vector stress_vec;      // generalized stress state
  // stress_vec(0) -> t_n     ... normal contact force
  // stress_vec(1) -> t_s(1)  ... tangential contact force in xi dir.
  // stress_vec(2) -> t_s(2) ... tangentail contact force in eta dir.
  // stress_vec(3) -> gap     ... current gap
  
  Matrix tangent_matrix;  // material tangent
  
};

#endif


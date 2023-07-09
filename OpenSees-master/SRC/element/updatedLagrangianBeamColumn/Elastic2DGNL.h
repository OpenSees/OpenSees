/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** See file 'COPYRIGHT'  in main directory for information on usage   **
** and redistribution of OpenSees, and for a DISCLAIMER OF ALL        **
** WARRANTIES.                                                        **
**                                                                    **
** Elastic2dGNL.h: interface for the Elastic2dGNL class               **
** Developed by:                                                      **
**    Rohit Kaul       (rkaul@stanford.edu)                           **
**    Greg Deierlein   (ggd@stanford.edu)                             **
**                                                                    **
**           John A. Blume Earthquake Engineering Center              **
**                    Stanford University                             **
** ****************************************************************** **/

// Written: rkaul
// Created: 7/30
//
// Description: This file contains the class definition for Elastic2dGNL.

// Elastic2dGNL is a subclass of UpdatedLagrangianBeam2D, that can be 
// used to model 2d beam column elements with large deformation effects. 
// Most of the virtual Element methods have been implemented by the parent 
// class (see UpdatedLagrangianBeam2D.h for details, 
// including the recorder arguments)

#ifndef Elastic2dGNL_H
#define Elastic2dGNL_H

// List of included files
// UpdatedLagrangianBeam2D - parent class of this class
#include "UpdatedLagrangianBeam2D.h"

class Elastic2dGNL : public UpdatedLagrangianBeam2D  
{
 public:	
  // Arguments passed to the constructor include - tag, a unique id in 
  // the domain, A - cross section area of the beam, E - modulus of
  // Elasticity, I - Izz, Nd1 and Nd2 are the nodal numbers of connected
  // nodes.  Last two parameters are optional, by default geometric
  // nonlinearity is turned on and rho, mass density is set to zero.
  // This element assumes lumped mass for mass matrix.
  Elastic2dGNL(int tag, double A, double E, double I, int Nd1, int Nd2, 
	       bool islinear = false, double rho = 0.0);
  virtual ~Elastic2dGNL();
  
  // Prints the element info to OPS_Stream
  void Print(OPS_Stream &s, int flag =0);
  
  // Methods for sending and receiving the object over a channel
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  
  
 protected:
  // Implementation of pure virtual subclass methods

  // Get the elastic stiffness in local coordinates,
  // stored into K
  void getLocalStiff(Matrix &K);
  // Get the mass matrix in local coordinate system,
  // stored in M
  void getLocalMass(Matrix &M);
  
 private:
  // Data declarations
  double A, E, Iz;
};

#endif // !defined Elastic2dGNL

/*
 *
 * WARNING/ERROR format/representation
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Warnings are generally issued if the standard implementation of a
 * method may not be possible.  Further analysis may lead to incorrect
 * solution, even though it may be possible to continue.
 *           
 * Warnings
 * --------
 *                
 * WARNING (W_Level_ID) - Class::method(..) [ElementTag]
 * Short description ...
 *
 * Errors
 * ------
 *  
 * ERROR (E_ID) - Class::method(..) [ElementTag]
 * Short description ...
 *           
 * Analysis may be halted if an error is encountered, further
 * analysis will definitely be erroneous.
 *
 * (see UpdatedLagrangianBeam2D.h for details)
 */



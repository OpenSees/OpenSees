// $Revision: 1.1 $
// $Date: 
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ContactMaterial3D.h,v $
                                                                        
// Written: Kathryn Petek
// Created: February 2004

// Description: This file contains the class prototype for ContactMaterial3D.
//
// What: "@(#) ContactMaterial3D.h, revA"

#ifndef ContactMaterial3D_h
#define ContactMaterial3D_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <Tensor.h>

class ContactMaterial3D : public NDMaterial
{
public:
    // Full constructor
    ContactMaterial3D (int tag, double mu, double G, double c, double t);

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
    double getcohesion ()  {return cohesion;}

	//Get tensile strength function for use in contact element
	double getTensileStrength() {return tensileStrength;}

    void ScaleCohesion (const double len) {cohesion *= len;}
	void ScaleTensileStrength (const double len) {tensileStrength *= len;}

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

	int updateParameter(int responseID, Information &eleInformation);

protected:

    //material parameters

	static int matCount;
	int matN;

    static double* frictionCoeffx;   // friction coefficient tan(phi)
    static double* stiffnessx;       // shear stiffness of the interface
    
	double frictionCoeff;
	double stiffness;
	double cohesion;        // interface cohesion (force)
	double tensileStrength;  // interface tensile strength (force)

    void zero();

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

	int MyTag;                  // material tag for debugging

};

#endif


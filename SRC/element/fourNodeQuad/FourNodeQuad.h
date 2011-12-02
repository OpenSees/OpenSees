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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/FourNodeQuad.h,v $
                                                                        
                                                                        
// File: ~/element/FourNodeQuad.h
//
// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class definition for FourNodeQuad.
//
// What: "@(#) FourNodeQuad.h, revA"

#ifndef FourNodeQuad_h
#define FourNodeQuad_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class NDMaterial;
class QuadRule1d;

class FourNodeQuad : public Element
{
  public:
    FourNodeQuad (int tag, int nd1, int nd2, int nd3, int nd4,
		  NDMaterial &m, const char *type,
		  double t, double pressure = 0.0, double rho = 0.0,
		  double b1 = 0.0, double b2 = 0.0);
    FourNodeQuad ();
    virtual ~FourNodeQuad();

    int getNumExternalNodes () const;
    const ID &getExternalNodes ();
    int getNumDOF ();	
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState ();
    int revertToLastCommit ();
    int revertToStart ();

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff ();
    const Matrix &getSecantStiff ();    
    const Matrix &getDamp ();    
    const Matrix &getMass ();    

    void zeroLoad ();
    int addLoad(const Vector &addLoad);
    int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector &getResistingForce ();
    const Vector &getResistingForceIncInertia ();            

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    int displaySelf (Renderer &theViewer, int displayMode, float fact);
    void Print (ostream &s, int flag =0);

    int setResponse (char **argv, int argc, Information &eleInformation);
    int getResponse (int responseID, Information &eleInformation);

    int setParameter(char **argv, int argc, Information &info);
    int updateParameter(int parameterID, Information &info);

  protected:
    
  private:
    // private attributes - a copy for each object of the class

    NDMaterial * **theMaterial; // pointer to the ND material objects
    
    ID  connectedExternalNodes; // Tags of quad nodes

    Node *nd1Ptr;		// Pointers to quad nodes
    Node *nd2Ptr;
    Node *nd3Ptr;
    Node *nd4Ptr;

    Matrix K;		// Element stiffness Matrix
    Matrix C;		// Element damping matrix
    Matrix M;		// Element mass matrix
    Vector P;		// Element resisting force vector
	Vector Q;		// Applied nodal loads
	Vector b;		// Body forces
	Vector pressureLoad;	// Pressure load at nodes

    double thickness;	// Element thickness
    double rho;			// Mass per unit volume
	double pressure;	// Normal surface traction (pressure) over entire element
						// Note: positive for outward normal

    QuadRule1d *theQuadRule;	// Quadrature rule
	int order;					// Order of the quadrature rule

    Matrix J;		// Jacobian of transformation
    Matrix L;		// Inverse of J
    Matrix B;		// Strain interpolation matrix
	
    // static data - single copy for all objects of the class	    
    static Matrix N;		// Displacement interpolation matrix
	
    // private member functions - only objects of this class can call these
    void setJacobian (double xi, double eta);
    double formDetJ (double xi, double eta);
    void formBMatrix (double xi, double eta);
    static void formNMatrix (double xi, double eta);
	void getMaterialIndices(int pointNum, int &i, int &j);
	void setPressureLoadAtNodes(void);
};


#endif


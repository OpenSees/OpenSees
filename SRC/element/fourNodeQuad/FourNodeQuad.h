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
// $Date: 2000-12-18 10:40:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/FourNodeQuad.h,v $
                                                                        
// Written: MHS
// Created: Feb 2000
// Revised: Dec 2000 for efficiency
//
// Description: This file contains the class definition for FourNodeQuad.

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
class Response;

class FourNodeQuad : public Element
{
  public:
    FourNodeQuad(int tag, int nd1, int nd2, int nd3, int nd4,
		  NDMaterial &m, const char *type,
		  double t, double pressure = 0.0, double rho = 0.0,
		  double b1 = 0.0, double b2 = 0.0);
    FourNodeQuad();
    virtual ~FourNodeQuad();

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getSecantStiff(void);
    const Matrix &getDamp(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(const Vector &addLoad);
    int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);
    void Print(ostream &s, int flag =0);

	Response *setResponse(char **argv, int argc, Information &eleInformation);
    int getResponse(int responseID, Information &eleInformation);

    int setParameter(char **argv, int argc, Information &info);
    int updateParameter(int parameterID, Information &info);

  protected:
    
  private:
    // private attributes - a copy for each object of the class

    NDMaterial **theMaterial; // pointer to the ND material objects
    
    ID connectedExternalNodes; // Tags of quad nodes

    Node *nd1Ptr;		// Pointers to quad nodes
    Node *nd2Ptr;
    Node *nd3Ptr;
    Node *nd4Ptr;

    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector
	Vector Q;		// Applied nodal loads
	double b[2];		// Body forces
	Vector pressureLoad;	// Pressure load at nodes

    double thickness;	// Element thickness
    double rho;			// Mass per unit volume
	double pressure;	// Normal surface traction (pressure) over entire element
						// Note: positive for outward normal

	static double shp[3][4];	// Stores shape functions and derivatives (overwritten)
	static double pts[4][2];	// Stores quadrature points
	static double wts[4];		// Stores quadrature weights

    // private member functions - only objects of this class can call these
    double shapeFunction(double xi, double eta);
	void setPressureLoadAtNodes(void);
};

#endif


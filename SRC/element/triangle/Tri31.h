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
                                                                        
// $Revision: 1.00 $
// $Date: 2010/09/08 20:01:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/triangular/Tri31.h,v $
                                                                        
// Written: Roozbeh Geraili Mikola (roozbehg@berkeley.edu)
// Created: Sep 2010
// Revised: --------
//
// Description: This file contains the class definition for Tri31.

#ifndef Tri31_h
#define Tri31_h

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

class Tri31 : public Element
{
  public:
    Tri31(int tag, int nd1, int nd2, int nd3,
	  NDMaterial &m, const char *type,
	  double t, double pressure = 0.0, 
	  double rho = 0.0,
	  double b1 = 0.0, double b2 = 0.0);
    Tri31();
    ~Tri31();

    const char *getClassType(void) const {return "Tri31";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    
    const Matrix &getMass(void);    

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.
    friend class PyLiq1;
    friend class TzLiq1;

  protected:
    
  private:
    // private attributes - a copy for each object of the class

    NDMaterial **theMaterial; // pointer to the ND material objects
    
    ID connectedExternalNodes; // Tags of Tri31 nodes

    Node *theNodes[3];

    static double matrixData[36];  // array data for matrix
    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector
    Vector Q;		        // Applied nodal loads
    double b[2];		// Body forces

	double appliedB[2]; // Body forces applied with load pattern
	int applyLoad;      // flag for body force in load

    Vector pressureLoad;	// Pressure load at nodes

    double thickness;	        // Element thickness
    double pressure;	        // Normal surface traction (pressure) over entire element
					 // Note: positive for outward normal
    double rho;
    static double shp[3][3];	// Stores shape functions and derivatives (overwritten)
    static double pts[1][2];	// Stores quadrature points
    static double wts[1];		// Stores quadrature weights

    // private member functions - only objects of this class can call these
    double shapeFunction(double xi, double eta);
    void setPressureLoadAtNodes(void);

    Matrix *Ki;

    static constexpr int numgp = 1; // number of gauss points
	static constexpr int numnodes = 3; // number of nodes
};

#endif

